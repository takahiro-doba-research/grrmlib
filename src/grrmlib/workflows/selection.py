from collections.abc import Iterable
from pathlib import Path

import polars as pl

from ..core import Molecules
from ..io import key_to_folder, read_xyz, read_xyzs, to_dataframe, write_xyzs
from ..operations import expand_derivatives, select_charges, select_lowest


HARTREE_IN_KCAL = 627.5095


def _read_results(
    folder: Path,
    basename: str,
    key_names: list[str],
    attrs: list[str],
    path_results: str | Path
) -> pl.DataFrame:
    # 計算の結果（タイトル行が JSON の xyz）を読み、全ての結果を path_results に CSV で書く
    mols = read_xyzs(folder, basename, ignore_errors=True)
    results = to_dataframe(mols, key_names, attrs).sort(key_names)
    results.write_csv(path_results)
    print(f"{basename}: {results.height} results, {results.filter(pl.col('success').ne_missing(True)).height} failed")
    return results


def select_sp(
    folder: str | Path,
    basename: str,
    path_results: str | Path,
    path_selected: str | Path,
    *,
    key_names: Iterable[str],
    group_by: Iterable[str]
) -> pl.DataFrame:
    """
    folder の下の SP の結果 basename を読み、group_by の列が同じものの中で最もエネルギーの低いものを選ぶ。

    - 結果は xyz で、タイトル行の JSON に charge, mult, scfenergy, success を書いたもの（read_xyz で読む）

    - path_results：全ての結果（key_names, scfenergy, success）
    - path_selected：選んだもの。計算に失敗したものは選ばない
    - 例：分離 EQ では group_by=["SG", "charge", "mult"]（SG ごとに SEQ を選ぶ）、
      置換基を導入した配座では角度を除いた列（置換基の組み合わせごとに角度を選ぶ）
    """
    key_names = list(key_names)
    sp = _read_results(Path(folder), basename, key_names, ["scfenergy", "success"], path_results)
    selected = select_lowest(sp.filter(pl.col("success")), list(group_by)).sort(key_names)
    selected.write_csv(path_selected)
    print(f"selected: {selected.height}")
    return selected


def make_opt(
    path_selected: str | Path,
    folder_sp: str | Path,
    folder_opt: str | Path,
    *,
    key_names: Iterable[str],
    basename: str = "init.xyz"
) -> int:
    """
    path_selected（select_sp の結果）の行の初期構造 basename を、folder_sp から folder_opt の同じ key のフォルダに写す。
    初期構造は key ごとに決まっているので、書き直しても中身は変わらない（overwrite=True）。書いた数を返す。
    """
    key_names = list(key_names)
    selected = pl.read_csv(path_selected)

    mols = Molecules()
    for key in selected.select(key_names).iter_rows():
        mols[key] = read_xyz(key_to_folder(folder_sp, key, key_names) / basename)

    write_xyzs(mols, folder_opt, basename, prefixes=key_names, record=True, overwrite=True)
    print(f"opt inputs: {len(mols)} written")
    return len(mols)


def select_opt(
    folder: str | Path,
    basename: str,
    path_table: str | Path,
    path_results: str | Path,
    path_selected: str | Path,
    *,
    key_names: Iterable[str],
    substituents: dict[str, list[int]] | None = None,
    total_charge: int = 0
) -> pl.DataFrame:
    """
    folder の下の最適化の結果 basename から、グループ（G）ごとに各分子の SG の電荷を選ぶ。
    分子の電荷の合計が total_charge になる組み合わせのうち、エネルギーの合計が最も低いもの（select_charges）。

    - 結果は xyz で、タイトル行の JSON に charge, mult, scfenergy, success, optsteps を書いたもの（read_xyz で読む）

    - path_table：separation の対応表（eq_list_to_xyzs の結果）。G の中で番号の最も小さい EQ の分子を使う
      （同じ G の EQ はどれも同じ SG の組に分かれる）
    - substituents を渡すと、置換基を導入した結果として、G と誘導体（"<名前>_derivative" の列）ごとに選ぶ
      （expand_derivatives。分子に含まれない置換基の結果はどの誘導体にも使う）
    - path_results：全ての結果（key_names, scfenergy, success, optsteps）。収束しなかったものは使わない
    - path_selected：選んだもの（G, 誘導体, fragment, key_names, scfenergy）
    - 電荷の選び方が 2 通り以上ある G は、選んだ電荷の組ごとに系の数と次の候補との差の最小値を表示する
    """
    key_names = list(key_names)
    opt = _read_results(Path(folder), basename, key_names, ["scfenergy", "success", "optsteps"], path_results)
    results = opt.filter(pl.col("success")).drop("success", "optsteps")

    derivative_keys = []
    if substituents is not None:
        results = expand_derivatives(results, substituents)
        derivative_keys = [f"{name}_derivative" for name in substituents]
    group_by = ["G", *derivative_keys]

    table = pl.read_csv(path_table)
    fragments = table.filter(pl.col("EQ") == pl.col("EQ").min().over("G")).select("G", "fragment", "SG")
    candidates = fragments.join(results, on="SG")

    # 全ての G と誘導体で、全ての分子に結果があることを確かめる
    n_fragments = candidates.group_by(group_by).agg(pl.col("fragment").n_unique()).join(
        fragments.group_by("G").len(), on="G"
    )
    n_expected = fragments["G"].n_unique() * (results.select(derivative_keys).n_unique() if derivative_keys else 1)
    if n_fragments.height != n_expected or n_fragments.filter(pl.col("fragment") != pl.col("len")).height:
        raise ValueError("some (G, derivative) have fragments without optimized structures")

    selected = select_charges(candidates, group_by)

    # 電荷の選び方が 2 通り以上ある G について、選んだ電荷の組を数える
    patterns = (
        selected
        .filter(pl.col("margin").is_not_null())
        .group_by(group_by, maintain_order=True)
        .agg(
            pl.when(pl.col("charge") > 0)
            .then(pl.format("SG={}(+{})", "SG", "charge"))
            .otherwise(pl.format("SG={}({})", "SG", "charge"))
            .filter(pl.col("charge") != 0).str.join(" ").alias("charged"),
            pl.col("margin").first(),
        )
        .group_by("G", "charged", maintain_order=True)
        .agg(pl.len(), (pl.col("margin").min() * HARTREE_IN_KCAL).alias("margin"))
    )
    for (g,), rows in patterns.group_by("G", maintain_order=True):
        text = ", ".join(
            f"{r['charged'] or 'neutral'} in {r['len']} ({r['margin']:.1f} kcal/mol or more)"
            for r in rows.iter_rows(named=True)
        )
        print(f"G={g}: {text}")

    selected = (
        selected
        .select(*group_by, "fragment", *key_names, "scfenergy")
        .sort(*group_by, "fragment")
    )
    selected.write_csv(path_selected)
    n_charged = selected.filter(pl.col("charge") != 0).select(group_by).n_unique()
    print(f"selected: {selected.height} rows in {selected.select(group_by).n_unique()} (G, derivative), {n_charged} with charged fragments")
    return selected
