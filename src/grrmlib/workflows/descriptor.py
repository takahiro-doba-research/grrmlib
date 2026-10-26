from collections.abc import Iterable
from pathlib import Path

import polars as pl


def sum_energies(
    path_selected: str | Path,
    path_long: str | Path,
    path_table: str | Path,
    *,
    substituents: Iterable[str]
) -> pl.DataFrame:
    """
    select_opt（substituents あり）で選んだ各分子のエネルギーを、G と誘導体ごとに合計する。

    - path_long：G と誘導体ごとのエネルギー（G, 置換基の番号, scfenergy）
    - path_table：誘導体を行、G を列（"G=0", "G=1", ...）にした表。値は scfenergy（Hartree。相対エネルギーにはしない）
    - 表を返す。値のない (G, 誘導体) があると ValueError
    """
    substituents = list(substituents)
    selected = pl.read_csv(path_selected)
    energy = (
        # 誘導体の列（*_derivative）を置換基の名前にする。分子ごとの置換基の列は使わない
        selected
        .select("G", *[pl.col(f"{name}_derivative").alias(name) for name in substituents], "scfenergy")
        .group_by("G", *substituents)
        .agg(pl.col("scfenergy").sum())
        .sort("G", *substituents)
    )
    energy.write_csv(path_long)

    table = (
        energy
        .with_columns(pl.format("G={}", "G").alias("column"))
        .pivot(on="column", index=substituents, values="scfenergy", sort_columns=False)
        .sort(substituents)
    )
    if table.null_count().sum_horizontal().item():
        raise ValueError("some (G, derivative) have no energy")
    table.write_csv(path_table)
    print(f"energy: {energy.height} (G, derivative), table {table.height} derivatives x {table.width - len(substituents)} G")
    return table
