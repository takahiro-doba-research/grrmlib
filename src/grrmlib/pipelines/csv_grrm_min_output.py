from pathlib import Path

import polars as pl

from ..io import read_grrm_min_outputs, to_dataframe


def collect_grrm_min_outputs(
    folder: str | Path,
    basename: str,
    path_csv: str | Path,
    *,
    key_names: list[str],
    version: str = "GRRM23"
) -> pl.DataFrame:
    """
    folder/<key>/basename の GRRM MIN の log から scfenergy と status を集め、
    path_csv に CSV で書き、返す。読めなかった log の行は null にする。
    """
    mols = read_grrm_min_outputs(folder, basename, version=version, ignore_errors=True)
    df = to_dataframe(mols, key_names, ["scfenergy", "status"])
    df.write_csv(path_csv)

    # log を読めなかったジョブは null なので、それも失敗に数える
    n_failed = df.filter(pl.col("status").ne_missing("Minimum point was found")).height
    print(f"{n_failed} / {df.height} {basename} failed.")
    return df
