from pathlib import Path

import polars as pl

from ..io import read_gaussian_outputs, to_dataframe


def collect_gaussian_outputs(
    folder: str | Path,
    basename: str,
    path_csv: str | Path,
    *,
    key_names: list[str]
) -> pl.DataFrame:
    """
    folder/<key>/basename の Gaussian の log から scfenergy と success を集め、
    path_csv に CSV で書き、返す。読めなかった log の行は null にする。
    """
    mols = read_gaussian_outputs(folder, basename, ignore_errors=True)
    df = to_dataframe(mols, key_names, ["scfenergy", "success"])
    df.write_csv(path_csv)

    # log を読めなかったジョブは null なので、それも失敗に数える
    n_failed = df.filter(pl.col("success").ne_missing(True)).height
    print(f"{n_failed} / {df.height} {basename} failed.")
    return df
