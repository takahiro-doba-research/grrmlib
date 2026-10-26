from collections.abc import Iterable
from pathlib import Path

import polars as pl

from ..core import Molecules
from ..io import key_to_folder, read_gaussian_input, write_grrm_inputs
from ..operations import select_lowest
from ._templates import read_grrm_template


def gaussian_sp_to_grrm_min(
    path_csv: str | Path,
    folder: str | Path,
    path_template: str | Path,
    path_selected: str | Path,
    *,
    key_names: list[str],
    group_by: Iterable[str] | str,
    basename_sp: str = "gaussian_sp.com",
    basename_min: str = "grrm_min.com",
    overwrite: bool = False
) -> pl.DataFrame:
    """
    Gaussian SP の結果の表（collect_gaussian_outputs の CSV）から、group_by ごとに
    scfenergy が最も低い構造を選び、その Gaussian の入力の構造で GRRM MIN の入力を書く。

    - 入力は folder/<key>/basename_sp から読み、同じフォルダの basename_min に書く
    - route, options は path_template の GRRM の入力から写す
    - 選んだ行を path_selected に CSV で書き、返す

    Example
    -------
    separation では key_names=SEPARATION_KEYS, group_by=["SG", "charge", "mult"]
    """
    df = select_lowest(pl.read_csv(path_csv), group_by)
    df.write_csv(path_selected)

    mols = Molecules()
    for row in df.select(key_names).iter_rows():
        key = row if len(row) > 1 else row[0]
        mols[key] = read_gaussian_input(key_to_folder(folder, key, key_names) / basename_sp)

    write_grrm_inputs(
        mols,
        folder,
        basename_min,
        prefixes=key_names,
        with_notes=False,
        overwrite=overwrite,
        **read_grrm_template(path_template)
    )
    return df
