from pathlib import Path

import polars as pl

from ..pipelines import eq_list_to_xyzs


def make_separation_sp(
    path_eq_list: str | Path,
    folder: str | Path,
    path_table: str | Path,
    *,
    charges: int | list[int],
    mults: int | list[int],
    basename: str = "init.xyz"
) -> pl.DataFrame:
    """
    EQ を分子ごとに分離し、SP の初期構造を folder/SG=*/SEQ=*/charge=*/mult=*/basename に書く（eq_list_to_xyzs）。
    対応表（G, EQ, SEQ, SG, fragment, labels）を path_table に書き、件数を表示して返す。
    """
    table = eq_list_to_xyzs(path_eq_list, folder, path_table, charges=charges, mults=mults, basename=basename)

    n_inputs = sum(1 for _ in Path(folder).rglob(basename))
    print(f"EQ: {table['EQ'].n_unique()}")
    print(f"G: {table['G'].n_unique()}")
    print(f"SEQ: {table['SEQ'].n_unique()}")
    print(f"SG: {table['SG'].n_unique()}")
    print(f"inputs: {n_inputs}")
    return table
