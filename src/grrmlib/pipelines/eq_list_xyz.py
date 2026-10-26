from pathlib import Path

import polars as pl

from ..io import read_eq_list, write_xyzs
from ..operations import product_charge_mult, separate_eqs


# folder/SG=1/SEQ=2/charge=0/mult=1/ のようにフォルダを分ける
SEPARATION_KEYS = ["SG", "SEQ", "charge", "mult"]


def eq_list_to_xyzs(
    path_eq_list: str | Path,
    folder: str | Path,
    path_table: str | Path,
    *,
    charges: int | list[int],
    mults: int | list[int],
    basename: str = "seq.xyz",
    overwrite: bool = False
) -> pl.DataFrame:
    """
    EQ_list の EQ をフラグメントに分け、同じものを除いて、
    有効な電荷・多重度の組み合わせごとに xyz を書く。

    - xyz は folder/SG=*/SEQ=*/charge=*/mult=*/basename に書く。
      タイトル行には charge, mult, labels を JSON で書く（write_xyz の record=True）
    - 計算の入力は xyz から別に作る（Gaussian なら xyzs_to_gaussian_inputs）
    - G, EQ, SEQ, SG の対応表を path_table に CSV で書き、返す
    """
    eqs = read_eq_list(path_eq_list)
    seqs, table = separate_eqs(eqs)
    seqs = seqs.expand(lambda seq: product_charge_mult(seq, charges, mults)).flatten()

    write_xyzs(
        seqs,
        folder,
        basename,
        prefixes=SEPARATION_KEYS,
        record=True,
        overwrite=overwrite
    )
    table.write_csv(path_table)
    return table
