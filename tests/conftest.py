"""
テスト全体で使う fixture。

tests/data には公開してよい実際の計算の入力・出力を置く（大きい log は先頭だけに縮めている）。
tests/data_private には未公開の計算データを置く。git の管理外なので、ない環境ではそのテストを skip する。
"""
from collections.abc import Callable
from pathlib import Path

import numpy as np
import pytest

import grrmlib as gl


@pytest.fixture
def data() -> Path:
    return Path(__file__).parent / "data"


@pytest.fixture
def data_private() -> Callable[[str], Path]:
    """
    data_private("gaussian/sp.log") のように使う。ファイルがなければそのテストを skip する。
    """
    folder = Path(__file__).parent / "data_private"

    def get(relpath: str) -> Path:
        path = folder / relpath
        if not path.exists():
            pytest.skip(f"tests/data_private/{relpath} がない")
        return path

    return get


@pytest.fixture
def water() -> gl.Molecule:
    return gl.Molecule(
        charge=0,
        mult=1,
        labels=np.array([1, 2, 3]),
        symbols=["O", "H", "H"],
        atomcoords=np.array([
            [0.0, 0.0, 0.0],
            [0.96, 0.0, 0.0],
            [-0.24, 0.93, 0.0],
        ]),
    )


@pytest.fixture
def eqs(data: Path) -> gl.Molecules:
    # SC-AFIR の EQ_list の先頭 3 件（Pd 錯体、53 原子）。
    # どれも 3 つのフラグメント（29 / 12 / 12 原子）に分かれ、EQ 0 と EQ 2 は結合が同じ
    return gl.read_eq_list(data / "grrm/td_c1_193_SCAFIR6_EQ_list.log")


@pytest.fixture
def pts(data: Path) -> gl.Molecules:
    # 同じ計算の PT_list の先頭 3 件。CONNECTION は 8 - 109, 49 - ??, 30 - 110
    return gl.read_pt_list(data / "grrm/td_c1_193_SCAFIR6_PT_list.log")
