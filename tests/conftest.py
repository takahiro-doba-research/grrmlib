"""
テスト全体で使う fixture。

tests/data のファイルは各形式に合わせて手で作ったもので、実際の計算の出力ではない。
"""
from pathlib import Path

import numpy as np
import pytest

import grrmlib as gl


@pytest.fixture
def data() -> Path:
    return Path(__file__).parent / "data"


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
    # 水 2 分子が 5 Å 離れて並んだ構造が 2 つ
    return gl.read_eq_list(data / "grrm/EQ_list.log")


@pytest.fixture
def pts(data: Path) -> gl.Molecules:
    # TS 0: EQ 0 - EQ 1, TS 1: EQ 1 - ??
    return gl.read_pt_list(data / "grrm/TS_list.log")
