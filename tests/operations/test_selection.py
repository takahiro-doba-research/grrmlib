import pytest
import polars as pl

import grrmlib as gl


def test_select_lowest():
    df = pl.DataFrame({
        "SG": [0, 0, 1, 1, 1, 2],
        "SEQ": [0, 6, 1, 4, 7, 3],
        "scfenergy": [-1.0, -2.0, -3.0, -1.0, -2.0, -1.0],
    })
    df_selected = gl.select_lowest(df, "SG")

    assert df_selected.columns == ["SG", "SEQ", "scfenergy"]
    assert df_selected.rows() == [(0, 6, -2.0), (1, 1, -3.0), (2, 3, -1.0)]


def test_select_lowest_keeps_column_order():
    df = pl.DataFrame({
        "SG": [0, 0],
        "SEQ": [0, 1],
        "charge": [0, 0],
        "scfenergy": [-1.0, -2.0],
    })
    # group_by の列が先頭に来ないようにする
    df_selected = gl.select_lowest(df, ["SG", "charge"])

    assert df_selected.columns == ["SG", "SEQ", "charge", "scfenergy"]
    assert df_selected.rows() == [(0, 1, 0, -2.0)]


def test_select_lowest_nulls_last():
    # 失敗したジョブ (null) は選ばない。全て失敗したグループは最初の行
    df = pl.DataFrame({
        "SG": [0, 0, 1, 1],
        "SEQ": [0, 1, 2, 3],
        "scfenergy": [None, -1.0, None, None],
    })
    assert gl.select_lowest(df, "SG").rows() == [(0, 1, -1.0), (1, 2, None)]


def test_select_charges():
    # G 0 は 2 つの分子に分かれ、どちらも ±1 になれる。G 1 は中性の分子 1 つだけ
    candidates = pl.DataFrame({
        "G": [0, 0, 0, 0, 1],
        "fragment": [0, 0, 1, 1, 0],
        "SG": [10, 10, 11, 11, 12],
        "charge": [1, -1, 1, -1, 0],
        "scfenergy": [-1.0, -1.2, -2.0, -2.5, -3.0],
    })
    selected = gl.select_charges(candidates, "G")

    # 合計 0 の組は (+1, -1) = -3.5 と (-1, +1) = -3.2。差は 0.3
    assert selected.select("G", "SG", "charge").rows() == [(0, 10, 1), (0, 11, -1), (1, 12, 0)]
    assert selected["margin"].to_list() == pytest.approx([0.3, 0.3, None], nan_ok=True)


def test_select_charges_no_combination():
    candidates = pl.DataFrame({
        "G": [0], "fragment": [0], "SG": [0], "charge": [1], "scfenergy": [-1.0],
    })
    with pytest.raises(ValueError, match="no combination"):
        gl.select_charges(candidates, "G")
