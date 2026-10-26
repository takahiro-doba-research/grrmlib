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
