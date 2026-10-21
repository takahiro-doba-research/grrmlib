import grrmlib as gl


def test_to_dataframe_tuple_keys(eqs):
    mols = gl.Molecules({(1, 2): eqs[0], (1, 3): eqs[1]})
    df = gl.to_dataframe(mols, ["SG", "SEQ"], ["scfenergy", "status"])
    
    assert df.columns == ["SG", "SEQ", "scfenergy", "status"]
    assert df["SEQ"].to_list() == [2, 3]
    assert df["scfenergy"].to_list() == [eqs[0].scfenergy, eqs[1].scfenergy]
    # 持っていない属性は null
    assert df["status"].to_list() == [None, None]


def test_to_dataframe_scalar_keys(eqs):
    df = gl.to_dataframe(eqs, "EQ", "scfenergy")
    
    assert df.columns == ["EQ", "scfenergy"]
    assert df["EQ"].to_list() == [0, 1, 2]
