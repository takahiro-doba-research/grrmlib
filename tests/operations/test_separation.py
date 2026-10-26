import grrmlib as gl


def test_separate_eqs_keys(eqs):
    seqs, _ = gl.separate_eqs(eqs)

    # 3 つの EQ から 3 つずつ、計 9 つのフラグメント (SEQ 0-8) ができ、
    # 結合が同じものが同じ SG にまとまる: SG 0 = {0, 6}, SG 1 = {1, 4, 7}, SG 2 = {2, 5, 8}, SG 3 = {3}
    assert list(seqs) == [(0, 0), (0, 6), (1, 1), (1, 4), (1, 7), (2, 2), (2, 5), (2, 8), (3, 3)]
    assert len(seqs[(1, 4)]) == 12
    assert len(seqs[(0, 6)]) == 29


def test_separate_eqs_table(eqs):
    seqs, df = gl.separate_eqs(eqs)

    assert df.columns == ["G", "EQ", "SEQ", "SG", "fragment", "labels"]
    # EQ 0 と EQ 2 は結合が同じなので同じ G にまとまる
    assert df["G"].to_list() == [0, 0, 0, 1, 1, 1, 0, 0, 0]
    assert df["EQ"].to_list() == [0, 0, 0, 1, 1, 1, 2, 2, 2]
    assert df["SEQ"].to_list() == list(range(9))
    assert df["SG"].to_list() == [0, 1, 2, 3, 1, 2, 0, 1, 2]
    assert df["fragment"].to_list() == [0, 1, 2] * 3

    # labels から元の EQ の原子をたどれる
    for eq, seq, sg, labels in df.select(["EQ", "SEQ", "SG", "labels"]).iter_rows():
        labels = [int(l) for l in labels.split("_")]
        assert seqs[(sg, seq)].symbols == eqs[eq].select_atoms(labels).symbols
    # 1 つの EQ のフラグメントで全原子を 1 回ずつ使う
    for _, df_eq in df.group_by("EQ"):
        labels = sorted(int(l) for s in df_eq["labels"] for l in s.split("_"))
        assert labels == list(range(1, 54))
