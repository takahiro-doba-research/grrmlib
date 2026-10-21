import networkx as nx

import grrmlib as gl


def test_reaction_path_network(eqs, pts):
    rpn = gl.ReactionPathNetwork(eqs=eqs, pts=pts)

    # 行き先が不明な PT は "??" という EQ につながる
    assert set(rpn.nodes) == {0, 1, "??"}
    rpn.remove_eq("??")
    assert set(rpn.nodes) == {0, 1}

    assert list(rpn.pts_between(0, 1)) == [0]
    assert len(rpn.pts_between(0, 0)) == 0
    assert nx.shortest_path(rpn, 0, 1) == [0, 1]


def test_cluster(eqs, pts):
    rpn = gl.ReactionPathNetwork(eqs=eqs, pts=pts).remove_eq("??")

    # すべての EQ を 1 つのグループにまとめる
    grpn = rpn.cluster(lambda eq0, eq1: True)
    assert list(grpn.nodes) == [0]
    assert list(grpn.eqs_of_group(0)) == [0, 1]
    assert grpn.group_of_eq(1) == 0
    assert grpn.group_of_eq(9) is None
    assert list(grpn.pts_between(0, 0)) == [0]

    # EQ ごとに別のグループにする
    grpn = rpn.cluster(lambda eq0, eq1: False)
    assert list(grpn.pts_between(0, 1)) == [0]
