import networkx as nx

import grrmlib as gl


def test_reaction_path_network(eqs, pts):
    rpn = gl.ReactionPathNetwork(eqs=eqs, pts=pts)

    # PT の CONNECTION にある EQ もノードになり、行き先が不明な PT は "??" につながる
    assert set(rpn.nodes) == {0, 1, 2, 8, 109, 49, "??", 30, 110}
    rpn.remove_eq("??")
    assert "??" not in rpn

    assert list(rpn.pts_between(8, 109)) == [0]
    assert list(rpn.pts_between(109, 8)) == [0]
    assert len(rpn.pts_between(0, 1)) == 0
    assert nx.shortest_path(rpn, 30, 110) == [30, 110]


def test_cluster(eqs, pts):
    # 先頭 3 件の PT は先頭 3 件の EQ をつながないので、つなぎ先を付け替えて使う
    pts = pts.copy()
    pts[0].connection = (0, 1)
    pts[1].connection = (1, "??")
    pts[2].connection = (0, 2)
    rpn = gl.ReactionPathNetwork(eqs=eqs, pts=pts).remove_eq("??")

    # EQ 0 と EQ 2 は結合が同じなので同じグループ
    grpn = rpn.cluster(gl.is_identical)
    assert list(grpn.eqs_of_group(0)) == [0, 2]
    assert list(grpn.eqs_of_group(1)) == [1]
    assert grpn.group_of_eq(2) == 0
    assert grpn.group_of_eq(9) is None
    assert list(grpn.pts_between(0, 1)) == [0]
    assert list(grpn.pts_between(0, 0)) == [2]
