import numpy as np
import pytest

import grrmlib as gl


def test_get_distance(eqs):
    eq = eqs[0]
    expected = np.linalg.norm(eq.atomcoords[0] - eq.atomcoords[1])

    # 原子はラベルで指定する
    assert gl.get_distance(eq, 1, 2) == pytest.approx(expected)
    # ベンゼン環の C-C
    assert 1.3 < gl.get_distance(eq, 1, 2) < 1.5


def test_get_distance_uses_labels(eqs):
    frag = eqs[0].select_atoms([25, 27])
    expected = np.linalg.norm(eqs[0].atomcoords[24] - eqs[0].atomcoords[26])

    assert gl.get_distance(frag, 25, 27) == pytest.approx(expected)


def test_get_dihedral_angle(eqs):
    # ベンゼン環はほぼ平面
    angle = gl.get_dihedral_angle(eqs[0], 1, 2, 3, 4, degrees=True)
    assert abs(angle) < 1.0


def test_get_distance_matrix(eqs):
    D = gl.get_distance_matrix(eqs[0])

    assert D.shape == (53, 53)
    assert np.allclose(D, D.T)
    assert D[0, 1] == pytest.approx(gl.get_distance(eqs[0], 1, 2))


def test_get_adj_matrix(eqs):
    A = gl.get_adj_matrix(eqs[0])
    eq = eqs[0]

    def neighbors(label):
        return [int(eq.labels[j]) for j in np.where(A[label - 1])[0]]

    assert neighbors(1) == [2, 6, 7, 25]
    assert neighbors(25) == [1, 27, 35]  # Pd
    assert np.array_equal(A, A.T)
    assert not A.diagonal().any()


def test_separate(eqs):
    frags = gl.separate(eqs[0])

    assert list(frags) == [0, 1, 2]
    assert [len(f) for f in frags.values()] == [29, 12, 12]
    # ラベルは元の分子のものが残る
    assert list(frags[1].labels) == list(range(13, 25))
    assert list(frags[2].labels) == list(range(42, 54))
    assert "Pd" in frags[0].symbols


def test_is_identical(eqs):
    assert gl.is_identical(eqs[0], eqs[2])
    assert not gl.is_identical(eqs[0], eqs[1])

    frags0 = gl.separate(eqs[0])
    frags1 = gl.separate(eqs[1])
    assert gl.is_identical(frags0[1], frags1[1])
    assert not gl.is_identical(frags0[0], frags1[0])


def test_cluster(eqs):
    frags = eqs.expand(gl.separate).flatten()
    groups = frags.cluster(gl.is_identical)

    assert {g: list(m) for g, m in groups.items()} == {
        0: [(0, 0), (2, 0)],
        1: [(0, 1), (1, 1), (2, 1)],
        2: [(0, 2), (1, 2), (2, 2)],
        3: [(1, 0)],
    }


def test_rotate(eqs):
    eq = eqs[0]
    coords = gl.rotate(eq, 1, 2, 120, degrees=True)

    assert coords.shape == (53, 3)
    # 軸の上の原子は動かない
    assert np.allclose(coords[[0, 1]], eq.atomcoords[[0, 1]])
    # 軸からの距離は変わらない
    axis = eq.atomcoords[1] - eq.atomcoords[0]
    axis /= np.linalg.norm(axis)

    def dist_from_axis(x):
        v = x - eq.atomcoords[0]
        return np.linalg.norm(v - np.outer(v @ axis, axis), axis=1)

    assert np.allclose(dist_from_axis(coords), dist_from_axis(eq.atomcoords))
    # 360° 回すと元に戻る
    assert np.allclose(gl.rotate(eq, 1, 2, 360, degrees=True), eq.atomcoords)


def test_overlay(eqs):
    # 同じ構造どうしなら座標は変わらない
    coords = gl.overlay(eqs[0], eqs[0], 1, 2, 3, 1, 2, 3)
    assert np.allclose(coords, eqs[0].atomcoords, atol=1e-6)

    # EQ 1 の原子 1 を EQ 0 の原子 1 に重ね、1-2 の向きをそろえる
    coords = gl.overlay(eqs[0], eqs[1], 1, 2, 3, 1, 2, 3)
    assert np.allclose(coords[0], eqs[0].atomcoords[0])
    v0 = eqs[0].atomcoords[1] - eqs[0].atomcoords[0]
    v1 = coords[1] - coords[0]
    assert np.allclose(v0 / np.linalg.norm(v0), v1 / np.linalg.norm(v1))
    # 剛体の移動なので原子間距離は変わらない
    assert np.allclose(
        gl.get_distance_matrix(eqs[1]),
        np.linalg.norm(coords[:, None] - coords[None], axis=-1),
    )
