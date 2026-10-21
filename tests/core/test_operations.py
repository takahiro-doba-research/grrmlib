import numpy as np
import pytest

import grrmlib as gl


def test_get_distance(eqs):
    # 原子はラベルで指定する
    assert gl.get_distance(eqs[0], 1, 2) == pytest.approx(0.96)
    assert gl.get_distance(eqs[0], 1, 4) == pytest.approx(5.0)


def test_get_dihedral_angle(eqs):
    # 水 2 分子は同じ平面にある
    angle = gl.get_dihedral_angle(eqs[0], 3, 1, 4, 6, degrees=True)
    assert angle == pytest.approx(0.0, abs=1e-6)


def test_get_distance_matrix(eqs):
    D = gl.get_distance_matrix(eqs[0])

    assert D.shape == (6, 6)
    assert np.allclose(D, D.T)
    assert D[0, 3] == pytest.approx(5.0)


def test_get_adj_matrix(eqs):
    A = gl.get_adj_matrix(eqs[0])

    assert A[0, 1] and A[0, 2]
    assert not A[0, 3]
    assert not A[1, 2]
    assert not A.diagonal().any()


def test_separate(eqs):
    frags = gl.separate(eqs[0])

    assert list(frags) == [0, 1]
    assert list(frags[1].labels) == [4, 5, 6]
    assert frags[1].symbols == ["O", "H", "H"]


def test_is_identical(eqs):
    frags = gl.separate(eqs[0])

    assert gl.is_identical(eqs[0], eqs[1])
    # ラベルが違うと別の分子とみなす
    assert not gl.is_identical(frags[0], frags[1])
    assert gl.is_identical(frags[0], frags[1].reset_labels())


def test_cluster(eqs):
    frags = eqs.expand(gl.separate).flatten()
    groups = frags.cluster(gl.is_identical)

    assert len(groups) == 2
    assert set(groups[0]) == {(0, 0), (1, 0)}


def test_rotate(eqs):
    # 1-4 軸まわりに 180° 回すと、軸の上にない原子は y が反転する
    coords = gl.rotate(eqs[0], 1, 4, 180, degrees=True)

    assert coords.shape == (6, 3)
    assert np.allclose(coords[[0, 1, 3, 4]], eqs[0].atomcoords[[0, 1, 3, 4]])
    assert coords[2, 1] == pytest.approx(-0.93)


def test_overlay(eqs):
    # 同じ構造どうしなら座標は変わらない
    coords = gl.overlay(eqs[0], eqs[1], 1, 2, 3, 1, 2, 3)
    assert np.allclose(coords, eqs[0].atomcoords, atol=1e-6)
