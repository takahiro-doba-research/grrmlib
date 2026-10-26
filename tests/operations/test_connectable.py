import numpy as np
import pytest

import grrmlib as gl


@pytest.fixture
def core(data):
    return gl.read_connectable(data / "connectable/core.com")


@pytest.fixture
def sub(data):
    return gl.read_connectable(data / "connectable/sub.com")


def test_read_connectivity(core):
    connectivity = gl.read_connectivity(core)

    assert connectivity[0]["labels"] == [1, 2, 3, 4]
    assert connectivity[0]["labels_ref"] == [1, 5, 2]
    assert connectivity[0]["target"] == 1
    assert connectivity[1]["labels"] == [5]
    assert "angle" not in connectivity[0]


def test_step_to_angles(sub):
    subs = gl.step_to_angles(sub)

    assert list(subs) == [0.0, 120.0, 240.0]
    assert subs[120.0].notes[0] == [0, 5, 2, 120.0]


def test_reset_connectable_notes(core):
    mol = gl.reset_connectable_notes(core.remove_atoms([3]), offset=10)

    assert list(mol.labels) == [11, 12, 13, 14]
    # notes の中のラベルも付け直す
    assert mol.notes[0] == [0, 14, 12]
    assert mol.notes[3] == [1, 11, 12]


def test_reset_connectable_notes_removed_ref(core):
    # 参照する原子 (H5) を除くと、その notes は fragment だけになる
    mol = gl.reset_connectable_notes(core.remove_atoms([5]))
    assert mol.notes[0] == [0]


def test_connect(core, sub):
    mol = gl.connect(core, gl.step_to_angles(sub)[0.0])

    # core の H5 を sub の H5 に置き換えるので、原子数と元素は変わらない
    assert mol.symbols == ["C", "H", "H", "H", "H"]
    assert list(mol.labels) == [1, 2, 3, 4, 5]
    # 新しい H は core の C1–H5 の方向に、sub の C1–H5 の距離で置かれる
    np.testing.assert_allclose(mol.atomcoords, core.atomcoords, atol=1e-8)


def test_connect_angle(core, sub):
    # 結合軸のまわりに回しても、軸の上の原子は動かない
    mol = gl.connect(core, gl.step_to_angles(sub)[120.0])
    np.testing.assert_allclose(mol.atomcoords, core.atomcoords, atol=1e-8)


def test_connect_all_iter(core, sub):
    subs = gl.Molecules({0: sub}).expand(gl.step_to_angles).flatten()
    out = list(gl.connect_all_iter(gl.Molecules({0: core}), [subs]))
    assert [names for names, _ in out] == [(0, 0, 0.0), (0, 0, 120.0), (0, 0, 240.0)]


def test_connect_all(core, sub):
    # 置換基の key は (置換基, 角度)
    subs = gl.Molecules({0: sub}).expand(gl.step_to_angles).flatten()
    mols = gl.connect_all(gl.Molecules({0: core}), [subs])
    # key は (seq, 置換基, 角度)。角度は int
    assert list(mols) == [(0, 0, 0), (0, 0, 120), (0, 0, 240)]
