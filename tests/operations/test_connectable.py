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


@pytest.mark.xfail(
    reason="connect が古い形 (座標の配列と添字) で overlay / rotate を呼んでいる",
    raises=AttributeError,
    strict=True,
)
def test_connect(core, sub):
    # TODO: 直したら、つないだ構造（原子数や結合距離）も確かめる。
    #       そのためには tests/data/connectable の notes の意味を実データに合わせる必要がある
    mol = gl.connect(core, gl.step_to_angles(sub)[0.0])
    assert isinstance(mol, gl.Molecule)


@pytest.mark.xfail(
    reason="connect が古い形 (座標の配列と添字) で overlay / rotate を呼んでいる",
    raises=AttributeError,
    strict=True,
)
def test_connect_all_iter(core, sub):
    subs = gl.step_to_angles(sub)
    out = list(gl.connect_all_iter(gl.Molecules({0: core}), [subs]))
    assert len(out) > 0
