import numpy as np

import grrmlib as gl


def test_read_connectable(data):
    mol = gl.read_connectable(data / "connectable/sub.com")

    assert (mol.charge, mol.mult) == (0, 1)
    assert mol.symbols == ["C", "H", "H", "H", "H"]
    assert mol.notes == [[0, 5, 2, 3], [0], [0], [0], [1, 1, 2]]


def test_write_and_read_connectable(tmp_path, water):
    water.charge = -1
    water.mult = 2
    water.notes = [[0, 2, 3, 120.0], [0], [1]]
    gl.write_connectable(water, tmp_path / "a.com", title=["test"])
    mol = gl.read_connectable(tmp_path / "a.com")

    assert mol.title == ["test"]
    assert (mol.charge, mol.mult) == (-1, 2)
    assert mol.notes == [[0, 2, 3, 120.0], [0], [1]]
    assert np.allclose(mol.atomcoords, water.atomcoords)


def test_write_and_read_connectables(tmp_path, water):
    water.notes = [[0], [0], [1]]
    gl.write_connectables(gl.Molecules({0: water, 1: water}), tmp_path, "c.com")

    assert set(gl.read_connectables(tmp_path, "c.com")) == {0, 1}
