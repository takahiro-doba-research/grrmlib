import numpy as np

import grrmlib as gl


def test_read_xyz(data):
    mol = gl.read_xyz(data / "xyz/water.xyz")

    assert mol.title == "water"
    assert mol.symbols == ["O", "H", "H"]
    assert list(mol.labels) == [1, 2, 3]
    assert np.allclose(mol.atomcoords[2], [-0.24, 0.93, 0.0])


def test_write_and_read_xyz(tmp_path, water):
    # xyz には notes を書かない
    water.notes = [[1], [1], [2]]
    gl.write_xyz(water, tmp_path / "a.xyz", title="test")
    mol = gl.read_xyz(tmp_path / "a.xyz")

    assert mol.title == "test"
    assert mol.symbols == water.symbols
    assert np.allclose(mol.atomcoords, water.atomcoords)
    assert (tmp_path / "a.xyz").read_text().splitlines()[2].split()[4:] == []


def test_write_and_read_xyzs(tmp_path, water):
    mols = gl.Molecules({(3, 4): water, 5: water})
    gl.write_xyzs(mols, tmp_path, "mol.xyz")

    assert (tmp_path / "3/4/mol.xyz").exists()
    assert (tmp_path / "5/mol.xyz").exists()
    assert set(gl.read_xyzs(tmp_path, "mol.xyz")) == {(3, 4), 5}
