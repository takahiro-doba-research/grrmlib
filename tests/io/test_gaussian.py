import numpy as np
import pytest

import grrmlib as gl


def test_read_gaussian_input(data):
    mol = gl.read_gaussian_input(data / "gaussian/sp.com")

    assert mol.link0 == ["%NProcShared=8", "%Mem=8GB", "%Chk=sp.chk"]
    assert mol.route == ["# B3LYP/6-31G(d) SP"]
    assert mol.title == ["water"]
    assert (mol.charge, mol.mult) == (0, 1)
    assert mol.symbols == ["O", "H", "H"]
    assert list(mol.labels) == [1, 2, 3]
    assert np.allclose(mol.atomcoords[1], [0.96, 0.0, 0.0])
    assert mol.notes is None


def test_write_and_read_gaussian_input(tmp_path, water):
    water.notes = [[1, 2, 3, 120.0], [1], [2]]
    gl.write_gaussian_input(
        water,
        tmp_path / "a.com",
        link0=["%Mem=1GB"],
        route=["# HF/STO-3G"],
        title=["test"],
    )
    mol = gl.read_gaussian_input(tmp_path / "a.com")

    assert mol.link0 == ["%Mem=1GB"]
    assert mol.route == ["# HF/STO-3G"]
    assert mol.title == ["test"]
    assert np.allclose(mol.atomcoords, water.atomcoords)
    assert mol.notes == [[1, 2, 3, 120.0], [1], [2]]


def test_write_gaussian_input_without_notes(tmp_path, water):
    water.notes = [[1], [1], [2]]
    gl.write_gaussian_input(water, tmp_path / "a.com", with_notes=False)

    assert gl.read_gaussian_input(tmp_path / "a.com").notes is None


def test_write_gaussian_input_defaults(tmp_path, water):
    water.charge = None
    water.mult = None
    gl.write_gaussian_input(water, tmp_path / "a.com")
    mol = gl.read_gaussian_input(tmp_path / "a.com")

    assert mol.link0 == []
    assert mol.route == ["#"]
    assert (mol.charge, mol.mult) == (0, 1)


def test_write_and_read_gaussian_inputs(tmp_path, water):
    mols = gl.Molecules({(0, 0): water, (0, 1): water})
    gl.write_gaussian_inputs(mols, tmp_path, "sp.com", route=["# HF/STO-3G"])

    mols2 = gl.read_gaussian_inputs(tmp_path, "sp.com")
    assert set(mols2) == {(0, 0), (0, 1)}
    assert mols2[(0, 1)].route == ["# HF/STO-3G"]


def test_read_gaussian_output_broken(data):
    with pytest.raises(Exception):
        gl.read_gaussian_output(data / "gaussian/broken.log")


def test_write_gaussian_scan(tmp_path, eqs):
    gl.write_gaussian_scan(eqs, tmp_path / "scan.log")
    text = (tmp_path / "scan.log").read_text()

    assert text.count("SCF Done") == len(eqs)
    assert "scan point     2 out of     2" in text
    assert text.endswith("Normal termination of Gaussian 16\n")
