import numpy as np

import grrmlib as gl


def test_read_eq_list(eqs):
    assert list(eqs) == [0, 1]

    eq = eqs[0]
    assert eq.name == 0
    assert eq.scfenergy == -152.8
    assert eq.afirenergy == -152.79
    assert eq.mult == 1
    assert eq.zpve == 0.04
    assert len(eq.nmeigen) == 12
    assert eq.symbols == ["O", "H", "H", "O", "H", "H"]
    assert list(eq.labels) == [1, 2, 3, 4, 5, 6]
    assert np.allclose(eqs[1].atomcoords[0], [0.1, 0.0, 0.0])


def test_read_pt_list(pts):
    assert list(pts) == [0, 1]
    assert pts[0].connection == (0, 1)
    assert pts[1].connection == (1, "??")
    assert len(pts[0].nmeigen) == 12


def test_read_grrm_input(data):
    mol = gl.read_grrm_input(data / "grrm/min.com")

    assert mol.infile == "grrm_min"
    assert mol.route == "# MIN/B3LYP/6-31G(d)"
    assert mol.title is None
    assert (mol.charge, mol.mult) == (0, 1)
    assert mol.symbols == ["O", "H", "H"]
    assert mol.options == ["KeepSymmetry", "Freq"]


def test_write_and_read_grrm_input(tmp_path, water):
    gl.write_grrm_input(
        water,
        tmp_path / "a.com",
        infile="x",
        route="# MIN/HF/STO-3G",
        title="test",
        options=["Freq"],
    )
    mol = gl.read_grrm_input(tmp_path / "a.com")

    assert mol.infile == "x"
    assert mol.route == "# MIN/HF/STO-3G"
    assert mol.title == "test"
    assert mol.options == ["Freq"]
    assert np.allclose(mol.atomcoords, water.atomcoords)


def test_write_and_read_grrm_input_without_options(tmp_path, water):
    water.notes = [[1], [1], [2]]
    gl.write_grrm_input(water, tmp_path / "a.com")
    mol = gl.read_grrm_input(tmp_path / "a.com")

    assert mol.infile is None
    assert mol.route == "#"
    assert mol.options is None
    assert mol.notes == [[1], [1], [2]]


def test_write_grrm_inputs_with_prefixes(tmp_path, water):
    mols = gl.Molecules({(1, 0): water})
    gl.write_grrm_inputs(mols, tmp_path, "min.com", prefixes=["SG", "SEQ"])

    assert (tmp_path / "SG=1/SEQ=0/min.com").exists()
    assert list(gl.read_grrm_inputs(tmp_path, "min.com")) == [(1, 0)]
