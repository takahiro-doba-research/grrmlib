import numpy as np
import pytest

import grrmlib as gl


def test_read_gaussian_input(data):
    mol = gl.read_gaussian_input(data / "gaussian/SP_connection.com")
    
    assert mol.link0 == ["%mem=2000MW"]
    assert mol.route == ["#p B3LYP/LanL2DZ SCF=SP"]
    assert mol.title == ["single point calculation with gaussian"]
    assert (mol.charge, mol.mult) == (0, 1)
    assert len(mol) == 53
    assert mol.symbols[24] == "Pd"
    assert list(mol.labels) == list(range(1, 54))
    assert np.allclose(mol.atomcoords[0], [1.29312168, 1.88612545, -0.82599387])
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
    assert "scan point     3 out of     3" in text
    assert text.endswith("Normal termination of Gaussian 16\n")


# --------------------------------------------
# 未公開のデータ（tests/data_private）を使うテスト
# --------------------------------------------
def test_read_gaussian_input_multiline_route(data_private):
    mol = gl.read_gaussian_input(data_private("gaussian/TDB-01-175-lower-OPT.gjf"))
    
    assert mol.link0 == ["%Chk=TDB-01-175-lower-OPT.chk", "%NProcShared=16", "%Mem=1600MW"]
    # 複数行の route はそのまま 1 行ずつ残す
    assert mol.route == [
        "#P M06L/Def2TZVP",
        "Opt=MaxCycle=100",
        "Freq=NoRaman",
        "SCF=MaxCycle=512",
        "Integral=Grid=UltraFine",
        "SCRF=(IEFPCM,Solvent=DiChloroEthane)",
    ]
    assert mol.title == ["Cp*Fe(CyJohnPhosS)+"]
    assert (mol.charge, mol.mult) == (1, 1)
    assert len(mol) == 83
    assert mol.symbols[:3] == ["Fe", "P", "S"]


def test_read_gaussian_output(data_private):
    mol = gl.read_gaussian_output(data_private("gaussian/TDB-01-175-lower-OPT.log"))
    mol_in = gl.read_gaussian_input(data_private("gaussian/TDB-01-175-lower-OPT.gjf"))
    
    assert mol.success is True
    assert (mol.charge, mol.mult) == (1, 1)
    assert mol.symbols == mol_in.symbols
    assert list(mol.labels) == list(range(1, 84))
    assert mol.atomcoords.shape == (83, 3)
    # 最後の "SCF Done:  E(RM06L) =  -3326.88401594" を Hartree で読む
    assert mol.scfenergy == pytest.approx(-3326.88401594, abs=1e-6)


def test_read_gaussian_outputs(tmp_path, data, data_private):
    import shutil
    
    (tmp_path / "1/0").mkdir(parents=True)
    shutil.copy(data_private("gaussian/TDB-01-175-lower-OPT.log"), tmp_path / "1/0/sp.log")
    (tmp_path / "2/0").mkdir(parents=True)
    shutil.copy(data / "gaussian/broken.log", tmp_path / "2/0/sp.log")
    
    with pytest.raises(Exception):
        gl.read_gaussian_outputs(tmp_path, "sp.log")
    
    mols = gl.read_gaussian_outputs(tmp_path, "sp.log", ignore_errors=True)
    assert mols[(1, 0)].success is True
    assert mols[(2, 0)].atomcoords is None
