import numpy as np
import pytest

import grrmlib as gl
from grrmlib.io.grrm._spin import spin2_to_mult


def test_read_eq_list(eqs):
    assert list(eqs) == [0, 1, 2]

    eq = eqs[0]
    assert eq.name == 0
    # Energy    = -1004.231330835170 (-1004.384498065987 :    0.000000000000)
    assert eq.scfenergy == -1004.384498065987
    assert eq.afirenergy == -1004.231330835170
    assert eq.mult == 1
    assert eq.zpve == 0.427857910044
    assert len(eq.nmeigen) == 153
    assert len(eq) == 53
    assert eq.symbols[24] == "Pd"
    assert list(eq.labels) == list(range(1, 54))
    assert np.allclose(eq.atomcoords[0], [1.293121677440, 1.886125453449, -0.825993870107])


def test_read_pt_list(pts):
    assert list(pts) == [0, 1, 2]
    assert [pt.connection for pt in pts.values()] == [(8, 109), (49, "??"), (30, 110)]
    assert pts[0].scfenergy == -1004.400669166796
    assert len(pts[0].nmeigen) == 153


@pytest.mark.parametrize(
    ("spin2", "mult"),
    [
        (0.0, 1),
        (0.75, 2),
        (2.0, 3),
        (3.75, 4),
        (6.0, 5),
        # スピン混入で <S^2> がずれていても最も近い多重度にする
        (0.7712, 2),
        (2.0315, 3),
    ],
)
def test_spin2_to_mult(spin2, mult):
    assert spin2_to_mult(spin2) == mult
    assert isinstance(spin2_to_mult(spin2), int)


def test_read_eq_list_mult(tmp_path, data):
    # 本物の EQ_list の Spin(**2) だけを二重項 (0.75) と三重項 (2.0) に書き換える
    text = (data / "grrm/td_c1_193_SCAFIR6_EQ_list.log").read_text()
    blocks = text.split("# Geometry of ")
    blocks[1] = blocks[1].replace("Spin(**2) =    0.000000000000", "Spin(**2) =    0.750000000000")
    blocks[2] = blocks[2].replace("Spin(**2) =    0.000000000000", "Spin(**2) =    2.000000000000")
    (tmp_path / "EQ_list.log").write_text("# Geometry of ".join(blocks))
    
    eqs = gl.read_eq_list(tmp_path / "EQ_list.log")
    assert [eq.mult for eq in eqs.values()] == [2, 3, 1]
    
    (tmp_path / "PT_list.log").write_text(
        (data / "grrm/td_c1_193_SCAFIR6_PT_list.log").read_text()
        .replace("Spin(**2) =    0.000000000000", "Spin(**2) =    0.750000000000")
    )
    assert {pt.mult for pt in gl.read_pt_list(tmp_path / "PT_list.log").values()} == {2}


def test_read_grrm_input(data):
    mol = gl.read_grrm_input(data / "grrm/MIN_connection.com")

    assert mol.infile is None
    assert mol.route == "# MIN/B3LYP/Def2SVP"
    assert mol.title is None
    assert (mol.charge, mol.mult) == (0, 1)
    assert len(mol) == 53
    assert mol.symbols[:3] == ["C", "C", "C"]
    assert mol.notes is None
    assert mol.options == ["GauProc=24", "GauMem=2000", "EigenCheck", "MinFreqValue=50"]


def test_read_grrm_input_oniom(data):
    mol = gl.read_grrm_input(data / "grrm/td_c1_193_SCAFIR6.com")

    assert mol.infile == "td_c1_193_SCAFIR5"
    assert mol.route == "# SC-AFIR/ONIOM(B3LYP/LanL2DZ:UFF) Int(Grid=FineGrid)"
    assert (mol.charge, mol.mult) == (0, 1)
    assert len(mol) == 53
    # ONIOM の原子の指定と層はそのまま残す
    assert mol.symbols[0] == "C-C_R--0.136243"
    assert mol.notes[0] == [1, "H"]
    assert mol.options[0] == "SaveDataInDirectory"
    assert mol.options[-1] == "grrm root=/home/GRRMdv2616"


def test_write_and_read_grrm_input_oniom(tmp_path, data):
    mol = gl.read_grrm_input(data / "grrm/td_c1_193_Repath1.com")
    gl.write_grrm_input(
        mol,
        tmp_path / "a.com",
        infile=mol.infile,
        route=mol.route,
        options=mol.options,
    )
    mol2 = gl.read_grrm_input(tmp_path / "a.com")

    assert mol2.route == "# Repath/ONIOM(B3LYP/LanL2DZ:UFF) Int(Grid=FineGrid)"
    assert mol2.symbols == mol.symbols
    assert mol2.notes == mol.notes
    assert mol2.options == mol.options
    assert np.allclose(mol2.atomcoords, mol.atomcoords)


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


# --------------------------------------------
# 未公開のデータ（tests/data_private）を使うテスト
# --------------------------------------------
def test_read_eq_list_full(data_private):
    eqs = gl.read_eq_list(data_private("grrm/td_c1_193_SCAFIR6_EQ_list.log"))
    assert len(eqs) == 3954


def test_read_pt_list_full(data_private):
    pts = gl.read_pt_list(data_private("grrm/td_c1_193_SCAFIR6_PT_list.log"))
    assert len(pts) == 7578


def test_read_lists_without_afir(data_private):
    # AFIR を使わない計算の list は Energy の行に括弧がない
    eqs = gl.read_eq_list(data_private("grrm/TDB-01-150-LUP_EQ_list.log"))
    tss = gl.read_pt_list(data_private("grrm/TDB-01-150-LUP_TS_list.log"))
    pts = gl.read_pt_list(data_private("grrm/TDB-01-150-LUP_PT_list.log"))

    assert len(eqs) == 3 and len(tss) == 1 and len(pts) == 2
    assert eqs[0].scfenergy == -2683.923924264470
    assert eqs[0].afirenergy is None
    assert tss[0].connection == (0, 1)
    assert pts[1].connection == (1, "??")


def test_read_grrm_min_output(data_private):
    itrs, opt, freq = gl.read_grrm_min_output(data_private("grrm/TDB-01-150-MIN.log"))

    assert len(itrs) == 59
    assert all(itr.mult == 1 for itr in itrs)
    assert opt.status == "Minimum point was found"
    assert opt.mult == 1
    assert opt.scfenergy == pytest.approx(-2683.92399803277)
    assert len(opt) == 55
    assert opt.hessian.shape == (165, 165)
    # 今は下三角だけに値が入り、上三角は NaN
    assert np.isfinite(opt.hessian[np.tril_indices(165)]).all()
    assert freq is not None and len(freq.vibfreqs) > 0


def test_read_grrm_min_outputs(tmp_path, data_private):
    import shutil

    (tmp_path / "1").mkdir()
    shutil.copy(data_private("grrm/TDB-01-150-MIN.log"), tmp_path / "1/min.log")
    (tmp_path / "2").mkdir()
    (tmp_path / "2/min.log").write_text("broken\n")

    mols = gl.read_grrm_min_outputs(tmp_path, "min.log", ignore_errors=True)
    assert mols[1].status == "Minimum point was found"
    assert mols[2].atomcoords is None


def test_read_lup_ts(data_private):
    from grrmlib.io.grrm.lup_ts import read_lup_ts

    lup = read_lup_ts(data_private("grrm/TDB-01-150-LUP_TS0.log"))
    assert len(lup.irc) == 217
    assert all(mol.mult == 1 for mol in lup.irc)
