import shutil
from pathlib import Path

import polars as pl
import pytest

import grrmlib as gl


@pytest.fixture
def workflow(tmp_path, data, monkeypatch):
    # workflow はカレントディレクトリにファイルを書く
    monkeypatch.chdir(tmp_path)
    shutil.copy(data / "grrm/td_c1_193_SCAFIR6_EQ_list.log", "EQ_list.log")
    shutil.copy(data / "gaussian/SP_connection.com", "template_sp.com")
    shutil.copy(data / "grrm/MIN_connection.com", "template_min.com")
    return gl.SeparationWorkflow(
        path_eq_list="EQ_list.log",
        path_gaussian_sp="template_sp.com",
        path_grrm_min="template_min.com",
    )


def _keys(paths: list[Path]) -> list[tuple[int, ...]]:
    # separation/SG=0/SEQ=0/charge=0/mult=1/xxx.com -> (0, 0, 0, 1)
    return sorted(
        tuple(int(part.split("=")[1]) for part in p.parent.parts[1:]) for p in paths
    )


def test_write_gaussian_sp(workflow, capsys):
    workflow.write_gaussian_sp(charges=[0, 1], mults=[1, 2])
    workflow.list_gaussian_sp()
    paths = [Path(p) for p in Path("gaussian_sp.txt").read_text().split()]

    # 3 つの EQ から 3 つずつ、計 9 つのフラグメント (SEQ 0-8) ができ、
    # 結合が同じものが同じ SG にまとまる: SG 0 = {0, 6}, SG 1 = {1, 4, 7}, SG 2 = {2, 5, 8}, SG 3 = {3}。
    # それぞれ (charge, mult) = (0, 1), (1, 2) の 2 通り
    seqs = [(0, 0), (0, 6), (1, 1), (1, 4), (1, 7), (2, 2), (2, 5), (2, 8), (3, 3)]
    assert _keys(paths) == sorted(
        (sg, seq, charge, mult) for sg, seq in seqs for charge, mult in [(0, 1), (1, 2)]
    )
    assert "18 gaussian_sp jobs listed." in capsys.readouterr().out

    mol = gl.read_gaussian_input("separation/SG=1/SEQ=4/charge=0/mult=1/gaussian_sp.com")
    assert mol.link0 == ["%mem=2000MW"]
    assert mol.route == ["#p B3LYP/LanL2DZ SCF=SP"]
    assert mol.title == ["single point calculation with gaussian"]
    assert len(mol) == 12
    assert mol.notes is None


def test_write_gaussian_sp_drops_chk(workflow):
    # template の link0 から %Chk だけを除く
    text = Path("template_sp.com").read_text()
    Path("template_sp.com").write_text("%Chk=sp.chk\n" + text)
    workflow.write_gaussian_sp(charges=0, mults=1)

    mol = gl.read_gaussian_input("separation/SG=3/SEQ=3/charge=0/mult=1/gaussian_sp.com")
    assert mol.link0 == ["%mem=2000MW"]


def test_write_grrm_min(workflow):
    workflow.write_gaussian_sp(charges=0, mults=1)

    # Gaussian の計算結果の代わりに gaussian_sp.csv を作る
    rows = [
        (0, 0, -1.0), (0, 6, -2.0),
        (1, 1, -3.0), (1, 4, -1.0), (1, 7, -2.0),
        (2, 2, -1.0), (2, 5, -1.5), (2, 8, -0.5),
        (3, 3, -1.0),
    ]
    pl.DataFrame(
        {
            "SG": [r[0] for r in rows],
            "SEQ": [r[1] for r in rows],
            "charge": [0] * len(rows),
            "mult": [1] * len(rows),
            "scfenergy": [r[2] for r in rows],
            "success": [True] * len(rows),
        }
    ).write_csv("gaussian_sp.csv")

    workflow.write_grrm_min()
    workflow.list_grrm_min()

    # SG ごとに最も安定な SEQ だけを選ぶ
    paths = [Path(p) for p in Path("grrm_min.txt").read_text().split()]
    assert _keys(paths) == [(0, 6, 0, 1), (1, 1, 0, 1), (2, 5, 0, 1), (3, 3, 0, 1)]

    mol = gl.read_grrm_input("separation/SG=0/SEQ=6/charge=0/mult=1/grrm_min.com")
    assert mol.route == "# MIN/B3LYP/Def2SVP"
    assert mol.options == ["GauProc=24", "GauMem=2000", "EigenCheck", "MinFreqValue=50"]
    assert len(mol) == 29


# --------------------------------------------
# 未公開のデータ（tests/data_private）を使うテスト
# --------------------------------------------
def test_analyze_gaussian_sp(workflow, data, data_private, capsys):
    # 計算が終わった log と、途中で落ちた log を 1 つずつ置く
    folder = Path("separation/SG=0/SEQ=0/charge=1/mult=1")
    folder.mkdir(parents=True)
    shutil.copy(data_private("gaussian/TDB-01-175-lower-OPT.log"), folder / "gaussian_sp.log")
    folder = Path("separation/SG=1/SEQ=1/charge=0/mult=1")
    folder.mkdir(parents=True)
    shutil.copy(data / "gaussian/broken.log", folder / "gaussian_sp.log")
    
    workflow.analyze_gaussian_sp()
    # 読めなかった log も失敗に数える
    assert "1 gaussian_sp jobs failed." in capsys.readouterr().out
    df = pl.read_csv("gaussian_sp.csv")
    
    assert df.columns == ["SG", "SEQ", "charge", "mult", "scfenergy", "success"]
    assert df.row(0) == (0, 0, 1, 1, pytest.approx(-3326.88401594, abs=1e-6), True)
    # 読めなかった log は scfenergy と success が null
    assert df.row(1) == (1, 1, 0, 1, None, None)


def test_analyze_grrm_min(workflow, data, data_private, capsys):
    folder = Path("separation/SG=0/SEQ=0/charge=0/mult=1")
    folder.mkdir(parents=True)
    shutil.copy(data_private("grrm/TDB-01-150-MIN.log"), folder / "grrm_min.log")
    folder = Path("separation/SG=1/SEQ=1/charge=0/mult=1")
    folder.mkdir(parents=True)
    (folder / "grrm_min.log").write_text("broken\n")
    
    workflow.analyze_grrm_min()
    assert "1 grrm_min jobs failed." in capsys.readouterr().out
    df = pl.read_csv("grrm_min.csv")
    
    assert df.columns == ["SG", "SEQ", "charge", "mult", "scfenergy", "status"]
    assert df.row(0) == (0, 0, 0, 1, pytest.approx(-2683.92399803277), "Minimum point was found")
    assert df.row(1) == (1, 1, 0, 1, None, None)
