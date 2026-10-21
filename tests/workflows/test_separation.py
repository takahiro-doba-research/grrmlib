import shutil
from pathlib import Path

import polars as pl
import pytest

import grrmlib as gl


@pytest.fixture
def workflow(tmp_path, data, monkeypatch):
    # workflow はカレントディレクトリにファイルを書く
    monkeypatch.chdir(tmp_path)
    shutil.copy(data / "grrm/EQ_list.log", "EQ_list.log")
    shutil.copy(data / "gaussian/sp.com", "template_sp.com")
    shutil.copy(data / "grrm/min.com", "template_min.com")
    return gl.SeparationWorkflow(
        path_eq_list="EQ_list.log",
        path_gaussian_sp="template_sp.com",
        path_grrm_min="template_min.com",
    )


def test_write_gaussian_sp(workflow, capsys):
    workflow.write_gaussian_sp(charges=[0, 1], mults=[1, 2])
    workflow.list_gaussian_sp()
    paths = [Path(p) for p in Path("gaussian_sp.txt").read_text().split()]

    # 2 つの EQ から水が 4 つ (SEQ 0-3) できる。is_identical はラベルも比べるので、
    # ラベル 1-3 の水 (SEQ 0, 2) とラベル 4-6 の水 (SEQ 1, 3) は別の SG になる。
    # それぞれ (charge, mult) = (0, 1), (1, 2) の 2 通り
    keys = sorted(
        tuple(int(part.split("=")[1]) for part in p.parent.parts[1:]) for p in paths
    )
    assert keys == [
        (0, 0, 0, 1), (0, 0, 1, 2), (0, 2, 0, 1), (0, 2, 1, 2),
        (1, 1, 0, 1), (1, 1, 1, 2), (1, 3, 0, 1), (1, 3, 1, 2),
    ]
    assert "8 gaussian_sp jobs listed." in capsys.readouterr().out

    mol = gl.read_gaussian_input(paths[0])
    # template の link0 から %Chk だけを除く
    assert mol.link0 == ["%NProcShared=8", "%Mem=8GB"]
    assert mol.route == ["# B3LYP/6-31G(d) SP"]
    assert mol.notes is None


def test_write_grrm_min(workflow):
    workflow.write_gaussian_sp(charges=0, mults=1)

    # Gaussian の計算結果の代わりに gaussian_sp.csv を作る
    pl.DataFrame(
        {
            "SG": [0, 0, 1, 1],
            "SEQ": [0, 2, 1, 3],
            "charge": [0, 0, 0, 0],
            "mult": [1, 1, 1, 1],
            "scfenergy": [-76.1, -76.2, -76.4, -76.3],
            "success": [True, True, True, True],
        }
    ).write_csv("gaussian_sp.csv")

    workflow.write_grrm_min()
    workflow.list_grrm_min()

    # SG ごとに最も安定な SEQ だけを選ぶ
    assert Path("grrm_min.txt").read_text().split() == [
        "separation/SG=0/SEQ=2/charge=0/mult=1/grrm_min.com",
        "separation/SG=1/SEQ=1/charge=0/mult=1/grrm_min.com",
    ]
    mol = gl.read_grrm_input("separation/SG=0/SEQ=2/charge=0/mult=1/grrm_min.com")
    assert mol.route == "# MIN/B3LYP/6-31G(d)"
    assert mol.options == ["KeepSymmetry", "Freq"]
