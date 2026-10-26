import shutil

import polars as pl
import pytest

import grrmlib as gl


@pytest.fixture
def paths(tmp_path, data):
    shutil.copy(data / "gaussian/SP_connection.com", tmp_path / "template_sp.com")
    shutil.copy(data / "grrm/MIN_connection.com", tmp_path / "template_min.com")
    return {
        "eq_list": data / "grrm/td_c1_193_SCAFIR6_EQ_list.log",
        "template_sp": tmp_path / "template_sp.com",
        "template_min": tmp_path / "template_min.com",
        "folder": tmp_path / "separation",
        "table": tmp_path / "separation.csv",
    }


def test_gaussian_sp_to_grrm_min(paths, tmp_path, folder_keys):
    gl.eq_list_to_xyzs(paths["eq_list"], paths["folder"], paths["table"], charges=0, mults=1)
    gl.xyzs_to_gaussian_inputs(paths["folder"], "seq.xyz", paths["template_sp"])

    # Gaussian の計算結果の代わりに collect_gaussian_outputs の CSV を作る。
    # SG 2 の SEQ 5 は失敗 (null) とする
    rows = [
        (0, 0, -1.0), (0, 6, -2.0),
        (1, 1, -3.0), (1, 4, -1.0), (1, 7, -2.0),
        (2, 2, -1.0), (2, 5, None), (2, 8, -0.5),
        (3, 3, -1.0),
    ]
    pl.DataFrame(
        {
            "SG": [r[0] for r in rows],
            "SEQ": [r[1] for r in rows],
            "charge": [0] * len(rows),
            "mult": [1] * len(rows),
            "scfenergy": [r[2] for r in rows],
            "success": [r[2] is not None for r in rows],
        }
    ).write_csv(tmp_path / "gaussian_sp.csv")

    df = gl.gaussian_sp_to_grrm_min(
        tmp_path / "gaussian_sp.csv",
        paths["folder"],
        paths["template_min"],
        tmp_path / "gaussian_sp_selected.csv",
        key_names=gl.SEPARATION_KEYS,
        group_by=["SG", "charge", "mult"],
    )

    # SG ごとに最も安定な SEQ だけを選ぶ。失敗したジョブは選ばない
    expected = [(0, 6, 0, 1), (1, 1, 0, 1), (2, 2, 0, 1), (3, 3, 0, 1)]
    assert folder_keys(paths["folder"], "grrm_min.com") == expected
    assert df.select(gl.SEPARATION_KEYS).rows() == expected
    assert pl.read_csv(tmp_path / "gaussian_sp_selected.csv").equals(df)

    mol = gl.read_grrm_input(paths["folder"] / "SG=0/SEQ=6/charge=0/mult=1/grrm_min.com")
    assert mol.route == "# MIN/B3LYP/Def2SVP"
    assert mol.options == ["GauProc=24", "GauMem=2000", "EigenCheck", "MinFreqValue=50"]
    assert len(mol) == 29
    assert mol.notes is None
