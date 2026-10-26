import shutil

import polars as pl
import pytest

import grrmlib as gl


# --------------------------------------------
# 未公開のデータ（tests/data_private）を使うテスト
# --------------------------------------------
def test_collect_grrm_min_outputs(tmp_path, data_private, capsys):
    folder = tmp_path / "separation/SG=0/SEQ=0/charge=0/mult=1"
    folder.mkdir(parents=True)
    shutil.copy(data_private("grrm/TDB-01-150-MIN.log"), folder / "grrm_min.log")
    folder = tmp_path / "separation/SG=1/SEQ=1/charge=0/mult=1"
    folder.mkdir(parents=True)
    (folder / "grrm_min.log").write_text("broken\n")

    df = gl.collect_grrm_min_outputs(
        tmp_path / "separation",
        "grrm_min.log",
        tmp_path / "grrm_min.csv",
        key_names=gl.SEPARATION_KEYS,
    )
    assert "1 / 2 grrm_min.log failed." in capsys.readouterr().out
    assert pl.read_csv(tmp_path / "grrm_min.csv").equals(df)

    assert df.columns == ["SG", "SEQ", "charge", "mult", "scfenergy", "status"]
    assert df.row(0) == (0, 0, 0, 1, pytest.approx(-2683.92399803277), "Minimum point was found")
    assert df.row(1) == (1, 1, 0, 1, None, None)
