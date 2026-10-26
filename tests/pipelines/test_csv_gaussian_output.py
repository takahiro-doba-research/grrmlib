import shutil

import polars as pl
import pytest

import grrmlib as gl


# --------------------------------------------
# 未公開のデータ（tests/data_private）を使うテスト
# --------------------------------------------
def test_collect_gaussian_outputs(tmp_path, data, data_private, capsys):
    # 計算が終わった log と、途中で落ちた log を 1 つずつ置く
    folder = tmp_path / "separation/SG=0/SEQ=0/charge=1/mult=1"
    folder.mkdir(parents=True)
    shutil.copy(data_private("gaussian/TDB-01-175-lower-OPT.log"), folder / "gaussian_sp.log")
    folder = tmp_path / "separation/SG=1/SEQ=1/charge=0/mult=1"
    folder.mkdir(parents=True)
    shutil.copy(data / "gaussian/broken.log", folder / "gaussian_sp.log")

    df = gl.collect_gaussian_outputs(
        tmp_path / "separation",
        "gaussian_sp.log",
        tmp_path / "gaussian_sp.csv",
        key_names=gl.SEPARATION_KEYS,
    )
    # 読めなかった log も失敗に数える
    assert "1 / 2 gaussian_sp.log failed." in capsys.readouterr().out
    assert pl.read_csv(tmp_path / "gaussian_sp.csv").equals(df)

    assert df.columns == ["SG", "SEQ", "charge", "mult", "scfenergy", "success"]
    assert df.row(0) == (0, 0, 1, 1, pytest.approx(-3326.88401594, abs=1e-6), True)
    # 読めなかった log は scfenergy と success が null
    assert df.row(1) == (1, 1, 0, 1, None, None)
