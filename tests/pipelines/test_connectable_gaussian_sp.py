import shutil

import pytest

import grrmlib as gl


def test_connection_keys():
    keys = gl.connection_keys("work/seq", ["work/arene", "alkene"])
    assert keys == ["seq", "arene", "arene_angle", "alkene", "alkene_angle"]

    keys = gl.connection_keys("work/seq", ["work/arene", "alkene"], with_angles=False)
    assert keys == ["seq", "arene", "alkene"]


def test_connectables_to_gaussian_sp_missing_folder(tmp_path, data):
    with pytest.raises(FileNotFoundError):
        gl.connectables_to_gaussian_sp(
            tmp_path / "seq",
            [tmp_path / "sub"],
            data / "gaussian/SP_connection.com",
            tmp_path / "connection",
        )


def test_connectables_to_gaussian_sp(tmp_path, data):
    (tmp_path / "seq").mkdir()
    (tmp_path / "sub").mkdir()
    shutil.copy(data / "connectable/core.com", tmp_path / "seq/0.com")
    shutil.copy(data / "connectable/sub.com", tmp_path / "sub/0.com")

    gl.connectables_to_gaussian_sp(
        tmp_path / "seq",
        [tmp_path / "sub"],
        data / "gaussian/SP_connection.com",
        tmp_path / "connection",
    )
    paths = sorted((tmp_path / "connection").rglob("gaussian_sp.com"))
    assert [p.parent.relative_to(tmp_path / "connection").as_posix() for p in paths] == [
        f"seq=0/sub=0/sub_angle={angle}" for angle in [0, 120, 240]
    ]
