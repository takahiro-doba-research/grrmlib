from pathlib import Path

import pytest

import grrmlib as gl
from grrmlib.io import folder_to_key, key_to_folder


@pytest.mark.parametrize(
    ("key", "prefixes", "expected"),
    [
        (5, None, "root/5"),
        ((1, 2), None, "root/1/2"),
        ((1, None, -1), None, "root/1/None/-1"),
        ((1, 2), ["SG", "SEQ"], "root/SG=1/SEQ=2"),
        (3, "SG", "root/SG=3"),
    ],
)
def test_key_to_folder(key, prefixes, expected):
    assert key_to_folder("root", key, prefixes) == Path(expected)


def test_key_to_folder_length_mismatch():
    with pytest.raises(ValueError):
        key_to_folder("root", (1, 2), ["SG"])


@pytest.mark.parametrize(
    ("path", "expected"),
    [
        ("root/5/a.com", 5),
        ("root/1/2/a.com", (1, 2)),
        ("root/SG=1/SEQ=None/charge=-1/a.com", (1, None, -1)),
        ("root/x/2.5/a.com", ("x", "2.5")),
    ],
)
def test_folder_to_key(path, expected):
    assert folder_to_key("root", path) == expected


def test_write_and_read_with_prefixes(tmp_path, water):
    mols = gl.Molecules({(1, 2, -1, 2): water, (1, None, 0, 1): water})
    gl.write_xyzs(mols, tmp_path, "mol.xyz", prefixes=["SG", "SEQ", "charge", "mult"])

    assert (tmp_path / "SG=1/SEQ=2/charge=-1/mult=2/mol.xyz").exists()
    assert (tmp_path / "SG=1/SEQ=None/charge=0/mult=1/mol.xyz").exists()
    assert set(gl.read_xyzs(tmp_path, "mol.xyz")) == set(mols)


def test_ignore_errors(tmp_path, water, capsys):
    gl.write_xyz(water, tmp_path / "1/mol.xyz")
    (tmp_path / "2").mkdir()
    (tmp_path / "2/mol.xyz").write_text("broken\n")

    with pytest.raises(ValueError):
        gl.read_xyzs(tmp_path, "mol.xyz")

    mols = gl.read_xyzs(tmp_path, "mol.xyz", ignore_errors=True)
    assert len(mols[1]) == 3
    # 読めなかったファイルは空の Molecule になり、パスを表示する
    assert mols[2].atomcoords is None
    assert "2/mol.xyz" in capsys.readouterr().out


def test_no_overwrite(tmp_path, water):
    gl.write_xyz(water, tmp_path / "mol.xyz")

    with pytest.raises(FileExistsError):
        gl.write_xyz(water, tmp_path / "mol.xyz")

    gl.write_xyz(water, tmp_path / "mol.xyz", overwrite=True)
