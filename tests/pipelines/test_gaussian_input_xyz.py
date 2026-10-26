import shutil

import pytest

import grrmlib as gl


@pytest.fixture
def paths(tmp_path, data, water):
    shutil.copy(data / "gaussian/SP_connection.com", tmp_path / "template_sp.com")
    water.charge, water.mult = 1, 2
    mols = gl.Molecules({(0, 1): water, (1, 1): water})
    gl.write_xyzs(mols, tmp_path / "mols", "mol.xyz", prefixes=["A", "B"], record=True)
    return {
        "template_sp": tmp_path / "template_sp.com",
        "folder": tmp_path / "mols",
    }


def test_xyzs_to_gaussian_inputs(paths):
    n = gl.xyzs_to_gaussian_inputs(paths["folder"], "mol.xyz", paths["template_sp"])
    assert n == 2

    # xyz と同じフォルダに書く
    mols = gl.read_gaussian_inputs(paths["folder"], "gaussian_sp.com")
    assert set(mols) == {(0, 1), (1, 1)}

    mol = mols[(0, 1)]
    assert mol.link0 == ["%mem=2000MW"]
    assert mol.route == ["#p B3LYP/LanL2DZ SCF=SP"]
    assert mol.title == ["single point calculation with gaussian"]
    assert (mol.charge, mol.mult) == (1, 2)
    assert mol.symbols == ["O", "H", "H"]
    assert mol.notes is None


def test_xyzs_to_gaussian_inputs_drops_chk(paths):
    # template の link0 から %Chk だけを除く
    text = paths["template_sp"].read_text()
    paths["template_sp"].write_text("%Chk=sp.chk\n" + text)
    gl.xyzs_to_gaussian_inputs(paths["folder"], "mol.xyz", paths["template_sp"])

    mol = gl.read_gaussian_input(paths["folder"] / "A=0/B=1/gaussian_sp.com")
    assert mol.link0 == ["%mem=2000MW"]


def test_xyzs_to_gaussian_inputs_overwrite(paths):
    gl.xyzs_to_gaussian_inputs(paths["folder"], "mol.xyz", paths["template_sp"])
    with pytest.raises(FileExistsError):
        gl.xyzs_to_gaussian_inputs(paths["folder"], "mol.xyz", paths["template_sp"])
    gl.xyzs_to_gaussian_inputs(
        paths["folder"], "mol.xyz", paths["template_sp"], overwrite=True
    )
