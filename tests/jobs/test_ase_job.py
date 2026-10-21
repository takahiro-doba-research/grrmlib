import pytest
from ase.calculators.emt import EMT

import grrmlib as gl
from grrmlib.calculators import HARTREE_IN_EV, optimize, single_point, to_atoms
from grrmlib.jobs.ase_job import ASEJobRunner, read_result, read_results, write_result


def test_to_atoms(water):
    atoms = to_atoms(water)

    assert atoms.get_chemical_formula() == "H2O"
    assert atoms.info == {"charge": 0, "spin": 1}


def test_to_atoms_requires_charge_mult(water):
    water.mult = None
    with pytest.raises(ValueError):
        to_atoms(water)


def test_single_point(water):
    mol = single_point(water, EMT())
    atoms = to_atoms(water)
    atoms.calc = EMT()

    assert mol.scfenergy == pytest.approx(atoms.get_potential_energy() / HARTREE_IN_EV)
    assert mol.success and mol.status == "SP"
    # 元の分子は変えない
    assert not hasattr(water, "status")


def test_optimize(water):
    mol = optimize(water, EMT(), fmax=0.05)

    assert mol.success and mol.status == "Converged"
    assert mol.optsteps > 0
    assert mol.scfenergy < single_point(water, EMT()).scfenergy


def test_write_and_read_result(tmp_path, water):
    water.labels = water.labels + 10
    water.notes = [[1, 2, 3, 120.0], [1], [2]]
    water.scfenergy = -1.0
    write_result(water, (1, "a"), tmp_path / "sp.xyz")
    key, mol = read_result(tmp_path / "sp.xyz")

    assert key == (1, "a")
    assert list(mol.labels) == [11, 12, 13]
    assert mol.notes == [[1, 2, 3, 120.0], [1], [2]]
    assert mol.scfenergy == -1.0
    assert not (tmp_path / "sp.xyz.tmp").exists()


def test_runner(tmp_path, water):
    mols = gl.Molecules({(0, 0): water, (0, 1): water, 1: water})
    runner = ASEJobRunner(EMT, nprocs=2)

    opts = runner.run(mols, tmp_path, "opt", fmax=0.05)
    assert set(opts) == set(mols)
    assert all(mol.success for mol in opts.values())
    assert (tmp_path / "0/1/opt.xyz").exists()

    sps = runner.run(opts, tmp_path, "sp")
    assert all(mol.status == "SP" for mol in sps.values())
    assert set(read_results(tmp_path, "opt.xyz")) == set(mols)


def test_runner_resumes(tmp_path, water, capsys):
    mols = gl.Molecules({0: water, 1: water})
    runner = ASEJobRunner(EMT)
    runner.run(gl.Molecules({0: water}), tmp_path, "sp")
    capsys.readouterr()

    runner.run(mols, tmp_path, "sp")
    assert "sp: 1 done, 1 to run" in capsys.readouterr().out


def test_runner_records_errors(tmp_path, water):
    water.mult = None  # to_atoms で ValueError になる
    results = ASEJobRunner(EMT).run(gl.Molecules({0: water}), tmp_path, "sp")

    assert results[0].success is False
    assert results[0].status.startswith("Error: ValueError")
