import pytest
from ase.calculators.emt import EMT

from grrmlib.calculators import HARTREE_IN_EV, optimize, single_point, to_atoms


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
