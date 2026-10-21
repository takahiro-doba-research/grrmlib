import numpy as np
import pytest

import grrmlib as gl


def test_unknown_attribute():
    with pytest.raises(TypeError):
        gl.Molecule(foo=1)


def test_allowed_extra_attribute():
    mol = gl.Molecule(scfenergy=-1.0, elapsed=2.0, optsteps=3)
    assert (mol.scfenergy, mol.elapsed, mol.optsteps) == (-1.0, 2.0, 3)


def test_validate(water):
    water.validate()

    water.symbols = ["O", "H"]
    with pytest.raises(ValueError):
        water.validate()


def test_iter_atoms(water):
    atoms = list(water.iter_atoms())
    assert [s for s, _, _ in atoms] == ["O", "H", "H"]
    assert all(n is None for _, _, n in atoms)

    water.notes = [[1], [2], [3]]
    assert [n for _, _, n in water.iter_atoms()] == [[1], [2], [3]]


def test_copy_is_deep(water):
    mol = water.copy()
    mol.atomcoords[0, 0] = 9.0
    assert water.atomcoords[0, 0] == 0.0


def test_select_and_remove_atoms_keep_labels(water):
    water.notes = [[1], [2], [3]]

    mol = water.select_atoms([1, 3])
    assert list(mol.labels) == [1, 3]
    assert mol.symbols == ["O", "H"]
    assert mol.notes == [[1], [3]]

    mol = water.remove_atoms([1])
    assert list(mol.labels) == [2, 3]
    assert mol.symbols == ["H", "H"]


def test_reset_labels(water):
    mol = water.select_atoms([2, 3]).reset_labels(offset=10)
    assert list(mol.labels) == [11, 12]


def test_labels_to_indices(water):
    mol = water.select_atoms([2, 3])
    assert mol.labels_to_indices([3, 2]) == [1, 0]

    with pytest.raises(ValueError):
        mol.labels_to_indices([1])


def test_join(water):
    other = water.reset_labels(offset=3)
    mol = water.join(other)

    assert len(mol) == 6
    assert list(mol.labels) == [1, 2, 3, 4, 5, 6]
    assert mol.symbols == ["O", "H", "H", "O", "H", "H"]


@pytest.mark.parametrize(
    ("charge", "mult", "expected"),
    [(0, 1, True), (0, 2, False), (1, 2, True), (-1, 2, True), (0, 3, True)],
)
def test_charge_mult_is_valid(water, charge, mult, expected):
    water.charge = charge
    water.mult = mult
    assert water.charge_mult_is_valid() is expected


def test_charge_mult_is_valid_rejects_non_integer(water):
    water.mult = 1.0
    with pytest.raises(gl.MultiplicityError):
        water.charge_mult_is_valid()


def test_molecules_methods(eqs):
    assert list(eqs.filter(lambda eq: eq.scfenergy < -152.75)) == [0]
    assert eqs.smallest("scfenergy").name == 0
    assert eqs.largest("scfenergy").name == 1
    assert eqs.map(lambda eq: eq.reset_labels(10))[0].labels[0] == 11
    assert list(gl.Molecules({"a": eqs[0], "b": eqs[1]}).reset_keys()) == [0, 1]


def test_grouped_molecules(eqs):
    frags = eqs.expand(gl.separate)

    assert set(frags) == {0, 1}
    assert set(frags.flatten()) == {(0, 0), (0, 1), (1, 0), (1, 1)}
    assert list(frags.reset_keys().flatten()) == [(0, 0), (0, 1), (1, 2), (1, 3)]
    assert list(frags.contract(lambda mols: mols[0])) == [0, 1]
    assert np.allclose(
        frags.map(lambda m: m.reset_labels())[0][1].labels, [1, 2, 3]
    )
