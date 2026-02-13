import copy

import numpy as np


def with_notes_from(mol0, mol1):
    mol0_new = mol0.copy()
    mol0_new.notes = copy.deepcopy(mol1.notes)
    return mol0_new


def reset_connectable_notes(mol, offset=0):
    mol_new = mol.copy()
    mol_new.labels = np.arange(offset + 1, len(mol.labels) + offset + 1)
    mapping = dict(zip(mol.labels, mol_new.labels))
    mol_new.notes = [
        [n[0], mapping[n[1]], mapping[n[2]], *n[3:]]
        if len(n) >= 3 else n
        for n in mol.notes
    ]
    return mol_new