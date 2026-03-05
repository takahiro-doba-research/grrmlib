import numpy as np

from .atomcoords import overlay, rotate
from .exceptions import FragmentNotFoundError, ManyFragmentsFoundError
from ..core import Molecules, GroupedMolecules


def read_connectivity(mol):
    """
    Format of notes
    ---------------
    note = [fragment]
    note = [fragment, label1, label2]
    note = [fragment, label1, label2, angle]
    
    - fragment : int
    - label1, label2 : int (for overlay and rotate functions)
    - angle : float (in degrees)
    
    Example
    -------
    connectivity = {
        0: {"labels": [1, 7], "labels_ref": [1, 2, 3], "target": 1, "angle": 180},
        1: {"labels": [2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]},
    }
    """
    connectivity = {}
    label_to_note = dict(zip(mol.labels, mol.notes))
    
    for label, note in label_to_note.items():
        fragment = note[0]
        entry = connectivity.setdefault(fragment, {"labels": []})
        entry["labels"].append(label)
        
        if len(note) >= 3:
            entry["labels_ref"] = [label, note[1], note[2]]
            entry["target"] = label_to_note[note[1]][0]
        
        if len(note) == 4:
            entry["angle"] = note[3]
    
    return connectivity


def step_to_angles(mols):
    grouped_arenes = GroupedMolecules()
    
    for name, mol in mols.items():
        connectivity = read_connectivity(mol)
        
        for fragment, data in connectivity.items():
            if "angle" in data:
                mols_new = Molecules()
                step = data["angle"]
                angles = [360 * i / step for i in range(step)]
                
                for angle in angles:
                    mol.notes = [[*n[:3], angle] if len(n) == 4 else n for n in mol.notes]
                    mols_new[angle] = mol
        
        grouped_arenes[name] = mols_new
    
    return grouped_arenes


def step_to_angles(mol):
    mols_new = Molecules()
    connectivity = read_connectivity(mol)
    
    for fragment, data in connectivity.items():
        if "angle" in data:
            mols_new = Molecules()
            step = data["angle"]
            angles = [360 * i / step for i in range(step)]
            
            for angle in angles:
                mol.notes = [[*n[:3], angle] if len(n) == 4 else n for n in mol.notes]
                mols_new[angle] = mol
    
    return mols_new


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


def connect(mol0, mol1):
    connectivity0 = read_connectivity(mol0)
    connectivity1 = read_connectivity(mol1)
    
    fragments0 = [f for f, d in connectivity1.items() if "angle" in d]
    if len(fragments0) != 1:
        raise ManyFragmentsFoundError(
            "exactly one connectable fragment must exist in mol1; "
            f"{len(fragments0)} fragments found"
        )
    fragment0 = fragments0[0]
    
    fragment1 = connectivity1[fragment0]["target"]
    if fragment1 not in connectivity0:
        raise FragmentNotFoundError(
            f"fragment {fragment1} is not connectable to mol0"
        )
    
    indices0 = mol0.labels_to_indices(connectivity0[fragment1]["labels_ref"])
    indices1 = mol1.labels_to_indices(connectivity1[fragment0]["labels_ref"])
    
    atomcoords1 = overlay(
        mol0.atomcoords,
        mol1.atomcoords,
        *[indices0[i] for i in [1, 0, 2]],
        *indices1
    )
    atomcoords1 = rotate(
        atomcoords1,
        *indices1[:2],
        connectivity1[fragment0]["angle"],
        degrees=True
    )
    mol1 = mol1.copy()
    mol1.atomcoords = atomcoords1
    
    mol0 = mol0.remove_atoms(connectivity0[fragment1]["labels"])
    mol1 = mol1.remove_atoms(connectivity1[fragment0]["labels"])
    mol0 = reset_connectable_notes(mol0)
    mol1 = reset_connectable_notes(mol1, offset=max(mol0.labels))
    return mol0.join(mol1)