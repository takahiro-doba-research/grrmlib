from typing import Iterator

import numpy as np

from .exceptions import FragmentNotFoundError, ManyFragmentsFoundError
from ..core import Molecule, Molecules, GroupedMolecules, overlay, rotate


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


def step_to_angles(mol):
    mols = Molecules()
    connectivity = read_connectivity(mol)
    
    for fragment, data in connectivity.items():
        if "angle" in data:
            step = data["angle"]
            angles = [360 * i / step for i in range(step)]
            
            for angle in angles:
                mol_copy = mol.copy()
                mol_copy.notes = [[*n[:3], angle] if len(n) == 4 else n for n in mol.notes]
                mols[angle] = mol_copy
    
    return mols


def reset_connectable_notes(mol, offset=0):
    mol_new = mol.copy()
    mol_new.labels = np.arange(offset + 1, len(mol.labels) + offset + 1)
    mapping = dict(zip(mol.labels, mol_new.labels))
    # 参照する原子が除かれた notes は使えないので、fragment だけにする
    mol_new.notes = [
        [n[0], mapping[n[1]], mapping[n[2]], *n[3:]]
        if len(n) >= 3 and n[1] in mapping and n[2] in mapping
        else n[:1]
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
    
    labels_ref0 = connectivity0[fragment1]["labels_ref"]
    labels_ref1 = connectivity1[fragment0]["labels_ref"]

    # mol0 で残す原子に mol1 で除く原子を重ね、mol1 で残す原子を mol0 で除く原子の方向に置く
    mol1 = mol1.copy()
    mol1.atomcoords = overlay(
        mol0,
        mol1,
        *[labels_ref0[i] for i in [1, 0, 2]],
        *labels_ref1
    )
    mol1.atomcoords = rotate(
        mol1,
        *labels_ref1[:2],
        connectivity1[fragment0]["angle"],
        degrees=True
    )
    
    mol0 = mol0.remove_atoms(connectivity0[fragment1]["labels"])
    mol1 = mol1.remove_atoms(connectivity1[fragment0]["labels"])
    mol0 = reset_connectable_notes(mol0)
    mol1 = reset_connectable_notes(mol1, offset=max(mol0.labels))
    return mol0.join(mol1)


def connect_all_iter(
    seqs: Molecules,
    subs_list: list[Molecules]
) -> Iterator[tuple[tuple, Molecule]]:
    """
    Example
    -------
    seqs = Molecules({0: Molecule(), 1: Molecule(), ...})
    subs_list = [
        Molecules({(0, 0.0): Molecule(), (0, 180.0): Molecule(), ...}),
        Molecules({(0, 0.0): Molecule(), (1, 0.0): Molecule(), ...}),
        Molecules({(0, 0.0): Molecule(), (1, 0.0): Molecule(), ...}),
        ...
    ]
    keys indicate (substituent number, angle)
    """
    
    def connect_recursive(names, seq, subs_list):
        if not subs_list:
            yield names, seq
            return
        
        connected = False
        
        for name, sub in subs_list[0].items():
            try:
                mol_new = connect(seq, sub)
                connected = True
                yield from connect_recursive(
                    names + name,
                    mol_new,
                    subs_list[1:]
                )
            except FragmentNotFoundError:
                continue  # or break
        
        if not connected:
            yield from connect_recursive(
                names + (None, None),
                seq,
                subs_list[1:]
            )
    
    for name, seq in seqs.items():
        yield from connect_recursive(
            (name, ),
            seq,
            subs_list
        )


def connect_all(
    seqs: Molecules,
    subs_list: list[Molecules]
) -> Molecules:
    """
    connect_all_iter の結果を Molecules にする。
    key は (seq, 置換基 1, 角度 1, 置換基 2, 角度 2, ...)。角度は int にし、
    つながらなかった置換基は (None, None)。
    """
    return Molecules({
        tuple(int(n) if n is not None else n for n in names): mol
        for names, mol in connect_all_iter(seqs, subs_list)
    })


def build_connection_folder(folder, names):
    """
    Example
    -------
    names = [
        ('SEQ', 7286),
        ('arene', 0), ('angle', 0.0),
        ('alkene', 0), ('angle', 240.0),
        ('pg', 0), ('angle', 0.0),
        ('backbone', 5), ('angle', 240.0),
        ('pyridone', 9), ('angle', 0.0)
    ]
    """
    folder_new = Path(folder)
    
    for key, value in names:
        value = round(value) if key == "angle" else value
        folder_new /= f"{key}{value}"
    
    return folder_new


def connect_all_and_write_as_gaussian_input(
    mols_tuple,
    gmols_tuples,
    mol_method,
    folder="connection",
    basename="gaussian_input.com"
):
    writer = GaussianInputWriter(
        header=mol_method.header,
        footer=mol_method.footer,
    )
    
    for mol, names in connection_factory(mols_tuple, gmols_tuples):
        folder_new = build_connection_folder(folder, names)
        folder_new.mkdir(parents=True, exist_ok=True)
        
        try:
            writer.write(mol, folder_new / basename)
        except FileExistsError as e:
            print(e)