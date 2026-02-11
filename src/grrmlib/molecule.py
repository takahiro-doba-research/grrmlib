import copy

import networkx as nx
import numpy as np
from scipy.spatial import distance

from .data import covalent_radius
from .geometry import overlay, rotate
from .molecules import Molecules


class Molecule:
    
    def __init__(
        self,
        name=None,
        functional=None,
        basis_set=None,
        comments=None,
        charge=None,
        mult=None,
        labels=None,
        symbols=None,
        atomcoords=None,
        notes=None,
        scfenergy=None,
        afirenergy=None,
        zpve=None,
        grads=None,
        hessian=None,
        nmeigen=None,
        status=None,
        **kwargs
    ):
        self.name = name
        self.functional = functional
        self.basis_set = basis_set
        self.comments = comments
        self.charge = charge
        self.mult = mult
        self.labels = labels
        self.symbols = symbols
        self.atomcoords = atomcoords
        self.notes = notes
        self.scfenergy = scfenergy
        self.afirenergy = afirenergy
        self.zpve = zpve
        self.grads = grads
        self.hessian = hessian
        self.nmeigen = nmeigen
        self.status = status
        
        for k, v in kwargs.items():
            setattr(self, k, v)
    
    def validate(self):
        if self.labels is None or self.symbols is None or self.atomcoords is None:
            raise ValueError("labels, symbols, and atomcoords must be set")
        
        n = len(self.labels)
        if len(self.symbols) != n:
            raise ValueError("symbols length mismatch")
        if len(self.atomcoords) != n:
            raise ValueError("atomcoords length mismatch")
        if self.atomcoords.ndim != 2 or self.atomcoords.shape[1] != 3:
            raise ValueError("atomcoords must be (N, 3)")
        if self.notes is not None:
            if len(self.notes) != n:
                raise ValueError("notes length mismatch")
    
    def copy(self):
        return copy.deepcopy(self)
    
    def reset_labels(self, offset=0):
        mol = self.copy()
        mol.labels = np.arange(offset + 1, len(mol.labels) + offset + 1)
        return mol
    
    def labels_to_indices(self, labels):
        label_to_index = {l: i for i, l in enumerate(self.labels)}
        try:
            return [label_to_index[l] for l in labels]
        except KeyError as e:
            raise ValueError(f"label not found: {e.args[0]}")
    
    def label_to_index(self, label):
        return self.labels_to_indices([label])[0]
    
    def _select_by_indices(self, indices):
        indices = list(indices)
        mol = self.copy()
        mol.labels = mol.labels[indices]
        mol.symbols = [self.symbols[i] for i in indices]
        mol.atomcoords = mol.atomcoords[indices]
        mol.notes = (
            [self.notes[i] for i in indices]
            if self.notes is not None
            else None
        )
        return mol
    
    def select_atoms(self, labels):
        labels = set(labels)
        indices = [i for i, l in enumerate(self.labels) if l in labels]
        return self._select_by_indices(indices)
    
    def remove_atoms(self, labels):
        labels = set(labels)
        indices = [i for i, l in enumerate(self.labels) if l not in labels]
        return self._select_by_indices(indices)
    
    def join(self, mol):
        self.validate()
        mol.validate()
        mol = mol.reset_labels(max(self.labels))
        return self.__class__(
            labels=np.concatenate([self.labels, mol.labels]),
            symbols=self.symbols + mol.symbols,
            atomcoords=np.vstack([self.atomcoords, mol.atomcoords]),
            notes=self.notes + mol.notes if self.notes and mol.notes else None
        )
    
    def get_adj_matrix(self, threshold=1.25):
        self.validate()
        D = distance.cdist(self.atomcoords, self.atomcoords)
        r = np.array([covalent_radius(s) for s in self.symbols])
        R = r[:, None] + r[None, :]
        A = D < R * threshold
        np.fill_diagonal(A, False)
        return A
    
    def separate(self):
        self.validate()
        A = self.get_adj_matrix()
        G = nx.from_numpy_array(A)
        components = sorted(nx.connected_components(G), key=len, reverse=True)
        
        mols = Molecules()
        for i, component in enumerate(components):
            indices = sorted(component)
            mol = self._select_by_indices(indices)
            mols[i] = mol
        
        return mols
    
    def to_connectable(self, notes):
        mol = self.copy()
        mol.__class__ = ConnectableMolecule
        mol.notes = notes
        mol.data = mol.read_notes()
        return mol
    
    def to_gv(self, path):
        self.validate()
        
        lines = [
            f"# {self.functional or 'B3LYP'}/{self.basis_set or '6-31G'}\n",
            "\n",
            f"{self.comments or 'title'}\n",
            "\n",
            f"{self.charge or 0} {self.mult or 1}\n"
        ]
        
        if self.notes is not None:
            lines += [
                f"{s:2s}  {x:17.12f} {y:17.12f} {z:17.12f} {' '.join(map(str, n))}\n"
                for s, (x, y, z), n in zip(self.symbols, self.atomcoords, self.notes)
            ]
        else:
            lines += [
                f"{s:2s}  {x:17.12f} {y:17.12f} {z:17.12f}\n"
                for s, (x, y, z) in zip(self.symbols, self.atomcoords)
            ]
        
        lines += ["\n"]
        
        with open(path, "w") as f:
            f.writelines(lines)
    
    def to_grrm(self, path):
        pass


class EQ(Molecule):
    pass


class PT(Molecule):
    
    def __init__(self, connection=None, **kwargs):
        super().__init__(**kwargs)
        self.connection = connection


class ConnectableMolecule(Molecule):
    """
    Format of notes
    ---------------
    note = [fragment]
    note = [fragment, label1, label2]
    note = [fragment, label1, label2, k]
    
    - fragment : int
    - label1, label2 : int (for overlay and rotate functions)
    - k : int (number of angular configurations)
    """
    
    def __init__(self, data=None, **kwargs):
        super().__init__(**kwargs)
        self.data = self.read_notes() if data is None else data
    
    def read_notes(self):
        """
        Example
        -------
        data = {
            0: {"labels": [1, 7], "labels_ref": [1, 2, 3], "target": 1, "angles": [0.0]},
            1: {"labels": [2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]},
        }
        """
        data = {}
        
        for label, note in zip(self.labels, self.notes):
            fragment = note[0]
            entry = data.setdefault(fragment, {"labels": []})
            entry["labels"].append(label)
            
            if len(note) >= 3:
                entry["labels_ref"] = [label, note[1], note[2]]
                entry["target"] = [n[0] for l, n in zip(self.labels, self.notes) if l == note[1]][0]
            
            if len(note) == 4:
                k = note[3]
                entry["angles"] = [360 * i / k for i in range(k)]
        
        return data
    
    def reset_labels(self, offset=0):
        mol = super().reset_labels(offset)
        mapping = dict(zip(self.labels, mol.labels))
        mol.notes = [
            [n[0], mapping[n[1]], mapping[n[2]], *n[3:]]
            if len(n) >= 3 else n
            for n in self.notes
        ]
        mol.data = mol.read_notes()
        return mol
        
    def connect(self, mol):
        self.validate()
        mol.validate()
        
        fragment_self = [f for f, d in mol.data.items() if "angles" in d]
        
        if len(fragment_self) != 1:
            raise ValueError("exactly one connectable fragment must exist in mol")
        
        fragment_self = fragment_self[0]
        fragment_mol = mol.data[fragment_self]["target"]
        
        if fragment_mol not in self.data:
            raise ValueError(f"fragment {fragment_mol} is not connectable to this molecule")
        
        indices_self = self.labels_to_indices(self.data[fragment_mol]["labels_ref"])
        indices_mol = mol.labels_to_indices(mol.data[fragment_self]["labels_ref"])
        
        atomcoords_overlay = overlay(
            self.atomcoords,
            mol.atomcoords,
            *[indices_self[i] for i in [1, 0, 2]],
            *indices_mol
        )
        
        mols = Molecules()
        angles = mol.data[fragment_self]["angles"]
        self_new = self.remove_atoms(self.data[fragment_mol]["labels"])
        
        for angle in angles:
            mol_new = mol.copy()
            mol_new.atomcoords = rotate(
                atomcoords_overlay,
                *indices_mol[:2],
                angle,
                degrees=True
            )
            mol_new = mol_new.remove_atoms(mol.data[fragment_self]["labels"])
            mols[angle] = self_new.join(mol_new)
        
        return mols