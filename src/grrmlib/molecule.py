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
    
    def reset_labels(self):
        mol = self.copy()
        mol.labels = np.arange(1, len(mol.labels) + 1)
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
        return Molecule(
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
    
    def overlay(self, labels, mol_ref, labels_ref):
        """
        Overlay this molecule onto mol_ref.
        labels, labels_ref: sequences of 3 atom labels (order matters)
        """
        self.validate()
        mol_ref.validate()
        indices = self.labels_to_indices(labels)
        indices_ref = mol_ref.labels_to_indices(labels_ref)
        mol = self.copy()
        mol.atomcoords = overlay(
            mol_ref.atomcoords,
            mol.atomcoords,
            *indices_ref,
            *indices
        )
        return mol
    
    def rotate(self, labels, angle, degrees=False):
        """
        Rotate molecule around the axis defined by two atom labels.
        Right-hand rule: positive angle rotates from labels[0] -> labels[1].
        labels: (label_a, label_b)
        """
        self.validate()
        indices = self.labels_to_indices(labels)
        mol = self.copy()
        mol.atomcoords = rotate(
            mol.atomcoords,
            *indices,
            angle,
            degrees
        )
        return mol
    
    def to_connectable(self, notes):
        return ConnectableMolecule(
            name=self.name,
            functional=self.functional,
            basis_set=self.basis_set,
            comments=self.comments,
            charge=self.charge,
            mult=self.mult,
            labels=self.labels,
            symbols=self.symbols,
            atomcoords=self.atomcoords,
            notes=notes,
            scfenergy=self.scfenergy,
            afirenergy=self.afirenergy,
            zpve=self.zpve,
            grads=self.grads,
            hessian=self.hessian,
            nmeigen=self.nmeigen,
            status=self.status,
        )
    
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
            0: {"labels": [1, 2, 3, 4, 5, 6]},
            1: {"labels": [7, 8, 9, 10], "labels_ref": [7, 1, 8], "angles": [0, 2pi/3, 4pi/3]},
            2: {"labels": [11, 12, 13, 14], "labels_ref": [12, 3, 14], "angles": [0, pi]},
        }
        """
        data = {}
        
        for label, note in zip(self.labels, self.notes):
            fragment = note[0]
            entry = data.setdefault(fragment, {"labels": []})
            entry["labels"].append(label)
            
            if len(note) >= 3:
                entry["labels_ref"] = [label, note[1], note[2]]
            
            if len(note) == 4:
                k = note[3]
                entry["angles"] = [2 * np.pi * i / k for i in range(k)]
        
        return data
    
    def reset_labels(self):
        mol = super().reset_labels()
        mapping = dict(zip(self.labels, mol.labels))
        mol.notes = [
            [n[0], mapping[n[1]], mapping[n[2]], *n[3:]]
            if len(n) >= 3 else n
            for n in self.notes
        ]
        mol.data = mol.read_notes()
        return mol
    
    def overlay(self):
        pass
    
    def rotate(self):
        pass