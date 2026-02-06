import copy

import networkx as nx
import numpy as np
from scipy.spatial import distance

from .data import covalent_radius
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
            mols[f"{self.name}F{i}"] = mol
        
        return mols
    
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


class SEQ(Molecule):
    pass


class ConnectableMolecule(Molecule):

    def connect(self, mol):
        pass