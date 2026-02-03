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
    
    def copy(self):
        return copy.deepcopy(self)
    
    def _apply_atom_mask(self, mask):
        mask = np.asarray(mask, dtype=bool)
        
        n = len(self.atomcoords)
        if len(mask) != n or len(self.symbols) != n:
            raise ValueError("mask / symbols / atomcoords length mismatch")
        if self.notes is not None and len(self.notes) != n:
            raise ValueError("notes length mismatch")
        
        mol = self.copy()
        mol.symbols = [s for s, m in zip(self.symbols, mask) if m]
        mol.atomcoords = self.atomcoords[mask]
        mol.notes = [n for n, m in zip(self.notes, mask) if m] if self.notes else None
        return mol
    
    def remove_atoms(self, labels):
        indices = np.asarray(labels, dtype=int) - 1
        mask = np.ones(len(self.atomcoords), dtype=bool)
        mask[indices] = False
        return self._apply_atom_mask(mask)
    
    def select_atoms(self, labels):
        indices = np.asarray(labels, dtype=int) - 1
        mask = np.zeros(len(self.atomcoords), dtype=bool)
        mask[indices] = True
        return self._apply_atom_mask(mask)
    
    def join(self, mol):
        return Molecule(
            symbols=self.symbols + mol.symbols,
            atomcoords=np.vstack([self.atomcoords, mol.atomcoords]),
            notes=self.notes + mol.notes if self.notes and mol.notes else None
        )
    
    def get_adj_matrix(self, threshold=1.25):
        D = distance.cdist(self.atomcoords, self.atomcoords)
        r = np.array([covalent_radius(s) for s in self.symbols])
        R = r[:, None] + r[None, :]
        A = D < R * threshold
        np.fill_diagonal(A, False)
        return A
    
    def separate(self):
        A = self.get_adj_matrix()
        G = nx.from_numpy_array(A)
        indices_list = sorted(nx.connected_components(G), key=len, reverse=True)
        mols = {
            f"F{i}": self.select_atoms([index + 1 for index in indices])
            for i, indices in enumerate(indices_list)
        }
        return Molecules(mols)
    
    def to_gv(self, path):
        lines = [
            f"# {self.functional or 'B3LYP'}/{self.basis_set or '6-31G'}\n",
            "\n",
            f"{self.comments or 'title'}\n",
            "\n",
            f"{self.charge or 0} {self.mult or 1}\n"
        ]
        lines += [
            f"{sym:2s}  {coord[0]:17.12f} {coord[1]:17.12f} {coord[2]:17.12f}\n"
            for sym, coord in zip(self.symbols, self.atomcoords)
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