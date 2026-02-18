from __future__ import annotations

import copy
from typing import TYPE_CHECKING

import networkx as nx
import numpy as np
from scipy.spatial import distance

from .data import covalent_radius

if TYPE_CHECKING:
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

    def copy(self) -> Molecule:
        return copy.deepcopy(self)
    
    def iter_atoms(self, with_notes=False):
        self.validate()
        if with_notes:
            for symbol, atomcoord, note in zip(self.symbols, self.atomcoords, self.notes):
                yield symbol, atomcoord, note
        else:
            for symbol, atomcoord in zip(self.symbols, self.atomcoords):
                yield symbol, atomcoord
    
    def reset_labels(self, offset=0) -> Molecule:
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
    
    def _select_by_indices(self, indices) -> Molecule:
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
    
    def select_atoms(self, labels) -> Molecule:
        labels = set(labels)
        indices = [i for i, l in enumerate(self.labels) if l in labels]
        return self._select_by_indices(indices)
    
    def remove_atoms(self, labels) -> Molecule:
        labels = set(labels)
        indices = [i for i, l in enumerate(self.labels) if l not in labels]
        return self._select_by_indices(indices)
    
    def join(self, mol) -> Molecule:
        self.validate()
        mol.validate()
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
    
    def separate(self) -> Molecules:
        from .molecules import Molecules
        
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