from __future__ import annotations

import copy
from typing import (
    TYPE_CHECKING,
    Any,
    Iterable,
    Iterator,
    Self,
    Sequence
)

import networkx as nx
import numpy as np
from periodictable import elements
from scipy.spatial import distance

from .data import covalent_radius

if TYPE_CHECKING:
    from .molecules import Molecules


class Molecule:
    """
    See cclib for the attribute names.
    (https://cclib.github.io/data.html)
    """
    _allowed_extra = {
        "afirenergy",
        "chk",
        "enthalpy",
        "enthalpy_correction",
        "enthalpy_correction_repl",
        "enthalpy_repl",
        "freeenergy",
        "freeenergy_correction",
        "freeenergy_correction_repl",
        "freeenergy_repl",
        "grads",
        "hessian",
        "infile",
        "mem",
        "nmeigen",
        "nprocshared",
        "options",
        "rotenergy",
        "rotentropy",
        "route",
        "scfenergy",
        "scfentropy",
        "scfzpvenergy",
        "scfzpvenergy_repl",
        "status",
        "success",
        "title",
        "transenergy",
        "transentropy",
        "vibenergy",
        "vibenergy_repl",
        "vibentropy",
        "vibentropy_repl",
        "vibfreqs",
        "zpve",
        "zpve_repl",
    }
    
    def __init__(
        self,
        name: str | None = None,
        charge: int | None = None,
        mult: int | None = None,
        labels: np.ndarray | None = None,
        symbols: Sequence[str] | None = None,
        atomcoords: np.ndarray | None = None,
        notes: Sequence[Sequence[int]] | None = None,
        **kwargs: Any
    ) -> None:
        
        self.name = name
        self.charge = charge
        self.mult = mult
        self.labels = labels
        self.symbols = symbols
        self.atomcoords = atomcoords
        self.notes = notes
        
        for k, v in kwargs.items():
            if not hasattr(self, k) and k not in self._allowed_extra:
                raise TypeError(f"Unknown attribute: {k}")
            setattr(self, k, v)
    
    def __len__(self) -> int:
        return len(self.atomcoords)
    
    def validate(self) -> None:
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
    
    def copy(self) -> Self:
        return copy.deepcopy(self)
    
    def with_attr_from(self, mol: Self, attr: str) -> Self:
        self_new = self.copy()
        setattr(self_new, attr, getattr(mol, attr))
        return self_new
    
    def iter_atoms(self, with_notes: bool = False) -> Iterator[
        tuple[str, np.ndarray] | tuple[str, np.ndarray, Sequence[int]]
    ]:
        self.validate()
        if with_notes:
            for symbol, atomcoord, note in zip(self.symbols, self.atomcoords, self.notes):
                yield symbol, atomcoord, note
        else:
            for symbol, atomcoord in zip(self.symbols, self.atomcoords):
                yield symbol, atomcoord
    
    def reset_labels(self, offset: int = 0) -> Self:
        mol = self.copy()
        mol.labels = np.arange(offset + 1, len(mol.labels) + offset + 1)
        return mol
    
    def labels_to_indices(self, labels: Iterable[int]) -> list[int]:
        label_to_index = {l: i for i, l in enumerate(self.labels)}
        try:
            return [label_to_index[l] for l in labels]
        except KeyError as e:
            raise ValueError(f"label not found: {e.args[0]}")
    
    def label_to_index(self, label: int) -> int:
        return self.labels_to_indices([label])[0]
    
    def _select_by_indices(self, indices: Iterable[int]) -> Self:
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
    
    def select_atoms(self, labels: Iterable[int]) -> Self:
        labels = set(labels)
        indices = [i for i, l in enumerate(self.labels) if l in labels]
        return self._select_by_indices(indices)
    
    def remove_atoms(self, labels: Iterable[int]) -> Self:
        labels = set(labels)
        indices = [i for i, l in enumerate(self.labels) if l not in labels]
        return self._select_by_indices(indices)
    
    def join(self, mol: Self) -> Self:
        self.validate()
        mol.validate()
        return self.__class__(
            name=self.name,
            charge=self.charge,
            mult=self.mult,
            labels=np.concatenate([self.labels, mol.labels]),
            symbols=self.symbols + mol.symbols,
            atomcoords=np.vstack([self.atomcoords, mol.atomcoords]),
            notes=self.notes + mol.notes if self.notes and mol.notes else None
        )
    
    def charge_mult_is_valid(self) -> bool:
        from ..operations import ChargeError, MultiplicityError
        
        if not isinstance(self.charge, int):
            raise ChargeError("The charge must be an integer.")
        
        if not isinstance(self.mult, int) or self.mult < 1:
            raise MultiplicityError("The multiplicity must be a positive integer.")
        
        total_Z = sum(elements.symbol(s).number for s in self.symbols)
        total_e = total_Z - self.charge
        
        return total_e % 2 != self.mult % 2
    
    def get_distance(self, label0: int, label1: int) -> float:
        from ..operations import get_distance
        indices = self.labels_to_indices([label0, label1])
        return get_distance(self.atomcoords, *indices)
    
    def get_adj_matrix(self, threshold: float = 1.25) -> np.ndarray:
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
        #components = sorted(nx.connected_components(G), key=len, reverse=True)
        components = sorted(nx.connected_components(G), key=min)
        
        mols = Molecules()
        for i, component in enumerate(components):
            indices = sorted(component)
            mol = self._select_by_indices(indices)
            mols[i] = mol
        
        return mols