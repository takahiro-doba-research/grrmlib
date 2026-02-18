from __future__ import annotations

import copy
from collections import UserDict
from collections.abc import Callable
from typing import TYPE_CHECKING, Hashable, Mapping, Self

if TYPE_CHECKING:
    from .molecule import Molecule
    from .molecules import Molecules


class GroupedMolecules(UserDict[Hashable, "Molecules"]):
    """
    Dictionary of Molecules.
    """
    
    def __init__(self, gmols: Mapping[Hashable, Molecules] | None = None):
        super().__init__(gmols or {})

    def copy(self) -> Self:
        return copy.deepcopy(self)
    
    def map(self, func: Callable[[Molecule], Molecule]) -> Self:
        gmols_new = self.__class__()
        for group, mols in self.items():
            mols_new = mols.__class__()
            for name, mol in mols.items():
                mols_new[name] = func(mol)
            gmols_new[group] = mols_new
        return gmols_new
    
    def map_groups(self, func: Callable[[Molecules], Molecules]) -> Self:
        return self.__class__(
            {group: func(mols) for group, mols in self.items()}
        )
    
    def flatten(self) -> Molecules:
        from .molecules import Molecules
        
        mols_new = Molecules()
        
        for group, mols in self.items():
            for name, mol in mols.items():
                mols_new[(group, name)] = mol.copy()
        
        return mols_new