from __future__ import annotations

from collections import UserDict
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .molecules import Molecules


class GroupedMolecules(UserDict):
    """
    Dictionary of Molecules.
    """
    
    def __init__(self, gmols=None):
        super().__init__(gmols or {})
    
    def map(self, func) -> Molecules:
        gmols_new = self.__class__()
        for group, mols in self.items():
            mols_new = mols.__class__()
            for name, mol in mols.items():
                mols_new[name] = func(mol)
            gmols_new[group] = mols_new
        return gmols_new
    
    def map_group(self, func):
        pass
    
    def flatten(self) -> Molecules:
        from .molecules import Molecules

        mols_new = Molecules()
        index = 0
        
        for _, mols in self.items():
            for _, mol in mols.items():
                mols_new[index] = mol
                index += 1
        
        return mols_new