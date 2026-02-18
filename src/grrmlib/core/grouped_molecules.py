from __future__ import annotations

from collections import UserDict
from collections.abc import Callable
from typing import TYPE_CHECKING, Mapping, Self

if TYPE_CHECKING:
    from .molecule import Molecule
    from .molecules import Molecules


class GroupedMolecules(UserDict[int, "Molecules"]):
    """
    Dictionary of Molecules.
    """
    
    def __init__(self, gmols: Mapping[int, Molecules] | None = None):
        super().__init__(gmols or {})
    
    def map(self, func: Callable[[Molecule], Molecule]) -> Self:
        gmols_new = self.__class__()
        for group, mols in self.items():
            mols_new = mols.__class__()
            for name, mol in mols.items():
                mols_new[name] = func(mol)
            gmols_new[group] = mols_new
        return gmols_new
    
    def map_groups(self, func: Callable[[Molecules], Molecules]) -> Self:
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