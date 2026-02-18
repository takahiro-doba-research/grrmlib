from __future__ import annotations

from collections import UserDict
from collections.abc import Callable
from typing import TYPE_CHECKING, Hashable, Mapping, Self

from .data import atomic_number

if TYPE_CHECKING:
    from .molecule import Molecule
    from .grouped_molecules import GroupedMolecules


class Molecules(UserDict[Hashable, "Molecule"]):
    """
    Dictionary of Molecule.
    """
    
    def __init__(self, mols: Mapping[Hashable, Molecule] | None = None) -> None:
        super().__init__(mols or {})
    
    def copy(self) -> Self:
        return copy.deepcopy(self)
    
    def map(self, func: Callable[[Molecule], Molecule]) -> Self:
        mols_new = self.__class__()
        for name, mol in self.items():
            mols_new[name] = func(mol)
        return mols_new
    
    def filter(self, predicate: Callable[[Molecule], bool]) -> Self:
        return self.__class__({k: v for k, v in self.items() if predicate(v)})
    
    def smallest(self, attr: str) -> Molecule:
        return min(self.values(), key=lambda m: getattr(m, attr))
    
    def largest(self, attr: str) -> Molecule:
        return max(self.values(), key=lambda m: getattr(m, attr))
    
    def separate(self) -> GroupedMolecules:
        from .grouped_molecules import GroupedMolecules

        gmols = GroupedMolecules()
        
        for name, mol in self.items():
            gmols[name] = mol.separate()
        
        return gmols
    
    def cluster(
        self,
        predicate: Callable[[Molecule, Molecule], bool]
    ) -> GroupedMolecules:
        from .grouped_molecules import GroupedMolecules

        gmols = GroupedMolecules()
        index = 0
        
        for name, mol in self.items():
            for group, mols in gmols.items():
                mol_rep = next(iter(mols.values()))
                if predicate(mol, mol_rep):
                    gmols[group][name] = mol
                    break
            else:
                gmols[index] = Molecules({name: mol})
                index += 1
        
        return gmols
    
    def reset_keys(self) -> Self:
        return self.__class__(
            {i: mol.copy() for i, mol in enumerate(self.values())}
        )
    
    def set_group(
        self,
        predicate: Callable[[Molecule, Molecule], bool]
    ) -> Self:
        mols = self.copy()
        group = 0
        mols_rep = []
        
        for mol in mols.values():
            for mol_rep in mols_rep:
                if predicate(mol, mol_rep):
                    mol.group = mol_rep.group
                    break
            else:
                mol.group = group
                mols_rep.append(mol)
                group += 1
        
        return mols