from collections import UserDict
from pathlib import Path


class GroupedMolecules(UserDict):
    """
    Dictionary of Molecules.
    """
    
    def __init__(self, grouped_mols=None):
        super().__init__(grouped_mols or {})
    
    def map(self, func):
        grouped_mols_new = self.__class__()
        for group, mols in self.items():
            mols_new = mols.__class__()
            for name, mol in mols.items():
                mols_new[name] = func(mol)
            grouped_mols_new[group] = mols_new
        return grouped_mols_new
    
    def map_group(self, func):
        pass
    
    def to_gv_folder(
        self,
        folder,
        name_group="group",
        name_mol="molecule",
        basename="gaussian.com"
    ):
        folder = Path(folder)
        
        for group, mols in self.items():
            for name, mol in mols.items():
                path = (
                    folder
                    / f"{name_group}{group}"
                    / f"{name_mol}{name}"
                )
                path.mkdir(parents=True, exist_ok=False)
                mol.to_gv(path / basename)