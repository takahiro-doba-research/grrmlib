from pathlib import Path

from ..core import Molecule, Molecules


class ConnectableWriter:
    
    def __init__(
        self,
        *,
        header: list[str] = None,
        footer: list[str] = None,
    ) -> None:
        self.header = header
        self.footer = footer
    
    def write(
        self,
        mol: Molecule,
        path: str | Path = "connectable.com"
    ) -> None:
        path = Path(path)
        
        lines = []
        
        if self.header:
            lines += self.header
        
        for s, (x, y, z), n in mol.iter_atoms(with_notes=True):
            lines.append(f"{s:2s}  {x:17.12f} {y:17.12f} {z:17.12f} {' '.join(map(str, n))}\n")
        
        if self.footer:
            lines += self.footer
        
        with path.open("x") as f:
            f.writelines(lines)
    
    def write_grouped(
        self,
        gmols: Molecules,
        folder: str | Path,
        prefix_group: str = "group",
        prefix_mol: str = "molecule",
        basename: str = "connectable.com",
    ) -> None:
        folder = Path(folder)
        
        for group, mols in gmols.items():
            for name, mol in mols.items():
                folder_new = (
                    folder
                    / f"{prefix_group}{group}"
                    / f"{prefix_mol}{name}"
                )
                folder_new.mkdir(parents=True, exist_ok=True)
                try:
                    self.write(mol, folder_new / basename)
                except FileExistsError as e:
                    print(e)