from pathlib import Path


class ConnectableWriter:
    
    def __init__(
        self,
        *,
        header=None,
        footer=None,
    ):
        self.header = header
        self.footer = footer
    
    def write(self, mol, path):
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
        gmols,
        folder,
        prefix_group="group",
        prefix_mol="molecule",
        basename="connectable.com",
    ):
        folder = Path(folder)
        
        for group, mols in gmols.items():
            for name, mol in mols.items():
                folder_new = (
                    folder
                    / f"{prefix_group}{group}"
                    / f"{prefix_mol}{name}"
                )
                folder_new.mkdir(parents=True, exist_ok=True)
                self.write(mol, folder_new / basename)