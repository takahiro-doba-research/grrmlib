from pathlib import Path


class GaussianInputWriter:
    
    def __init__(
        self,
        *,
        nprocshared=None,
        mem=None,
        chk=None,
        header=None,
        footer=None,
        with_notes=False,
    ):
        self.nprocshared = nprocshared
        self.mem = mem
        self.chk = chk
        self.header = header
        self.footer = footer
        self.with_notes = with_notes
    
    def write(self, mol, path):
        path = Path(path)
        
        lines = []
        
        if self.nprocshared is not None:
            lines.append(f"%nprocshared={self.nprocshared}\n")
        
        if self.mem is not None:
            lines.append(f"%mem={self.mem}\n")
        
        if self.chk is not None:
            lines.append(f"%chk={self.chk}\n")
        
        if self.header:
            lines += self.header
        
        if self.with_notes:
            for s, (x, y, z), n in mol.iter_atoms(with_notes=True):
                lines.append(f"{s:2s}  {x:17.12f} {y:17.12f} {z:17.12f} {' '.join(map(str, n))}\n")
        else:
            for s, (x, y, z) in mol.iter_atoms():
                lines.append(f"{s:2s}  {x:17.12f} {y:17.12f} {z:17.12f}\n")
        
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
        basename="gaussian_input.com",
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


class GaussianOutputWriter:
    pass