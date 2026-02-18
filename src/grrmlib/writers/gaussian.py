from pathlib import Path

from ..core import (
    Molecule,
    Molecules,
    GroupedMolecules,
    atomic_number
)


class GaussianInputWriter:
    
    def __init__(
        self,
        *,
        nprocshared: int | None = None,
        mem: str | None = None,
        chk: str | None = None,
        header: list[str] | None = None,
        footer: list[str] | None = None,
        with_notes: bool = False,
    ) -> None:
        self.nprocshared = nprocshared
        self.mem = mem
        self.chk = chk
        self.header = header
        self.footer = footer
        self.with_notes = with_notes
    
    def write(
        self,
        mol: Molecule,
        path: str = "gaussian_input.com"
    ) -> None:
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
        gmols: GroupedMolecules,
        folder: str | Path,
        prefix_group: str = "group",
        prefix_mol: str = "molecule",
        basename: str = "gaussian_input.com",
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


class GaussianOutputWriter:
    
    def write_scan(
        self,
        mols: Molecules,
        path: str | Path = "gaussian_scan.log"
    ) -> None:
        path = Path(path)
        num = len(mols)
        lines = [" #p\n", " \n"]
        
        for i, mol in enumerate(mols.values(), start=1):
            lines += [
                " GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad\n",
                "                          Input orientation:                          \n",
                " ---------------------------------------------------------------------\n",
                " Center     Atomic      Atomic             Coordinates (Angstroms)    \n",
                " Number     Number       Type             X           Y           Z   \n",
                " ---------------------------------------------------------------------\n",
                " ---------------------------------------------------------------------\n",
                " \n",
                " GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad\n",
                f" Step number   1 out of a maximum of   2 on scan point {i:5d} out of {num:5d}\n",
                " \n",
                " GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad\n",
                "                          Input orientation:                          \n",
                " ---------------------------------------------------------------------\n",
                " Center     Atomic      Atomic             Coordinates (Angstroms)    \n",
                " Number     Number       Type             X           Y           Z   \n",
                " ---------------------------------------------------------------------\n",
                *[
                    f"{j:7d} {atomic_number(s):10d}           0     {x:11.6f} {y:11.6f} {z:11.6f}\n"
                    for j, (s, (x, y, z)) in enumerate(mol.iter_atoms(), start=1)
                ],
                " ---------------------------------------------------------------------\n",
                f" SCF Done:  E({mol.functional or 'B3LYP'}) = {mol.scfenergy:15.12f}     A.U.\n",
                " \n",
                " GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad\n",
                f" Step number   2 out of a maximum of   2 on scan point {i:5d} out of {num:5d}\n",
                " \n",
            ]
        
        lines += [
            " GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad\n",
            " Normal termination of Gaussian 16\n"
        ]
        
        with path.open("x") as f:
            f.writelines(lines)