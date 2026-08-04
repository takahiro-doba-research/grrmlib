from pathlib import Path

import cclib
import numpy as np
from periodictable import elements

from ..core import Molecule, Molecules


# ----------------------
# Functions for reading.
# ----------------------
def read_gaussian_output(path: str | Path) -> Molecule:
    path = Path(path)
    parser = cclib.io.ccopen(path)
    data = parser.parse()
    success = data.metadata["success"]
    
    if success:
        symbols = [str(elements[n]) for n in data.atomnos]
        labels = np.arange(1, len(symbols) + 1)
        return Molecule(
            charge=data.charge,
            mult=data.mult,
            labels=labels,
            symbols=symbols,
            atomcoords=data.atomcoords[-1],
            scfenergy=data.scfenergies[-1],
            success=success,
        )
    
    else:
        return Molecule(success=success)


def read_gaussian_outputs(
    folder: str | Path,
    basename: str | Path
) -> Molecules:
    folder = Path(folder)
    basename = str(basename)
    mols = Molecules()
    
    for path in sorted(folder.rglob(basename)):
        key = tuple(map(int, path.parent.relative_to(folder).parts))
        key = key[0] if len(key) == 1 else key
        mols[key] = read_gaussian_output(path)
    
    return mols


# ----------------------
# Functions for writing.
# ----------------------
def write_gaussian_scan(
    mols: Molecules,
    path: str | Path,
    *,
    overwrite: bool = False
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    
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
                f"{j:7d} {elements.symbol(s).number:10d}           0     {x:11.6f} {y:11.6f} {z:11.6f}\n"
                for j, (s, (x, y, z), _) in enumerate(mol.iter_atoms(), start=1)
            ],
            " ---------------------------------------------------------------------\n",
            f" SCF Done:  E = {mol.scfenergy:15.12f}     A.U.\n",
            " \n",
            " GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad\n",
            f" Step number   2 out of a maximum of   2 on scan point {i:5d} out of {num:5d}\n",
            " \n",
        ]
    
    lines += [
        " GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad\n",
        " Normal termination of Gaussian 16\n"
    ]
    
    with path.open("w" if overwrite else "x") as f:
        f.writelines(lines)