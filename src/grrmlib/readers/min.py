from ..core.data import atomic_number
from ._optopt import _read_optopt


class MIN:

    def __init__(self, name=None, itrs=None, optimized=None):
        self.name = name
        self.itrs = itrs
        self.optimized = optimized

    def to_gv(self, path):
        num = len(self.itrs)
        lines = [" #p\n", " \n"]
        
        for i, molecule in enumerate(self.itrs):
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
                f" Step number   1 out of a maximum of   2 on scan point {i+1:5d} out of {num:5d}\n",
                " \n",
                " GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad\n",
                "                          Input orientation:                          \n",
                " ---------------------------------------------------------------------\n",
                " Center     Atomic      Atomic             Coordinates (Angstroms)    \n",
                " Number     Number       Type             X           Y           Z   \n",
                " ---------------------------------------------------------------------\n",
                *[
                    f"{i+1:7d} {atomic_number(sym):10d}           0     {coord[0]:11.6f} {coord[1]:11.6f} {coord[2]:11.6f}\n"
                    for i, (sym, coord) in enumerate(zip(molecule.symbols, molecule.atomcoords))
                ],
                " ---------------------------------------------------------------------\n",
                f" SCF Done:  E({molecule.functional}) = {molecule.scfenergy:15.12f}     A.U.\n",
                " \n",
                " GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad\n",
                f" Step number   2 out of a maximum of   2 on scan point {i+1:5d} out of {num:5d}\n",
                " \n",
            ]
        
        lines += [
            " GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad\n",
            " Normal termination of Gaussian 16\n"
        ]
        
        with open(path, "w") as f:
            f.writelines(lines)


def read_min(path):
    
    with open(path, "r") as f:
        lines = f.readlines()
    
    indices_optopt = [i for i, line in enumerate(lines) if line.startswith("OPTOPT")]
    lines_optopt = lines[indices_optopt[0]:indices_optopt[1] + 1]
    itrs, optimized = _read_optopt(lines_optopt)
    return MIN(
        itrs=itrs,
        optimized=optimized
    )