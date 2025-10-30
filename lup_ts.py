from .data import atomic_number
from ._ircirc import _read_ircirc


class LUPTS:

    def __init__(self, name=None, irc=None):
        self.name = name
        self.irc = irc
        
    def to_gv(self, path, reverse=False):
        irc = self.irc[::-1] if reverse else self.irc
        num = len(irc)
        lines = [" #p\n", " \n"]
        
        for i, molecule in enumerate(irc):
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


def read_lup_ts(path):
    
    with open(path, "r") as f:
        lines = f.readlines()
    
    indices_ircirc = [i for i, line in enumerate(lines) if line.startswith("IRCIRC")]
    lines_ircirc = lines[indices_ircirc[0]:indices_ircirc[1] + 1]
    forward_step, _, forward_optimized, backward_step, _, backward_optimized = _read_ircirc(lines_ircirc)
    irc = [backward_optimized] + backward_step[::-1] + forward_step + [forward_optimized]
    return LUPTS(irc=irc)