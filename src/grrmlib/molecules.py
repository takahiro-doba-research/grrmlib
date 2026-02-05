import copy
from collections import UserDict

import numpy as np

from .data import atomic_number
from .geometry import get_distance


class Molecules(UserDict):
    
    def __init__(self, mols=None):
        super().__init__(mols or {})
    
    def copy(self):
        return copy.deepcopy(self)
    
    def set_group(self, predicate):
        mols = self.copy()
        group = 0
        mols_rep = []
        
        for mol in mols.values():
            for mol_rep in mols_rep:
                if predicate(mol, mol_rep):
                    mol.group = mol_rep.group
                    break
            else:
                mol.group = f"G{group}"
                mols_rep.append(mol)
                group += 1
        
        return mols
    
    def separate(self):
        index = 0
        seqs_all = SEQs()
        
        for eq in self.values():
            seqs = eq.separate()
            for seq in seqs.values():
                seqs_all[f"SEQ{index}"] = seq
                index += 1
        
        return seqs_all
        
    def distance_longer(self, label0, label1, distance):
        mols_ = {
            k: mol for k, mol in self.items()
            if distance < get_distance(mol.atomcoords, label0, label1)
        }
        return Molecules(mols_)
    
    def distance_shorter(self, label0, label1, distance):
        mols_ = {
            k: mol for k, mol in self.items()
            if get_distance(mol.atomcoords, label0, label1) < distance
        }
        return Molecules(mols_)
    
    def distance_between(self, label0, label1, distance0, distance1):
        mols_ = {
            k: mol for k, mol in self.items()
            if distance0 < get_distance(mol.atomcoords, label0, label1) < distance1
        }
        return Molecules(mols_)
    
    def filter(self, predicate):
        return Molecules({k: v for k, v in self.items() if predicate(v)})
    
    def smallest(self, attr):
        return min(self.values(), key=lambda m: getattr(m, attr))
    
    def largest(self, attr):
        return max(self.values(), key=lambda m: getattr(m, attr))
    
    def to_gv(self, path):
        num = len(self)
        lines = [" #p\n", " \n"]
        
        for i, mol in enumerate(self.values()):
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
                    for i, (sym, coord) in enumerate(zip(mol.symbols, mol.atomcoords))
                ],
                " ---------------------------------------------------------------------\n",
                f" SCF Done:  E({mol.functional or 'B3LYP'}) = {mol.scfenergy:15.12f}     A.U.\n",
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


class EQs(Molecules):
    
    def __init__(self, eqs=None):
        super().__init__(mols=eqs)


class PTs(Molecules):
    
    def __init__(self, pts=None):
        super().__init__(mols=pts)


class SEQs(Molecules):
    
    def __init__(self, seqs=None):
        super().__init__(mols=seqs)

    def set_group(self, predicate):
        raise NotImplementedError(
            "SEQs does not support set_group(); use set_sgroup() instead"
        )
    
    def set_sgroup(self, predicate):
        mols = self.copy()
        sgroup = 0
        mols_rep = []
        
        for mol in mols.values():
            for mol_rep in mols_rep:
                if predicate(mol, mol_rep):
                    mol.sgroup = mol_rep.sgroup
                    break
            else:
                mol.sgroup = f"SG{sgroup}"
                mols_rep.append(mol)
                sgroup += 1
        
        return mols