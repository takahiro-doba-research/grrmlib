import copy

import numpy as np


class Molecule:
    
    def __init__(
        self,
        name=None,
        functional=None,
        basis_set=None,
        comments=None,
        charge=None,
        mult=None,
        symbols=None,
        atomcoords=None,
        notes=None,
        scfenergy=None,
        afirenergy=None,
        zpve=None,
        grads=None,
        hessian=None,
        nmeigen=None,
        status=None,
        **kwargs
    ):
        self.name = name
        self.functional = functional
        self.basis_set = basis_set
        self.comments = comments
        self.charge = charge
        self.mult = mult
        self.symbols = symbols
        self.atomcoords = atomcoords
        self.notes = notes
        self.scfenergy = scfenergy
        self.afirenergy = afirenergy
        self.zpve = zpve
        self.grads = grads
        self.hessian = hessian
        self.nmeigen = nmeigen
        self.status = status
        
        for k, v in kwargs.items():
            setattr(self, k, v)
    
    def copy(self):
        return copy.deepcopy(self)
    
    def remove_atoms(self, labels):
        indices = np.asarray(labels, dtype=int) - 1
        
        if indices.size:
            if indices.min() < 0 or indices.max() >= len(self.atomcoords):
                raise IndexError("Atom label out of range")
        
        mask = np.ones(len(self.atomcoords), dtype=bool)
        mask[indices] = False
        mol = self.copy()
        mol.symbols = [s for i, s in enumerate(self.symbols) if i + 1 not in labels]
        mask = np.ones(self.atomcoords.shape[0], dtype=bool)
        mask[np.array(sorted(labels)) - 1] = False
        mol.atomcoords = self.atomcoords[mask]
        mol.notes = (
            [n for i, n in enumerate(self.notes) if i + 1 not in labels]
            if self.notes else None
        )
        return mol
    
    def join(self, mol):
        return Molecule(
            symbols=self.symbols + mol.symbols,
            atomcoords=np.vstack([self.atomcoords, mol.atomcoords]),
            notes=self.notes + mol.notes if self.notes and mol.notes else None
        )
    
    def to_gv(self, path):
        lines = [
            f"# {self.functional or 'B3LYP'}/{self.basis_set or '6-31G'}\n",
            "\n",
            f"{self.comments or 'title'}\n",
            "\n",
            f"{self.charge or 0} {self.mult or 1}\n"
        ]
        lines += [
            f"{sym:2s}  {coord[0]:17.12f} {coord[1]:17.12f} {coord[2]:17.12f}\n"
            for sym, coord in zip(self.symbols, self.atomcoords)
        ]
        lines += ["\n"]
        
        with open(path, "w") as f:
            f.writelines(lines)
    
    def to_grrm(self, path):
        pass


class EQ(Molecule):
    pass


class PT(Molecule):
    
    def __init__(self, connection=None, **kwargs):
        super().__init__(**kwargs)
        self.connection = connection