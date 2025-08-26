class Molecule:
    
    def __init__(
        self,
        name=None,
        functional="B3LYP",
        basis_set="6-31G",
        comments="title",
        charge=0,
        mult=1,
        symbols=None,
        atomcoords=None,
        scfenergy=None,
        afirenergy=None,
        zpve=None,
        grads=None,
        hessian=None,
        nmeigen=None,
        connection=None,
        status=None,
    ):
        self.name = name
        self.functional = functional
        self.basis_set = basis_set
        self.comments = comments
        self.charge = charge
        self.mult = mult
        self.symbols = symbols
        self.atomcoords = atomcoords
        self.scfenergy = scfenergy
        self.afirenergy = afirenergy
        self.zpve = zpve
        self.grads = grads
        self.hessian = hessian
        self.nmeigen = nmeigen
        self.connection = connection
        self.status = status
    
    def to_gv(self, path):
        lines = [
            f"# {self.functional}/{self.basis_set}\n",
            "\n",
            f"{self.comments}\n",
            "\n",
            f"{self.charge} {self.mult}\n"
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