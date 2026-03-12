class Trajectory(list["Molecule"]):
    
    def __init__(self, mols=None):
        super().__init__(mols or [])