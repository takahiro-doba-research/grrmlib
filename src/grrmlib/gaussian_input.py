import numpy as np

from .molecule import Molecule


def read_gaussian_input(path):
    
    with open(path, "r") as f:
        lines = f.readlines()
    
    indices_blank = [i for i, line in enumerate(lines) if line == "\n"]
    header = lines[:indices_blank[1] + 2]
    footer = lines[indices_blank[2]:]
    lines_coord = lines[indices_blank[1] + 2:indices_blank[2]]
    labels = np.arange(1, len(lines_coord) + 1)
    symbols = [line.split()[0] for line in lines_coord]
    atomcoords = np.array([list(map(float, line.split()[1:4])) for line in lines_coord])
    
    notes = [list(map(int, line.split()[4:])) for line in lines_coord]
    if all(not note for note in notes):
        notes = None
    
    return Molecule(
        labels=labels,
        symbols=symbols,
        atomcoords=atomcoords,
        notes=notes,
        header=header,
        footer=footer,
    )