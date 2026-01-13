import re
import numpy as np

from .molecule import Molecule
from ._optopt import _read_optopt


def _read_ircirc(lines):
    """
    Read IRCIRC section.
    """
    index_forward = [i for i, line in enumerate(lines) if line.startswith("IRC FOLLOWING (FORWARD)")][0]
    index_backward = [i for i, line in enumerate(lines) if line.startswith("IRC FOLLOWING (BACKWARD)")][0]
    index_profile = [i for i, line in enumerate(lines) if line.startswith("Energy profile along IRC")][0]
    
    """
    Read forward IRC
    """
    lines_forward = lines[index_forward:index_backward]
    indices_step = [i for i, line in enumerate(lines_forward) if line.startswith("# STEP")]
    indices_optopt = [i for i, line in enumerate(lines_forward) if line.startswith("OPTOPT")]
    indices_step.append(indices_optopt[0] - 2)
    molecules_forward_step = []
    
    for i, j in zip(indices_step, indices_step[1:]):
        lines_step = lines_forward[i:j]
        name = "STEP" + re.search(r"STEP (\d+)", lines_step[0]).group(1)
        index_energy = [i for i, line in enumerate(lines_step) if "ENERGY" in line][0]
        lines_coord = lines_step[1:index_energy]
        symbols = [line.split()[0] for line in lines_coord]
        atomcoords = np.array([list(map(float, line.split()[1:4])) for line in lines_coord])
        scfenergy = float(re.search(r"ENERGY\s*=\s*(-?\d+\.?\d+)", lines_step[index_energy]).group(1))
        mult = float(re.search(r"Spin\(\*\*2\)\s*=\s*(-?\d+\.?\d+)", lines_step[index_energy + 1]).group(1))
        molecule = Molecule(
            name=name,
            mult=mult,
            symbols=symbols,
            atomcoords=atomcoords,
            scfenergy=scfenergy,
        )
        molecules_forward_step.append(molecule)
    
    molecules_forward_itr, molecule_forward_optimized = _read_optopt(lines_forward[indices_optopt[0]:indices_optopt[1] + 1])
    
    #indices_freqfreq = [i for i, line in enumerate(lines_forward) if line.startswith("FREQFREQ")]
    
    """
    Read backward IRC
    """
    lines_backward = lines[index_backward:index_profile - 1]
    indices_step = [i for i, line in enumerate(lines_backward) if line.startswith("# STEP")]
    indices_optopt = [i for i, line in enumerate(lines_backward) if line.startswith("OPTOPT")]
    indices_step.append(indices_optopt[0] - 2)
    molecules_backward_step = []
    
    for i, j in zip(indices_step, indices_step[1:]):
        lines_step = lines_backward[i:j]
        name = "STEP" + re.search(r"STEP (\d+)", lines_step[0]).group(1)
        index_energy = [i for i, line in enumerate(lines_step) if "ENERGY" in line][0]
        lines_coord = lines_step[1:index_energy]
        symbols = [line.split()[0] for line in lines_coord]
        atomcoords = np.array([list(map(float, line.split()[1:4])) for line in lines_coord])
        scfenergy = float(re.search(r"ENERGY\s*=\s*(-?\d+\.?\d+)", lines_step[index_energy]).group(1))
        mult = float(re.search(r"Spin\(\*\*2\)\s*=\s*(-?\d+\.?\d+)", lines_step[index_energy + 1]).group(1))
        molecule = Molecule(
            name=name,
            mult=mult,
            symbols=symbols,
            atomcoords=atomcoords,
            scfenergy=scfenergy,
        )
        molecules_backward_step.append(molecule)
    
    molecules_backward_itr, molecule_backward_optimized = _read_optopt(lines_backward[indices_optopt[0]:indices_optopt[1] + 1])
    
    #indices_freqfreq = [i for i, line in enumerate(lines_backward) if line.startswith("FREQFREQ")]
    
    return (
        molecules_forward_step,
        molecules_forward_itr,
        molecule_forward_optimized,
        molecules_backward_step,
        molecules_backward_itr,
        molecule_backward_optimized,
    )