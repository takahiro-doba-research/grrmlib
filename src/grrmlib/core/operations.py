import networkx as nx
import numpy as np

from .atomcoords import (
    _get_distance,
    _get_dihedral_angle,
    _overlay,
    _rotate,
    _get_distance_matrix,
)
from .data import covalent_radius
from .molecule import Molecule
from .molecules import Molecules


def get_distance(
    mol: Molecule,
    label0: int,
    label1: int
) -> float:
    return _get_distance(
        mol.atomcoords,
        *mol.labels_to_indices([label0, label1])
    )


def get_dihedral_angle(
    mol: Molecule,
    label0: int,
    label1: int,
    label2: int,
    label3: int,
    degrees=False
) -> float:
    return _get_dihedral_angle(
        mol.atomcoords,
        *mol.labels_to_indices([label0, label1, label2, label3]),
        degrees
    )


def overlay(
    mol0: Molecule,
    mol1: Molecule,
    label0: int,
    label1: int,
    label2: int,
    label3: int,
    label4: int,
    label5: int
) -> np.ndarray:
    return _overlay(
        mol0.atomcoords,
        mol1.atomcoords,
        *mol0.labels_to_indices([label0, label1, label2]),
        *mol1.labels_to_indices([label3, label4, label5])
    )


def rotate(
    mol: Molecule,
    label0: int,
    label1: int,
    angle: float,
    degrees=False
) -> np.ndarray:
    return _rotate(
        mol.atomcoords,
        *mol.labels_to_indices([label0, label1]),
        angle,
        degrees
    )


def get_distance_matrix(mol: Molecule) -> np.ndarray:
    return _get_distance_matrix(mol.atomcoords)


def get_adj_matrix(
    mol: Molecule,
    threshold: float = 1.25
) -> np.ndarray:
    D = get_distance_matrix(mol)
    r = np.array([covalent_radius(s) for s in mol.symbols])
    R = r[:, None] + r[None, :]
    A = D < R * threshold
    np.fill_diagonal(A, False)
    return A


def separate(
    mol: Molecule,
    threshold: float = 1.25,
    old: bool = False
) -> Molecules:
    A = get_adj_matrix(mol, threshold)
    G = nx.from_numpy_array(A)
    
    if old:
        components = sorted(nx.connected_components(G), key=len, reverse=True)
    else:
        components = sorted(nx.connected_components(G), key=min)
    
    mols = Molecules()
    
    for i, component in enumerate(components):
        indices = sorted(component)
        fragment = mol._select_by_indices(indices)
        mols[i] = fragment
    
    return mols