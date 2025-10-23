import numpy as np
from scipy.spatial import distance

from .data import covalent_radius


def get_adj_matrix(symbols, atomcoords, threshold=1.25):
    arr_distance = distance.cdist(atomcoords, atomcoords)
    arr_radius = np.array([covalent_radius(s) for s in symbols])
    arr_radius = arr_radius[:, None] + arr_radius[None, :]
    arr_adj = (arr_distance < arr_radius * threshold).astype(int)
    return arr_adj


def get_distance(atomcoords, label0, label1):
    d = np.sqrt(np.sum((atomcoords[label0 - 1] - atomcoords[label1 - 1]) ** 2))
    return d


def get_dihedral_angle(atomcoords, label0, label1, label2, label3):
    """
    Returns dihedral angle in °
    Assumes A-B-C-D connection and planes of ABC and BCD.
    """
    A = atomcoords[label0 - 1]
    B = atomcoords[label1 - 1]
    C = atomcoords[label2 - 1]
    D = atomcoords[label3 - 1]
    
    n0 = np.cross(B - A, C - A)
    n1 = np.cross(C - B, D - B)
    
    n0 /= np.linalg.norm(n0)
    n1 /= np.linalg.norm(n1)
    cos_theta = np.dot(n0, n1)
    cos_theta = np.clip(cos_theta, -1.0, 1.0)  # 数値誤差対策
    theta = np.arccos(cos_theta)
    
    sign = np.sign(np.dot(np.cross(n0, n1), C - B))
    theta *= sign
    angle = np.degrees(theta)
    return angle