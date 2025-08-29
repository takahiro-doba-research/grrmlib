import numpy as np
from scipy.spatial import distance

from .data import covalent_radius


def get_adj_matrix(symbols, atomcoords, threshold=1.25):
    arr_distance = distance.cdist(atomcoords, atomcoords)
    arr_radius = np.array([covalent_radius(s) for s in symbols])
    arr_radius = arr_radius[:, None] + arr_radius[None, :]
    arr_adj = (arr_distance < arr_radius * threshold).astype(int)
    return arr_adj