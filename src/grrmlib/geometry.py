import numpy as np
from scipy.spatial import distance
from scipy.spatial.transform import Rotation as R


def get_distance(atomcoords, index_a, index_b):
    return np.linalg.norm(atomcoords[index_a] - atomcoords[index_b])


def get_dihedral_angle(atomcoords, index_a, index_b, index_c, index_d, degrees=False):
    """
    Returns dihedral angle.
    a-b-c-d dihedral, right-hand rule, range (-pi, pi].
    If degrees=True, returns degrees instead of radians.
    """
    a = atomcoords[index_a]
    b = atomcoords[index_b]
    c = atomcoords[index_c]
    d = atomcoords[index_d]
    
    ab = b - a
    bc = c - b
    cd = d - c
    
    n1 = np.cross(ab, bc)
    n2 = np.cross(bc, cd)
    
    eps = 1e-12
    n1_norm = np.linalg.norm(n1)
    n2_norm = np.linalg.norm(n2)
    bc_norm = np.linalg.norm(bc)
    
    if n1_norm < eps or n2_norm < eps or bc_norm < eps:
        return np.nan
    
    n1 /= n1_norm
    n2 /= n2_norm
    bc /= bc_norm
    
    angle = np.arctan2(
        np.dot(np.cross(n1, n2), bc),
        np.dot(n1, n2)
    )
    
    if degrees:
        angle = np.degrees(angle)
    
    return angle


def overlay(
    atomcoords0,
    atomcoords1,
    index_a,
    index_b,
    index_c,
    index_d,
    index_e,
    index_f,
):
    """
    Overlays atomcoords1 onto atomcoords0 such that:
      1) d is translated to a
      2) vector d–e is aligned to a–b
      3) dihedral angle f–e–d–c is set to 0
    
    atomcoords0   atomcoords1
          c             f                 (c)f
    a – b    +    d – e     ->    a,d – (b)e
    """
    # 1) d is translated to a
    a = atomcoords0[index_a]
    d = atomcoords1[index_d]
    atomcoords1 = atomcoords1 - d + a
    
    # 2) vector d–e is aligned to a–b
    b = atomcoords0[index_b]
    d = atomcoords1[index_d]
    e = atomcoords1[index_e]
    ab = b - a
    de = e - d
    ab /= np.linalg.norm(ab)
    de /= np.linalg.norm(de)
    axis = np.cross(de, ab)
    sin_theta = np.linalg.norm(axis)
    cos_theta = np.dot(de, ab)
    
    if sin_theta < 1e-12:
        if cos_theta < 0:
            # antiparallel
            f = atomcoords1[index_f]
            ef = f - e
            ef /= np.linalg.norm(ef)
            axis_tmp = np.cross(de, ef)
            axis_tmp /= np.linalg.norm(axis_tmp)
            rot = R.from_rotvec(axis_tmp * np.pi)
        else:
            # parallel
            rot = R.identity()
    else:
        axis /= sin_theta
        rot = R.from_rotvec(axis * np.arctan2(sin_theta, cos_theta))
    
    atomcoords1 = rot.apply(atomcoords1 - d) + d
    
    # 3) dihedral angle f–e–d–c is set to 0
    c = atomcoords0[index_c]
    d = atomcoords1[index_d]
    e = atomcoords1[index_e]
    f = atomcoords1[index_f]
    fe = e - f
    ed = d - e
    dc = c - d
    n0 = np.cross(fe, ed)
    n1 = np.cross(ed, dc)
    n0 /= np.linalg.norm(n0)
    n1 /= np.linalg.norm(n1)
    ed /= np.linalg.norm(ed)
    angle = np.arctan2(
        np.dot(np.cross(n0, n1), ed),
        np.dot(n0, n1)
    )
    rot = R.from_rotvec(ed * angle)
    atomcoords1 = rot.apply(atomcoords1 - d) + d
    return atomcoords1


def rotate(atomcoords, index_a, index_b, angle, degrees=False):
    """
    Rotates atomcoords around the a – b axis by `angle` (radians).
    Right-hand rule: positive angle rotates according to a –> b direction.
    """
    if degrees:
        angle = np.deg2rad(angle)
    a = atomcoords[index_a]
    b = atomcoords[index_b]
    axis = b - a
    axis /= np.linalg.norm(axis)
    rot = R.from_rotvec(axis * angle)
    atomcoords = rot.apply(atomcoords - a) + a
    return atomcoords