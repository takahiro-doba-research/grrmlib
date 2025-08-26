def covalent_radii(atom):
    """
    Chem. Eur. J. 2009, 15, 186-197.
    """
    if atom == 'H':
        return 0.384  #0.32 * 1.2
    elif atom == 'C':
        return 0.75
    elif atom == 'N':
        return 0.71
    elif atom == 'O':
        return 0.63
    elif atom == 'He':
        return 0.46 * 1.2
    elif atom == 'Li':
        return 1.33 * 0.8  #1.064
    elif atom == 'Be':
        return 1.02 * 0.9  #0.918
    elif atom == 'B':
        return 0.85
    elif atom == 'F':
        return 0.64
    elif atom == 'Ne':
        return 0.67
    elif atom == 'Na':
        return 1.55 * 0.8  #1.240
    elif atom == 'Mg':
        return 1.39 * 0.9  #1.251
    elif atom == 'Al':
        return 1.26
    elif atom == 'Si':
        return 1.16
    elif atom == 'P':
        return 1.11
    elif atom == 'S':
        return 1.03
    elif atom == 'Cl':
        return 0.99
    elif atom == 'Ar':
        return 0.96
    elif atom == 'K':
        return 1.96 * 0.8  #1.568
    elif atom == 'Ca':
        return 1.71 * 0.9  #1.539
    elif atom == 'Sc':
        return 1.48
    elif atom == 'Ti':
        return 1.36
    elif atom == 'V':
        return 1.34
    elif atom == 'Cr':
        return 1.22
    elif atom == 'Mn':
        return 1.19
    elif atom == 'Fe':
        return 1.16
    elif atom == 'Co':
        return 1.11
    elif atom == 'Ni':
        return 1.10
    elif atom == 'Cu':
        return 1.12
    elif atom == 'Zn':
        return 1.18
    elif atom == 'Ga':
        return 1.24
    elif atom == 'Ge':
        return 1.21
    elif atom == 'As':
        return 1.21
    elif atom == 'Se':
        return 1.16
    elif atom == 'Br':
        return 1.14
    elif atom == 'Kr':
        return 1.17
    elif atom == 'Rb':
        return 2.10 * 0.8  #1.680
    elif atom == 'Sr':
        return 1.85 * 0.9  #1.665
    elif atom == 'Y':
        return 1.63
    elif atom == 'Zr':
        return 1.54
    elif atom == 'Nb':
        return 1.47
    elif atom == 'Mo':
        return 1.38
    elif atom == 'Tc':
        return 1.28
    elif atom == 'Ru':
        return 1.25
    elif atom == 'Rh':
        return 1.25
    elif atom == 'Pd':
        return 1.20
    elif atom == 'Ag':
        return 1.28
    elif atom == 'Cd':
        return 1.36
    elif atom == 'In':
        return 1.42
    elif atom == 'Sn':
        return 1.40
    elif atom == 'Sb':
        return 1.40
    elif atom == 'Te':
        return 1.36
    elif atom == 'I':
        return 1.33
    elif atom == 'Xe':
        return 1.31
    elif atom == 'Cs':
        return 2.32 * 0.8  #1.856
    elif atom == 'Ba':
        return 1.96 * 0.9  #1.764
    elif atom == 'La':
        return 1.80
    elif atom == 'Ce':
        return 1.63
    elif atom == 'Pr':
        return 1.76
    elif atom == 'Nd':
        return 1.74
    elif atom == 'Pm':
        return 1.73
    elif atom == 'Sm':
        return 1.72
    elif atom == 'Eu':
        return 1.68
    elif atom == 'Gd':
        return 1.69
    elif atom == 'Tb':
        return 1.68
    elif atom == 'Dy':
        return 1.67
    elif atom == 'Ho':
        return 1.66
    elif atom == 'Er':
        return 1.65
    elif atom == 'Tm':
        return 1.64
    elif atom == 'Yb':
        return 1.70
    elif atom == 'Lu':
        return 1.62
    elif atom == 'Hf':
        return 1.52
    elif atom == 'Ta':
        return 1.46
    elif atom == 'W':
        return 1.37
    elif atom == 'Re':
        return 1.31
    elif atom == 'Os':
        return 1.29
    elif atom == 'Ir':
        return 1.22
    elif atom == 'Pt':
        return 1.23
    elif atom == 'Au':
        return 1.24
    elif atom == 'Hg':
        return 1.33
    elif atom == 'Tl':
        return 1.44
    elif atom == 'Pb':
        return 1.44
    elif atom == 'Bi':
        return 1.51
    elif atom == 'Po':
        return 1.45
    elif atom == 'At':
        return 1.47
    elif atom == 'Rn':
        return 1.42
    elif atom == 'Fr':
        return 2.23 * 0.8  #1.784
    elif atom == 'Ra':
        return 2.01 * 0.9  #1.809
    elif atom == 'Ac':
        return 1.86
    elif atom == 'Th':
        return 1.75
    elif atom == 'Pa':
        return 1.69
    elif atom == 'U':
        return 1.70
    elif atom == 'Np':
        return 1.71
    elif atom == 'Pu':
        return 1.72
    elif atom == 'Am':
        return 1.66
    elif atom == 'Cm':
        return 1.66
    elif atom == 'Bk':
        return 1.68
    elif atom == 'Cf':
        return 1.68
    elif atom == 'Es':
        return 1.65
    elif atom == 'Fm':
        return 1.67
    elif atom == 'Md':
        return 1.73
    elif atom == 'No':
        return 1.76
    elif atom == 'Lr':
        return 1.61
    elif atom == 'Rf':
        return 1.57
    elif atom == 'Db':
        return 1.49
    elif atom == 'Sg':
        return 1.43
    elif atom == 'Bh':
        return 1.41
    elif atom == 'Hs':
        return 1.34
    elif atom == 'Mt':
        return 1.29
    elif atom == 'Ds':
        return 1.28
    elif atom == 'Rg':
        return 1.21
    elif atom == 'Cn':
        return 1.22
    elif atom == 'Nh':
        return 1.36
    elif atom == 'Fl':
        return 1.43
    elif atom == 'Mc':
        return 1.62
    elif atom == 'Lv':
        return 1.75
    elif atom == 'Ts':
        return 1.65
    elif atom == 'Og':
        return 1.57
    elif atom == 'TV':
        return 1.50
    else:
        return None


def atomic_number(atom):
    if atom == 'H':
        return 1
    elif atom == 'C':
        return 6
    elif atom == 'N':
        return 7
    elif atom == 'O':
        return 8
    elif atom == 'He':
        return 2
    elif atom == 'Li':
        return 3
    elif atom == 'Be':
        return 4
    elif atom == 'B':
        return 5
    elif atom == 'F':
        return 9
    elif atom == 'Ne':
        return 10
    elif atom == 'Na':
        return 11
    elif atom == 'Mg':
        return 12
    elif atom == 'Al':
        return 13
    elif atom == 'Si':
        return 14
    elif atom == 'P':
        return 15
    elif atom == 'S':
        return 16
    elif atom == 'Cl':
        return 17
    elif atom == 'Ar':
        return 18
    elif atom == 'K':
        return 19
    elif atom == 'Ca':
        return 20
    elif atom == 'Sc':
        return 21
    elif atom == 'Ti':
        return 22
    elif atom == 'V':
        return 23
    elif atom == 'Cr':
        return 24
    elif atom == 'Mn':
        return 25
    elif atom == 'Fe':
        return 26
    elif atom == 'Co':
        return 27
    elif atom == 'Ni':
        return 28
    elif atom == 'Cu':
        return 29
    elif atom == 'Zn':
        return 30
    elif atom == 'Ga':
        return 31
    elif atom == 'Ge':
        return 32
    elif atom == 'As':
        return 33
    elif atom == 'Se':
        return 34
    elif atom == 'Br':
        return 35
    elif atom == 'Kr':
        return 36
    elif atom == 'Rb':
        return 37
    elif atom == 'Sr':
        return 38
    elif atom == 'Y':
        return 39
    elif atom == 'Zr':
        return 40
    elif atom == 'Nb':
        return 41
    elif atom == 'Mo':
        return 42
    elif atom == 'Tc':
        return 43
    elif atom == 'Ru':
        return 44
    elif atom == 'Rh':
        return 45
    elif atom == 'Pd':
        return 46
    elif atom == 'Ag':
        return 47
    elif atom == 'Cd':
        return 48
    elif atom == 'In':
        return 49
    elif atom == 'Sn':
        return 50
    elif atom == 'Sb':
        return 51
    elif atom == 'Te':
        return 52
    elif atom == 'I':
        return 53
    elif atom == 'Xe':
        return 54
    elif atom == 'Cs':
        return 55
    elif atom == 'Ba':
        return 56
    elif atom == 'La':
        return 57
    elif atom == 'Ce':
        return 58
    elif atom == 'Pr':
        return 59
    elif atom == 'Nd':
        return 60
    elif atom == 'Pm':
        return 61
    elif atom == 'Sm':
        return 62
    elif atom == 'Eu':
        return 63
    elif atom == 'Gd':
        return 64
    elif atom == 'Tb':
        return 65
    elif atom == 'Dy':
        return 66
    elif atom == 'Ho':
        return 67
    elif atom == 'Er':
        return 68
    elif atom == 'Tm':
        return 69
    elif atom == 'Yb':
        return 70
    elif atom == 'Lu':
        return 71
    elif atom == 'Hf':
        return 72
    elif atom == 'Ta':
        return 73
    elif atom == 'W':
        return 74
    elif atom == 'Re':
        return 75
    elif atom == 'Os':
        return 76
    elif atom == 'Ir':
        return 77
    elif atom == 'Pt':
        return 78
    elif atom == 'Au':
        return 79
    elif atom == 'Hg':
        return 80
    elif atom == 'Tl':
        return 81
    elif atom == 'Pb':
        return 82
    elif atom == 'Bi':
        return 83
    elif atom == 'Po':
        return 84
    elif atom == 'At':
        return 85
    elif atom == 'Rn':
        return 86
    elif atom == 'Fr':
        return 87
    elif atom == 'Ra':
        return 88
    elif atom == 'Ac':
        return 89
    elif atom == 'Th':
        return 90
    elif atom == 'Pa':
        return 91
    elif atom == 'U':
        return 92
    elif atom == 'Np':
        return 93
    elif atom == 'Pu':
        return 94
    elif atom == 'Am':
        return 95
    elif atom == 'Cm':
        return 96
    elif atom == 'Bk':
        return 97
    elif atom == 'Cf':
        return 98
    elif atom == 'Es':
        return 99
    elif atom == 'Fm':
        return 100
    elif atom == 'Md':
        return 101
    elif atom == 'No':
        return 102
    elif atom == 'Lr':
        return 103
    elif atom == 'Rf':
        return 104
    elif atom == 'Db':
        return 105
    elif atom == 'Sg':
        return 106
    elif atom == 'Bh':
        return 107
    elif atom == 'Hs':
        return 108
    elif atom == 'Mt':
        return 109
    elif atom == 'Ds':
        return 110
    elif atom == 'Rg':
        return 111
    elif atom == 'Cn':
        return 112
    elif atom == 'Nh':
        return 113
    elif atom == 'Fl':
        return 114
    elif atom == 'Mc':
        return 115
    elif atom == 'Lv':
        return 116
    elif atom == 'Ts':
        return 117
    elif atom == 'Og':
        return 118
    else:
        return None
        