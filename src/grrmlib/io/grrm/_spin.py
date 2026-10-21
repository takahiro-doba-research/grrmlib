import math


def spin2_to_mult(spin2: float) -> int:
    """
    GRRM が出力する Spin(**2) = <S^2> = S(S+1) から多重度 2S+1 を求める。
    2S+1 = sqrt(1 + 4 S(S+1))。スピン混入で <S^2> が少しずれていても、四捨五入して整数にする。
    """
    return round(math.sqrt(1 + 4 * spin2))
