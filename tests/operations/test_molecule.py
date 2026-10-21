import grrmlib as gl


def test_product_charge_mult(water):
    mols = gl.product_charge_mult(water, [-1, 0, 1], [1, 2])

    # 電子数と多重度の偶奇が合う組み合わせだけを残す
    assert set(mols) == {(-1, 2), (0, 1), (1, 2)}
    assert (mols[(1, 2)].charge, mols[(1, 2)].mult) == (1, 2)


def test_product_charge_mult_scalar(water):
    assert list(gl.product_charge_mult(water, 0, 1)) == [(0, 1)]
