import polars as pl

from ..core import Molecules, is_identical, separate


def separate_eqs(eqs: Molecules) -> tuple[Molecules, pl.DataFrame]:
    """
    EQ をフラグメント (SEQ) に分け、is_identical で同じものを SG にまとめる。
    SEQ には EQ, fragment の順に番号を振る。

    Returns
    -------
    seqs : key が (SG, SEQ) の Molecules
    table : G, EQ, SEQ, SG の対応表。1 行が 1 つの SEQ。
        G は EQ を、SG は SEQ を is_identical でまとめたグループ。
        labels は SEQ の原子の元の EQ でのラベルを "_" でつないだもの。
    """
    # キーは (EQ, fragment)
    frags = eqs.expand(lambda eq: separate(eq)).flatten()
    sgs = frags.reset_keys().cluster(is_identical)

    eq_to_g = {
        eq: g
        for g, mols in eqs.cluster(is_identical).items()
        for eq in mols
    }
    seq_to_sg = {
        seq: sg
        for sg, mols in sgs.items()
        for seq in mols
    }
    rows = [
        {
            "G": eq_to_g[eq],
            "EQ": eq,
            "SEQ": seq,
            "SG": seq_to_sg[seq],
            "fragment": fragment,
            "labels": "_".join(str(l) for l in frag.labels),
        }
        for seq, ((eq, fragment), frag) in enumerate(frags.items())
    ]
    return sgs.flatten(), pl.DataFrame(rows)
