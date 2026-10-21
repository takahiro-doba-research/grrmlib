from collections.abc import Iterable

import polars as pl

from ..core import Molecules


def to_dataframe(
    mols: Molecules,
    key_names: Iterable[str] | str,
    attrs: Iterable[str] | str
) -> pl.DataFrame:
    """
    key を key_names の列に、各分子の attrs を残りの列にした DataFrame を作る。
    分子が持っていない属性は null にする。
    """
    key_names = [key_names] if isinstance(key_names, str) else list(key_names)
    attrs = [attrs] if isinstance(attrs, str) else list(attrs)
    rows = []

    for key, mol in mols.items():
        row = list(key) if isinstance(key, tuple) else [key]
        row += [getattr(mol, attr, None) for attr in attrs]
        rows.append(row)

    return pl.DataFrame(
        rows,
        schema=key_names + attrs,
        orient="row"
    )
