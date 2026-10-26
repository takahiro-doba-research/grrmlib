from collections.abc import Iterable

import polars as pl


def select_lowest(
    df: pl.DataFrame,
    group_by: Iterable[str] | str,
    by: str = "scfenergy"
) -> pl.DataFrame:
    """
    group_by の列が同じ行のうち、by の値が最も小さい行を 1 つずつ選ぶ。
    by が null の行（計算に失敗したジョブなど）は後回しにし、
    グループの全ての行が null のときは最初の行を選ぶ。列の順番は df のまま。
    """
    group_by = [group_by] if isinstance(group_by, str) else list(group_by)
    return (
        df
        .group_by(group_by, maintain_order=True)
        .agg(pl.all().sort_by(by, nulls_last=True, maintain_order=True).first())
        .select(df.columns)
    )
