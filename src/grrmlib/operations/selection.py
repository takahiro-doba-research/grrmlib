import itertools
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


def select_charges(
    candidates: pl.DataFrame,
    group_by: Iterable[str] | str,
    *,
    part: str = "fragment",
    charge: str = "charge",
    by: str = "scfenergy",
    total_charge: int = 0
) -> pl.DataFrame:
    """
    group_by の列が同じ行（例：1 つの G）を 1 つの系とし、part の列（例：分子）ごとに 1 行ずつ選ぶ。
    part ごとの電荷の候補の組み合わせのうち、電荷の合計が total_charge になり、by の合計が最も小さいものを選ぶ。

    - 戻り値は選んだ行（系ごとに part の数だけ）で、列 margin に次の候補との by の合計の差を加える
      （候補が 1 通りの系では null）
    - 電荷の合計が total_charge になる組み合わせがない系があると ValueError
    """
    group_by = [group_by] if isinstance(group_by, str) else list(group_by)
    rows = []

    for key, frags in candidates.sort(*group_by, part, charge).group_by(group_by, maintain_order=True):
        choices = [f.to_dicts() for _, f in frags.group_by(part, maintain_order=True)]
        combos = [c for c in itertools.product(*choices) if sum(r[charge] for r in c) == total_charge]
        if not combos:
            raise ValueError(f"{dict(zip(group_by, key))}: no combination of charges sums to {total_charge}")
        combos.sort(key=lambda c: sum(r[by] for r in c))
        margin = sum(r[by] for r in combos[1]) - sum(r[by] for r in combos[0]) if len(combos) > 1 else None
        rows.extend({**r, "margin": margin} for r in combos[0])

    return pl.DataFrame(rows, schema={**candidates.schema, "margin": pl.Float64})
