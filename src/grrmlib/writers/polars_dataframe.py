from typing import Iterable

import polars as pl

from ..core import Molecules


class PolarsExporter:
    
    def export(
        self,
        mols: Molecules,
        prefixes: Iterable[str] | str,
        attrs: Iterable[str] | str
    ) -> pl.DataFrame:
        prefixes = [prefixes] if isinstance(prefixes, str) else prefixes
        attrs = [attrs] if isinstance(attrs, str) else attrs
        rows = []
        
        for name, mol in mols.items():
            row = list(name)
            for attr in attrs:
                if hasattr(mol, attr):
                    row.append(getattr(mol, attr))
                else:
                    row.append(None)
            rows.append(row)
        
        return pl.DataFrame(
            rows,
            schema=prefixes + attrs,
            orient="row"
        )