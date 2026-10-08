"""
Shared reader for the all_hits.tsv / merged_all_hits.tsv tables.

Those tables list every panel gene x sample, with a `hit` column (TRUE/FALSE). Rows with
hit == FALSE carry '.' in every hit-related column, which would change how pandas types those
columns (e.g. a True/False column becomes text, and "False" is truthy). Scripts that only care
about hits should read the table through read_hit_rows(): it drops the non-hit rows first and
then parses what is left exactly as pd.read_csv(path, **kwargs) would have parsed the old
hits-only table, so downstream code sees the same values and dtypes as before.
"""
import io

import pandas as pd

HIT_COL = "hit"


def is_hit(values):
    """Boolean mask for a `hit` column read as text (TRUE/FALSE, any case)."""
    return values.astype(str).str.strip().str.upper() == "TRUE"


def read_hit_rows(path, **read_csv_kwargs):
    raw = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if HIT_COL in raw.columns:
        raw = raw[is_hit(raw[HIT_COL])]
    buf = io.StringIO(raw.to_csv(sep="\t", index=False))
    read_csv_kwargs.setdefault("sep", "\t")
    return pd.read_csv(buf, **read_csv_kwargs)
