"""The matched primary record, in one place.

The corpus has two parts.  Temperature-resolved series live in
data/trinomial_benchmark.csv, one analysis unit per (family, variant, step).
Records measured at a single condition live in their own files, and each such
record is its own analysis unit.  Until 2026-09-14 every script that needed the
single-condition records listed their files itself; when a record was added,
each list had to be found and extended by hand.  They now read this module.

Run from analysis/:   python corpus.py
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

DATA = Path(__file__).resolve().parent.parent / "data"

SERIES_FILE = DATA / "trinomial_benchmark.csv"

# (file, filter on the note column or None, group label used by the figure)
SINGLE_CONDITION = (
    ("ladh_adh_primary.csv", None, "ladh"),
    ("bsao_grant1989.csv", None, "bsao"),
    # Cha et al. report three determinations and their average; the average is
    # the record, the determinations enter only the robustness checks.
    ("cha1989_yadh.csv", "average", "yadh"),
    # Agrawal et al. 2004: the one temperature whose effects are stated as
    # numbers (20 C); the others are plotted only and are not digitized.
    ("agrawal2004_ectsase.csv", None, "ectsase"),
)


_SPECIES = {"Escherichia coli": "ec", "Homo sapiens": "hs"}


def display_family(family: str, variant) -> str:
    """Family label naming the species where one family spans two.

    Thymidylate synthase appears as the human enzyme (wild type) and the E. coli
    enzyme (Y209W in the series, wild type as a single record), so "TSase WT"
    would name two enzymes.  The species is read from the benchmark's system
    column, never assumed.  DHFR families already carry it (ecDHFR, hsDHFR).
    """
    if family != "TSase":
        return family
    s = series()
    sysname = s[(s.family == family) & (s.variant.astype(str) == str(variant))].system.iloc[0]
    return next(v for k, v in _SPECIES.items() if sysname.startswith(k)) + family


def series() -> pd.DataFrame:
    return pd.read_csv(SERIES_FILE)


def single_condition() -> pd.DataFrame:
    """Every single-condition record, one row per analysis unit."""
    frames = []
    for fname, note, grp in SINGLE_CONDITION:
        d = pd.read_csv(DATA / fname)
        if note is not None:
            d = d[d.note.str.contains(note)]
        frames.append(d.assign(grp=grp))
    return pd.concat(frames, ignore_index=True)


def counts() -> dict:
    s, u = series(), single_condition()
    n_series = s.groupby(["family", "variant", "step"]).ngroups
    return {"records": len(s) + len(u), "series_records": len(s),
            "series": n_series, "single": len(u), "units": n_series + len(u)}


if __name__ == "__main__":
    c = counts()
    print(f"records {c['records']} = {c['series_records']} in {c['series']} "
          f"temperature series + {c['single']} single-condition records")
    print(f"analysis units {c['units']}")
