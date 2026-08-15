"""Sweep every derived number in the manuscript against the code that produces it.

Motivation.  A number written into prose while composing a sentence is
indistinguishable, to the author, from a number recalled from a computation.
Whitelist verification -- checking the numbers one thinks to check -- cannot
catch a fabricated value, because a fabricated value is exactly the one that
does not come to mind as needing a check.  This module inverts that: it
recomputes every derived quantity from the code and asserts that the manuscript
contains it, and flags any place where the manuscript instead contains a
*different* value for the same quantity (the stale-value failure).

Run from analysis/:   python audit_numbers.py
Exit status is non-zero if any derived quantity is missing or contradicted.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import masses as M
from completion import (bsao_homogeneity, bsao_masking_by_temperature,
                        yadh_completion)
from network_geometry import BENCH, bypass_to_destroy, endpoint, endpoint_closed
from partial_id import F_min_exact

GSC = M.gamma_sc("C")
F0 = M.offset_F0("C")
ROOT = Path("/home/dong/Workspace/WritePaper/MDPI_Entropy_QuantumB/v3_quantum/manuscript")
DOCS = {"main": ROOT / "pnas_manuscript.tex", "si": ROOT / "si/si_body.tex"}


def derived():
    """(label, value, decimals, which documents must contain it)."""
    out = [
        ("gamma_SC", GSC, 5, ("main", "si")),
        ("F0", F0, 6, ("main", "si")),
    ]
    for a, lab in ((1.3, "hi"), (7.3, "lo")):
        r = yadh_completion(a)
        out.append((f"YADH F_int {lab}", r["F"], 3, ("main", "si")))
        out.append((f"YADH gamma_int {lab}", r["gamma"], 2, ("main", "si")))
    T, m, me = bsao_masking_by_temperature()
    out.append(("BSAO m(25C)", m[3], 3, ("main", "si")))
    out.append(("BSAO sd m(25C)", me[3], 3, ("main", "si")))
    mu, se, Q, dof = bsao_homogeneity()
    out.append(("BSAO Q all six", Q, 1, ("main", "si")))
    out.append(("bypass phi* r=1.31 vs F0",
                bypass_to_destroy(7.13, 1.73, F0, 1.31), 3, ("si",)))
    out.append(("bypass phi* r=1.31 vs 0",
                bypass_to_destroy(7.13, 1.73, 0.0, 1.31), 3, ("si",)))
    out.append(("E_r(0.15) r=1.31", endpoint_closed(7.13, 1.73, 0.15, 1.31), 4, ("si",)))
    # sensitivity illustration, now carrying r
    from completion import commitment_precision_needed, intrinsic_from_commitment
    for m, dp in ((3.0, 3), (5.0, 3)):
        c = m * (7.13 - 1)
        xh = intrinsic_from_commitment(7.13, c)
        xd = intrinsic_from_commitment(1.73, 1.3107 * c)
        Fo = float(np.log(7.13) - GSC * np.log(1.73))
        gain = (float(np.log(xh) - GSC * np.log(xd)) - Fo) / abs(F0)
        out.append((f"YADH gain at m={m:.0f}", gain, 1, ("si",)))
        out.append((f"YADH c precision m={m:.0f}",
                    100 * commitment_precision_needed(7.13, 1.73, c, 1.3107), 1,
                    ("main", "si")))
    out.append(("E(0.35) r=1", endpoint_closed(7.13, 1.73, 0.35, 1.0), 4, ("si",)))
    for lab, KH, KD, r in BENCH:
        out.append((f"endpoint {lab}", endpoint(KH, KD, 0.0, r), 4, ("si",)))
    ya = pd.read_csv("../data/cha1989_yadh.csv")
    a = ya[ya.note.str.contains("average")].iloc[0]
    out.append(("YADH F_obs", F_min_exact(a.K_HT, a.K_DT)[0], 3, ("main", "si")))
    return out


def _texts():
    return {k: p.read_text() for k, p in DOCS.items()}


def run():
    texts = _texts()
    bad = 0
    print(f"{'quantity':32s} {'value':>11s}  documents")
    for lab, val, dp, docs in derived():
        s = f"{abs(val):.{dp}f}"
        for d in docs:
            t = texts[d]
            if s in t:
                verdict = "found"
            else:
                # a near-miss at the same precision is a STALE value, not a gap
                near = re.findall(rf"{re.escape(s[:dp - 1])}\d", t)
                verdict = f"MISSING (near: {sorted(set(near))[:3]})" if near else "MISSING"
                bad += 1
            print(f"  {lab:30s} {val:+11.6f}  {d}: {verdict}")
    print(f"\n{'FAIL' if bad else 'PASS'}: {bad} derived quantities missing or contradicted")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(run())
