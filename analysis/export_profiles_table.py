"""Profile table for the SI (Table `tab:profiles`).

Regenerates manuscript/sm_table_profiles.tex from results/offset_profiles.csv.
Previously hand-maintained, which let it fall out of step with the benchmark
when new series were added.
"""
from __future__ import annotations

import pathlib

import numpy as np
import pandas as pd
from scipy import stats

import masses as M

F0 = abs(M.offset_F0("C"))
# The manuscript tree is not part of the public package; fall back to
# results/ so the table is still reproducible from the shipped code.
_MS = pathlib.Path("../manuscript")
OUT = str(_MS / "sm_table_profiles.tex") if _MS.is_dir() else "../results/sm_table_profiles.tex"


def main() -> None:
    d = pd.read_csv("../results/offset_profiles.csv")
    d["dof"] = 2 * d.n_T - 5
    d["p_fit"] = stats.chi2.sf(d.chi2_min, d.dof.clip(lower=1))
    # adequate fits first, then by interval width: the ordering the text uses
    d["adequate"] = d.p_fit > 0.05
    d = d.sort_values(["adequate", "width"], ascending=[False, True])
    lines = []
    for _, r in d.iterrows():
        lines.append(
            f"{r.family} & {r.variant} & {r.step} & {int(r.n_T)} & "
            f"{r.chi2_min:.1f}/{int(r.dof)} & {r.p_fit:.3f} & "
            f"{r.F_hat:+.2f} & [{r.F_lo:+.2f}, {r.F_hi:+.2f}] & "
            f"{r.width:.3f} ({r.width/F0:.1f})\\\\"
        )
    open(OUT, "w").write("\n".join(lines) + "\n")
    print(f"wrote {OUT}: {len(d)} series, "
          f"{int(d.adequate.sum())} adequate fits, {int((~d.adequate).sum())} rejected")


if __name__ == "__main__":
    main()
