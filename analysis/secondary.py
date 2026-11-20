"""Compare secondary H/T and D/T effects with the exponents in the same papers.

These mixed-label measurements use different primary backgrounds. The
comparison checks published numbers without interpreting the secondary
exponent as evidence against a kinetic mechanism.
"""
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    package = Path(__file__).resolve().parents[1]
    records = pd.read_csv(package / "data" / "secondary_adh.csv")
    records["gamma_obs"] = np.log(records.K_HT) / np.log(records.K_DT)
    records["abs_difference"] = (records.gamma_obs - records.exp_pub).abs()
    records["within_quoted_error"] = records.abs_difference <= records.exp_pub_se
    columns = ["family", "variant", "T_C", "source", "gamma_obs", "exp_pub",
               "abs_difference", "within_quoted_error"]
    report = (
        "SECONDARY COMPETITIVE ISOTOPE EFFECTS\n\n"
        + records[columns].to_string(index=False, float_format=lambda x: f"{x:.4f}")
        + f"\n\n{len(records)} comparisons with a published exponent. "
        + f"Maximum absolute difference: {records.abs_difference.max():.4f}.\n"
        + f"Within quoted source uncertainty: {records.within_quoted_error.sum()}.\n"
        + "Different primary tracer backgrounds and unresolved reverse commitments "
        + "limit interpretation of these secondary measurements.\n"
    )
    (package / "results" / "secondary_adh.txt").write_text(report)
    records.to_csv(package / "results" / "secondary_adh.csv", index=False)
    print(report)


if __name__ == "__main__":
    main()
