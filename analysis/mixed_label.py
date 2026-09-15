"""Forward masking in the mixed-labeling secondary experiment is ordered.

In the mixed-labeling design (Kohen & Jensen, J. Am. Chem. Soc. 124, 3858, 2002) the
secondary H/T effect is measured for C-H cleavage and the secondary D/T effect
for C-D cleavage.  Index a molecule by the isotope at the transferred position,
then at the secondary one.  With isotope-independent dissociation each secondary
ratio is the series map

    K_i = x_i (1 + c_i) / (x_i + c_i),    c_i = k_off / k_iT,

whose reference molecules are HT and DT.  So c_D / c_H = k_HT / k_DT, the primary
H/D effect on the tritiated reference: above one for a normal primary effect.
That is the ordering of Proposition S4, and its proof uses no property of the
reference exponent beyond exceeding one, so forward masking can only lower an
observed secondary exponent below its intrinsic value.

Until 2026-09-15 the supplement presented c_H = 5, c_D = 1 (observed exponent
7.85 from an intrinsic 4.8) as what masking does in this design.  That ratio
requires an inverse primary effect of five on the reference.  Kohen and Jensen's
inflation comes instead from a reverse commitment acting with a secondary
equilibrium isotope effect (their eqs 13-14), outside forward masking.

Run from analysis/:   python mixed_label.py
"""
from __future__ import annotations

import numpy as np


def observed(x, c):
    """The series map for one secondary ratio."""
    return x * (1.0 + c) / (x + c)


def observed_exponent(x_D, gamma_int, c_H, c_D):
    """Observed secondary exponent for an intrinsic pair on the gamma_int locus."""
    x_H = x_D ** gamma_int
    return float(np.log(observed(x_H, c_H)) / np.log(observed(x_D, c_D)))


def check_ordered(n=400000, seed=20260915, g_lo=1.05, g_hi=8.0):
    """F_obs <= F_int for ordered maps (c_D >= c_H), at arbitrary reference exponents.

    Draws the reference exponent g, an intrinsic pair with F_int >= 0 against g,
    c_H over twenty-four decades and r = c_D/c_H >= 1.  Returns (draws,
    violations, largest excess of F_obs over F_int).  For contrast also returns
    how often r < 1 lifts a pair lying ON the g-ray above it.
    """
    rng = np.random.default_rng(seed)
    g = rng.uniform(g_lo, g_hi, n)
    x_D = 1.0 + rng.uniform(1e-3, 3.0, n)
    x_H = x_D ** rng.uniform(g, g + 3.0)
    c_H = np.exp(rng.uniform(-8.0, 16.0, n))
    r = np.exp(rng.uniform(0.0, 6.0, n))
    F_int = np.log(x_H) - g * np.log(x_D)
    F_obs = np.log(observed(x_H, c_H)) - g * np.log(observed(x_D, r * c_H))
    excess = F_obs - F_int
    bad = int(np.sum(excess > 1e-12))
    # anti-ordered contrast: pairs exactly on the ray, r in (0.05, 1)
    r_bad = rng.uniform(0.05, 1.0, n)
    x_Hr = x_D ** g
    F_anti = np.log(observed(x_Hr, c_H)) - g * np.log(observed(x_D, r_bad * c_H))
    lifted = int(np.sum(F_anti > 1e-12))
    return n, bad, float(excess.max()), lifted


def main():
    print("mixed-labeling secondary experiment, intrinsic pair on the 4.8 locus, x_D = 1.10")
    for c_H, c_D, tag in ((1.0, 5.0, "ordered (normal primary effect)"),
                          (5.0, 5.0, "shared"),
                          (5.0, 1.0, "anti-ordered (inverse primary effect of 5)")):
        print(f"  c_H={c_H:g}  c_D={c_D:g}  exponent {observed_exponent(1.10, 4.8, c_H, c_D):.4f}   {tag}")
    n, bad, worst, lifted = check_ordered()
    print(f"\nordered maps, reference exponent 1.05-8: {bad} of {n} with F_obs > F_int "
          f"(largest excess {worst:.1e}); anti-ordered maps lift a ray pair in {lifted} of {n}")
    return bad


if __name__ == "__main__":
    raise SystemExit(1 if main() else 0)
