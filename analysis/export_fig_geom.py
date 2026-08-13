"""Data for the geometry figure: why the offset, and why the set is one-sided.

Panel A is the mass-scaling plane.  A point is (X_D, X_H) = (ln K_DT, ln K_HT);
the semiclassical locus is the ray X_H = gamma_SC X_D, and the offset F is the
VERTICAL distance from it.  As the commitment falls, the inferred intrinsic pair
traces a curve, and that curve only ever moves upward: dF/dc < 0 in the regime
that matters.  So an observation above the locus can never be brought onto it by
any commitment, while one below it can be.  That is Theorem 1, drawn.

Panel B is the ridge.  The exponent is the SLOPE of the ray from the origin and
converges to gamma_SC; the offset is the vertical displacement and does not move.
"""
from __future__ import annotations

import numpy as np

import masses as M

G = M.gamma_sc("C"); F0 = M.offset_F0("C")
PH, PD = M.ridge_P("C"); CH, CD = M.ridge_C("C")
OUT = "../figures/tikz/data/"


def traj(kht, kdt, cmax=4000.0, n=400):
    """(X_D, X_H) as the shared commitment runs from just above K_HT-1 upward."""
    a = kht - 1.0
    c = a * (1.0 + np.logspace(-2.1, 3.4, n))
    xh = kht * c / (1.0 + c - kht)
    xd = kdt * c / (1.0 + c - kdt)
    keep = (xh > 0) & (xd > 0) & (np.log(xh) < 3.6)
    return np.log(xd[keep]), np.log(xh[keep])


CASES = [("yadh", 7.13, 1.73), ("light", 3.10, 1.51)]
for name, kh, kd in CASES:
    XD, XH = traj(kh, kd)
    with open(OUT + f"fg_{name}.dat", "w") as f:
        f.write("xd xh\n")
        for u, v in zip(XD, XH):
            f.write(f"{u:.6f} {v:.6f}\n")
    with open(OUT + f"fg_{name}_obs.dat", "w") as f:
        f.write("xd xh\n")
        f.write(f"{np.log(kd):.6f} {np.log(kh):.6f}\n")
    Fo = np.log(kh) - G * np.log(kd)
    print(f"  {name:6s} observed ({np.log(kd):.3f},{np.log(kh):.3f})  "
          f"F_obs {Fo:+.4f}  crosses locus: {'no' if Fo > 0 else 'yes'}")

# where the negative-offset trajectory crosses the locus (F(c) = 0)
from scipy.optimize import brentq
kh, kd = 3.10, 1.51
F = lambda c: np.log(kh*c/(1+c-kh)) - G*np.log(kd*c/(1+c-kd))
a = kh - 1.0
c0 = brentq(F, a*1.0001, a*1e4)
xh0, xd0 = kh*c0/(1+c0-kh), kd*c0/(1+c0-kd)
with open(OUT + "fg_cross.dat", "w") as f:
    f.write("xd xh\n")
    f.write(f"{np.log(xd0):.6f} {np.log(xh0):.6f}\n")
print(f"  ecDHFR trajectory crosses the locus at c = {c0:.2f}, "
      f"(X_D,X_H) = ({np.log(xd0):.3f},{np.log(xh0):.3f})")

# the semiclassical locus
with open(OUT + "fg_locus.dat", "w") as f:
    f.write("xd xh\n")
    for x in np.linspace(0.25, 0.80, 60):
        f.write(f"{x:.6f} {G*x:.6f}\n")

# panel B: the ridge, exponent and offset against the ray parameter
with open(OUT + "fg_ridge.dat", "w") as f:
    f.write("s expo off\n")
    for s in np.logspace(-2, 4.2, 300):
        XH, XD = PH + s * CH, PD + s * CD
        f.write(f"{s:.6g} {XH/XD:.6f} {XH-G*XD:.8f}\n")
print(f"  ridge: exponent {PH/PD:.4f} -> {G:.4f}, offset constant at {F0:+.6f}")
