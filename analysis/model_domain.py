"""How far the ground-channel model reaches when its parameters are bounded.

F0 is the supremum of the model's offset over ALL positive action A and gating
width w (Theorem S1).  It is approached only as w grows without bound with A
proportional to w^2, a limit in which the rate itself vanishes.  The comparison
of the census with F0 is therefore a comparison with the most permissive
version of the model: an exclusion drawn from F0 holds under any restriction of
the parameters, but a unit that cannot exclude the unrestricted model may
exclude a restricted one.  This module quantifies that dependence.

At fixed intrinsic H/T effect K the model's offset increases with w (the proof
of Theorem S1), from the rigid value at w -> 0 to its supremum, and the action
that keeps X_H = ln K increases with w for the effects met here.  With w <= w_max
and A <= A_max the attainable offsets at that K are therefore the interval
between the rigid value and the offset at the largest admissible width.  An
observed pair can be reproduced by the bounded model exactly when its allowed
curve (one intrinsic pair per commitment) passes through that band.

Physical scales: w = mu_H omega sigma^2 / hbar and A = mu_H omega R0^2 / (2 hbar)
for wells of angular frequency omega, sigma the r.m.s. fluctuation of the
donor-acceptor distance and R0 the mean separation of the two minima.

Run from analysis/:   python model_domain.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import corpus
from partial_id import F0, GSC, MUD, MUT

T_, D_ = MUT, MUD                      # square-root reduced-mass ratios t and d

# the domains reported in the supplement: (largest w, largest A)
DOMAINS = (("unrestricted", np.inf, np.inf),
           ("w <= 20, A <= 1000", 20.0, 1000.0),
           ("w <= 5, A <= 100", 5.0, 100.0),
           ("w <= 2, A <= 50", 2.0, 50.0))


def _p(w, m):
    return 0.5 * np.log((1 + w * T_) / (1 + w * m))


def _q(w, m):
    return (T_ - m) / ((1 + w * m) * (1 + w * T_))


def action(XH, w):
    """Action that gives the intrinsic log H/T effect XH at gating width w."""
    return (XH - _p(w, 1.0)) / _q(w, 1.0)


def offset(XH, w):
    """Model offset at intrinsic log H/T effect XH and gating width w."""
    return XH - GSC * (_p(w, D_) + action(XH, w) * _q(w, D_))


def rigid_offset(XH):
    """The w -> 0 limit: ln k_L linear in sqrt(mu_L)."""
    return XH * (1.0 - GSC * (T_ - D_) / (T_ - 1.0))


def _largest_width(XH, w_max, A_max, iters=80):
    """Largest w <= w_max at which the action stays positive and <= A_max."""
    XH = np.asarray(XH, dtype=float)
    hi = np.full(XH.shape, min(w_max, 1e12))
    ok = (action(XH, hi) <= A_max) & (action(XH, hi) > 0)
    lo = np.full(XH.shape, 1e-12)
    # bisect where the upper end is not admissible (action increases with w here)
    for _ in range(iters):
        mid = np.sqrt(lo * hi)
        good = (action(XH, mid) <= A_max) & (action(XH, mid) > 0)
        lo = np.where(~ok & good, mid, lo)
        hi = np.where(~ok & ~good, mid, hi)
    return np.where(ok, hi, lo)


def sup_offset(K, w_max=np.inf, A_max=np.inf):
    """Largest offset of the bounded model at intrinsic H/T effect K."""
    XH = np.log(np.asarray(K, dtype=float))
    if not np.isfinite(w_max) and not np.isfinite(A_max):
        return np.where(np.exp(XH) >= np.sqrt(T_), F0, np.nan) + 0 * XH
    return offset(XH, _largest_width(XH, w_max, A_max))


def reproducible(kht, kdt, w_max=np.inf, A_max=np.inf, n=40001):
    """Can the bounded model, at some commitment, reproduce the observed pair?"""
    if not kht > kdt > 1:
        return False
    a = kht - 1.0
    c = a * (1.0 + np.logspace(-9, 9, n))
    xH = kht * c / (c - a)
    xD = kdt * c / (c - (kdt - 1.0))
    F = np.log(xH) - GSC * np.log(xD)
    XH = np.log(xH)
    top = sup_offset(xH, w_max, A_max)
    if not np.isfinite(w_max) and not np.isfinite(A_max):
        return bool(np.any((F < top) & (F > rigid_offset(XH))))
    return bool(np.any((F <= top) & (F > rigid_offset(XH))))


def units_not_excluding(w_max=np.inf, A_max=np.inf):
    """Number of analysis units whose central values the bounded model can
    reproduce at every admissible temperature (so the unit cannot exclude it
    even without measurement error), and the labels of the others."""
    s = corpus.series()
    groups = [(f"{f} {v} ({st})", g) for (f, v, st), g in
              s.groupby(["family", "variant", "step"], sort=False)]
    groups += [(f"{r.family} {r.variant}, {r.T_C:g} C ({r.source})", pd.DataFrame([r]))
               for _, r in corpus.single_condition().iterrows()]
    n, excl = 0, []
    for label, g in groups:
        adm = g[(g.K_HT > g.K_DT) & (g.K_DT > 1)]
        if all(reproducible(kh, kd, w_max, A_max) for kh, kd in zip(adm.K_HT, adm.K_DT)):
            n += 1
        else:
            excl.append(label)
    return n, excl


def scales(w, A, wavenumber=3000.0):
    """(sigma, R0) in angstrom for a well of the given wavenumber (cm^-1)."""
    import masses as M
    hbar, amu, c_cm = 1.054571817e-34, 1.66053906660e-27, 2.99792458e10
    mu = M.reduced(M.M_H_ATOMIC, "C") * amu
    omega = 2 * np.pi * c_cm * wavenumber
    b2 = mu * omega / hbar                         # = w / sigma^2, in m^-2
    return float(np.sqrt(w / b2) * 1e10), float(np.sqrt(2 * A / b2) * 1e10)


def width_for_gap(K, gap):
    """Natural logarithm of the gating width at which the model comes within
    `gap` of F0 at effect K (exponentiate for w)."""
    from scipy.optimize import brentq
    XH = float(np.log(K))
    return float(brentq(lambda lw: F0 - offset(XH, np.exp(lw)) - gap, np.log(1e-3), np.log(1e9)))


def table():
    rows = []
    for name, w_max, A_max in DOMAINS:
        n, _ = units_not_excluding(w_max, A_max)
        sig, R0 = scales(w_max, A_max) if np.isfinite(w_max) else (np.inf, np.inf)
        rows.append(dict(domain=name, w_max=w_max, A_max=A_max, sigma_A=sig, R0_A=R0,
                         sup_at_7=float(sup_offset(7.0, w_max, A_max)), not_excluding=n))
    return pd.DataFrame(rows)


def main():
    t = table()
    for r in t.itertuples():
        print(f"{r.domain:20s} sigma <= {r.sigma_A:5.2f} A  R0 <= {r.R0_A:5.2f} A   "
              f"sup F at K = 7: {r.sup_at_7:+.4f}   units that cannot exclude: {r.not_excluding}")
    w = np.exp(width_for_gap(7.0, 0.01))
    print(f"within 0.01 of F0 at K = 7 needs w = {w:.0f} (sigma = {scales(w, 1.0)[0]:.2f} A)")
    t.to_csv("../results/model_domain.csv", index=False, float_format="%.6f")
    print("wrote ../results/model_domain.csv")


if __name__ == "__main__":
    main()
