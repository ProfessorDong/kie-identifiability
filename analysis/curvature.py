"""Curvature and the one-sidedness of the identified set.

Masking acts on log-effects as X_L -> h(X_L) with the SAME h for both isotope
pairs.  Mass scaling is homogeneity (a ray through the origin); masking is
curvature.  A concave h fixing the origin is subhomogeneous, so it drags every
mass-scaled pair below the ray.  This script verifies:

  1. the scheme of Eq.(1) satisfies the four axioms;
  2. the conclusion holds for arbitrary shared concave-through-origin maps;
  3. curvature fixes the direction (parallel/convex reverses it);
  4. the conclusion survives composition of series bottlenecks;
  5. it fails when h is not shared, and when any single axiom is dropped.
"""
from __future__ import annotations

import numpy as np

import masses as M

GSC = M.gamma_sc("C")
SEED = 20260813


def h_series(c):
    """Bottleneck in SERIES: 1/(V/K) = 1/k_on + k_off/(k_on k_i).  Concave."""
    return lambda t: np.log(np.exp(t) * (1.0 + c) / (np.exp(t) + c))


def h_parallel(c):
    """Isotope-blind route in PARALLEL: conductances add.  Convex."""
    return lambda t: np.log((np.exp(t) + c) / (1.0 + c))


def random_concave(rng, k=6):
    """h(t) = sum a_j (1 - exp(-b_j t)); concave, h(0)=0, normalised to h'(0)=1."""
    a = rng.uniform(0.1, 1.5, k)
    b = rng.uniform(0.1, 3.0, k)
    a = a / np.sum(a * b)
    return lambda t: np.sum(
        a[:, None] * (1.0 - np.exp(-b[:, None] * np.atleast_1d(t))), axis=0
    )


def offset(h, tH, tD, hD=None):
    """F_obs for intrinsic log-effects (tH, tD); hD lets the maps differ."""
    return np.asarray(h(tH)).ravel()[0] - GSC * np.asarray((hD or h)(tD)).ravel()[0]


def _draw(rng):
    """An intrinsic pair with F_int >= 0."""
    tD = rng.uniform(0.05, 1.2)
    Fi = rng.uniform(0.0, 0.8)
    return GSC * tD + Fi, tD, Fi


def check_axioms(tol=1e-9):
    """The scheme of Eq.(1): h(0)=0, 0 <= h' <= 1, h'' < 0."""
    out = []
    for c in (1e-2, 1e-1, 1.0, 1e1, 1e3):
        h = h_series(c)
        t = np.linspace(1e-6, 3.0, 4001)
        d1 = np.gradient(h(t), t)
        d2 = np.gradient(d1, t)
        out.append(
            dict(c=c, h0=float(h(0.0)), dmin=d1.min(), dmax=d1.max(), d2max=d2[5:-5].max())
        )
    return out


def run(n=20000):
    rng = np.random.default_rng(SEED)
    say = print
    say("Axioms for the scheme of Eq.(1):  h(0)=0, 0<h'<1, h''<0")
    for r in check_axioms():
        say(f"   c={r['c']:8.0e}  h(0)={r['h0']:+.1e}  h' in "
            f"[{r['dmin']:.4f},{r['dmax']:.4f}]  max h''={r['d2max']:+.2e}")

    say("\nClaim: shared h with (A1)-(A4) and F_int>=0  =>  F_obs <= F_int")
    trials = [
        ("series masking, the scheme of Eq.(1)",
         lambda: h_series(10 ** rng.uniform(-2, 3)), None),
        ("arbitrary shared concave-through-origin h",
         lambda: random_concave(rng), None),
        ("chain of 2-5 composed series bottlenecks",
         lambda: _compose([h_series(10 ** rng.uniform(-2, 3))
                           for _ in range(rng.integers(2, 6))]), None),
    ]
    for name, mk, _ in trials:
        bad = 0
        for _ in range(n):
            tH, tD, Fi = _draw(rng)
            if offset(mk(), tH, tD) > Fi + 1e-9:
                bad += 1
        say(f"   {name:44s} violations {bad}/{n}")

    say("\nNegative controls: drop one axiom, or drop sharing")
    neg = [
        ("(A4) convex instead of concave (parallel route)",
         lambda: (h_parallel(10 ** rng.uniform(-2, 3)), None)),
        ("(A3) slope > 1 (masking that amplifies)",
         lambda: ((lambda a: (lambda t: a * t))(rng.uniform(1.01, 2.0)), None)),
        ("sharing: h_H != h_D, H masked less than D",
         lambda: (h_series(10 ** rng.uniform(1, 3)),
                  h_series(10 ** rng.uniform(-2, 1)))),
    ]
    for name, mk in neg:
        bad = 0
        for _ in range(n):
            tH, tD, Fi = _draw(rng)
            hH, hD = mk()
            if offset(hH, tH, tD, hD) > Fi + 1e-9:
                bad += 1
        say(f"   {name:44s} violations {bad}/{n}")

    say("\nUnequal tritium references (the primary protocol)")
    for label, (bad, tot, worst) in unequal_reference().items():
        say(f"   {label:34s} violations {bad}/{tot}, max excess {worst:+.2e}")

    say("\nDirection on the ray (F_int = 0): sign of F_obs follows -sign(h'')")
    for label, mk in (("concave (series)", h_series), ("convex (parallel)", h_parallel)):
        up = 0
        for _ in range(n):
            tD = rng.uniform(0.05, 1.2)
            if offset(mk(10 ** rng.uniform(-2, 3)), GSC * tD, tD) > 0:
                up += 1
        say(f"   {label:20s} F_obs > 0 in {up:6d}/{n}")



def unequal_reference(n=200000):
    """The primary protocol shares h only approximately: the H/T and D/T tritium
    references differ at the non-transferred position, so c_D = r*c_H with r the
    secondary H/D effect on the reference.  Normal secondary effects give r>1,
    which masks H more and pushes F_obs down.  An inverse effect (r<1) would
    reverse this and break the argument."""
    rng = np.random.default_rng(SEED + 1)
    K = lambda x, c: x * (1 + c) / (x + c)
    out = {}
    for lo, hi, label in ((1.0, 2.0, "r >= 1 (normal, as observed)"),
                          (0.2, 1.0, "r < 1 (inverse, hypothetical)")):
        bad = 0
        worst = -np.inf
        for _ in range(n):
            tD = rng.uniform(0.05, 1.5)
            Fi = rng.uniform(0.0, 0.8) if lo >= 1.0 else 0.0
            cH = 10 ** rng.uniform(-3, 4)
            r = rng.uniform(lo, hi)
            F = np.log(K(np.exp(GSC * tD + Fi), cH)) - GSC * np.log(K(np.exp(tD), r * cH))
            if F > Fi + 1e-9:
                bad += 1
            worst = max(worst, F - Fi)
        out[label] = (bad, n, worst)
    return out


def _compose(hs):
    def f(t):
        for h in hs:
            t = h(t)
        return t
    return f


if __name__ == "__main__":
    run()
