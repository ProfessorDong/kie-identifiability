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


def series_steps(n_draw=400000, seed=20261022):
    """Chains of isotope-sensitive steps in series (Supplementary Note 3).

    Each step i masks its own intrinsic pair through the series map at its own
    reciprocal commitment c_i, and the observed effects are the products, so
    the observed offset is the sum of the per-step observed offsets.

    Returns (example, (violations, n, worst excess, total offset at the worst),
    violations with every step's own offset nonnegative):
      * example: two steps with k_off/k_T = 1e3, per-step D/T effect 2 and
        per-step offset 0.1 -> observed aggregate offset (aggregate 0.2);
      * separate pairs, one step above the ray and one below, total offset
        >= 0: how often, and by how much, F_obs exceeds F_int.  Draws: per-step
        D/T effect log-uniform on [1.05, 5], offsets uniform on [0, 0.5] and
        [-0.5, 0], c_i log-uniform on [1e-2, 1e3].
    Until 2026-10-22 the supplement quoted these numbers without code.
    """
    ser = lambda x, c: x * (1 + c) / (x + c)
    xd = 2.0
    xh = np.exp(0.1 + GSC * np.log(xd))
    ex = 2 * (np.log(ser(xh, 1e3)) - GSC * np.log(ser(xd, 1e3)))
    rng = np.random.default_rng(seed)
    n = n_draw
    xd = np.exp(rng.uniform(np.log(1.05), np.log(5.0), (n, 2)))
    f1 = rng.uniform(0.0, 0.5, n)
    f2 = rng.uniform(-0.5, 0.0, n)
    keep = f1 + f2 >= 0
    f = np.c_[f1, f2][keep]
    xd = xd[keep]
    c = np.exp(rng.uniform(np.log(1e-2), np.log(1e3), (len(f), 2)))
    xh = np.exp(f + GSC * np.log(xd))
    fobs = (np.log(ser(xh, c)) - GSC * np.log(ser(xd, c))).sum(axis=1)
    fint = f.sum(axis=1)
    exc = fobs - fint
    i = int(np.argmax(exc))
    viol = (int((exc > 1e-12).sum()), len(f), float(exc[i]), float(fint[i]))
    # the same draws with each step's own offset nonnegative
    f_pos = np.abs(f)
    xh2 = np.exp(f_pos + GSC * np.log(xd))
    exc2 = (np.log(ser(xh2, c)) - GSC * np.log(ser(xd, c))).sum(axis=1) - f_pos.sum(axis=1)
    return float(ex), viol, int((exc2 > 1e-12).sum())


def sign_and_binding_draws(n=400000, seed=20261022):
    """Random-draw checks quoted in Supplementary Note 2 and Note 4.

    (1) Corollary S4.1: with c_D = r c_H and r >= 1, a non-positive intrinsic
        offset never gives a positive observed one.  Draws: x_D log-uniform on
        [1.01, 10], F_int uniform on [-1, 0], c_H log-uniform on [1e-2, 1e3],
        r uniform on [1, 3]; kept when x_H > x_D.
    (2) Proposition S7: with binding effects alpha_i, beta_i uniform on
        [0.8, 1.2] and the corrected pair in normal order, F_int > F_obs - F_bind
        whenever the corrected criterion L~_H >= gamma holds, and how often it
        fails when that criterion does not hold.  Draws: x_D as above, F_int
        uniform on [-0.5, 1], c log-uniform on [1e-2, 1e3].
    Returns ((kept, violations), (kept with criterion, violations,
    kept without criterion, failures without it)).  Until 2026-10-22 the
    supplement quoted draw counts for these checks without code.
    """
    M = lambda x, c: x * (1 + c) / (x + c)
    rng = np.random.default_rng(seed)
    xd = np.exp(rng.uniform(np.log(1.01), np.log(10.0), n))
    f = rng.uniform(-1.0, 0.0, n)
    xh = np.exp(f + GSC * np.log(xd))
    c = np.exp(rng.uniform(np.log(1e-2), np.log(1e3), n))
    r = rng.uniform(1.0, 3.0, n)
    ok = xh > xd
    fobs = np.log(M(xh, c)) - GSC * np.log(M(xd, r * c))
    s41 = (int(ok.sum()), int((fobs[ok] > 1e-12).sum()))
    xd = np.exp(rng.uniform(np.log(1.01), np.log(10.0), n))
    f = rng.uniform(-0.5, 1.0, n)
    xh = np.exp(f + GSC * np.log(xd))
    c = np.exp(rng.uniform(np.log(1e-2), np.log(1e3), n))
    aH, bH, aD, bD = (rng.uniform(0.8, 1.2, n) for _ in range(4))
    KH = aH * xh * (1 + c) / (xh + bH * c)
    KD = aD * xd * (1 + c) / (xd + bD * c)
    tH, tD = KH / aH, KD / aD
    norm = (tH > tD) & (tD > 1)
    crit = norm & ((tH - 1) / (tD - 1) >= GSC)
    fbind = np.log(aH / bH) - GSC * np.log(aD / bD)
    lhs = f
    rhs = np.log(KH) - GSC * np.log(KD) - fbind
    fail = lhs <= rhs + 1e-12
    nocrit = norm & ~crit
    s7 = (int(crit.sum()), int(fail[crit].sum()), int(nocrit.sum()), int(fail[nocrit].sum()))
    return s41, s7


def remaining_draw_checks(seed=20261022):
    """Three random-draw statements of the supplement, with code since 2026-10-22.

    (1) Supplementary Note 5: with separable rates k_ab = p_a s_b, normal
        primary order and an observed normal secondary H/T effect, the observed
        ratio R satisfies R <= r when r >= 1 and R < 1 when r < 1 (10^6 draws).
    (2) Proposition S5: for the mixed class, the offset along the mass-scaling
        ray takes both signs exactly when x* > 1 (2 x 10^4 random maps).
    (3) Supplementary Note 8: the two-viscosity inversion returns the intrinsic
        effect exactly (2 x 10^4 random admissible triples; worst relative error).
    Returns (exceptions to (1) among the draws that qualify, qualifying draws),
    (exceptions to (2), maps), worst relative error of (3).
    """
    rng = np.random.default_rng(seed)
    n = 1_000_000
    pH, pD = np.sort(np.exp(rng.uniform(0, np.log(50.0), (n, 2))), axis=1)[:, ::-1].T
    sT = np.ones(n)
    sH, sD = np.exp(rng.uniform(-0.5, 0.5, (2, n)))
    d = np.exp(rng.uniform(np.log(1e-2), np.log(1e3), n))
    v = lambda p, s: p * s / (d + p * s)
    SHT = v(pH, sH) / v(pH, sT)
    SDT = v(pD, sD) / v(pD, sT)
    R = SHT / SDT
    r = sH / sD
    q = SHT > 1
    exc1 = int((q & (r >= 1) & (R > r * (1 + 1e-12))).sum() + (q & (r < 1) & (R >= 1)).sum())
    # (2) mixed maps: K(x) = (Ax+B)/(Cx+D), A+B = C+D, AD > BC, B, C > 0
    m = 20000
    exc2 = 0
    xd = np.geomspace(1 + 1e-6, 1e4, 4000)
    for _ in range(m):
        A, C = np.exp(rng.uniform(-3, 3, 2))
        B = np.exp(rng.uniform(-3, 3))
        D = A + B - C
        if D <= 0 or A * D <= B * C:
            continue
        K = lambda x: (A * x + B) / (C * x + D)
        F = np.log(K(xd ** GSC)) - GSC * np.log(K(xd))
        both = F.max() > 1e-12 and F.min() < -1e-12
        xstar = np.sqrt(B * D / (A * C))
        if both != (xstar > 1) and xstar > 1:
            # A positive lobe that double precision cannot resolve (x* within a
            # fraction of a percent of 1): look for it at 50 digits close to the
            # origin, where the ray offset is (gamma(gamma-1)/2) h''(0) t^2 with
            # the sign of x* - 1.
            import mpmath as mp
            mp.mp.dps = 50
            Am, Bm, Cm_, Dm = (mp.mpf(float(v)) for v in (A, B, C, D))
            Km = lambda z: (Am * z + Bm) / (Cm_ * z + Dm)
            g = mp.mpf(float(GSC))
            both = any(mp.log(Km(z ** g)) - g * mp.log(Km(z)) > 0
                       for z in (mp.mpf(1) + mp.mpf(10) ** -k for k in range(1, 12)))
        exc2 += both != (xstar > 1)
    # (3) two-viscosity inversion
    import completion as Cm
    x = np.exp(rng.uniform(np.log(1.05), np.log(30.0), m))
    c = np.exp(rng.uniform(np.log(1e-2), np.log(1e3), m))
    rr = np.exp(rng.uniform(np.log(1.1), np.log(10.0), m))
    Kf = lambda x_, c_: x_ * (1 + c_) / (x_ + c_)
    rec = Cm.intrinsic_from_two_viscosities(Kf(x, c), Kf(x, c / rr), rr)
    err3 = float(np.max(np.abs(rec / x - 1)))
    return (exc1, int(q.sum())), (int(exc2), m), err3


if __name__ == "__main__":
    run()
