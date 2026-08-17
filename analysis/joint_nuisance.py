"""Joint identified set under binding isotope effects and an isotope-blind bypass.

The supplement bounds each of these nuisances on its own: it reports how large a
binding offset F_bind the exclusion tolerates, and how large a bypass ratio phi
it tolerates, each with the other absent.  Those two numbers do not compose.
Both mechanisms act on the same observation map,

    K_i = alpha_i * (x_i + phi_i)/(1 + phi_i)
                  * (q_i + 1 + phi_i)/(beta_i q_i + x_i + phi_i),

which is the first-passage result for a complex that reaches product either
through the isotope-sensitive step or through a blind route, and returns to free
enzyme at an isotope-dependent rate, entered at an isotope-dependent rate.  With
unequal tritium references, phi_D = r phi_H and q_D = r q_H.

Setting phi = 0 recovers Proposition S8; setting alpha = beta = 1 recovers
Proposition S6.  Holding one at its null while bounding the other is exactly
the marginal calculation the supplement reports, and this module computes the
joint region instead.

Run from analysis/:   python joint_nuisance.py
"""
from __future__ import annotations

import numpy as np

import masses as M

G = M.gamma_sc("C")
F0 = M.offset_F0("C")

# the decisive observation
KH_OBS, KD_OBS, R_YADH = 7.130, 1.730, 1.31


def x_from_obs(K, phi, q, alpha, beta):
    """Invert the joint map for the intrinsic effect on the sensitive step.

    K = alpha (x+phi)(1+phi+q) / [(1+phi)(beta q + x + phi)] gives
    x = A beta q/(1+phi+q-A) - phi with A = K(1+phi)/alpha.
    """
    A = K * (1.0 + phi) / alpha
    den = 1.0 + phi + q - A
    with np.errstate(divide="ignore", invalid="ignore"):
        x = A * beta * q / den - phi
    return np.where(den > 0, x, np.nan)


def endpoint(KH=KH_OBS, KD=KD_OBS, F_bind=0.0, phi_H=0.0, r=R_YADH,
             put_binding_in="beta", n=240000):
    """Lower endpoint of the identified set for F_int, given the nuisances.

    Sweeps the reciprocal commitment over a wide range and returns the infimum
    of F_int over admissible commitments.  The exclusion survives exactly when
    this exceeds F0.

    F_bind = ln(alpha_H/beta_H) - gamma ln(alpha_D/beta_D) does not by itself
    fix the four binding constants, so `put_binding_in` selects which one
    carries it; the two choices are compared in `binding_placement()`.
    """
    aH = aD = bH = bD = 1.0
    if put_binding_in == "beta":
        bD = np.exp(F_bind / G)          # ln(aD/bD) = -F_bind/gamma
    elif put_binding_in == "alpha":
        aH = np.exp(F_bind)              # ln(aH/bH) = F_bind
    else:
        raise ValueError(put_binding_in)

    phi_D = r * phi_H
    qH = np.geomspace(1e-3, 1e9, n)
    qD = r * qH
    xH = x_from_obs(KH, phi_H, qH, aH, bH)
    xD = x_from_obs(KD, phi_D, qD, aD, bD)
    ok = np.isfinite(xH) & np.isfinite(xD) & (xH > 1.0) & (xD > 1.0)
    if not ok.any():
        return np.nan
    F = np.log(xH[ok]) - G * np.log(xD[ok])
    return float(np.min(F))


def marginal_tolerances(target=F0):
    """The two one-at-a-time tolerances the supplement reports."""
    def bisect(f, lo, hi, tol=1e-7):
        if f(lo) <= 0:
            return lo
        while hi - lo > tol:
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if f(mid) > 0 else (lo, mid)
        return 0.5 * (lo + hi)

    fb = bisect(lambda v: endpoint(F_bind=v) - target, 0.0, 2.0)
    ph = bisect(lambda v: endpoint(phi_H=v) - target, 0.0, 50.0)
    return fb, ph


def joint_boundary(target=F0, nphi=61):
    """For each bypass ratio, the largest binding offset the exclusion survives."""
    def bisect(phi, lo=0.0, hi=2.0, tol=1e-6):
        if endpoint(F_bind=lo, phi_H=phi) - target <= 0:
            return 0.0
        while hi - lo > tol:
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if endpoint(F_bind=mid, phi_H=phi) - target > 0 else (lo, mid)
        return 0.5 * (lo + hi)

    phis = np.linspace(0.0, 0.30, nphi)
    return phis, np.array([bisect(p) for p in phis])


def binding_placement(F_bind=0.09, phi_H=0.08):
    """The endpoint depends on where the binding effect sits, not only on F_bind."""
    return (endpoint(F_bind=F_bind, phi_H=phi_H, put_binding_in="beta"),
            endpoint(F_bind=F_bind, phi_H=phi_H, put_binding_in="alpha"))


def main():
    print("validation against the supplement's own endpoint routine")
    from network_geometry import endpoint_closed
    for phi in (0.0, 0.15, 0.35):
        mine = endpoint(F_bind=0.0, phi_H=phi)
        theirs = endpoint_closed(KH_OBS, KD_OBS, phi, R_YADH)
        print(f"   phi={phi:.2f}:  joint module {mine:+.5f}   "
              f"supplement {theirs:+.5f}   diff {abs(mine-theirs):.2e}")

    fb, ph = marginal_tolerances()
    print(f"\nmarginal tolerances against F0 (each with the other absent)")
    print(f"   binding only : F_bind up to {fb:.4f}")
    print(f"   bypass  only : phi_H  up to {ph:.4f}")

    print(f"\nthe reviewer's joint point: F_bind=0.090, phi_H=0.080")
    e = endpoint(F_bind=0.090, phi_H=0.080)
    print(f"   both individually inside their marginal tolerances")
    print(f"   joint endpoint = {e:+.5f}   F0 = {F0:+.5f}   "
          f"exclusion survives: {e > F0}")

    phis, fbs = joint_boundary()
    print(f"\njoint tolerance boundary, binding offset admissible at each bypass")
    for p in (0.00, 0.02, 0.05, 0.08, 0.10, 0.15, 0.20):
        i = int(np.argmin(np.abs(phis - p)))
        print(f"   phi_H={phis[i]:.3f}   F_bind up to {fbs[i]:.4f}")

    b, a = binding_placement()
    print(f"\nplacement of the binding effect at F_bind=0.09, phi_H=0.08")
    print(f"   carried by beta_D : endpoint {b:+.5f}")
    print(f"   carried by alpha_H: endpoint {a:+.5f}")
    return phis, fbs


if __name__ == "__main__":
    main()


# --------------------------------------------------------- profiled boundary
# F_bind = ln(alpha_H/beta_H) - gamma ln(alpha_D/beta_D) fixes only a
# combination of the four binding constants, and the endpoint depends on more
# than that combination once a bypass is present.  A boundary drawn in
# (F_bind, phi) is therefore a worst case only if it is minimized over the
# binding constants consistent with F_bind, inside stated individual bounds.
BIND_BOX = (0.90, 1.10)      # individual association / dissociation effects


def endpoint_profiled(F_bind, phi_H, r=R_YADH, box=BIND_BOX, n=13):
    """Worst-case endpoint over binding constants consistent with F_bind.

    Free parameters are alpha_H, beta_H, alpha_D; beta_D follows from the
    constraint.  All four are held inside `box`, which must be stated as an
    assumption because F_bind alone does not bound them.
    """
    lo, hi = box
    # F_bind is attainable inside the box only up to (1+gamma)ln(hi/lo); beyond
    # that no binding constants realize it and the point is infeasible, not safe
    if abs(F_bind) > (1.0 + G) * np.log(hi / lo) + 1e-12:
        return np.nan
    grid = np.linspace(lo, hi, n)
    worst = np.inf
    for aH in grid:
        for bH in grid:
            # ln(aD/bD) = [ln(aH/bH) - F_bind]/gamma
            lr = (np.log(aH / bH) - F_bind) / G
            for aD in grid:
                bD = aD * np.exp(-lr)
                if not (lo <= bD <= hi):
                    continue
                e = _endpoint_ab(aH, bH, aD, bD, phi_H, r)
                if np.isfinite(e):
                    worst = min(worst, e)
    return worst if np.isfinite(worst) else np.nan


def _endpoint_ab(aH, bH, aD, bD, phi_H, r=R_YADH, n=120000):
    """Endpoint with the four binding constants given explicitly."""
    phi_D = r * phi_H
    qH = np.geomspace(1e-3, 1e9, n)
    qD = r * qH
    xH = x_from_obs(KH_OBS, phi_H, qH, aH, bH)
    xD = x_from_obs(KD_OBS, phi_D, qD, aD, bD)
    ok = np.isfinite(xH) & np.isfinite(xD) & (xH > 1.0) & (xD > 1.0)
    if not ok.any():
        return np.nan
    return float(np.min(np.log(xH[ok]) - G * np.log(xD[ok])))


def joint_boundary_profiled(target=F0, nphi=13, box=BIND_BOX):
    """Largest F_bind surviving at each bypass, worst case over the box."""
    def bisect(phi, lo=0.0, hi=0.30, tol=2e-3):
        if endpoint_profiled(lo, phi, box=box) - target <= 0:
            return 0.0
        while hi - lo > tol:
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if endpoint_profiled(mid, phi, box=box) - target > 0 else (lo, mid)
        return 0.5 * (lo + hi)
    phis = np.linspace(0.0, 0.18, nphi)
    return phis, np.array([bisect(p) for p in phis])


def corner_sufficient(F_bind, phi_H, box, n=9):
    """Check the worst case really sits on the box boundary before using corners."""
    lo, hi = box
    full = endpoint_profiled(F_bind, phi_H, box=box, n=n)
    corners = np.inf
    for aH in (lo, hi):
        for bH in (lo, hi):
            lr = (np.log(aH / bH) - F_bind) / G
            for aD in (lo, hi):
                bD = aD * np.exp(-lr)
                if lo <= bD <= hi:
                    e = _endpoint_ab(aH, bH, aD, bD, phi_H)
                    if np.isfinite(e):
                        corners = min(corners, e)
    return full, corners, abs(full - corners) < 1e-6


def box_sensitivity(halfwidths=(0.02, 0.05, 0.10, 0.15, 0.20), nphi=8):
    """How the profiled joint boundary depends on the assumed binding range.

    The box is an assumption, not a measurement, so the honest presentation is
    the boundary as a function of how wide it is taken to be.
    """
    out = []
    for h in halfwidths:
        box = (1.0 - h, 1.0 + h)
        phis = np.linspace(0.0, 0.16, nphi)
        fbs = []
        for p in phis:
            lo, hi = 0.0, 0.30
            def ok(v):
                e = endpoint_profiled(v, p, box=box, n=7)
                return np.isfinite(e) and e - F0 > 0
            if not ok(lo):
                fbs.append(0.0); continue
            while hi - lo > 2e-3:
                mid = 0.5 * (lo + hi)
                lo, hi = (mid, hi) if ok(mid) else (lo, mid)
            fbs.append(0.5 * (lo + hi))
        fbs = np.array(fbs)
        m = fbs > 0
        if m.sum() > 1:
            A = np.vstack([phis[m], np.ones(m.sum())]).T
            sl, ic = np.linalg.lstsq(A, fbs[m], rcond=None)[0]
        else:
            sl = ic = np.nan
        out.append((h, ic, sl, float(fbs[0]), float(phis[m].max()) if m.any() else 0.0))
    return out
