"""Data for the atlas figure: topology fixes the geometry.

One representative of each class the enumeration produces.  Nothing here is
drawn by hand: every coefficient is the symbolic Moebius result for a named
mechanism at stated rate constants, taken from network_atlas.

  concave   committed step, no blind route      set opens above
  mixed     the same step with a blind bypass   direction turns on x* against 1
  blind     step below the last commitment      nothing is identified
  inverse   blind route plus a late exit        AD < BC, directions reverse

The inverse case is the smallest network in the enumeration that shows it:
speeding the labeled step moves flux out of the state that can still reach
product by the blind route, and into one that leaks back to free enzyme, so a
normal intrinsic effect is observed as an inverse one.

Panel B plots the observation map in log coordinates, where mass scaling is
homogeneity and masking is curvature.  Panel C plots the offset along the
mass-scaling ray, which is the geometry the map implies and is exactly the
quantity the verification profiles.
"""
from __future__ import annotations

import sympy as sp
import numpy as np

import masses as M
from network_atlas import (classify, curvature_point, draws, enumerate_family,
                           moebius, profile_geometry)

G = M.gamma_sc("C")
OUT = "../figures/tikz/data/"

FAM = {m.name: m for m in enumerate_family(3)}


def _sub(**kw):
    return {sp.Symbol(k, positive=True): v for k, v in kw.items()}


def coeffs(name, sub):
    return moebius(FAM[name]), sub


CLASSES = [
    ("concave", "n=1; iso on chemistry", _sub(k_off=3.0, k_T=1.0)),
    ("mixed", "n=1; iso on chemistry; blind bypass",
     _sub(k_off=12.5, k_byp=0.6, k_T=1.0)),
    ("blind", "n=2; iso on chemistry", _sub(k_off=3.0, k_f1=1.0, k_T=1.0)),
    ("inverse", "n=2; iso on step 1; blind bypass; late exit",
     _sub(k_off=0.6, k_byp=1.0, k_cat=2.5, k_off2=3.5, k_T=5.0)),
]


def resolve(name, sub):
    """(A,B,C,D) for a named mechanism; the inverse case uses an atlas draw."""
    co = moebius(FAM[name])
    if sub is None:
        sub = draws(co, name, k=12)[0]
    return co, sub, [float(sp.N(z.subs(sub))) for z in co]


def main():
    print(f"{'class':9s} {'A':>9s} {'B':>9s} {'C':>9s} {'D':>9s} "
          f"{'det':>10s} {'x*':>9s}  geometry")
    for tag, name, sub in CLASSES:
        co, sub, (A, B, C, D) = resolve(name, sub)
        det = A * D - B * C
        xs = curvature_point(co, sub)
        _, geo, _ = classify(co, sub)
        pgeo, mag = profile_geometry(co, sub)
        assert pgeo == geo, f"{tag}: {geo} vs profiled {pgeo}"
        print(f"{tag:9s} {A:9.4f} {B:9.4f} {C:9.4f} {D:9.4f} "
              f"{det:10.4f} {xs:9.4f}  {geo}")

        K = lambda x: (A * x + B) / (C * x + D)

        # panel B: the observation map in log coordinates
        x = np.geomspace(1.0, 40.0, 220)
        with open(OUT + f"fn_map_{tag}.dat", "w") as f:
            f.write("lx lk\n")
            for u in x:
                f.write(f"{np.log(u):.6f} {np.log(K(u)):.6f}\n")

        # panel C: the offset along the mass-scaling ray
        xd = np.geomspace(1.0, 2.2, 260)
        F = np.log(K(xd ** G)) - G * np.log(K(xd))
        with open(OUT + f"fn_off_{tag}.dat", "w") as f:
            f.write("xd F\n")
            for u, v in zip(xd, F):
                f.write(f"{u:.6f} {v:.6f}\n")

        # the curvature switch, where it falls inside the plotted range
        if np.isfinite(xs) and 1.0 < xs < 40.0:
            with open(OUT + f"fn_xstar_{tag}.dat", "w") as f:
                f.write("lx lk\n")
                f.write(f"{np.log(xs):.6f} {np.log(K(xs)):.6f}\n")
            print(f"          curvature switch at x* = {xs:.4f}")

    # census over the full enumeration, for the panel label
    fam4 = enumerate_family(4)
    from network_atlas import blind_by_topology
    nb = sum(1 for m in fam4 if blind_by_topology(m))
    nby = sum(1 for m in fam4 if not blind_by_topology(m) and "bypass" in m.note)
    with open(OUT + "fn_census.dat", "w") as f:
        f.write("class n\n")
        f.write(f"concave {len(fam4) - nb - nby}\n")
        f.write(f"mixed {nby}\n")
        f.write(f"blind {nb}\n")
    print(f"\ncensus over {len(fam4)} mechanisms: "
          f"concave {len(fam4) - nb - nby}, conditional {nby}, blind {nb}")


if __name__ == "__main__":
    main()
