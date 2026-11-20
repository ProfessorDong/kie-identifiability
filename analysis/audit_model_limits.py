"""Checks of the limiting behavior of the tunneling model.

Each statement is tested against the model as implemented in qtunnel.py: which
constraints cap the exponent, where the semiclassical value is approached, the
ridge limit, the finite-scale bound, the dimensions of the parameters and the
mass convention.
"""
from __future__ import annotations

import numpy as np
import sympy as sp

import masses as M
import qtunnel as Q

MU = {"H": Q.MU_H, "D": Q.MU_D, "T": Q.MU_T}
GSC_ATOMIC = Q.gamma_rigid() * Q.MU_D / Q.MU_H


def hdr(t):
    print("\n" + "=" * 76 + f"\n{t}\n" + "=" * 76)


# ---------------------------------------------------------------- CHECK 1
hdr("CHECK 1: an UPPER bound sigma/R0 <= eps does NOT cap gamma_TUN")
print("counterexample: w = n, A = n^3  =>  r -> 0 but A/w^2 -> inf")
print(f"{'n':>8}{'r=sigma/R0':>14}{'A/w^2':>12}{'gamma_TUN':>12}")
for n in (10, 1e2, 1e3, 1e4, 1e5, 1e6):
    w, A = n, n ** 3
    r = np.sqrt(w / (2 * A))
    print(f"{n:8.0e}{r:14.3e}{A/w**2:12.3g}{Q.gamma_tunneling(A, w):12.6f}")
print(f"\nsemiclassical gamma_SC (reduced masses, as in the model) = {GSC_ATOMIC:.6f}")
print("=> sup over {r <= eps} equals gamma_SC for EVERY eps > 0.")

print("\nOptimizing only on the boundary w = 2*A*eps^2 gives smaller numbers.")
print("Set-theoretic check: {r<=eps} grows with eps, so its max cannot decrease.")
print(f"{'eps':>8}{'max on boundary r=eps':>26}{'true sup over r<=eps':>24}")
for eps in (0.01, 0.05, 0.10, 0.20, 1 / 3, 1.0):
    onbdry = max(Q.gamma_tunneling(A, 2 * A * eps ** 2)
                 for A in np.logspace(-3, 14, 6000))
    # true sup over the region: free to take A >> w
    region = max(Q.gamma_tunneling(A, w)
                 for w in np.logspace(-3, 6, 300)
                 for A in [w / (2 * eps ** 2) * f for f in (1, 10, 1e3, 1e6, 1e9)])
    print(f"{eps:8.3f}{onbdry:26.4f}{region:24.4f}")
print("The boundary maxima are NOT the constrained sup.")

# ---------------------------------------------------------------- CHECK 1b
hdr("CHECK 1b: a LOWER bound sigma/R0 >= eps DOES cap gamma_TUN")
print(f"{'eps':>8}{'sup gamma_TUN over r >= eps':>30}")
for eps in (0.01, 0.05, 0.10, 0.20, 1 / 3, 1.0):
    best = -np.inf
    for w in np.logspace(-4, 8, 900):
        Amax = w / (2 * eps ** 2)          # r >= eps  <=>  A <= w/(2 eps^2)
        for A in np.logspace(-6, np.log10(Amax), 60):
            best = max(best, Q.gamma_tunneling(A, w))
    print(f"{eps:8.3f}{best:30.4f}")
print("Caps are strictly below gamma_SC and DECREASE with eps:")
print("the discriminating constraint is a LOWER bound on the gating amplitude.")

# ---------------------------------------------------------------- CHECK 2
hdr("CHECK 2: gamma_SC is reached only as the intrinsic KIEs DIVERGE")
print(f"{'A/w^2':>10}{'gamma_TUN':>12}{'ln K_HT^int':>14}{'K_HT^int':>14}")
for s in (1, 10, 1e2, 1e3, 1e4, 1e5):
    w = 1e6
    A = s * w ** 2
    lh, ld = Q.log_kie_intrinsic(A, w)
    print(f"{s:10.0e}{lh/ld:12.6f}{lh:14.4f}{np.exp(lh):14.3e}")
print("=> the supremum is an asymptotic slope at divergent KIE, not a finite")
print("   boundary point in observable space.")

# ---------------------------------------------------------------- RIDGE
hdr("THE RIDGE LIMIT")
P = {k: 0.5 * np.log(Q.MU_T / v) for k, v in MU.items()}
C = {k: (Q.MU_T - v) / (v * Q.MU_T) for k, v in MU.items()}
print(f"P_H={P['H']:.6f}  P_D={P['D']:.6f}   C_H={C['H']:.6f}  C_D={C['D']:.6f}")
print(f"C_H/C_D = {C['H']/C['D']:.10f}   gamma_SC = {GSC_ATOMIC:.10f}")

print("\nverify ln K_LT^int -> P_L + s C_L  along A = s w^2, w -> inf:")
for s in (0.5, 5.0, 50.0):
    for w in (1e3, 1e5, 1e7):
        A = s * w ** 2
        lh, ld = Q.log_kie_intrinsic(A, w)
        print(f"  s={s:5.1f} w={w:8.0e}: lnK_HT {lh:10.6f} vs {P['H']+s*C['H']:10.6f}"
              f"   lnK_DT {ld:10.6f} vs {P['D']+s*C['D']:10.6f}")

print("\n=> in the ridge limit the attainable (lnK_DT, lnK_HT) set is the RAY")
print(f"   base (P_D,P_H) = ({P['D']:.4f},{P['H']:.4f})  direction (C_D,C_H)")
print(f"   whose slope C_H/C_D = gamma_SC exactly.")
print("   The semiclassical reference is the ray from the ORIGIN with the same slope.")
print("   The two are PARALLEL LINES, offset forever.")

perp = abs(P["D"] * GSC_ATOMIC - P["H"]) / np.sqrt(1 + GSC_ATOMIC ** 2)
sig = P["H"] - GSC_ATOMIC * P["D"]
print(f"\n   signed mass-scaling offset  F = lnK_HT - gamma_SC lnK_DT = {sig:.6f}")
print(f"   perpendicular separation                                   = {perp:.6f}")
print("   This is a FIXED number set by the masses alone: it does not vanish,")
print("   and it does not depend on the KIE magnitude or on A, w separately.")

# symbolic confirmation that C_H/C_D == gamma_SC
a, b, c = sp.symbols("a b c", positive=True)
CH = (c - a) / (a * c); CD = (c - b) / (b * c)
gsc = (1 / a - 1 / c) / (1 / b - 1 / c)
print(f"\n   symbolic: simplify(C_H/C_D - gamma_SC) = {sp.simplify(CH/CD - gsc)}")

# ---------------------------------------------------------------- CHECK 2b
hdr("CHECK 2b: finite-scale bound gamma_inf(L_H) < gamma_SC")
print(f"{'K_HT^int':>10}{'gamma_inf':>14}")
for k in (1.5, 2.0, 3.0, 5.0, 8.0, 11.0, 15.0, 100.0, 1e4):
    LH = np.log(k)
    s = (LH - P["H"]) / C["H"]
    g = LH / (P["D"] + s * C["D"]) if s > 0 else np.nan
    print(f"{k:10.1f}{g:14.6f}")
print(f"gamma_SC = {GSC_ATOMIC:.6f}   (approached only as K_HT^int -> inf)")

# ---------------------------------------------------------------- CHECK 3
hdr("CHECK 3: dimensions of A and w")
print("kappa = sqrt(k_f)/(2 hbar):  [k_f]=M T^-2, [hbar]=M L^2 T^-1")
print("  => [kappa] = M^(1/2) T^-1 / (M L^2 T^-1) = M^(-1/2) L^-2")
print("  so A = kappa R0^2 carries M^(-1/2): NOT dimensionless.")
print("  correct: A = kappa sqrt(m_H) R0^2,  w = 2 kappa sqrt(m_H) sigma^2.")

# ---------------------------------------------------------------- CHECK 4
hdr("CHECK 4: mass-convention consistency")
def gsc_of(mH, mD, mT):
    return (mH ** -0.5 - mT ** -0.5) / (mD ** -0.5 - mT ** -0.5)
print(f"mass numbers 1,2,3          : gamma_SC = {gsc_of(1,2,3):.6f}")
print(f"atomic masses               : gamma_SC = "
      f"{gsc_of(M.M_H_ATOMIC, M.M_D_ATOMIC, M.M_T_ATOMIC):.6f}")
print(f"C-H/C-D/C-T reduced masses  : gamma_SC = {gsc_of(*M.triple('C')):.6f}")
print("The supremum identity holds in whatever masses are used, so the tunneling")
print("supremum EQUALS gamma_SC in the same convention.  Comparing a model")
print("evaluated with one set of masses against gamma_SC from another mixes")
print("conventions and manufactures a gap.")
