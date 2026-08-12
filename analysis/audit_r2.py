"""Independent check of the second review's four central claims."""
from __future__ import annotations

import numpy as np
import sympy as sp

import masses as M

GSC = M.gamma_sc("C")
F0 = M.offset_F0("C")
MUH, MUD, MUT = M.mu_ratios("C")          # 1, d, t


def X(A, w):
    f = lambda mu: -0.5 * np.log1p(w * mu) - A * mu / (1.0 + w * mu)
    lt = f(MUT)
    return f(MUH) - lt, f(MUD) - lt


def F(A, w):
    xh, xd = X(A, w)
    return xh - GSC * xd


def hdr(t):
    print("\n" + "=" * 76 + f"\n{t}\n" + "=" * 76)


# ---------------------------------------------------------------- CLAIM I
hdr("CLAIM I: Theorem 1 is FALSE -- F exceeds F0 at small A")
print(f"F0 = {F0:+.6f}   mu ratios: 1, {MUD:.6f}, {MUT:.6f}   sqrt(t) = {np.sqrt(MUT):.6f}")
print(f"\n{'A':>10}{'w':>10}{'F':>12}{'> F0 ?':>9}{'K_HT^int':>11}")
for A in (0.01, 0.05, 0.2, 1.0):
    for w in (1e-6, 1e-2):
        xh, _ = X(A, w)
        print(f"{A:10.3f}{w:10.0e}{F(A,w):12.6f}{'YES' if F(A,w) > F0 else 'no':>9}"
              f"{np.exp(xh):11.4f}")

print("\n  analytic: dF/dA = q_H - gamma_SC q_D")
d, t, g, wv = sp.symbols("d t gamma w", positive=True)
qH = (t - 1) / ((1 + wv) * (1 + wv * t))
qD = (t - d) / ((1 + wv * d) * (1 + wv * t))
gsc_sym = d * (t - 1) / (t - d)
expr = sp.simplify(qH - gsc_sym * qD)
print(f"  symbolic simplification: {expr}")
tgt = -(t - 1) * (d - 1) / ((1 + wv * t) * (1 + wv) * (1 + wv * d))
print(f"  equals -(t-1)(d-1)/[(1+wt)(1+w)(1+wd)] : "
      f"{sp.simplify(expr - tgt) == 0}")
print("  => dF/dA < 0 strictly, so the supremum in A is at A -> 0, not A -> inf.")

print(f"\n  global sup over (A,w): F(A->0, w->0) = {F(1e-9, 1e-9):+.8f}")
print("  => sup F = 0, NOT F0.  Theorem 1 as stated is FALSE.  CLAIM CONFIRMED.")

# the scale-qualified version
hdr("CLAIM I(b): the corrected scale-dependent envelope B(K)")


def supF_at_K(K, n=4000):
    """max F over (A,w) subject to ln K_HT^int = ln K."""
    L = np.log(K)
    best = -np.inf
    for w in np.logspace(-8, 8, n):
        f = lambda mu: -0.5 * np.log1p(w * mu)
        a_coef = -(MUH / (1 + w * MUH) - MUT / (1 + w * MUT))
        const = f(MUH) - f(MUT)
        if abs(a_coef) < 1e-300:
            continue
        A = (L - const) / a_coef
        if A <= 0:
            continue
        best = max(best, F(A, w))
    return best


def B_closed(K):
    """Reviewer's closed form for 1 < K < sqrt(t)."""
    return np.log(K) - 0.5 * GSC * np.log(
        K ** 2 * (MUT - 1) / (MUT - MUD + K ** 2 * (MUD - 1)))


print(f"threshold sqrt(t) = {np.sqrt(MUT):.6f}")
print(f"{'K_HT^int':>10}{'numeric sup':>14}{'closed form':>14}{'F0':>12}")
for K in (1.01, 1.05, 1.10, 1.20, 1.2689, 1.5, 3.0, 11.0):
    num = supF_at_K(K)
    cf = B_closed(K) if K < np.sqrt(MUT) else F0
    print(f"{K:10.4f}{num:14.6f}{cf:14.6f}{F0:12.6f}")
print("  => Theorem 1 holds only for K_HT^int >= sqrt(mu_T/mu_H) = 1.2689.")
print("     Below that the envelope rises toward 0.  My numerical check started")
print("     at K = 1.5, which sat entirely inside the region where it holds.")

# --------------------------------------------------------------- CLAIM III
hdr("CLAIM III: the experimental recommendation has the commitment reversed")
x = 8.0
print("observation map K_obs = x(1+c)/(x+c), c = k_off/k_chem")
print(f"{'c':>10}{'K_obs':>10}{'masking':>12}")
for c in (0.01, 0.1, 1, 10, 100, 1e4):
    kobs = x * (1 + c) / (x + c)
    print(f"{c:10.2f}{kobs:10.4f}{x-kobs:12.4f}")
print(f"  intrinsic x = {x}")
print("  K_obs -> x as c -> INFINITY.  Masking is reduced by INCREASING")
print("  c = k_off/k_chem.  The manuscript recommends REDUCING it.")
print("  CLAIM CONFIRMED: sign error in the design recommendation.")

# ---------------------------------------------------------------- CLAIM IV
hdr("CLAIM IV: the model-free bounds ignore measurement uncertainty")
kht, sh, kdt, sd = 5.04, 0.03, 1.65, 0.02          # ecDHFR W133F at 25 C
Fobs = np.log(kht) - GSC * np.log(kdt)
se = np.sqrt((sh / kht) ** 2 + (GSC * sd / kdt) ** 2)
print(f"ecDHFR W133F, 25 C:  K_HT = {kht} +- {sh},  K_DT = {kdt} +- {sd}")
print(f"  F_obs        = {Fobs:+.5f}")
print(f"  delta-method SE (independence) = {se:.5f}")
print(f"  distance to F0 = {F0 - Fobs:+.5f}  ->  {abs(F0-Fobs)/se:.2f} SE")
print("  CLAIM CONFIRMED: the reported shortfall is well inside one standard")
print("  error, so 'falls short by 0.018' is not an inferential statement.")
