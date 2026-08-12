"""Adversarial audit of the reconstructed results, before any manuscript is written.

Each section tries to BREAK a claim rather than confirm it.  Two earlier errors
in this project (optimizing on a constraint boundary instead of a region, and
asserting a monotonicity that failed on 21/73 records) were both of the form
"checked the easy case, generalized without testing".  So every claim here is
probed off the easy case.
"""
from __future__ import annotations

import numpy as np
import sympy as sp

import masses as M
import qtunnel as Q
import ridge as Rg

GSC = M.gamma_sc("C")
F0 = M.offset_F0("C")
MUH, MUD, MUT = M.mu_ratios("C")
FAIL = []


def hdr(t):
    print("\n" + "=" * 78 + f"\n{t}\n" + "=" * 78)


def check(name, ok, detail=""):
    print(f"  {'PASS' if ok else 'FAIL'}  {name:56s} {detail}")
    if not ok:
        FAIL.append(name)


# --- rebuild the model in the C-H reduced-mass convention -------------------
def log_kie(A, w):
    """ln K_HT^int, ln K_DT^int for the reduced-mass triple."""
    f = lambda mu: -0.5 * np.log1p(w * mu) - A * mu / (1.0 + w * mu)
    lt = f(MUT)
    return f(MUH) - lt, f(MUD) - lt


def F_of(A, w):
    lh, ld = log_kie(A, w)
    return lh - GSC * ld


hdr("1. DOES THE MODEL ACTUALLY PREDICT F = F0?  (the dangerous claim)")
print("  F0 is the RIDGE value. Scan the whole (A,w) plane for the true range.")
As = np.logspace(-4, 8, 260)
ws = np.logspace(-6, 8, 260)
FF = np.array([[F_of(A, w) for w in ws] for A in As])
print(f"  F over the scanned plane: [{FF.min():+.4f}, {FF.max():+.4f}]")
print(f"  F0 = {F0:+.6f}   semiclassical F = 0")
check("F0 is NOT the model's only prediction", FF.min() < F0 - 0.1,
      f"model reaches F = {FF.min():+.3f}")
# rigid limit analytic
brack = (MUT - MUH) - GSC * (MUT - MUD)
print(f"\n  rigid limit (w->0): F -> A*[(mu_T-mu_H) - gamma_SC (mu_T-mu_D)]")
print(f"                        = A * ({brack:+.6f})   so F -> "
      f"{'-inf' if brack < 0 else '+inf'} as A grows")
check("rigid-limit slope is negative", brack < 0, f"bracket = {brack:.6f}")
for A in (0.1, 1.0, 10.0, 100.0):
    print(f"    A={A:6.1f}, w=1e-6 : F = {F_of(A, 1e-6):+.4f}"
          f"   (analytic {A*brack:+.4f})")

print("\n  => the gated model spans a WIDE range of F. F0 is the value on the")
print("     ridge only. Any claim must be phrased about the ridge, not 'the")
print("     tunneling model predicts F0'.")

hdr("2. IS THE RIDGE PHYSICALLY REACHABLE?")
print("  A = kappa sqrt(m_H) R0^2,  w = 2 kappa sqrt(m_H) sigma^2")
print("  s = A/w^2 fixed with w -> inf requires sigma^2 -> inf or R0 -> inf.")
print(f"  {'s':>8}{'w':>10}{'sigma/R0':>12}{'F':>10}")
for s in (1.0, 10.0):
    for w in (1e1, 1e3, 1e5):
        A = s * w ** 2
        print(f"  {s:8.1f}{w:10.0e}{np.sqrt(w/(2*A)):12.3e}{F_of(A, w):+10.4f}")
print("  sigma/R0 -> 0 along the ridge, i.e. the ridge is the NARROW-barrier")
print("  corner, not a large-gating one.  It is an asymptotic boundary of the")
print("  attainable set, and F -> F0 there.")
check("F -> F0 along the ridge", abs(F_of(10*1e10**2, 1e10) - F0) < 1e-3,
      f"F = {F_of(10*1e10**2, 1e10):+.6f} vs F0 = {F0:+.6f}")

hdr("3. WHAT IS THE ATTAINABLE F AT PHYSICALLY PLAUSIBLE PARAMETERS?")
print("  Constrain to observed-scale intrinsic KIEs (2 <= K_HT^int <= 30) and")
print("  a gating ratio sigma/R0 in [0.01, 0.5].")
rows = []
for A in np.logspace(-3, 10, 700):
    for w in np.logspace(-6, 8, 700):
        r = np.sqrt(w / (2 * A))
        if not (0.01 <= r <= 0.5):
            continue
        lh, ld = log_kie(A, w)
        if not (np.log(2) <= lh <= np.log(30)):
            continue
        rows.append((lh - GSC * ld, lh, r))
rows = np.array(rows)
print(f"  admissible parameter points: {len(rows)}")
if len(rows):
    print(f"  F range there : [{rows[:,0].min():+.4f}, {rows[:,0].max():+.4f}]")
    print(f"  F0 = {F0:+.6f}, semiclassical 0")
    inside = ((rows[:, 0] > F0 - 0.01) & (rows[:, 0] < F0 + 0.01)).mean()
    print(f"  fraction within 0.01 of F0: {inside:.3f}")
check("physically constrained F still spans a wide range",
      len(rows) > 0 and (rows[:, 0].max() - rows[:, 0].min()) > 0.1,
      f"span {rows[:,0].max()-rows[:,0].min():.3f}" if len(rows) else "")

hdr("4. SYMBOLIC CHECK OF THE dF/dc EXPANSION")
c, kh, kd, g = sp.symbols("c K_H K_D gamma", positive=True)
xh = kh * c / (1 + c - kh)
xd = kd * c / (1 + c - kd)
F = sp.log(xh) - g * sp.log(xd)
dF = sp.simplify(sp.diff(F, c))
lead = sp.simplify(sp.limit(dF * c ** 2, c, sp.oo))
print(f"  lim c^2 dF/dc = {sp.simplify(lead)}")
target = g * (kd - 1) - (kh - 1)
check("leading coefficient equals gamma(K_D-1)-(K_H-1)",
      sp.simplify(lead - target) == 0, f"difference {sp.simplify(lead-target)}")

hdr("5. IDENTIFIABLE-SET PROPOSITION, STRESS-TESTED OFF THE DATA")
print("  random (K_HT, K_DT) with K_HT > K_DT > 1, including extreme values")
rng = np.random.default_rng(11)
bad = 0
for _ in range(4000):
    kdt = 1 + 10 ** rng.uniform(-3, 1.0)
    kht = kdt * (1 + 10 ** rng.uniform(-3, 1.5))
    lo = kht - 1.0
    cc = lo * (1 + np.logspace(-10, 10, 4000))
    x = lambda K: K * cc / (1.0 + cc - K)
    Fc = np.log(x(kht)) - GSC * np.log(x(kdt))
    interior = 0 < int(np.argmin(Fc)) < len(cc) - 1
    predicted = GSC * (kdt - 1) - (kht - 1) > 0
    if interior != predicted:
        bad += 1
check("L_H criterion predicts interior minimum on random draws", bad == 0,
      f"{bad}/4000 mismatches")

print("\n  also: is x(c) > K always (intrinsic exceeds observed)?")
viol = 0
for _ in range(2000):
    K = 1 + 10 ** rng.uniform(-3, 1.5)
    cc = (K - 1) * (1 + np.logspace(-6, 6, 500))
    if np.any(K * cc / (1 + cc - K) < K - 1e-12):
        viol += 1
check("commitment map always masks (x > K)", viol == 0, f"{viol} violations")

hdr("6. MASS-POWER INVERSION, IN THE REDUCED-MASS CONVENTION")
rig = (MUH - MUT) / (MUD - MUT)
CH, CD = M.ridge_C("C")
print(f"  rigid   exponent uses sqrt(mu) differences : {rig:.5f}")
print(f"  ridge   exponent uses 1/sqrt(mu) diffs     : {CH/CD:.5f}")
print(f"  gamma_SC (zero-point energy)               : {GSC:.5f}")
check("rigid exponent equals M.gamma_rigid", abs(rig - M.gamma_rigid("C")) < 1e-12)
check("ridge exponent equals gamma_SC", abs(CH / CD - GSC) < 1e-12)
check("C_L equals the ZPE mass factor",
      abs(CH - (1 / MUH - 1 / MUT)) < 1e-12 and abs(CD - (1 / MUD - 1 / MUT)) < 1e-12)

hdr("7. COMMITMENT MASKING VERSUS THE SIGNAL")
print("  observed offset for a purely semiclassical intrinsic pair")
worst = 0.0
for u in (1.2, 1.5, 2.0, 3.0, 5.0, 10.0):
    for cc in (0.03, 0.1, 0.3, 1, 3, 10, 100, 1e4):
        worst = min(worst, Rg.observed_offset_semiclassical(u, cc, "C"))
print(f"  most negative observed offset found: {worst:+.4f}")
print(f"  |F0| = {abs(F0):.4f}   ratio = {abs(worst)/abs(F0):.1f}x")
check("masking exceeds the signal by an order of magnitude",
      abs(worst) / abs(F0) > 10, f"{abs(worst)/abs(F0):.1f}x")

print("\n" + "=" * 78)
print(f"FAILURES: {len(FAIL)}" + ("" if not FAIL else "  -> " + "; ".join(FAIL)))
print("=" * 78)
raise SystemExit(1 if FAIL else 0)
