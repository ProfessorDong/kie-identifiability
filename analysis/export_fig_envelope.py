"""Data for the envelope figure: what the gated model can attain, and how it moves.

Panel A is the comparison the paper makes, drawn in the plane of the intrinsic
H/T log effect X_H and the offset F.  The ground-channel model of main-text
Eq. (XL) fills the region between two edges:

  * above, the envelope B(K) of partial_id.envelope_B, a supremum approached on
    the ridge and flat at F0 once K exceeds (mu_T/mu_H)^(1/4);
  * below, the rigid limit w -> 0, where X_L = A (t - mu_L) and every pair lies
    on the ray of exponent gamma_rigid, so F = (1 - gamma_SC/gamma_rigid) X_H.

Both edges are checked against a dense (A, w) grid before anything is written.
An observation enters as the curve of intrinsic pairs its compatible set
contains, (X_H(c), F(c)) under the base map with a shared reference.  Masking
only inflates, so every point of that curve has X_H >= ln K_HT; where that
exceeds the envelope threshold, comparing the curve with the region is the same
as comparing its lower endpoint with F0.

Panels B and C are the promoting-mode readout of completion.py: F(T) over
250-350 K at A = 5, w0 = 0.5 for the four wavenumbers the text quotes, and the
span over that interval against wavenumber, with the peak locations over the
text's sweep of A and w0.

Run from analysis/:   python export_fig_envelope.py
"""
from __future__ import annotations

import sys

import numpy as np

import masses as M
from corpus import series_record
import qtunnel as Q
from completion import CM_TO_K, offset_vs_T, sensitivity_peak
from partial_id import envelope_B, scale_threshold

G = M.gamma_sc("C"); F0 = M.offset_F0("C")
OUT = "../figures/tikz/data/"
X_MAX = 3.0                         # plotted range of X_H
SLOPE_RIGID = 1.0 - G / M.gamma_rigid("C")
CASES = [("yadh", 7.13, 1.73), ("light", *series_record("light enzyme", 5.0))]
NUS = (50, 200, 400, 800)
SWEEP_A = (1.0, 2.0, 5.0, 10.0, 20.0, 50.0)       # audit_numbers.check_promoting_mode
SWEEP_W0 = (0.1, 0.25, 0.5, 1.0, 2.0)


def check_region():
    """Every sampled (A, w) lies between the rigid edge and B; both edges are reached."""
    # the envelope is a supremum approached as 1/w (verify_derivation section 12),
    # so the grid must reach very wide gating to come close to it at large X_H
    A = np.geomspace(1e-3, 1e9, 1500)
    w = np.geomspace(1e-5, 1e9, 1500)
    AA, WW = np.meshgrid(A, w)
    lh, ld = Q.log_kie_intrinsic(AA, WW)
    F = lh - G * ld
    inside = (lh > 0) & (lh <= X_MAX)
    lh, F = lh[inside], F[inside]
    below_rigid = float(np.max(SLOPE_RIGID * lh - F))
    above_env = float(np.max(F - envelope_B(np.exp(lh))))
    # how closely the grid approaches each edge, binned in X_H
    bins = np.linspace(0.05, X_MAX, 30)
    gap_top = gap_bot = 0.0
    for lo, hi in zip(bins[:-1], bins[1:]):
        s = (lh >= lo) & (lh < hi)
        gap_top = max(gap_top, float(np.min(envelope_B(np.exp(lh[s])) - F[s])))
        gap_bot = max(gap_bot, float(np.min(F[s] - SLOPE_RIGID * lh[s])))
    return below_rigid, above_env, gap_top, gap_bot


def main():
    bad = 0
    b_r, a_e, g_t, g_b = check_region()
    print(f"  grid: below rigid edge by at most {b_r:.1e}, above envelope by at most "
          f"{a_e:.1e}; approaches top within {g_t:.1e}, bottom within {g_b:.1e}")
    if b_r > 1e-12 or a_e > 1e-12 or g_t > 1e-4 or g_b > 1e-4:
        print("  FAIL envelope: the drawn edges do not bound the sampled model")
        bad += 1

    X = np.r_[np.linspace(1e-4, 0.5, 400), np.linspace(0.5, X_MAX, 200)[1:]]
    with open(OUT + "fe_env.dat", "w") as f:
        f.write("xh F rigid\n")
        for x in X:
            f.write(f"{x:.6f} {float(envelope_B(np.exp(x))):.7f} {SLOPE_RIGID * x:.7f}\n")

    # observations as curves of compatible intrinsic pairs
    for name, kh, kd in CASES:
        a = kh - 1.0
        c = a * (1.0 + np.geomspace(1e-4, 1e5, 1500))[::-1]
        xh = kh * c / (1.0 + c - kh)
        xd = kd * c / (1.0 + c - kd)
        XH, FF = np.log(xh), np.log(xh) - G * np.log(xd)
        if np.any(XH < np.log(kh) - 1e-12):
            print(f"  FAIL envelope: {name} curve has X_H below ln K_HT")
            bad += 1
        keep = (XH <= X_MAX) & (FF <= 0.55)
        with open(OUT + f"fe_obs_{name}.dat", "w") as f:
            f.write("xh F\n")
            f.write(f"{np.log(kh):.6f} {np.log(kh) - G * np.log(kd):.6f}\n")
            for u, v in zip(XH[keep], FF[keep]):
                f.write(f"{u:.6f} {v:.6f}\n")
        if np.log(kh) < np.log(scale_threshold()):
            print(f"  FAIL envelope: {name} lies below the envelope threshold")
            bad += 1
    print(f"  envelope threshold X_H = {np.log(scale_threshold()):.4f} "
          f"(K = {scale_threshold():.4f}); rigid slope {SLOPE_RIGID:.4f}")

    # promoting mode: F(T) relative to 250 K, in units of |F0|
    T = np.linspace(250.0, 350.0, 201)
    spans = {}
    for nu in NUS:
        F = offset_vs_T(5.0, 0.5, nu, T)
        spans[nu] = abs(F[-1] - F[0]) / abs(F0)
        with open(OUT + f"fe_T_{nu}.dat", "w") as f:
            f.write("T dF\n")
            for t, v in zip(T, (F - F[0]) / abs(F0)):
                f.write(f"{t:.2f} {v:.6f}\n")
    print("  spans over 250-350 K at A=5, w0=0.5: " +
          ", ".join(f"{nu}: {s:.2f}|F0|" for nu, s in spans.items()))
    for nu, want in zip(NUS, ("0.07", "0.85", "1.00", "0.35")):
        if f"{spans[nu]:.2f}" != want:
            print(f"  FAIL envelope: span at {nu} cm-1 is {spans[nu]:.2f}, text {want}")
            bad += 1

    nus = np.arange(30.0, 1205.0, 5.0)
    span = np.array([abs(offset_vs_T(5.0, 0.5, nu, [350.0])[0]
                         - offset_vs_T(5.0, 0.5, nu, [250.0])[0]) for nu in nus]) / abs(F0)
    with open(OUT + "fe_span.dat", "w") as f:
        f.write("nu span\n")
        for u, v in zip(nus, span):
            f.write(f"{u:.0f} {v:.6f}\n")
    with open(OUT + "fe_span_pts.dat", "w") as f:
        f.write("nu span\n")
        for nu in NUS:
            f.write(f"{nu} {spans[nu]:.6f}\n")
    peaks = [sensitivity_peak(A, w0)[0] for A in SWEEP_A for w0 in SWEEP_W0]
    lo, hi = min(peaks), max(peaks)
    cross = 2 * 300 / CM_TO_K
    print(f"  sensitivity peak over the sweep: {lo:.0f} to {hi:.0f} cm-1; "
          f"crossover {cross:.0f} cm-1; A=5, w0=0.5 peak at "
          f"{nus[int(np.argmax(span))]:.0f} cm-1")
    if (f"{lo:.0f}", f"{hi:.0f}", f"{cross:.0f}") != ("140", "470", "417"):
        print("  FAIL envelope: peak range or crossover disagrees with the text")
        bad += 1
    with open(OUT + "fe_levels.tex", "w") as f:
        f.write("% Generated by analysis/export_fig_envelope.py; do not edit.\n")
        f.write(f"\\def\\feFzero{{{F0:.6f}}}\n")
        f.write(f"\\def\\feXthr{{{np.log(scale_threshold()):.6f}}}\n")
        f.write(f"\\def\\fePeakLo{{{lo:.0f}}}\n\\def\\fePeakHi{{{hi:.0f}}}\n")
        f.write(f"\\def\\feCross{{{cross:.1f}}}\n")
        for name, kh, kd in CASES:          # the observations, as open circles
            tag = name.capitalize()
            f.write(f"\\def\\fe{tag}X{{{np.log(kh):.6f}}}\n")
            f.write(f"\\def\\fe{tag}F{{{np.log(kh) - G * np.log(kd):.6f}}}\n")
    print("envelope figure data:", "OK" if not bad else f"{bad} FAILURES")
    return bad


if __name__ == "__main__":
    sys.exit(1 if main() else 0)
