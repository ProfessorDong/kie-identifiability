"""What experiment would separate gated tunneling from semiclassical kinetics?

Theory: the two classes share a boundary (sup gamma_TUN = gamma_SC), so they
separate only where the gating amplitude is bounded.  Here is the operational
version.

The truth must be SPECIFIED BY PHYSICS, not obtained by fitting.  Fitting the
tunneling model to the data drives it into the degenerate corner
(sigma/R0 -> 0, gamma_TUN -> gamma_SC) where it mimics semiclassics exactly --
that is the theorem in action, and simulating from such a fit merely asks
whether a semiclassical truth looks semiclassical.  Instead we fix the intrinsic
exponent at physically distinct values spanning the tunneling range, match the
intrinsic effect magnitude to the published scale, and let the commitment adjust
so the OBSERVED effects reproduce the measured series.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import least_squares, brentq
from scipy import stats

import qtunnel as Q
from discriminate import predict_sc, predict_qm, series, chi2_of, MODELS, R_KCAL

RNG = np.random.default_rng(31415)
OUT = []


def say(s=""):
    print(s)
    OUT.append(s)


# ------------------------------------------------------- specify a truth
def A_for(w, target):
    """A giving ln K_HT^int = target at gating w, or NaN if unreachable.

    At large w the exponential term dies as A/w^2 while the prefactor saturates
    at -1/2 ln(mu_H/mu_T) = 0.2747, so intrinsic effects above K_HT ~ 1.32
    simply cannot be produced by a strongly gated barrier.  That ceiling is
    physics, not a solver failure.
    """
    f = lambda A: Q.log_kie_intrinsic(A, w)[0] - target
    lo, hi = 1e-9, 1e9
    if not np.isfinite(f(lo)) or not np.isfinite(f(hi)) or f(lo) * f(hi) > 0:
        return np.nan
    return brentq(f, lo, hi, xtol=1e-12, rtol=1e-14)


def solve_Aw(gamma_target, lnKHT_target):
    """Find (A, w) giving the requested intrinsic exponent and H/T magnitude."""
    ws = np.logspace(-8, 4, 2000)
    gs, As = [], []
    for w in ws:
        A = A_for(w, lnKHT_target)
        if np.isfinite(A):
            As.append(A)
            gs.append(Q.gamma_tunneling(A, w))
        else:
            As.append(np.nan)
            gs.append(np.nan)
    gs, As = np.array(gs), np.array(As)
    ok = np.isfinite(gs)
    if not ok.any():
        return None
    wv, gv = ws[ok], gs[ok]
    # find a sign change of gamma(w) - target along the reachable branch
    d = gv - gamma_target
    idx = np.where(np.sign(d[:-1]) * np.sign(d[1:]) < 0)[0]
    if len(idx) == 0:
        return None
    i = idx[0]
    w = brentq(lambda x: Q.gamma_tunneling(A_for(x, lnKHT_target), x) - gamma_target,
               wv[i], wv[i + 1], xtol=1e-14, rtol=1e-14)
    return A_for(w, lnKHT_target), w


def reachable_gamma_range(lnKHT_target):
    """Range of intrinsic exponents attainable at a given H/T magnitude."""
    ws = np.logspace(-8, 4, 3000)
    gs = []
    for w in ws:
        A = A_for(w, lnKHT_target)
        if np.isfinite(A):
            gs.append(Q.gamma_tunneling(A, w))
    return (min(gs), max(gs)) if gs else (np.nan, np.nan)


def make_truth(gamma_target, lnKHT_int, T_obs, kht, sh, kdt, sd, theta_K=979.0):
    """(A, w0, theta, c0, Ec) reproducing the measured series at fixed gamma."""
    A, w_ref = solve_Aw(gamma_target, lnKHT_int)
    if A is None or not np.isfinite(A):
        return None
    # w(T) = w0 coth(theta/2T); pin w(298.15) = w_ref
    w0 = w_ref * np.tanh(theta_K / (2 * 298.15))

    def resid(cp):
        lnc0, Ec = cp
        w = Q.gating_w(T_obs, w0, theta_K)
        lh, ld = Q.log_kie_intrinsic(A, w)
        c = np.exp(lnc0 - Ec / (R_KCAL * T_obs))
        ph = Q.observed_from_intrinsic(np.exp(lh), c)
        pd_ = Q.observed_from_intrinsic(np.exp(ld), c)
        out = np.concatenate([(ph - kht) / sh, (pd_ - kdt) / sd])
        return np.where(np.isfinite(out), out, 1e6)

    best = None
    for _ in range(400):
        r = least_squares(resid, [RNG.uniform(-3, 3), RNG.uniform(-6, 6)],
                          max_nfev=4000)
        if best is None or r.cost < best.cost:
            best = r
    th = np.array([np.log(A), np.log(w0), np.log(theta_K), best.x[0], best.x[1]])
    return th, A, w_ref, chi2_of(best)


def fit_sc(T, kht, sh, kdt, sd, n_restart=40):
    _, _, init = MODELS["SC  semiclassical + commitment"]

    def resid(th):
        with np.errstate(all="ignore"):
            ph, pd_ = predict_sc(th, T)
        out = np.concatenate([(ph - kht) / sh, (pd_ - kdt) / sd])
        return np.where(np.isfinite(out), out, 1e6)

    best = None
    for _ in range(n_restart):
        try:
            r = least_squares(resid, init(RNG), max_nfev=4000)
        except Exception:
            continue
        if best is None or r.cost < best.cost:
            best = r
    return best


def main():
    T, kht, sh, kdt, sd = series("hyd")
    rel_h, rel_d = float(np.mean(sh / kht)), float(np.mean(sd / kdt))

    say("=" * 78)
    say("TRUTHS SPECIFIED BY PHYSICS, SPANNING THE TUNNELING RANGE")
    say("=" * 78)
    say(f"tunneling range of the intrinsic exponent: "
        f"[{Q.gamma_rigid():.3f}, {Q.gamma_rigid()*Q.MU_D/Q.MU_H:.3f})")
    say(f"semiclassical value (reduced mass): 3.340")
    say("")
    say("A large intrinsic effect and a large exponent cannot be had together:")
    say(f"{'K_HT^int':>10s}{'reachable gamma_TUN':>26s}")
    for k in (1.5, 2.0, 3.0, 5.0, 8.0, 11.0, 15.0):
        lo, hi = reachable_gamma_range(np.log(k))
        say(f"{k:10.1f}{f'[{lo:.3f}, {hi:.3f}]':>26s}")
    say("")
    say(f"{'gamma_TUN':>10s}{'A':>10s}{'w(298K)':>10s}{'sigma/R0':>10s}"
        f"{'chi2 to data':>14s}")
    truths = {}
    for g in (2.35, 2.50, 2.70, 2.90):
        out = make_truth(g, np.log(11.0), T, kht, sh, kdt, sd)
        if out is None:
            say(f"{g:10.2f}   unreachable at this intrinsic magnitude")
            continue
        th, A, w_ref, x2 = out
        sor = np.sqrt(w_ref / (2 * A))
        say(f"{g:10.2f}{A:10.3f}{w_ref:10.4f}{sor:10.3f}{x2:14.2f}")
        truths[g] = th
    say("")
    say("chi2 is against the 8 measured hydride observations with 2 free")
    say("commitment parameters (6 d.o.f.); a large value means that intrinsic")
    say("exponent cannot reproduce the data no matter how the commitment is set.")
    say("")

    say("=" * 78)
    say("POWER TO REJECT THE SEMICLASSICAL MODEL WHEN GATED TUNNELING IS TRUE")
    say("=" * 78)
    say(f"measured relative precision: H/T {100*rel_h:.2f}%, D/T {100*rel_d:.2f}%")
    say("temperatures evenly spaced over 5-35 C; alpha = 0.05")
    say("")
    say(f"{'gamma_TUN':>10s}{'n_T':>5s}{'precision':>13s}{'power':>9s}")
    rows, NSIM = [], 300
    for g, th_true in truths.items():
        for nT in (4, 8, 12):
            for scale, tag in ((1.0, "as measured"), (0.5, "2x better"),
                               (0.25, "4x better")):
                Tg = np.linspace(5, 35, nT) + 273.15
                ph, pd_ = predict_qm(th_true, Tg)
                s1, s2 = scale * rel_h * ph, scale * rel_d * pd_
                rej = 0
                for _ in range(NSIM):
                    y1 = ph + RNG.normal(0, s1)
                    y2 = pd_ + RNG.normal(0, s2)
                    r = fit_sc(Tg, y1, s1, y2, s2, n_restart=25)
                    if r is None:
                        continue
                    if 1 - stats.chi2.cdf(chi2_of(r), 2 * nT - 4) < 0.05:
                        rej += 1
                pw = rej / NSIM
                say(f"{g:10.2f}{nT:5d}{tag:>13s}{pw:9.1%}")
                rows.append(dict(gamma_tun=g, n_temperatures=nT,
                                 precision=tag, power=pw))
    pd.DataFrame(rows).to_csv("../results/design_power.csv", index=False)
    with open("../results/design_report.txt", "w") as f:
        f.write("\n".join(OUT) + "\n")
    print("\n[written] ../results/design_power.csv, design_report.txt")


if __name__ == "__main__":
    main()
