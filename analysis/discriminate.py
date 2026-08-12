"""Can a competitive multiple-isotope experiment tell tunneling from semiclassics?

Two model classes are fitted to the same temperature-resolved competitive KIEs,
each composed with the same commitment map:

  SC   semiclassical zero-point-energy scaling, ln K_HT^int = gamma ln K_DT^int
       with gamma fixed at the reduced-mass value; the intrinsic D/T effect is
       Arrhenius in 1/T.                                   4 parameters
  QM   gated vibronic tunneling, ground channel, with the exact thermal
       variance of a harmonic promoting mode.              5 parameters

Both are then masked by an Arrhenius reciprocal commitment c(T) = c0 exp(-Ec/RT),
shared between the H/T and D/T pairs because both are measured against the same
tritium reference.

The question is not which fits better.  It is whether the data can distinguish
them at all -- and if not, what experiment could.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy import stats

import qtunnel as Q

RNG = np.random.default_rng(20260810)
R_KCAL = 1.98720425e-3
GAMMA_SC = 3.34            # reduced-mass convention, as in the classical analysis
OUT = []


def say(s=""):
    print(s)
    OUT.append(s)


def head(t):
    say("=" * 78)
    say(t)
    say("=" * 78)


# ------------------------------------------------------------------ data
def series(name):
    d = pd.read_csv(f"../data/tsase_series/g_{name}_obs.dat", sep=r"\s+")
    return (d["T"].to_numpy() + 273.15, d["kht"].to_numpy(), d["sh"].to_numpy(),
            d["kdt"].to_numpy(), d["sd"].to_numpy())


# ------------------------------------------------------------- model SC
def predict_sc(th, T):
    """th = [ln A_u, E_u, ln c0, E_c]; intrinsic D/T Arrhenius, H/T by scaling."""
    lnAu, Eu, lnc0, Ec = th
    u = 1.0 + np.exp(lnAu - Eu / (R_KCAL * T))          # K_DT^int > 1 enforced
    kht_int = u ** GAMMA_SC
    c = np.exp(lnc0 - Ec / (R_KCAL * T))
    return Q.observed_from_intrinsic(kht_int, c), Q.observed_from_intrinsic(u, c)


# ------------------------------------------------------------- model QM
def predict_qm(th, T):
    """th = [ln A, ln w0, ln theta, ln c0, E_c]; gated vibronic tunneling."""
    lnA, lnw0, lnth, lnc0, Ec = th
    A, w0, theta = np.exp(lnA), np.exp(lnw0), np.exp(lnth)
    w = Q.gating_w(T, w0, theta)
    lh, ld = Q.log_kie_intrinsic(A, w)
    c = np.exp(lnc0 - Ec / (R_KCAL * T))
    return (Q.observed_from_intrinsic(np.exp(lh), c),
            Q.observed_from_intrinsic(np.exp(ld), c))


MODELS = {
    "SC  semiclassical + commitment": (predict_sc, 4,
                                       lambda r: np.array([r.uniform(-2, 4), r.uniform(-4, 6),
                                                           r.uniform(-3, 3), r.uniform(-6, 6)])),
    "QM  gated tunneling + commitment": (predict_qm, 5,
                                         lambda r: np.array([r.uniform(-2, 5), r.uniform(-6, 2),
                                                             r.uniform(4, 8), r.uniform(-3, 3),
                                                             r.uniform(-6, 6)])),
}


def fit(pred, init, T, kht, sh, kdt, sd, n_restart=600, rng=None):
    rng = rng or RNG

    def resid(th):
        with np.errstate(all="ignore"):
            ph, pd_ = pred(th, T)
        out = np.concatenate([(ph - kht) / sh, (pd_ - kdt) / sd])
        return np.where(np.isfinite(out), out, 1e6)

    best = None
    for _ in range(n_restart):
        try:
            r = least_squares(resid, init(rng), max_nfev=4000)
        except Exception:
            continue
        if best is None or r.cost < best.cost:
            best = r
    return best


def chi2_of(res):
    return 2.0 * res.cost


# ------------------------------------------------------- bootstrap calibration
def boot_p(pred, init, th_hat, T, sh, sd, chi2_obs, B=400):
    """Parametric bootstrap from the fitted null, refitted identically."""
    ge = 0
    for _ in range(B):
        ph, pd_ = pred(th_hat, T)
        y1 = ph + RNG.normal(0, sh)
        y2 = pd_ + RNG.normal(0, sd)
        r = fit(pred, init, T, y1, sh, y2, sd, n_restart=40)
        if r is not None and chi2_of(r) >= chi2_obs:
            ge += 1
    return (1 + ge) / (B + 1)


def main():
    head("1. THE TWO MODEL CLASSES ARE NOT SEPARATED BY THE SWAIN-SCHAAD EXPONENT")
    say(f"rigid-limit tunneling exponent      gamma_TUN(w->0)   = {Q.gamma_rigid():.4f}")
    say(f"fully gated, finite barrier width   gamma_TUN(w->inf) = {Q.gamma_fully_gated():.4f}")
    sup = Q.gamma_rigid() * Q.MU_D / Q.MU_H
    say(f"supremum over ALL (A, w)                              = {sup:.6f}")
    say(f"semiclassical exponent, bare atomic masses            = "
        f"{(1 - 3 ** -0.5) / (2 ** -0.5 - 3 ** -0.5):.6f}")
    say("These last two are equal identically, for any isotope triple.")
    say("The semiclassical locus is therefore the BOUNDARY of the tunneling")
    say("model's range, approached only as sigma/R0 -> infinity.")
    say("")
    for frac in (0.10, 0.20, 1 / 3):
        best = max(Q.gamma_tunneling(A, 2 * A * frac ** 2)
                   for A in np.logspace(-2, 4, 600))
        say(f"  physically admissible sigma/R0 <= {frac:.2f}:  "
            f"gamma_TUN <= {best:.4f}")
    say("")

    rows = []
    for name in ("hyd", "pro"):
        T, kht, sh, kdt, sd = series(name)
        label = {"hyd": "hsTSase hydride transfer",
                 "pro": "hsTSase proton abstraction"}[name]
        head(f"2. {label.upper()}")
        say(f"{'model':34s}{'chi2':>9s}{'dof':>5s}{'p_asym':>10s}"
            f"{'p_boot':>9s}{'AIC':>9s}")
        fits = {}
        for mname, (pred, npar, init) in MODELS.items():
            r = fit(pred, init, T, kht, sh, kdt, sd)
            x2 = chi2_of(r)
            dof = 2 * len(T) - npar
            pa = 1 - stats.chi2.cdf(x2, dof) if dof > 0 else np.nan
            pb = boot_p(pred, init, r.x, T, sh, sd, x2)
            aic = x2 + 2 * npar
            fits[mname] = (r, x2, dof, pa, pb, aic)
            say(f"{mname:34s}{x2:9.2f}{dof:5d}{pa:10.4f}{pb:9.3f}{aic:9.1f}")
            rows.append(dict(series=label, model=mname, chi2=x2, dof=dof,
                             p_asym=pa, p_boot=pb, aic=aic))
        say("")
        # what the tunneling fit says about the gating, if it is believable
        rq = fits["QM  gated tunneling + commitment"][0]
        A, w0, theta = np.exp(rq.x[0]), np.exp(rq.x[1]), np.exp(rq.x[2])
        say(f"  gated-tunneling parameters: A = {A:.3f}, w0 = {w0:.4f}, "
            f"theta = {theta:.0f} K")
        say(f"  implied sigma/R0 at 298 K   = "
            f"{np.sqrt(Q.gating_w(298.15, w0, theta) / (2 * A)):.3f}")
        gt = [Q.gamma_tunneling(A, Q.gating_w(t, w0, theta)) for t in T]
        say("  implied intrinsic exponent gamma_TUN(T) = "
            + ", ".join(f"{g:.3f}" for g in gt))
        say(f"  observed gamma_obs(T)                   = "
            + ", ".join(f"{np.log(a) / np.log(b):.3f}" for a, b in zip(kht, kdt)))
        say("")

    pd.DataFrame(rows).to_csv("../results/model_comparison.csv", index=False)

    head("3. DISCRIMINATION: CAN FOUR TEMPERATURES TELL THEM APART?")
    T, kht, sh, kdt, sd = series("hyd")
    pred_sc, _, init_sc = MODELS["SC  semiclassical + commitment"]
    pred_qm, _, init_qm = MODELS["QM  gated tunneling + commitment"]
    r_sc = fit(pred_sc, init_sc, T, kht, sh, kdt, sd)
    say("Simulating from the fitted SEMICLASSICAL model and asking how often the")
    say("tunneling model is rejected, and vice versa.  If neither is rejected,")
    say("the experiment does not discriminate, whatever it reports.")
    say("")
    n_sim, rej_qm, rej_sc = 200, 0, 0
    for _ in range(n_sim):
        ph, pd_ = pred_sc(r_sc.x, T)
        y1, y2 = ph + RNG.normal(0, sh), pd_ + RNG.normal(0, sd)
        rq = fit(pred_qm, init_qm, T, y1, sh, y2, sd, n_restart=40)
        rs = fit(pred_sc, init_sc, T, y1, sh, y2, sd, n_restart=40)
        if rq is not None and 1 - stats.chi2.cdf(chi2_of(rq), 3) < 0.05:
            rej_qm += 1
        if rs is not None and 1 - stats.chi2.cdf(chi2_of(rs), 4) < 0.05:
            rej_sc += 1
    say(f"data generated by SEMICLASSICAL truth, {n_sim} replicates at measured precision:")
    say(f"   tunneling model rejected in    {rej_qm/n_sim:6.1%} of replicates")
    say(f"   semiclassical model rejected in {rej_sc/n_sim:6.1%} of replicates (size check)")
    say("")

    with open("../results/discrimination_report.txt", "w") as f:
        f.write("\n".join(OUT) + "\n")
    print("\n[written] ../results/discrimination_report.txt")


if __name__ == "__main__":
    main()
