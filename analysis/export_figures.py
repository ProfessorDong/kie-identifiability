"""Export every plotted quantity as a pgfplots table.

Same discipline as the classical paper: figures are drawn from these tables, so
the numbers in the figures and the numbers in the text are one computation.
"""
from __future__ import annotations

import pathlib

import numpy as np
import pandas as pd

import qtunnel as Q
from discriminate import predict_sc, predict_qm, series, chi2_of, MODELS, fit

OUT = pathlib.Path("../figures/tikz/data")
OUT.mkdir(parents=True, exist_ok=True)

GAMMA_SC_ATOMIC = Q.gamma_rigid() * Q.MU_D / Q.MU_H     # = 3.2572525595
GAMMA_SC_REDUCED = 3.34


def w(name, df):
    df.to_csv(OUT / name, sep=" ", index=False, float_format="%.10g")
    print(f"  {name:26s} {len(df):5d} rows")


# ------------------------------------------------- fig 1: exponent landscape
def fig_landscape():
    ws = np.logspace(-4, 6, 400)
    for i, A in enumerate([0.5, 5.0, 50.0, 5e3, 5e6]):
        w(f"f1_A{i}.dat", pd.DataFrame({
            "w": ws, "g": [Q.gamma_tunneling(A, x) for x in ws]}))
    # the narrow-barrier Mobius envelope, which carries the supremum
    w("f1_envelope.dat", pd.DataFrame({
        "w": ws,
        "g": Q.gamma_rigid() * (1 + ws * Q.MU_D) / (1 + ws * Q.MU_H)}))
    w("f1_levels.dat", pd.DataFrame({
        "name": ["rigid", "gated", "supremum", "reduced"],
        "g": [Q.gamma_rigid(), Q.gamma_fully_gated(),
              GAMMA_SC_ATOMIC, GAMMA_SC_REDUCED]}))


# ------------------------------------------------------ fig 2: gating caps
def fig_caps():
    eps = np.logspace(np.log10(0.002), 0.0, 80)
    grid = np.logspace(-3, 14, 4000)
    caps = [max(Q.gamma_tunneling(A, 2 * A * e ** 2) for A in grid) for e in eps]
    w("f2_caps.dat", pd.DataFrame({"eps": eps, "cap": caps}))


# --------------------------------------------- fig 3: the two fits to data
def fig_fits():
    rows_obs, rows_fit, rows_res = [], [], []
    for tag, label in (("hyd", "hydride"), ("pro", "proton")):
        T, kht, sh, kdt, sd = series(tag)
        fits = {}
        for mname, (pred, npar, init) in MODELS.items():
            r = fit(pred, init, T, kht, sh, kdt, sd, n_restart=400)
            fits[mname] = (pred, r)
        Tf = np.linspace(T.min() - 2, T.max() + 2, 160)
        for mname, (pred, r) in fits.items():
            key = "sc" if mname.startswith("SC") else "qm"
            fh, fd = pred(r.x, Tf)
            w(f"f3_{tag}_{key}_fit.dat",
              pd.DataFrame({"T": Tf - 273.15, "kht": fh, "kdt": fd}))
            ph, pd_ = pred(r.x, T)
            for j in range(len(T)):
                rows_res.append(dict(series=label, model=key, T=T[j] - 273.15,
                                     zh=(ph[j] - kht[j]) / sh[j],
                                     zd=(pd_[j] - kdt[j]) / sd[j]))
        w(f"f3_{tag}_obs.dat", pd.DataFrame({
            "T": T - 273.15, "kht": kht, "sh": sh, "kdt": kdt, "sd": sd,
            "gobs": np.log(kht) / np.log(kdt)}))
    res = pd.DataFrame(rows_res)
    for tag in ("hydride", "proton"):
        for key in ("sc", "qm"):
            s = res[(res.series == tag) & (res.model == key)]
            w(f"f3_res_{tag[:3]}_{key}.dat",
              s[["T", "zh", "zd"]].reset_index(drop=True))


# ---------------------------------------------------- fig 4: design power
def fig_power():
    p = pathlib.Path("../results/design_power.csv")
    if not p.exists():
        print("  design_power.csv not present yet -- skipping fig 4")
        return
    d = pd.read_csv(p)
    if "gamma_tun" not in d.columns:
        print("  design_power.csv is from the superseded run -- skipping fig 4")
        return
    order = {"as measured": 1.0, "2x better": 0.5, "4x better": 0.25}
    d["relprec"] = d.precision.map(order)
    for g, sub in d.groupby("gamma_tun"):
        for nT, s2 in sub.groupby("n_temperatures"):
            s2 = s2.sort_values("relprec", ascending=False)
            w(f"f4_g{str(g).replace('.','')}_n{nT}.dat",
              s2[["relprec", "power"]].reset_index(drop=True))
    d.to_csv(OUT / "f4_all.dat", sep=" ", index=False)


if __name__ == "__main__":
    print("pgfplots tables:")
    fig_landscape()
    fig_caps()
    fig_fits()
    fig_power()
