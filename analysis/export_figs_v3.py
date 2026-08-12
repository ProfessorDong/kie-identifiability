"""Export every plotted quantity for the reconstructed manuscript."""
from __future__ import annotations

import pathlib

import numpy as np
import pandas as pd

import masses as M

OUT = pathlib.Path("../figures/tikz/data")
OUT.mkdir(parents=True, exist_ok=True)
GSC = M.gamma_sc("C")
F0 = M.offset_F0("C")
MUH, MUD, MUT = M.mu_ratios("C")


def w_(name, df):
    df.to_csv(OUT / name, sep=" ", index=False, float_format="%.10g")
    print(f"  {name:24s} {len(df):5d} rows")


def log_kie(A, w):
    f = lambda mu: -0.5 * np.log1p(w * mu) - A * mu / (1.0 + w * mu)
    lt = f(MUT)
    return f(MUH) - lt, f(MUD) - lt


# ---- fig 1a: the two loci in log-KIE space, and the benchmark -------------
x = np.linspace(0.0, 1.15, 60)
w_("g1_sc.dat", pd.DataFrame({"x": x, "y": GSC * x}))
w_("g1_tun.dat", pd.DataFrame({"x": x, "y": GSC * x + F0}))

# attainable set of the gated model, sampled
pts = []
for A in np.logspace(-3, 3, 240):
    for wv in np.logspace(-4, 6, 240):
        lh, ld = log_kie(A, wv)
        if 0.02 < ld < 1.15 and lh < 4.0:
            pts.append((ld, lh))
pts = np.array(pts)
w_("g1_attain.dat", pd.DataFrame({"x": pts[:, 0], "y": pts[:, 1]}))

d = pd.read_csv("../data/trinomial_benchmark.csv")
w_("g1_obs.dat", pd.DataFrame({"x": np.log(d.K_DT), "y": np.log(d.K_HT)}))

# ---- fig 1b: F versus gating, several barrier widths ----------------------
ws = np.logspace(-4, 8, 400)
for i, A in enumerate([0.05, 0.5, 5.0, 5e2, 5e5]):
    lh, ld = log_kie(A, ws)
    w_(f"g1b_A{i}.dat", pd.DataFrame({"w": ws, "F": lh - GSC * ld}))

# ---- fig 2a: F(c), the identifiable set ----------------------------------
def Fc(kht, kdt, c):
    xx = lambda K: K * c / (1.0 + c - K)
    return np.log(xx(kht)) - GSC * np.log(xx(kdt))


d["L_H"] = (d.K_HT - 1) / (d.K_DT - 1)
mono = d[d.L_H >= GSC].iloc[0]
inter = d[d.L_H < GSC].iloc[0]
for tag, r in (("mono", mono), ("inter", inter)):
    c = (r.K_HT - 1) * (1 + np.logspace(-3, 6, 400))
    w_(f"g2_{tag}.dat", pd.DataFrame({"c": c, "F": Fc(r.K_HT, r.K_DT, c)}))
    print(f"     {tag}: K_HT={r.K_HT} K_DT={r.K_DT} L_H={r.L_H:.3f}")

# ---- fig 2b: model-free bounds per series --------------------------------
b = pd.read_csv("../results/modelfree_bounds.csv").sort_values("bound")
b = b.reset_index(drop=True)
b["idx"] = np.arange(len(b))
w_("g2b_bounds.dat", b[["idx", "bound"]])
(OUT / "g2b_labels.tex").write_text("".join(
    f"\\node[font=\\fontsize{{6.4}}{{7.6}}\\selectfont,anchor=east] "
    f"at (axis cs:-1.02,{r.idx:.0f}) {{{r.family} {r.variant}}};\n"
    for _, r in b.iterrows()))
print(f"  g2b_labels.tex          {len(b)} labels")
print(f"\n  gamma_SC={GSC:.5f}  F0={F0:+.6f}")


if __name__ == "__main__":
    pass
