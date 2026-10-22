"""Independent determinations of one condition: are the published errors honest?

Every confidence bound in this analysis reads the published uncertainties as
standard errors.  Whether they are realistic can be checked wherever the same
enzyme form was measured twice, by different experiments, at the same
condition.  The census of 2026-10-20 supplied three such conditions:

  horse liver ADH F93W, benzyl alcohol, 25 C, pH 7
      Bahnson et al. 1993 (standard errors of the mean, 11 to 36 determinations)
      Chin & Klinman 2000 (standard errors of the mean, 38 and 41 determinations)
  E. coli DHFR wild type, 25 C, pH 9
      Sikorski et al. 2004 (at least three experiments; error type not stated)
      Sen et al. 2011 (standard deviations; temperature of these runs not printed)
      Sen et al. 2012 ([carbonyl-14C] reference; error type not stated)
  E. coli DHFR wild type, pH 7.0, 5 and 45 C
      Francis et al. 2013 (main-text Table 1)
      Liu et al. 2014 (SI Table S8; same laboratory, possibly the same data)

For each pair the difference of each effect, and of the observed offset
F_obs = ln K_HT - gamma ln K_DT, is divided by its combined standard error (the
two determinations are independent, so their errors add in quadrature; within a
determination the H/T-D/T correlation is unknown and the offset's error is
bracketed by rho in {-1, 0}).  If the errors are standard errors of independent
determinations, these z values are standard normal; their root mean square
estimates the factor by which the stated errors understate the between-study
scatter.  Quoted values (a later paper reprinting an earlier one) are not
determinations and are excluded.

Run from analysis/:   python repeat_determinations.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import corpus
import masses as M

G = M.gamma_sc("C")


def _series(variant, T):
    s = corpus.series()
    r = s[(s.variant == variant) & (s.T_C == T)].iloc[0]
    return (r.K_HT, r.K_HT_se, r.K_DT, r.K_DT_se)


def _single(fname, variant):
    d = pd.read_csv(corpus.DATA / fname)
    r = d[d.variant == variant].iloc[0]
    return (r.K_HT, r.K_HT_se, r.K_DT, r.K_DT_se)


# Francis et al. 2013, J Biol Chem 288, 35961, Table 1 (pH 7.0 rows; not in the
# corpus, which admits the fuller Liu et al. 2014 series of that condition)
FRANCIS2013 = {5.0: (1.64, 0.02, 1.23, 0.06), 45.0: (2.68, 0.01, 1.43, 0.01)}


def pairs():
    """(condition, label A, record A, label B, record B, note)."""
    f93_b = _single("ladh_adh_primary.csv", "F93W")
    f93_c = _single("ladh_chin2000.csv", "F93W (Chin 2000)")
    sik = _series("light enzyme", 25.0)
    sen11 = _single("dhfr_sen.csv", "WT (Sen 2011)")
    sen12 = _single("dhfr_sen.csv", "WT (Sen 2012)")
    out = [("LADH F93W, 25 C", "Bahnson 1993", f93_b, "Chin 2000", f93_c, "independent"),
           ("ecDHFR WT, 25 C, pH 9", "Sikorski 2004", sik, "Sen 2011", sen11, "independent"),
           ("ecDHFR WT, 25 C, pH 9", "Sikorski 2004", sik, "Sen 2012", sen12, "independent"),
           ("ecDHFR WT, 25 C, pH 9", "Sen 2011", sen11, "Sen 2012", sen12, "independent")]
    for T, fr in FRANCIS2013.items():
        out.append((f"ecDHFR WT, pH 7.0, {T:g} C", "Francis 2013", fr, "Liu 2014",
                    _series("WT (pH 7.0)", T), "same laboratory; possibly shared data"))
    return out


def z_scores(a, b, rho=0.0):
    """z for H/T, D/T and F_obs between two independent determinations."""
    (h1, sh1, d1, sd1), (h2, sh2, d2, sd2) = a, b
    zh = (h1 - h2) / np.hypot(sh1, sh2)
    zd = (d1 - d2) / np.hypot(sd1, sd2)

    def fvar(h, sh, d, sd):
        x, y = sh / h, sd / d
        return x * x + (G * y) ** 2 - 2 * rho * G * x * y
    f1 = np.log(h1) - G * np.log(d1)
    f2 = np.log(h2) - G * np.log(d2)
    zf = (f1 - f2) / np.sqrt(fvar(h1, sh1, d1, sd1) + fvar(h2, sh2, d2, sd2))
    return zh, zd, zf, f1 - f2


def summary():
    rows = []
    for cond, la, a, lb, b, note in pairs():
        zh, zd, zf0, df = z_scores(a, b, 0.0)
        _, _, zf1, _ = z_scores(a, b, -1.0)
        rows.append(dict(condition=cond, a=la, b=lb, z_HT=zh, z_DT=zd, dF=df,
                         z_F_rho0=zf0, z_F_rhom1=zf1, note=note))
    return pd.DataFrame(rows)


def rms_factor(d, independent_only=True):
    """Root-mean-square z over the effects (H/T and D/T), the estimated factor by
    which stated errors understate between-study scatter (1 if they do not)."""
    x = d[d.note == "independent"] if independent_only else d
    z = np.concatenate([x.z_HT.values, x.z_DT.values])
    return float(np.sqrt(np.mean(z ** 2))), len(z)


def main():
    d = summary()
    pd.set_option("display.width", 160)
    print(d.to_string(index=False, float_format=lambda v: f"{v:+.3f}"))
    s, n = rms_factor(d)
    s_all, n_all = rms_factor(d, False)
    print(f"\nRMS z over {n} independent effect comparisons: {s:.2f}")
    print(f"RMS z including the possibly shared pH 7.0 data ({n_all}): {s_all:.2f}")
    worst = d.loc[d[["z_HT", "z_DT"]].abs().max(axis=1).idxmax()]
    print(f"largest |z|: {worst.condition} {worst.a} vs {worst.b}: "
          f"H/T {worst.z_HT:+.2f}, D/T {worst.z_DT:+.2f}")
    d.to_csv("../results/repeat_determinations.csv", index=False)
    with open("../results/sm_table_repeat.tex", "w") as fh:
        fh.write(table(d))
    print("[written] ../results/repeat_determinations.csv, sm_table_repeat.tex")


def table(d):
    """Body of the supplement's table of repeat determinations (tab:repeat)."""
    lines = ["% Generated by analysis/repeat_determinations.py -- do not edit; edit the generator."]
    for r in d.itertuples():
        cond = r.condition.replace(" C", "\\,$^{\\circ}$C")
        mark = "" if r.note == "independent" else "$^{a}$"
        # Sen et al. 2011 print 25 C and pH 9 only for their H/D runs; the
        # match of their H/T and D/T runs to that condition is assumed
        sa = "$^{b}$" if r.a == "Sen 2011" else ""
        sb = "$^{b}$" if r.b == "Sen 2011" else ""
        lines.append(f"{cond} & {r.a}{sa} & {r.b}{sb}{mark} & ${r.z_HT:+.2f}$ & ${r.z_DT:+.2f}$ & "
                     f"${r.dF:+.4f}$ & ${r.z_F_rho0:+.2f}$ & ${r.z_F_rhom1:+.2f}$\\\\")
    return "\n".join(lines + ["\\bottomrule"]) + "\n"


if __name__ == "__main__":
    main()
