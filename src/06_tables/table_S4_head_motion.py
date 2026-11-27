#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Supplementary Table (head motion) - heritability of mean framewise displacement (FD) and FD-adjusted
cross-twin correlations of the 12 brain-state metrics. Cited in Results Section 3.5.

Rows
  Head motion       mean FD itself: ICC(1) in MZ and DZ pairs, Falconer h2, one-sided Fisher z test of r_MZ > r_DZ
  FO / DT / AR      for each of the four states:
                      rho(FD)      Spearman correlation between mean FD and the metric across all 240 participants,
                                   with a 95% CI from a bootstrap of twin pairs (pairs kept intact)
                      Unadjusted   r_MZ, r_DZ, h2, exactly as in the cross-twin correlation table
                      FD-adjusted  the same statistics after regressing the metric on mean FD
                      FD + age + sex   the same statistics after regressing the metric on mean FD, age and sex
                    Residuals come from one ordinary least-squares fit across all 240 participants.

Methods are identical to "Supplementary Table S5.py" (ICC(1) with F-based CIs, h2 = 2(r_MZ - r_DZ) with
95% percentile CIs from 10,000 bootstrap resamples of twin pairs within zygosity, one-sided Fisher r-to-z
with var(z) = 1/(n - 3/2)); the functions are imported from that script and the same bootstrap pair
indices are reused, so the unadjusted columns reproduce that table exactly.

INPUT   inputs/2026-11-13-master-sheet.xlsx, sheet "Master" (adds MeanFD_from_HCP_regressors, Age_in_Yrs, Gender)
OUTPUT  outputs/Supp-Table-FD-heritability.docx, outputs/Supp-Table-FD-heritability-numbers.xlsx
Run     python "Supplementary Table FD.py"   (from the Figure_Codes folder)
"""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("s5", HERE / "table_S3.py")
s5 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s5)

FD, AGE, SEX = "MeanFD_from_HCP_regressors", "Age_in_Yrs", "Gender"
ADJ = {"Unadjusted": [], "FD-adjusted": [FD], "FD + age + sex": [FD, AGE, SEX]}
METRIC_NAMES = {"FO": "Fractional occupancy (FO)", "DT": "Dwell time (DT, s)",
                "AR": "Appearance rate (AR, occurrences min⁻¹)"}
OUT_DOCX = Path(s5.OUTPUT_FOLDER) / "Supp-Table-FD-heritability.docx"
OUT_XLSX = Path(s5.OUTPUT_FOLDER) / "Supp-Table-FD-heritability-numbers.xlsx"
s5.TITLE = ("**Supplementary Table S4.** Heritability of head motion (mean framewise displacement, FD) and "
            "cross-twin correlations of the brain-state metrics before and after adjustment for FD, age and sex "
            "in monozygotic (MZ) and dizygotic (DZ) twin pairs")


def residualise(df, col, covs):
    if not covs:
        return df[col].astype(float)
    X = np.column_stack([np.ones(len(df))] + [df[c].astype(float) for c in covs])
    y = df[col].astype(float).to_numpy()
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return pd.Series(y - X @ beta + y.mean(), index=df.index)   # mean added back; ICC is unaffected


def twin_stats(df, col, idx_mz, idx_dz):
    w = s5.make_pairs(df, None, None, col=col)
    mz, dz = w[w.zyg == "MZ"], w[w.zyg == "DZ"]
    T1m, T2m, T1d, T2d = (x.to_numpy(float) for x in (mz.t1, mz.t2, dz.t1, dz.t2))
    cm, cd = s5.icc_inference(T1m, T2m), s5.icc_inference(T1d, T2d)
    bh = 2 * (s5.icc1(T1m[idx_mz], T2m[idx_mz])[0] - s5.icc1(T1d[idx_dz], T2d[idx_dz])[0])
    lo, hi = np.quantile(bh, [s5.ALPHA / 2, 1 - s5.ALPHA / 2])
    zi = (np.arctanh(cm["r"]) - np.arctanh(cd["r"])) / np.sqrt(1 / (cm["n"] - 1.5) + 1 / (cd["n"] - 1.5))
    return dict(r_MZ=cm["r"], r_MZ_lo=cm["ci_lo"], r_MZ_hi=cm["ci_hi"], r_MZ_p=cm["p"],
                r_DZ=cd["r"], r_DZ_lo=cd["ci_lo"], r_DZ_hi=cd["ci_hi"], r_DZ_p=cd["p"],
                h2=2 * (cm["r"] - cd["r"]), h2_lo=lo, h2_hi=hi, p_MZ_gt_DZ=stats.norm.sf(zi),
                n_MZ=cm["n"], n_DZ=cd["n"])


def rho_fd(df, col, rng, B=s5.N_BOOT):
    """Spearman rho(FD, metric) over participants; 95% CI from resampling whole twin pairs."""
    rho = stats.spearmanr(df[FD], df[col]).statistic
    pairs = df["TwinPairID"].to_numpy()
    upairs = np.unique(pairs)
    rows_of = {p: np.flatnonzero(pairs == p) for p in upairs}
    x, y = df[FD].to_numpy(float), df[col].to_numpy(float)
    boots = np.empty(B)
    for b in range(B):
        idx = np.concatenate([rows_of[p] for p in rng.choice(upairs, len(upairs))])
        boots[b] = stats.spearmanr(x[idx], y[idx]).statistic
    lo, hi = np.quantile(boots, [0.025, 0.975])
    return rho, lo, hi


def build_numbers(df):
    idx_mz, idx_dz = s5.bootstrap_pair_indices(df)
    rng = np.random.default_rng(s5.SEED + 21)
    rows = [dict(metric="Head motion", state="Mean FD", adjustment="Unadjusted",
                 **twin_stats(df, FD, idx_mz, idx_dz))]
    for m in s5.METRICS:
        for s in s5.STATE_COLS:
            col = f"{m}_{s}"
            rho, rlo, rhi = rho_fd(df, col, rng)
            for adj, covs in ADJ.items():
                tmp = df.copy()
                tmp["_y"] = residualise(df, col, covs)
                rows.append(dict(metric=m, state=s, adjustment=adj, rho_FD=rho, rho_FD_lo=rlo, rho_FD_hi=rhi,
                                 **twin_stats(tmp, "_y", idx_mz, idx_dz)))
    return pd.DataFrame(rows)


def p_cell(p):
    return s5.fmt_p_word(p)


def r2(x):
    return s5.minus(f"{0.0 if round(x, 2) == 0 else x:.2f}")


def write_docx(N, df):
    key = N.set_index(["metric", "state", "adjustment"])
    fd = key.loc[("Head motion", "Mean FD", "Unadjusted")]
    mz_fd = df.loc[df.Zygosity == "MZ", FD]; dz_fd = df.loc[df.Zygosity == "DZ", FD]
    h2 = lambda r: (s5.fmt_ci(r.h2, r.h2_lo, r.h2_hi, 2) + ("^{c}" if not 0 <= r.h2 <= 1 else ""))
    rows = [("Head motion", [["Mean FD", "—", r2(fd.r_MZ), r2(fd.r_DZ), h2(fd), p_cell(fd.p_MZ_gt_DZ),
                              "", "", "", "", "", "", "", ""]])]
    for m in s5.METRICS:
        g = []
        for s in s5.STATE_COLS:
            u, f, a = (key.loc[(m, s, k)] for k in ADJ)
            rho_txt = s5.fmt_ci(u.rho_FD, u.rho_FD_lo, u.rho_FD_hi, 2)
            if u.rho_FD_lo > 0 or u.rho_FD_hi < 0:
                rho_txt = f"**{rho_txt}**"
            g.append([s5.STATE_LABEL[s], rho_txt,
                      r2(u.r_MZ), r2(u.r_DZ), h2(u), p_cell(u.p_MZ_gt_DZ),
                      r2(f.r_MZ), r2(f.r_DZ), h2(f), p_cell(f.p_MZ_gt_DZ),
                      r2(a.r_MZ), r2(a.r_DZ), h2(a), p_cell(a.p_MZ_gt_DZ)])
        rows.append((METRIC_NAMES[m], g))
    doc = s5.new_document()
    s5.build(doc, widths=[1950, 1750, 700, 700, 1500, 780, 700, 700, 1500, 780, 700, 700, 1500, 880],
             head1=[("State", 1, True), ("*ρ* (FD)^{a} [95% CI]", 1, True), ("Unadjusted", 4, False),
                    ("FD-adjusted^{b}", 4, False), ("FD + age + sex adjusted^{b}", 4, False)],
             head2=["*r*_{MZ}", "*r*_{DZ}", "h^{2} [95% CI]", "*p*^{d}",
                    "*r*_{MZ}", "*r*_{DZ}", "h^{2} [95% CI]", "*p*^{d}",
                    "*r*_{MZ}", "*r*_{DZ}", "h^{2} [95% CI]", "*p*^{d}"],
             rows=rows,
             notes=[
                 "*Note.* 80 MZ and 40 DZ complete twin pairs (240 participants). Mean FD per participant from the "
                 f"HCP movement regressors (MZ {mz_fd.mean():.3f} ± {mz_fd.std():.3f} mm; DZ {dz_fd.mean():.3f} ± "
                 f"{dz_fd.std():.3f} mm); mean FD *r*_{{MZ}} = {fd.r_MZ:.2f} [{fd.r_MZ_lo:.2f}, {fd.r_MZ_hi:.2f}], "
                 f"*r*_{{DZ}} = {fd.r_DZ:.2f} [{s5.minus(f'{fd.r_DZ_lo:.2f}')}, {fd.r_DZ_hi:.2f}] (exact *F*-based 95% CIs). "
                 "*r*_{MZ} and *r*_{DZ} are one-way random-effects intraclass correlations, ICC(1); "
                 "h^{2} = 2(*r*_{MZ} − *r*_{DZ}) (Falconer’s formula), reported without truncation, with 95% "
                 "percentile CIs from 10,000 bootstrap resamples of complete twin pairs within each zygosity group. "
                 "Unadjusted values are identical to those in Supplementary Table S3. *p*-values are uncorrected; "
                 "bold indicates *p* < 0.05.",
                 "^{a}Spearman correlation between mean FD and the metric across all 240 participants; the 95% CI is "
                 "from 10,000 bootstrap resamples of whole twin pairs, which accounts for the non-independence of co-twins; "
                 "bold indicates a CI that excludes 0.",
                 "^{b}Each metric was regressed on mean FD (FD-adjusted) or on mean FD, age and sex (FD + age + sex "
                 "adjusted) by ordinary least squares across all 240 participants, and the residuals were used to "
                 "recompute the twin statistics.",
                 "^{c}Exceeds the theoretical maximum of 1; Falconer’s estimator is unconstrained, so this reflects "
                 "sampling variability.",
                 "^{d}One-sided Fisher *r*-to-*z* test of *r*_{MZ} > *r*_{DZ} [var(*z*) = 1/(*n* − 3/2)], equivalent "
                 "to a test of h^{2} > 0."])
    doc.save(OUT_DOCX)


def main():
    df = s5.load_master()
    for c in (FD, AGE, SEX):
        assert df[c].notna().all(), f"missing values in {c}"
    N = build_numbers(df)
    Path(s5.OUTPUT_FOLDER).mkdir(parents=True, exist_ok=True)
    N.replace({"State 1": s5.STATE_LABEL["State 1"], "State 2": s5.STATE_LABEL["State 2"],
               "State 3": s5.STATE_LABEL["State 3"], "State 4": s5.STATE_LABEL["State 4"]}).to_excel(
        OUT_XLSX, index=False, sheet_name="FD_table_numbers")
    write_docx(N, df)
    pd.set_option("display.width", 250)
    print(N[["metric", "state", "adjustment", "rho_FD", "r_MZ", "r_DZ", "h2", "h2_lo", "h2_hi", "p_MZ_gt_DZ"]]
          .round(3).to_string(index=False))
    print("Saved:", OUT_DOCX, OUT_XLSX, sep="\n  ")


if __name__ == "__main__":
    main()
