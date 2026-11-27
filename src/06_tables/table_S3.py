#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Supplementary Table S3 (corrected analysis) - journal-format Word table

Within-pair distances (means with 95% BCa bootstrap CIs; two-sided Mann-Whitney U), cross-twin ICC(1)
with F-based 95% CIs, the one-sided Fisher r-to-z test of r_MZ > r_DZ (BH-FDR marked) and Falconer
h2 = 2(r_MZ - r_DZ) with 95% pair-bootstrap CIs, for FO, DT and AR in States 1-4.

HOW TO RUN
  1. Edit the two lines under "EDIT THESE TWO LINES" (full path to the master Excel file and an output folder).
  2. Press Run (Spyder, VS Code, IDLE, Jupyter) or type:  python "2026-11-13-Table_S1_corrected.py"

INPUT   the "Master" sheet of 2026-11-13-DMN-C-twin-master.xlsx
        (one row per participant: Subject, TwinPairID, Zygosity, FO_State 1-4, DT_State 1-4, AR_State 1-4)
OUTPUT  <prefix>-table-S1-journal-format-corrected.docx, <prefix>-TableS1-corrected-numbers.xlsx

This script is self-contained: it does not import or need any other script.
Requires Python 3.10+ with numpy, scipy (>= 1.11), pandas, openpyxl, python-docx.
Random seeds are fixed, so every run gives identical numbers.
"""

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from docx import Document
from docx.enum.section import WD_ORIENT
from docx.enum.table import WD_TABLE_ALIGNMENT, WD_CELL_VERTICAL_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Pt, Twips
from docx.table import _Cell

# =============================================================================
# EDIT THESE TWO LINES
# =============================================================================
MASTER_FILE = r"inputs/2026-11-13-master-sheet.xlsx"   # full path, e.g. r"C:\Users\me\Documents\twin\2026-11-13-DMN-C-twin-master.xlsx"
OUTPUT_FOLDER = r"outputs"                            # folder for the results (created if it does not exist)
# =============================================================================

# Fixed settings (do not change unless you intend to change the analysis)
MASTER_SHEET = "Master"
METRICS = ["FO", "DT", "AR"]
STATE_COLS = ["State 1", "State 2", "State 3", "State 4"]
STATE_LABEL = {"State 1": "State 1: Visual Network-A", "State 2": "State 2: DMN",
               "State 3": "State 3: Visual Network-B", "State 4": "State 4: Hippocampal–Temporoparietal"}
ALPHA = 0.05
SEED = 20261110            # master random seed; each procedure uses SEED + a fixed offset
N_BOOT = 10_000            # bootstrap resamples
FILE_PREFIX = "2026-11-13"
FIG_DPI = 600
Z975 = stats.norm.ppf(1 - ALPHA / 2)

# Word table layout
FONT = "Arial"
TABLE_PT, NOTE_PT, TITLE_PT = 8, 8, 10
GROUPS = [("FO", "Fractional occupancy (FO)"), ("DT", "Dwell time (DT, s)"),
          ("AR", "Appearance rate (AR, occurrences min^{−1})")]   # DT/AR units inferred from FO = AR x DT / 60
DIST_DECIMALS = {"FO": 4, "DT": 3, "AR": 3}
# Alcohol and odor variables: (group heading, master-sheet column, row label). Pairs with a missing value are
# left out; their bootstrap uses its own random stream (SEED + 7), so the brain-state numbers are unchanged.
EXTRA = [("Alcohol consumption", "SSAGA_Alc_12_Drinks_Per_Day", "Drinks per drinking occasion"),
         ("Alcohol consumption", "Total_Drinks_7days", "Drinks per week"),
         ("Odor identification", "Odor_AgeAdj", "Age-adjusted score"),
         ("Odor identification", "Odor_Unadj", "Unadjusted score")]
EXTRA_DECIMALS = 2
TITLE = ("**Table S1.** Within-pair distances, cross-twin correlations and Falconer heritability estimates for "
         "brain-state dynamic metrics in monozygotic (MZ) and dizygotic (DZ) twin pairs")


# -----------------------------------------------------------------------------
# Data
# -----------------------------------------------------------------------------

def load_master(path=None):
    """Read the Master sheet and check its structure (stops with a clear message if anything is wrong)."""
    path = Path(path or MASTER_FILE)
    if not path.exists():
        sys.exit(f"Master file not found: {path}\nEdit MASTER_FILE at the top of this script.")
    df = pd.read_excel(path, sheet_name=MASTER_SHEET)
    df.columns = [re.sub(r"^(FO|DT|AR)_State(\d)$", r"\1_State \2", str(c).strip()) for c in df.columns]
    required = ["Subject", "TwinPairID", "Zygosity"] + [f"{m}_{s}" for m in METRICS for s in STATE_COLS]
    missing = [c for c in required if c not in df.columns]
    if missing:
        sys.exit(f"Columns missing from the Master sheet: {missing}")
    df["Subject"] = df["Subject"].astype(str).str.strip().str.replace(r"\.0$", "", regex=True)
    df["TwinPairID"] = df["TwinPairID"].astype(str).str.strip()
    df["Zygosity"] = df["Zygosity"].astype(str).str.strip()
    if df[required].isna().any().any():
        sys.exit("The Master sheet has empty cells; fill or remove them before running.")
    if df["Subject"].duplicated().any():
        sys.exit(f"Duplicated subject IDs: {df.loc[df.Subject.duplicated(), 'Subject'].tolist()}")
    if not set(df["Zygosity"]) <= {"MZ", "DZ"}:
        sys.exit(f"Zygosity must be MZ or DZ; found {sorted(set(df['Zygosity']))}")
    size = df.groupby("TwinPairID").size()
    if (size != 2).any():
        sys.exit(f"Every TwinPairID must have exactly 2 participants: {size[size != 2].to_dict()}")
    if (df.groupby("TwinPairID")["Zygosity"].nunique() > 1).any():
        sys.exit("A twin pair has mixed zygosity labels.")
    n = df.groupby("Zygosity")["TwinPairID"].nunique()
    print(f"Master sheet: {len(df)} participants, {n.get('MZ', 0)} MZ pairs, {n.get('DZ', 0)} DZ pairs")
    return df


def make_pairs(df, metric, state, col=None):
    """One row per twin pair for one metric and state (or for one column, col). Co-twins are matched by
    TwinPairID; twin 1 = the co-twin with the lower subject ID. Rows: MZ pairs then DZ pairs, by pair number."""
    col = col or f"{metric}_{state}"
    d = df[["Subject", "Zygosity", "TwinPairID", col]].copy()
    d["_sid"] = pd.to_numeric(d["Subject"], errors="coerce")
    d = d.sort_values(["TwinPairID", "_sid", "Subject"])
    rows = []
    for pid, g in d.groupby("TwinPairID", sort=False):
        rows.append(dict(pair=pid, zyg=g["Zygosity"].iloc[0], subj_t1=g["Subject"].iloc[0],
                         subj_t2=g["Subject"].iloc[1], t1=g[col].iloc[0], t2=g[col].iloc[1]))
    w = pd.DataFrame(rows)
    w["complete"] = w[["t1", "t2"]].notna().all(axis=1)
    w["pair_num"] = pd.to_numeric(w["pair"].str.extract(r"(\d+)")[0], errors="coerce")
    return w.sort_values(["zyg", "pair_num", "pair"], ascending=[False, True, True]).reset_index(drop=True)


# -----------------------------------------------------------------------------
# Statistics
# -----------------------------------------------------------------------------
def icc1(t1, t2):
    """One-way random-effects ICC(1) for pairs (Shrout & Fleiss ICC(1,1)), k = 2.
    Invariant to which member is labelled twin 1. Works on (..., n) arrays."""
    t1, t2 = np.asarray(t1, float), np.asarray(t2, float)
    n = t1.shape[-1]
    m = (t1 + t2) / 2
    M = m.mean(axis=-1, keepdims=True)
    msb = 2 * ((m - M) ** 2).sum(axis=-1) / (n - 1)
    msw = ((t1 - t2) ** 2 / 2).sum(axis=-1) / n
    return (msb - msw) / (msb + msw), msb / msw


def icc_inference(t1, t2):
    """ICC(1) with exact F test (two-sided, H0: rho = 0) and F-based 95% CI."""
    n = len(t1)
    icc, F = icc1(t1, t2)
    df1, df2 = n - 1, n
    p_upper = stats.f.sf(F, df1, df2)
    p_two = min(1.0, 2 * min(p_upper, stats.f.cdf(F, df1, df2)))
    FL = F / stats.f.ppf(1 - ALPHA / 2, df1, df2)
    FU = F * stats.f.ppf(1 - ALPHA / 2, df2, df1)
    return dict(r=float(icc), F=float(F), df1=df1, df2=df2, p=float(p_two), p_one_sided_gt0=float(p_upper),
                ci_lo=float((FL - 1) / (FL + 1)), ci_hi=float((FU - 1) / (FU + 1)), n=n)


def bh_fdr(p):
    p = np.asarray(p, float)
    q = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    if ok.any():
        q[ok] = stats.false_discovery_control(p[ok], method="bh")
    return q


def hodges_lehmann(a, b):
    """HL shift (a - b) with Moses distribution-free CI (normal approximation)."""
    D = np.sort((a[:, None] - b[None, :]).ravel())
    n1, n2 = len(a), len(b)
    N = n1 * n2
    k = int(np.floor(N / 2 - Z975 * np.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)))
    k = max(k, 1)
    return float(np.median(D)), float(D[k - 1]), float(D[N - k])


def rank_biserial(a, b):
    """r_rb = P(a > b) - P(a < b); positive = first group (MZ) has LARGER distances."""
    U = stats.mannwhitneyu(a, b, method="asymptotic").statistic
    return 2 * U / (len(a) * len(b)) - 1


def rank_biserial_boot(a, b, rng, B):
    ia = rng.integers(0, len(a), (B, len(a)))
    ib = rng.integers(0, len(b), (B, len(b)))
    A, Bm = a[ia], b[ib]
    # P(A>B) - P(A<B) via sign of all pairwise differences, chunked
    out = np.empty(B)
    for s in range(0, B, 500):
        d = A[s:s + 500, :, None] - Bm[s:s + 500, None, :]
        out[s:s + 500] = np.sign(d).mean(axis=(1, 2))
    return out


def bootstrap_mean_ci(d, rng):
    """95% BCa bootstrap CI of a mean (scipy >= 1.15 takes rng=, older versions random_state=)."""
    kw = dict(n_resamples=N_BOOT, method="BCa", confidence_level=1 - ALPHA)
    try:
        return stats.bootstrap((d,), np.mean, rng=rng, **kw).confidence_interval
    except TypeError:
        return stats.bootstrap((d,), np.mean, random_state=rng, **kw).confidence_interval


def distance_analysis(w, rng):
    """Within-pair absolute differences, descriptives and MZ-vs-DZ tests."""
    c = w[w.complete]
    a = np.abs(c[c.zyg == "MZ"].t1 - c[c.zyg == "MZ"].t2).to_numpy()
    b = np.abs(c[c.zyg == "DZ"].t1 - c[c.zyg == "DZ"].t2).to_numpy()
    out = {}
    for z, d in (("MZ", a), ("DZ", b)):
        q1, q3 = np.percentile(d, [25, 75])
        ci = bootstrap_mean_ci(d, rng)
        out.update({f"{z}_n_pairs": len(d), f"{z}_mean": d.mean(), f"{z}_sd": d.std(ddof=1),
                    f"{z}_mean_ci_lo": ci.low, f"{z}_mean_ci_hi": ci.high,
                    f"{z}_median": np.median(d), f"{z}_q1": q1, f"{z}_q3": q3, f"{z}_iqr": q3 - q1,
                    f"{z}_min": d.min(), f"{z}_max": d.max(),
                    f"{z}_skew": stats.skew(d, bias=False), f"{z}_shapiro_W": stats.shapiro(d).statistic,
                    f"{z}_shapiro_p": stats.shapiro(d).pvalue, f"{z}_n_zero": int((d == 0).sum())})
    out["n_tied_values"] = int(len(np.r_[a, b]) - len(np.unique(np.r_[a, b])))
    # Primary: Mann-Whitney U, two-sided, normal approximation with tie & continuity correction
    mw = stats.mannwhitneyu(a, b, alternative="two-sided", method="asymptotic", use_continuity=True)
    n1, n2 = len(a), len(b)
    out["MWU_U_MZ"] = mw.statistic
    out["MWU_z"] = float(stats.norm.isf(mw.pvalue / 2) * np.sign(mw.statistic - n1 * n2 / 2))
    out["MWU_p_two_sided"] = mw.pvalue
    out["MWU_p_one_sided_MZ_less"] = stats.mannwhitneyu(a, b, alternative="less", method="asymptotic",
                                                        use_continuity=True).pvalue
    out["rank_biserial"] = rank_biserial(a, b)
    rb = rank_biserial_boot(a, b, rng, N_BOOT)
    out["rank_biserial_ci_lo"], out["rank_biserial_ci_hi"] = np.quantile(rb, [ALPHA / 2, 1 - ALPHA / 2])
    out["P_MZ_more_similar"] = (1 - out["rank_biserial"]) / 2   # P(d_MZ < d_DZ) + .5 ties
    hl, lo, hi = hodges_lehmann(a, b)
    out.update(HL_shift_MZ_minus_DZ=hl, HL_ci_lo=lo, HL_ci_hi=hi)
    # Sensitivity: Welch t-test
    wt = stats.ttest_ind(a, b, equal_var=False)
    ci = wt.confidence_interval(confidence_level=1 - ALPHA)
    out.update(Welch_t=wt.statistic, Welch_df=wt.df, Welch_p=wt.pvalue,
               mean_diff_MZ_minus_DZ=a.mean() - b.mean(), mean_diff_ci_lo=ci.low, mean_diff_ci_hi=ci.high)
    sp = np.sqrt(((n1 - 1) * a.var(ddof=1) + (n2 - 1) * b.var(ddof=1)) / (n1 + n2 - 2))
    J = 1 - 3 / (4 * (n1 + n2) - 9)
    g = J * (a.mean() - b.mean()) / sp
    se_g = np.sqrt((n1 + n2) / (n1 * n2) + g ** 2 / (2 * (n1 + n2)))
    out.update(hedges_g=g, hedges_g_ci_lo=g - Z975 * se_g, hedges_g_ci_hi=g + Z975 * se_g)
    out["Levene_BF_p"] = stats.levene(a, b, center="median").pvalue
    return out, a, b
def run_distances(df):
    """Within-pair distances for all 12 metric-state combinations (fixed order FO, DT, AR x States 1-4,
    which keeps the random-number stream - and therefore every bootstrap CI - identical across scripts)."""
    rng_dist = np.random.default_rng(SEED + 5)
    pairs, rows = [], []
    for m in METRICS:
        for s in STATE_COLS:
            w = make_pairs(df, m, s)
            w["metric"], w["state"] = m, s
            w["abs_diff"] = (w.t1 - w.t2).abs()
            pairs.append(w)
            d, _, _ = distance_analysis(w, rng_dist)
            rows.append(dict(metric=m, state=s, **d))
    for grp, col, lab in EXTRA:                       # after the brain states: their random stream is unchanged
        w = make_pairs(df, None, None, col=col)
        w["metric"], w["state"] = grp, lab
        w["abs_diff"] = (w.t1 - w.t2).abs()
        pairs.append(w)
        d, _, _ = distance_analysis(w, rng_dist)
        rows.append(dict(metric=grp, state=lab, **d))
    return pd.concat(pairs, ignore_index=True), pd.DataFrame(rows)

def bootstrap_pair_indices(df):
    """Bootstrap resamples of complete twin pairs within each zygosity group (drawn once, reused for all
    12 metric-state combinations)."""
    tmpl = make_pairs(df, METRICS[0], STATE_COLS[0])
    n_mz, n_dz = int((tmpl.zyg == "MZ").sum()), int((tmpl.zyg == "DZ").sum())
    rng_boot = np.random.default_rng(SEED + 1)
    idx_mz = rng_boot.integers(0, n_mz, (N_BOOT, n_mz))
    idx_dz = rng_boot.integers(0, n_dz, (N_BOOT, n_dz))
    return idx_mz, idx_dz

def build_table_s1(df):
    """All numbers of Table S1: distances (Mann-Whitney), ICCs, the MZ > DZ test and Falconer h2 with bootstrap CIs."""
    pairs, dist = run_distances(df)
    idx_mz, idx_dz = bootstrap_pair_indices(df)
    D = dist.set_index(["metric", "state"])
    rows = []
    for m in METRICS:
        for s in STATE_COLS:
            w = make_pairs(df, m, s)
            mz, dz = w[w.zyg == "MZ"], w[w.zyg == "DZ"]
            T1m, T2m, T1d, T2d = mz.t1.to_numpy(), mz.t2.to_numpy(), dz.t1.to_numpy(), dz.t2.to_numpy()
            cm, cd = icc_inference(T1m, T2m), icc_inference(T1d, T2d)
            h2 = 2 * (cm["r"] - cd["r"])
            bh = 2 * (icc1(T1m[idx_mz], T2m[idx_mz])[0] - icc1(T1d[idx_dz], T2d[idx_dz])[0])
            lo, hi = np.quantile(bh, [ALPHA / 2, 1 - ALPHA / 2])
            zi = (np.arctanh(cm["r"]) - np.arctanh(cd["r"])) / np.sqrt(1 / (cm["n"] - 1.5) + 1 / (cd["n"] - 1.5))
            d = D.loc[(m, s)]
            rows.append({
                "Metric": m, "State": s,
                "MZ mean distance": d.MZ_mean, "MZ mean distance 95% CI low": d.MZ_mean_ci_lo,
                "MZ mean distance 95% CI high": d.MZ_mean_ci_hi,
                "DZ mean distance": d.DZ_mean, "DZ mean distance 95% CI low": d.DZ_mean_ci_lo,
                "DZ mean distance 95% CI high": d.DZ_mean_ci_hi,
                "Distance shift MZ-DZ (Hodges-Lehmann)": d.HL_shift_MZ_minus_DZ,
                "Distance shift 95% CI low": d.HL_ci_lo, "Distance shift 95% CI high": d.HL_ci_hi,
                "Distance p (Mann-Whitney U, two-sided)": d.MWU_p_two_sided,
                "r_MZ (ICC)": cm["r"], "r_MZ 95% CI low": cm["ci_lo"], "r_MZ 95% CI high": cm["ci_hi"],
                "r_DZ (ICC)": cd["r"], "r_DZ 95% CI low": cd["ci_lo"], "r_DZ 95% CI high": cd["ci_hi"],
                "r_MZ p-value": cm["p"], "r_DZ p-value": cd["p"],
                "Heritability h2 (Falconer)": h2,
                "h2 95% CI low (bootstrap percentile)": lo, "h2 95% CI high (bootstrap percentile)": hi,
                "p (r_MZ > r_DZ, Fisher z, one-sided)": stats.norm.sf(zi),
                "N_MZ pairs": int(d.MZ_n_pairs), "N_DZ pairs": int(d.DZ_n_pairs)})
    T = pd.DataFrame(rows)
    T.insert(T.columns.get_loc("N_MZ pairs"), "BH-FDR q (r_MZ > r_DZ)",
             bh_fdr(T["p (r_MZ > r_DZ, Fisher z, one-sided)"].to_numpy()))
    # alcohol and odor variables: complete pairs only, own bootstrap stream, no FDR (outside the 12-test family)
    rng_x = np.random.default_rng(SEED + 7)
    xrows = []
    for grp, col, lab in EXTRA:
        w = make_pairs(df, None, None, col=col)
        w = w[w.complete]
        mz, dz = w[w.zyg == "MZ"], w[w.zyg == "DZ"]
        T1m, T2m, T1d, T2d = mz.t1.to_numpy(float), mz.t2.to_numpy(float), dz.t1.to_numpy(float), dz.t2.to_numpy(float)
        cm, cd = icc_inference(T1m, T2m), icc_inference(T1d, T2d)
        im = rng_x.integers(0, len(T1m), (N_BOOT, len(T1m)))
        idd = rng_x.integers(0, len(T1d), (N_BOOT, len(T1d)))
        bh = 2 * (icc1(T1m[im], T2m[im])[0] - icc1(T1d[idd], T2d[idd])[0])
        lo, hi = np.quantile(bh, [ALPHA / 2, 1 - ALPHA / 2])
        zi = (np.arctanh(cm["r"]) - np.arctanh(cd["r"])) / np.sqrt(1 / (cm["n"] - 1.5) + 1 / (cd["n"] - 1.5))
        d = D.loc[(grp, lab)]
        xrows.append({
            "Metric": grp, "State": lab,
            "MZ mean distance": d.MZ_mean, "MZ mean distance 95% CI low": d.MZ_mean_ci_lo,
            "MZ mean distance 95% CI high": d.MZ_mean_ci_hi,
            "DZ mean distance": d.DZ_mean, "DZ mean distance 95% CI low": d.DZ_mean_ci_lo,
            "DZ mean distance 95% CI high": d.DZ_mean_ci_hi,
            "Distance shift MZ-DZ (Hodges-Lehmann)": d.HL_shift_MZ_minus_DZ,
            "Distance shift 95% CI low": d.HL_ci_lo, "Distance shift 95% CI high": d.HL_ci_hi,
            "Distance p (Mann-Whitney U, two-sided)": d.MWU_p_two_sided,
            "r_MZ (ICC)": cm["r"], "r_MZ 95% CI low": cm["ci_lo"], "r_MZ 95% CI high": cm["ci_hi"],
            "r_DZ (ICC)": cd["r"], "r_DZ 95% CI low": cd["ci_lo"], "r_DZ 95% CI high": cd["ci_hi"],
            "r_MZ p-value": cm["p"], "r_DZ p-value": cd["p"],
            "Heritability h2 (Falconer)": 2 * (cm["r"] - cd["r"]),
            "h2 95% CI low (bootstrap percentile)": lo, "h2 95% CI high (bootstrap percentile)": hi,
            "p (r_MZ > r_DZ, Fisher z, one-sided)": stats.norm.sf(zi), "BH-FDR q (r_MZ > r_DZ)": np.nan,
            "N_MZ pairs": len(T1m), "N_DZ pairs": len(T1d)})
    return pd.concat([T, pd.DataFrame(xrows)], ignore_index=True)


# -----------------------------------------------------------------------------
# Word table
# -----------------------------------------------------------------------------
def minus(s):
    return s.replace("-", "−")


def fmt_ci(est, lo, hi, nd):
    z = lambda x: 0.0 if round(x, nd) == 0 else x            # avoid "-0.000"
    return minus(f"{z(est):.{nd}f} [{z(lo):.{nd}f}, {z(hi):.{nd}f}]")


def fmt_p_word(p, dagger=False):
    """Returns (text, superscript, bold, dagger). p < 0.001 -> 'm × 10' + superscript exponent."""
    p = float(p)
    if p < 0.001:
        m, e = f"{p:.1e}".split("e")
        return f"{m} × 10", minus(str(int(e))), True, dagger
    txt = f"{p:.3f}"
    if p < ALPHA <= float(txt):                               # 0.0498 would otherwise show as 0.050
        txt = f"{p:.4f}"
    return txt, "", p < ALPHA, dagger


TOKEN = re.compile(r"(\*\*[^*]+\*\*|\*[^*]+\*|_\{[^}]+\}|\^\{[^}]+\})")


def add_runs(par, text, size=TABLE_PT, bold=False, italic=False):
    """Inline markup: **bold**, *italic*, _{subscript}, ^{superscript}."""
    pos = 0
    for m in TOKEN.finditer(text):
        if m.start() > pos:
            _run(par, text[pos:m.start()], size, bold, italic)
        t = m.group(0)
        if t.startswith("**"):
            _run(par, t[2:-2], size, True, italic)
        elif t.startswith("_{"):
            _run(par, t[2:-1], size, bold, italic, sub=True)
        elif t.startswith("^{"):
            _run(par, t[2:-1], size, bold, italic, sup=True)
        else:
            _run(par, t[1:-1], size, bold, True)
        pos = m.end()
    if pos < len(text):
        _run(par, text[pos:], size, bold, italic)


def _run(par, text, size, bold=False, italic=False, sub=False, sup=False):
    r = par.add_run(text)
    r.font.name = FONT
    r._element.rPr.rFonts.set(qn("w:eastAsia"), FONT)
    r.font.size = Pt(size)
    r.bold, r.italic = bold, italic
    r.font.subscript, r.font.superscript = sub, sup
    return r


def set_borders(cell, top=None, bottom=None):
    """top/bottom: None (no line) or width in eighths of a point (12 = 1.5 pt, 4 = 0.5 pt)."""
    tcPr = cell._tc.get_or_add_tcPr()
    old = tcPr.find(qn("w:tcBorders"))
    if old is not None:
        tcPr.remove(old)
    b = OxmlElement("w:tcBorders")
    tcPr.append(b)
    for edge, w in (("top", top), ("left", None), ("bottom", bottom), ("right", None)):
        el = OxmlElement(f"w:{edge}")
        if w:
            el.set(qn("w:val"), "single"); el.set(qn("w:sz"), str(w)); el.set(qn("w:color"), "000000")
        else:
            el.set(qn("w:val"), "nil")
        b.append(el)


def set_margins(cell, top=40, bottom=40, left=60, right=60):
    tcPr = cell._tc.get_or_add_tcPr()
    mar = OxmlElement("w:tcMar")
    for k, v in (("top", top), ("left", left), ("bottom", bottom), ("right", right)):
        el = OxmlElement(f"w:{k}"); el.set(qn("w:w"), str(v)); el.set(qn("w:type"), "dxa"); mar.append(el)
    tcPr.append(mar)


def write_cell(cell, content, align="center", italic=False, top=None, bottom=None, valign="center"):
    """content: markup string, or a (text, sup, bold, dagger) tuple from fmt_p."""
    for extra in cell.paragraphs[1:]:                         # merged cells keep one paragraph only
        extra._element.getparent().remove(extra._element)
    par = cell.paragraphs[0]
    par.alignment = WD_ALIGN_PARAGRAPH.LEFT if align == "left" else WD_ALIGN_PARAGRAPH.CENTER
    par.paragraph_format.space_after = Pt(0); par.paragraph_format.space_before = Pt(0)
    if isinstance(content, tuple):
        txt, sup, bold, dagger = content
        _run(par, txt, TABLE_PT, bold)
        if sup:
            _run(par, sup, TABLE_PT, bold, sup=True)
        if dagger:
            _run(par, "†", TABLE_PT, bold)
    else:
        add_runs(par, content, italic=italic)
    cell.vertical_alignment = (WD_CELL_VERTICAL_ALIGNMENT.BOTTOM if valign == "bottom"
                               else WD_CELL_VERTICAL_ALIGNMENT.CENTER)
    set_borders(cell, top, bottom); set_margins(cell)


def new_document():
    doc = Document()
    sec = doc.sections[0]
    sec.orientation = WD_ORIENT.LANDSCAPE
    sec.page_width, sec.page_height = Twips(16838), Twips(11906)
    for side in ("left_margin", "right_margin", "top_margin", "bottom_margin"):
        setattr(sec, side, Twips(1000))
    st = doc.styles["Normal"]; st.font.name = FONT; st.font.size = Pt(10)
    st.element.rPr.rFonts.set(qn("w:eastAsia"), FONT)
    return doc


def build(doc, widths, head1, head2, rows, notes):
    """head1: list of (text, n_cols, rowspan2); head2: texts for columns not covered by a rowspan."""
    p = doc.add_paragraph(); p.paragraph_format.space_after = Pt(8)
    add_runs(p, TITLE, size=TITLE_PT)
    ncol = len(widths)
    t = doc.add_table(rows=2, cols=ncol)
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    t.autofit = False
    tblPr = t._tbl.tblPr
    lay = OxmlElement("w:tblLayout"); lay.set(qn("w:type"), "fixed"); tblPr.append(lay)
    for i, w in enumerate(widths):
        t.columns[i].width = Twips(w)
        for rr in range(2):
            t.cell(rr, i).width = Twips(w)                    # set before merging so spans get summed widths
    # header row 1 (spanners / row-spanning cells) and row 2
    col, spanned_rows = 0, set()
    for text, n, rowspan in head1:
        a = t.cell(0, col)
        if n > 1:
            a = a.merge(t.cell(0, col + n - 1))
        if rowspan:
            a = a.merge(t.cell(1, col)); spanned_rows.add(col)
        write_cell(a, text, align="left" if col == 0 else "center", top=12, bottom=4,
                   valign="bottom" if rowspan else "center")
        col += n
    it = iter(head2)
    for c in range(ncol):
        if c not in spanned_rows:
            write_cell(t.cell(1, c), next(it), bottom=4)
    # bottom rule under row-spanning header cells is drawn on their continuation cell in row 2
    grid = 0
    for tc in t.rows[1]._tr.tc_lst:
        if grid in spanned_rows:
            set_borders(_Cell(tc, t), top=None, bottom=4)
        grid += tc.grid_span
    # body
    for gi, (gname, grows) in enumerate(rows):
        r = t.add_row()
        for c, w in enumerate(widths):
            r.cells[c].width = Twips(w)
        g = r.cells[0].merge(r.cells[ncol - 1])
        write_cell(g, gname, align="left", italic=True)
        for ri, vals in enumerate(grows):
            last = gi == len(rows) - 1 and ri == len(grows) - 1
            r = t.add_row()
            for c, w in enumerate(widths):
                r.cells[c].width = Twips(w)
            for c, v in enumerate(vals):
                write_cell(r.cells[c], v, align="left" if c == 0 else "center", bottom=12 if last else None)
    for row in t.rows:                                        # keep each row on one page
        trPr = row._tr.get_or_add_trPr(); cs = OxmlElement("w:cantSplit"); trPr.append(cs)
    hdr = t.rows[0]._tr.get_or_add_trPr(); th = OxmlElement("w:tblHeader"); hdr.append(th)
    doc.add_paragraph().paragraph_format.space_after = Pt(2)
    for n in notes:
        p = doc.add_paragraph(); p.paragraph_format.space_after = Pt(3)
        add_runs(p, n, size=NOTE_PT)


def word_table_corrected(T, out_path):
    T = T.set_index(["Metric", "State"])
    rows = []
    for m, gname in GROUPS:
        g = []
        for s in STATE_COLS:
            r = T.loc[(m, s)]
            h2 = fmt_ci(r["Heritability h2 (Falconer)"], r["h2 95% CI low (bootstrap percentile)"],
                        r["h2 95% CI high (bootstrap percentile)"], 2)
            if not 0 <= r["Heritability h2 (Falconer)"] <= 1:
                h2 += "^{c}"
            g.append([STATE_LABEL[s],
                      fmt_ci(r["MZ mean distance"], r["MZ mean distance 95% CI low"], r["MZ mean distance 95% CI high"], DIST_DECIMALS[m]),
                      fmt_ci(r["DZ mean distance"], r["DZ mean distance 95% CI low"], r["DZ mean distance 95% CI high"], DIST_DECIMALS[m]),
                      fmt_p_word(r["Distance p (Mann-Whitney U, two-sided)"]),
                      fmt_ci(r["r_MZ (ICC)"], r["r_MZ 95% CI low"], r["r_MZ 95% CI high"], 2),
                      fmt_ci(r["r_DZ (ICC)"], r["r_DZ 95% CI low"], r["r_DZ 95% CI high"], 2),
                      fmt_p_word(r["r_MZ p-value"]), fmt_p_word(r["r_DZ p-value"]),
                      fmt_p_word(r["p (r_MZ > r_DZ, Fisher z, one-sided)"], dagger=r["BH-FDR q (r_MZ > r_DZ)"] < ALPHA),
                      h2, str(int(r["N_MZ pairs"])), str(int(r["N_DZ pairs"]))])
        rows.append((gname, g))
    for grp in dict.fromkeys(e[0] for e in EXTRA):
        g = []
        for _, col, lab in [e for e in EXTRA if e[0] == grp]:
            r = T.loc[(grp, lab)]
            h2 = fmt_ci(r["Heritability h2 (Falconer)"], r["h2 95% CI low (bootstrap percentile)"],
                        r["h2 95% CI high (bootstrap percentile)"], 2)
            if not 0 <= r["Heritability h2 (Falconer)"] <= 1:
                h2 += "^{c}"
            nd = EXTRA_DECIMALS
            g.append([lab,
                      fmt_ci(r["MZ mean distance"], r["MZ mean distance 95% CI low"], r["MZ mean distance 95% CI high"], nd),
                      fmt_ci(r["DZ mean distance"], r["DZ mean distance 95% CI low"], r["DZ mean distance 95% CI high"], nd),
                      fmt_p_word(r["Distance p (Mann-Whitney U, two-sided)"]),
                      fmt_ci(r["r_MZ (ICC)"], r["r_MZ 95% CI low"], r["r_MZ 95% CI high"], 2),
                      fmt_ci(r["r_DZ (ICC)"], r["r_DZ 95% CI low"], r["r_DZ 95% CI high"], 2),
                      fmt_p_word(r["r_MZ p-value"]), fmt_p_word(r["r_DZ p-value"]),
                      fmt_p_word(r["p (r_MZ > r_DZ, Fisher z, one-sided)"]),
                      h2, str(int(r["N_MZ pairs"])), str(int(r["N_DZ pairs"]))])
        rows.append((grp, g))
    doc = new_document()
    build(doc, widths=[780, 1960, 1960, 700, 1460, 1580, 1100, 760, 1060, 1600, 500, 500],
          head1=[("State", 1, True), ("Within-pair distance", 3, False), ("Cross-twin correlation (ICC)", 4, False),
                 ("Falconer heritability", 2, False), ("Pairs, *n*", 2, False)],
          head2=["MZ, mean [95% CI]", "DZ, mean [95% CI]", "*p*^{a}", "*r*_{MZ} [95% CI]", "*r*_{DZ} [95% CI]",
                 "*p* (*r*_{MZ})", "*p* (*r*_{DZ})", "*p* (*r*_{MZ} > *r*_{DZ})^{b}", "h^{2} [95% CI]", "MZ", "DZ"],
          rows=rows,
          notes=["*Note.* Within-pair distance = |twin 1 − twin 2|, with group means and 95% bias-corrected and accelerated "
                 "(BCa) bootstrap confidence intervals (CIs; 10,000 resamples). *r*_{MZ} and *r*_{DZ} are one-way random-effects "
                 "intraclass correlations, ICC(1), which do not depend on which co-twin is labelled twin 1. Their 95% CIs are "
                 "exact *F*-based intervals, and *p* (*r*_{MZ}) and *p* (*r*_{DZ}) are two-sided *F* tests of ρ = 0; they show "
                 "twin resemblance but do not test MZ–DZ differences. h^{2} = 2(*r*_{MZ} − *r*_{DZ}) (Falconer’s formula), "
                 "reported without truncation, with 95% percentile CIs from 10,000 bootstrap resamples of complete twin pairs "
                 "within each zygosity group (random seed 20261111). *p*-values are uncorrected for multiple comparisons; bold "
                 "indicates *p* < 0.05. †Remains significant after Benjamini–Hochberg false-discovery-rate correction across "
                 "the 12 metric–state combinations (*q* = 0.025).",
                 "^{a}Two-sided Mann–Whitney *U* test comparing MZ and DZ within-pair distances.",
                 "^{b}One-sided Fisher *r*-to-*z* test of *r*_{MZ} > *r*_{DZ} for independent intraclass correlations "
                 "[var(*z*) = 1/(*n* − 3/2)]. Because h^{2} = 2(*r*_{MZ} − *r*_{DZ}), this is also the test of h^{2} > 0.",
                 "^{c}Exceeds the theoretical maximum of 1. Falconer’s estimator is unconstrained, so this reflects sampling "
                 "variability. Falconer estimates are not equivalent to ACE-model variance components."])
    doc.save(out_path)


def word_table_original(T, out_path):
    T = T.set_index(["Metric", "State"])
    rows = []
    for m, gname in GROUPS:
        g = []
        for s in STATE_COLS:
            r = T.loc[(m, s)]
            h2 = fmt_ci(r["Heritability h2 (Falconer)"], r["h2 95% CI low (bootstrap percentile)"],
                        r["h2 95% CI high (bootstrap percentile)"], 2)
            if not 0 <= r["Heritability h2 (Falconer)"] <= 1:
                h2 += "^{b}"
            g.append([STATE_LABEL[s],
                      fmt_ci(r["MZ mean distance"], r["MZ mean distance 95% CI low"], r["MZ mean distance 95% CI high"], DIST_DECIMALS[m]),
                      fmt_ci(r["DZ mean distance"], r["DZ mean distance 95% CI low"], r["DZ mean distance 95% CI high"], DIST_DECIMALS[m]),
                      fmt_p_word(r["Distance p (Mann-Whitney U, one-sided H1: MZ < DZ)"]),
                      fmt_ci(r["r_MZ (Pearson)"], r["r_MZ 95% CI low"], r["r_MZ 95% CI high"], 2),
                      fmt_ci(r["r_DZ (Pearson)"], r["r_DZ 95% CI low"], r["r_DZ 95% CI high"], 2),
                      fmt_p_word(r["r_MZ p-value"]), fmt_p_word(r["r_DZ p-value"]),
                      h2, str(int(r["N_MZ pairs"])), str(int(r["N_DZ pairs"]))])
        rows.append((gname, g))
    doc = new_document()
    build(doc, widths=[900, 1820, 1820, 820, 1560, 1640, 1180, 880, 1700, 560, 560],
          head1=[("State", 1, True), ("Within-pair distance", 3, False), ("Cross-twin correlation", 4, False),
                 ("Falconer h^{2} [95% CI]", 1, True), ("Pairs, *n*", 2, False)],
          head2=["MZ, mean [95% CI]", "DZ, mean [95% CI]", "*p*^{a}", "*r*_{MZ} [95% CI]", "*r*_{DZ} [95% CI]",
                 "*p* (*r*_{MZ})", "*p* (*r*_{DZ})", "MZ", "DZ"],
          rows=rows,
          notes=["*Note.* Within-pair distance = |twin 1 − twin 2|, with group means and 95% bias-corrected and accelerated "
                 "(BCa) bootstrap confidence intervals (CIs; 10,000 resamples). *r*_{MZ} and *r*_{DZ} are Pearson cross-twin "
                 "correlations, with twin 1 defined as the co-twin with the lower subject ID. Their 95% CIs are from Fisher’s "
                 "*z* transformation, and the *p*-values test *r* = 0 (two-sided); they do not test the MZ–DZ difference. "
                 "h^{2} = 2(*r*_{MZ} − *r*_{DZ}) (Falconer’s formula), reported without truncation, with 95% percentile CIs "
                 "from 10,000 bootstrap resamples of complete twin pairs within each zygosity group (random seed 20261111). "
                 "*p*-values are uncorrected for multiple comparisons; bold indicates *p* < 0.05.",
                 "^{a}One-sided Mann–Whitney *U* test (H_{1}: MZ distances < DZ distances); values above 0.5 indicate that MZ "
                 "distances tended to be larger than DZ distances.",
                 "^{b}Exceeds the theoretical maximum of 1. Falconer’s estimator is unconstrained, so this reflects sampling "
                 "variability. Falconer estimates are not equivalent to ACE-model variance components."])
    doc.save(out_path)

# -----------------------------------------------------------------------------
# Excel output
# -----------------------------------------------------------------------------

def save_excel(sheets, path):
    """Write {sheet name: table} to one Excel file (Arial 10, bold header row, Yes/No for true/false)."""
    from openpyxl import Workbook
    from openpyxl.styles import Font, Alignment
    from openpyxl.utils import get_column_letter
    wb = Workbook(); wb.remove(wb.active)
    for name, df in sheets.items():
        ws = wb.create_sheet(name[:31])
        for j, col in enumerate(df.columns, 1):
            c = ws.cell(1, j, str(col)); c.font = Font(name="Arial", size=10, bold=True)
            c.alignment = Alignment(wrap_text=True, vertical="top")
            ws.column_dimensions[get_column_letter(j)].width = min(max(len(str(col)) * 0.9, 10), 34)
        for i, row in enumerate(df.itertuples(index=False), 2):
            for j, v in enumerate(row, 1):
                if isinstance(v, (bool, np.bool_)):
                    v = "Yes" if v else "No"
                elif isinstance(v, np.integer):
                    v = int(v)
                elif isinstance(v, (float, np.floating)):
                    v = None if not np.isfinite(v) else float(v)
                ws.cell(i, j, v).font = Font(name="Arial", size=10)
        ws.freeze_panes = "A2"
    wb.save(path)

def main():
    df = load_master()
    out = Path(OUTPUT_FOLDER); out.mkdir(parents=True, exist_ok=True)
    T = build_table_s1(df)
    docx = out / f"{FILE_PREFIX}-table-S1-journal-format-corrected.docx"
    word_table_corrected(T, docx)
    xlsx = out / f"{FILE_PREFIX}-TableS1-corrected-numbers.xlsx"
    save_excel({"Table_S1_corrected": T}, xlsx)
    print("Saved:", str(docx), str(xlsx), sep="\n  ")
    return T


if __name__ == "__main__":
    main()
