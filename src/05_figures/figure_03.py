"""
Figure 3, version 11: six panels with large text, in the order of the Results text.

Row 1, brain-state dynamics and alcohol consumption
  A  Pair similarity (-|twin 1 - twin 2|) of all 12 brain-state measures with drinks per drinking occasion and drinks/week
  B  Absolute within-pair differences: Visual Network-A dwell time vs drinks per drinking occasion
  C  Visual Network-A dwell time and drinks per drinking occasion in the pair-similarity, signed-difference, within-pair and
     between-pair models
Row 2, brain-state dynamics and odor identification
  D  Within-pair differences (twin 2 - twin 1): FO DMN vs odor identification (mixed model, within-pair term)
  E  Within-pair differences (twin 2 - twin 1): FO Visual Network-B vs odor identification (same model)
  F  Pair similarity (-|twin 1 - twin 2|): DMN appearance rate vs odor identification (pair-similarity model)

Statistics printed in the panels come from the PI's workbook (Brain_Dynamics_Alcohol_Odor_Results_Final.xlsx,
from "GLM code_alcohol_ordinal.R"): beta, 95% CI, p and q_12 exactly as in the workbook, so the figure and the
results tables agree (P_SOURCE = "wald" restores the earlier Wald p from beta / SE with BH q across the 12 brain
measures). Drinks per drinking occasion is modelled as an ordinal outcome (HCP 0-6 code): its pair-similarity and within/between
estimates are log-odds (cumulative-link models); signed differences are linear (raw code units). Drinks per week
estimates are linear. Points are the raw pair-level data from the master sheet (twin 1 = lower Subject ID); lines
are simple least-squares fits with a shaded 95% band (through the origin for signed differences).

Usage:
  python figure3.py PI_RESULTS.xlsx MASTER.xlsx OUT_DIR
Requires: pandas, numpy, scipy, matplotlib, statsmodels, openpyxl, python-calamine
"""
import os
import sys
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyBboxPatch, FancyArrowPatch, Polygon
from scipy import stats
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")
PI_X, MASTER, OUT = sys.argv[1:4]
os.makedirs(OUT, exist_ok=True)
STEM = os.path.join(OUT, "2026-11-16-figure3-v12")
P_SOURCE = "workbook"          # "workbook": p and q_12 from the workbook; "wald": recomputed from beta / SE

# ------------------------------------------------------------------ style: large, high-contrast text
INK, INK2, MUTED, BAND = "#000000", "#000000", "#000000", "#ffffff"   # no grey: all lines and text black
BLUE, ORANGE = "#2a78d6", "#eb6834"
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 15, "axes.labelsize": 16, "xtick.labelsize": 15, "ytick.labelsize": 15,
    "axes.edgecolor": INK, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK, "text.color": INK,
    "axes.linewidth": 1.4, "xtick.major.width": 1.4, "ytick.major.width": 1.4, "xtick.major.size": 6,
    "ytick.major.size": 6, "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42,
    "legend.frameon": False, "legend.fontsize": 14, "axes.labelpad": 6, "axes.labelweight": "bold",
})
STAT_BOX = dict(boxstyle="round,pad=0.35", fc="white", ec=MUTED, lw=1.0, alpha=0.95)

# ------------------------------------------------------------------ PI results, Wald p and q
PI = pd.read_excel(PI_X, sheet_name=None, engine="calamine")
for k, d in PI.items():
    if "beta" in d:
        d["p_wald"] = 2 * stats.norm.sf(np.abs(d.beta / d.SE))
for k, grp in [("Alcohol_Pair_Similarity", ["outcome"]), ("Alcohol_Signed_Differences", ["outcome"]),
               ("Alcohol_Mixed_Effects", ["outcome", "effect"]), ("Odor_Brain_Mixed", ["effect"]),
               ("Odor_Brain_Pair", None)]:
    d = PI[k]
    if grp is None:
        d["q_wald"] = multipletests(d.p_wald, method="fdr_bh")[1]
    else:
        d["q_wald"] = d.groupby(grp).p_wald.transform(lambda s: pd.Series(multipletests(s, method="fdr_bh")[1], index=s.index))
    if P_SOURCE == "workbook":                    # the values printed in the results tables
        d["p_wald"], d["q_wald"] = d["p"], d["q_12"]
DAY, WEEK = "SSAGA_Alc_12_Drinks_Per_Day", "Total_Drinks_7days"


def pick(sheet, **kw):
    d = PI[sheet]
    for k, v in kw.items():
        d = d[d[k] == v]
    assert len(d) == 1, (sheet, kw)
    return d.iloc[0]


def f2(x):
    s = f"{x:.2f}".replace("-", "−")
    return s


def f3(x):
    return f"{x:.3f}".replace("-", "\u2212")


def fp(p):
    return "p < 0.001" if p < .001 else f"p = {p:.3f}"


def stat_text(r, q=True):
    s = f"β = {f3(r.beta)}\n95% CI [{f3(r.CI_low)}, {f3(r.CI_high)}]\n{fp(r.p_wald)}"
    if q:
        s += f", q = {r.q_wald:.3f}"
    return s


# ------------------------------------------------------------------ raw pair data
ms = pd.read_excel(MASTER, sheet_name="Master")
ms = ms.rename(columns={f"{m}_State{k}": f"{m}{k}" for m in ("FO", "DT", "AR") for k in range(1, 5)})
rows = []
for pid, g in ms.groupby("TwinPairID"):
    g = g.sort_values("Subject")
    if len(g) != 2:
        continue
    t1, t2 = g.iloc[0], g.iloc[1]
    r = {"pair": pid, "zyg": t1.Zygosity}
    for v in ["DT1", "FO2", "FO3", "AR2"]:
        r[v + "_sd"] = t2[v] - t1[v]
    for v, col in [("odor", "Odor_Unadj"), ("day", DAY)]:
        a, b = pd.to_numeric(t1[col], errors="coerce"), pd.to_numeric(t2[col], errors="coerce")
        r[v + "_sd"] = b - a
    rows.append(r)
P = pd.DataFrame(rows)
for v in ["DT1", "AR2", "odor", "day"]:
    P[v + "_ad"] = P[v + "_sd"].abs()

# ------------------------------------------------------------------ layout: row 1 alcohol (A-C), row 2 odor (D-F)
fig = plt.figure(figsize=(15.5, 13.3))
fig.patch.set_facecolor("white")
gs = fig.add_gridspec(2, 3, left=0.16, right=0.985, top=0.862, bottom=0.07, wspace=0.50, hspace=1.02,
                      width_ratios=[1.05, 1, 1])
axA, axB, axC = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[0, 2])
axD, axE, axF = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1]), fig.add_subplot(gs[1, 2])
SRC = {}
ABS = "|twin 1 − twin 2|"
SGN = "twin 2 − twin 1"


def label(ax, letter, title, dx=0.0):
    bb = ax.get_position()
    x = bb.x0 - 0.062 + dx
    fig.text(x, bb.y1 + 0.052, letter, fontsize=26, fontweight="bold", va="bottom")
    fig.text(x + 0.42 / fig.get_figwidth(), bb.y1 + 0.054, title, fontsize=17, fontweight="bold", va="bottom")


def box(ax, r, q=False):
    """Statistics as two lines just above the plotting area, so the data can fill the axes."""
    s = f"\u03b2 = {f3(r.beta)}, 95% CI [{f3(r.CI_low)}, {f3(r.CI_high)}]\n{fp(r.p_wald)}"
    if q:
        s += f", q = {r.q_wald:.3f}"
    ax.text(1.0, 1.025, s, transform=ax.transAxes, ha="right", va="bottom", fontsize=13.5, linespacing=1.2)


# (scatter_fit follows)

def scatter_fit(ax, x, y, zyg, through_origin=False, jitter_y=0.0, seed=1):
    rng = np.random.default_rng(seed)
    yj = y + (rng.uniform(-jitter_y, jitter_y, len(y)) if jitter_y else 0)
    for z, col, mk in [("MZ", BLUE, "o"), ("DZ", ORANGE, "s")]:
        m = zyg == z
        ax.scatter(x[m], yj[m], s=60 if mk == "o" else 52, marker=mk, facecolor=col, edgecolor="white",
                   linewidth=0.8, alpha=0.85, zorder=3, label={"MZ": "Monozygotic", "DZ": "Dizygotic"}[z])
    xs = np.linspace(x.min(), x.max(), 100)
    if through_origin:
        b = (x @ y) / (x @ x); res = y - b * x; s2 = res @ res / (len(x) - 1); se = np.sqrt(s2 / (x @ x))
        tcrit = stats.t.ppf(.975, len(x) - 1); yhat = b * xs; band = tcrit * se * np.abs(xs)
    else:
        X = np.column_stack([np.ones_like(x), x]); beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        res = y - X @ beta; s2 = res @ res / (len(x) - 2); cov = s2 * np.linalg.inv(X.T @ X)
        Xs = np.column_stack([np.ones_like(xs), xs]); yhat = Xs @ beta
        band = stats.t.ppf(.975, len(x) - 2) * np.sqrt(np.einsum("ij,jk,ik->i", Xs, cov, Xs))
    ax.fill_between(xs, yhat - band, yhat + band, color=INK, alpha=0.12, lw=0, zorder=1)   # 95% band (shaded)
    ax.plot(xs, yhat, color=INK, lw=2.6, zorder=4)




# ================================================================== A: pair similarity, all 12 brain measures
STATES = [("1", "Visual Network-A"), ("2", "DMN"), ("3", "Visual Network-B"), ("4", "Hippocampal–Temporoparietal")]
order = [(m, s, lab) for m in ["FO", "DT", "AR"] for s, lab in STATES]
names, recA = [], []
rKey = pick("Alcohol_Pair_Similarity", brain="DT_State 1", outcome=DAY)
for i, (m, s, lab) in enumerate(order):
    names.append("Hippocampal–\nTemporoparietal" if s == "4" else lab)   # two lines so it fits
    for oc, col, mk, dy in [(DAY, BLUE, "o", -0.17), (WEEK, ORANGE, "s", 0.17)]:
        r = pick("Alcohol_Pair_Similarity", brain=f"{m}_State {s}", outcome=oc)
        sig = r.p_wald < .05
        axA.plot([r.CI_low, r.CI_high], [i + dy, i + dy], color=col, lw=2.4, solid_capstyle="butt", zorder=3)
        axA.plot(r.beta, i + dy, marker=mk, ms=10 if mk == "o" else 9, mfc=col if sig else "white", mec=col, mew=2.0, zorder=4)
        recA.append({"brain": f"{m}_State {s}", "outcome": oc, "beta": r.beta, "CI_low": r.CI_low, "CI_high": r.CI_high,
                     "p_wald": r.p_wald, "q_wald_12": r.q_wald})
hi = 4   # dwell time, Visual Network-A
axA.axhspan(hi - 0.5, hi + 0.5, color="#fff1c2", zorder=0)
axA.axvline(0, color=INK2, lw=1.4, zorder=2)
for k in (3.5, 7.5):
    axA.axhline(k, color=MUTED, lw=0.9, ls=(0, (3, 3)), zorder=1)
axA.set_ylim(len(order) - 0.5, -0.5)
axA.set_yticks(range(len(order))); axA.set_yticklabels(names, fontsize=14)
axA.get_yticklabels()[hi].set_fontweight("bold")
for t in axA.get_yticklabels():
    if "\n" in t.get_text():
        t.set_fontsize(12.5); t.set_linespacing(0.9)
axA.tick_params(axis="y", length=0); axA.spines["left"].set_visible(False)
XA = 1.32
axA.set_xlim(-XA, XA); axA.set_xticks([-1, -0.5, 0, 0.5, 1])
axA.set_xticklabels(["−1.0", "−0.5", "0", "0.5", "1.0"])
axA.set_xlabel("β (95% CI)")
for g, gname in enumerate(["Fractional\noccupancy\n(proportion)", "Dwell time\n(s)", "Appearance\nrate\n(per min)"]):
    y0, y1 = g * 4 - 0.35, g * 4 + 3.35
    axA.annotate("", xy=(-0.54, y0), xytext=(-0.54, y1), xycoords=("axes fraction", "data"),
                 arrowprops=dict(arrowstyle="-", color=INK2, lw=1.4), annotation_clip=False)
    axA.text(-0.655, g * 4 + 1.5, gname, transform=axA.get_yaxis_transform(), rotation=90, ha="center", va="center",
             fontsize=12.5, fontweight="bold", linespacing=1.0)
axA.legend(handles=[Line2D([], [], color=BLUE, marker="o", ms=9, lw=2.4, label="Drinks per drinking occasion (log-odds)"),
                    Line2D([], [], color=ORANGE, marker="s", ms=8, lw=2.4, label="Weekly alcohol consumption (drinks)")],
           loc="upper center", bbox_to_anchor=(0.45, -0.20), ncol=1, handlelength=1.6, borderaxespad=0, fontsize=13.5)
axA.text(XA - 0.02, hi, fp(rKey.p_wald) + f"\nq = {rKey.q_wald:.3f}", va="center", ha="right", fontsize=13,
         fontweight="bold", color=BLUE, linespacing=1.0)
label(axA, "A", "Pair similarity: all brain-state measures", dx=-0.085)
SRC["A"] = pd.DataFrame(recA)

# ================================================================== B: absolute differences, Visual-A DT and drinks per drinking occasion
d = P.dropna(subset=["DT1_ad", "day_ad"])
scatter_fit(axB, d.DT1_ad.values, d.day_ad.values, d.zyg.values, jitter_y=0.13, seed=7)
box(axB, rKey, q=True)
axB.set_xlabel("|Δ dwell time| (s)\nVisual Network-A state")
axB.set_ylabel("|Δ drinks per drinking occasion|\n(SSAGA code, 0–6)")
axB.set_ylim(-0.45, 4.45); axB.set_yticks([0, 1, 2, 3, 4])
label(axB, "B", "Absolute difference", dx=-0.03)
SRC["B_points"] = d[["pair", "zyg", "DT1_ad", "day_ad"]]

# ================================================================== C: Visual-A DT and drinks per drinking occasion in each model
MODELS = [("Pair similarity", "Alcohol_Pair_Similarity", None),
          ("Signed difference", "Alcohol_Signed_Differences", None),
          ("Within pair", "Alcohol_Mixed_Effects", "within"),
          ("Between pair", "Alcohol_Mixed_Effects", "between")]
recC = []
for i, (lab, sh, eff) in enumerate(MODELS):
    r = pick(sh, brain="DT_State 1", outcome=DAY, **({"effect": eff} if eff else {}))
    sig = r.p_wald < .05
    axC.plot([r.CI_low, r.CI_high], [i, i], color=BLUE, lw=3.0, solid_capstyle="butt", zorder=3)
    axC.plot(r.beta, i, marker="o", ms=13, mfc=BLUE if sig else "white", mec=BLUE, mew=2.4, zorder=4)
    axC.text(1.28, i, fp(r.p_wald), va="center", ha="right", fontsize=14.5,
             fontweight="bold" if sig else "normal", color=BLUE if sig else INK2)
    recC.append({"model": lab.replace("\n", " "), "beta": r.beta, "CI_low": r.CI_low, "CI_high": r.CI_high, "p_wald": r.p_wald})
axC.axvline(0, color=INK2, lw=1.4, zorder=1)
axC.set_ylim(len(MODELS) - 0.5, -0.6)
axC.set_yticks(range(len(MODELS))); axC.set_yticklabels([m[0] for m in MODELS], fontsize=15)
axC.tick_params(axis="y", length=0); axC.spines["left"].set_visible(False)
axC.set_xlim(-1.3, 1.3); axC.set_xticks([-1.0, -0.5, 0, 0.5])
axC.set_xticklabels(["−1.0", "−0.5", "0", "0.5"])
axC.set_xlabel("β (95% CI, log-odds)\ndrinks per drinking occasion")
label(axC, "C", "Visual Network-A dwell time", dx=-0.06)
SRC["C"] = pd.DataFrame(recC)

# ================================================================== D, E: within-pair (signed differences), FO and odor ID
for ax, key, brain, state, letter, seed in [(axD, "FO2", "FO_State 2", "DMN", "D", 3),
                                            (axE, "FO3", "FO_State 3", "Visual Network-B", "E", 4)]:
    d = P.dropna(subset=[key + "_sd", "odor_sd"])
    scatter_fit(ax, d[key + "_sd"].values, d.odor_sd.values, d.zyg.values, through_origin=True, seed=seed)
    r = pick("Odor_Brain_Mixed", brain=brain, effect="within")
    box(ax, r, q=True)
    ax.axhline(0, color=MUTED, lw=1.0, zorder=0); ax.axvline(0, color=MUTED, lw=1.0, zorder=0)
    ax.set_xlabel(f"Δ fractional occupancy (proportion)\n(twin 2 − twin 1), {state} state")
    ax.set_ylabel("Δ odor identification score (points)")
    lim = np.nanmax(np.abs(d[key + "_sd"])) * 1.08; ax.set_xlim(-lim, lim)
    ax.set_ylim(-33, 33); ax.set_yticks([-30, -15, 0, 15, 30])
    label(ax, letter, "Within-pair difference", dx=-0.03)
    SRC[letter + "_points"] = d[["pair", "zyg", key + "_sd", "odor_sd"]]
rDf = pick("Odor_Brain_Mixed", brain="FO_State 2", effect="within")
rEf = pick("Odor_Brain_Mixed", brain="FO_State 3", effect="within")

# ================================================================== F: pair similarity, DMN AR and odor ID
# similarity = -|twin 1 - twin 2| (0 = identical twins), the quantity the pair-level model uses
d = P.dropna(subset=["AR2_ad", "odor_ad"]).copy()
d["AR2_sim"], d["odor_sim"] = -d.AR2_ad, -d.odor_ad
scatter_fit(axF, d.AR2_sim.values, d.odor_sim.values, d.zyg.values, seed=5)
rF = pick("Odor_Brain_Pair", brain="AR_State 2")
box(axF, rF, q=True)
axF.set_xlabel("Appearance-rate similarity,\n−|Δ| (per min), DMN state")
axF.set_ylabel("Odor similarity, −|Δ score| (points)")
axF.set_ylim(-31.5, 1.8); axF.set_yticks([-30, -20, -10, 0])
label(axF, "F", "Pair similarity", dx=-0.03)
SRC["F_points"] = d[["pair", "zyg", "AR2_sim", "odor_sim"]]

# row headings
for ax, txt in [(axA, "Brain-state dynamics and alcohol consumption"), (axD, "Brain-state dynamics and odor identification")]:
    y = ax.get_position().y1 + 0.112
    fig.text(0.012, y, txt, fontsize=20, fontweight="bold", color=INK2, va="bottom")
    fig.add_artist(Line2D([0.012, 0.985], [y - 0.007, y - 0.007], color=MUTED, lw=1.2))
    # one shared key for the twin-pair markers, at the right end of each row heading
    fig.legend(handles=[Line2D([], [], ls="", marker="o", ms=10, mfc=BLUE, mec="white", label="Monozygotic pair"),
                        Line2D([], [], ls="", marker="s", ms=9, mfc=ORANGE, mec="white", label="Dizygotic pair")],
               loc="lower right", bbox_to_anchor=(0.985, y - 0.002), ncol=2, fontsize=15, handletextpad=0.3,
               columnspacing=1.4, borderaxespad=0, frameon=False)


fig.savefig(STEM + ".pdf")
fig.savefig(STEM + ".png", dpi=500)   # under 8,000 px wide; the TIFF keeps 600 dpi
fig.savefig(STEM + ".tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
SRC["stats_printed"] = pd.DataFrame([
    {"panel": "A/B", "model": "Alcohol_Pair_Similarity DT_State 1 drinks per drinking occasion", **rKey[["beta", "CI_low", "CI_high", "p_wald", "q_wald"]].to_dict()},
    {"panel": "D", "model": "Odor_Brain_Mixed FO_State 2 within", **rDf[["beta", "CI_low", "CI_high", "p_wald", "q_wald"]].to_dict()},
    {"panel": "E", "model": "Odor_Brain_Mixed FO_State 3 within", **rEf[["beta", "CI_low", "CI_high", "p_wald", "q_wald"]].to_dict()},
    {"panel": "F", "model": "Odor_Brain_Pair AR_State 2", **rF[["beta", "CI_low", "CI_high", "p_wald", "q_wald"]].to_dict()}])
with pd.ExcelWriter(STEM + "-source-data.xlsx", engine="openpyxl") as xw:
    for k, v in SRC.items():
        v.to_excel(xw, sheet_name=k[:31], index=False)
print("saved", STEM)
print(SRC["stats_printed"].round(3).to_string(index=False))
