"""Figure 2 - twin similarity in brain-state dynamics and alcohol consumption (panels A-G).

Look of the lab's Figure 2 (metric colours; MZ filled / DZ open markers; MZ-DZ joining lines; brackets),
with: a labelled y-axis with units on every graph, p and q on every test, all statistics at 3 decimals,
significant values (p or q < 0.05) in red, no grey, and larger text.

Statistics are read from the twin-correlation workbook written by "Supplementary Table S5.py"
(outputs/2026-11-13-TableS1-corrected-numbers.xlsx); alcohol z-scaling uses participant SDs from the master sheet.
  A, D, E  mean |twin 1 - twin 2| with 95% BCa CI; p two-sided Mann-Whitney U (MZ vs DZ)
  B, F     ICC(1) for MZ and DZ; p one-sided Fisher r-to-z test of r_MZ > r_DZ
  C, G     Falconer h2 = 2(r_MZ - r_DZ) with 95% pair-bootstrap CI (test as in B / F)
  q        Benjamini-Hochberg FDR: across the 12 metric-state tests (A-C); across the 2 alcohol measures (D-G)

Run:  python "Figure 2.py" [NUMBERS.xlsx] [MASTER.xlsx] [OUT_DIR]
"""
import sys
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.legend_handler import HandlerBase
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from statsmodels.stats.multitest import multipletests

NUMBERS = sys.argv[1] if len(sys.argv) > 1 else r"C:\Users\Tanya\Downloads\Mean_diff_MZ_DZ.xlsx"
MASTER = sys.argv[2] if len(sys.argv) > 2 else r"C:\Users\Tanya\Desktop\Figures and Tables\Excels\Master-sheet.xlsx"
OUT = Path(sys.argv[3]) if len(sys.argv) > 3 else Path(r"C:\Users\Tanya\Desktop\Figures and Tables")

# =============================================================================
# DATA
# =============================================================================
STATE_TICK = ["Visual-A", "DMN", "Visual-B", "Hippocampal–\nTemporoparietal"]
METRICS = ["FO", "DT", "AR"]
METRIC_TITLES = {"FO": "Fractional occupancy", "DT": "Dwell time", "AR": "Appearance rate"}
A_YLABEL = {"FO": "|Δ fractional occupancy|\n(proportion)", "DT": "|Δ dwell time| (s)",
            "AR": "|Δ appearance rate|\n(per min)"}
ALC = ["Weekly alcohol consumption", "Drinks per drinking occasion"]
ALC_KEY = {"Weekly alcohol consumption": ("Drinks per week", "Total_Drinks_7days"),
           "Drinks per drinking occasion": ("Drinks per drinking occasion", "SSAGA_Alc_12_Drinks_Per_Day")}

N = pd.read_excel(NUMBERS)
P_DIST, P_ICC = "Distance p (Mann-Whitney U, two-sided)", "p (r_MZ > r_DZ, Fisher z, one-sided)"
brain = N[N.Metric.isin(METRICS)].copy()
brain["q_dist"] = multipletests(brain[P_DIST], method="fdr_bh")[1]
brain["q_icc"] = multipletests(brain[P_ICC], method="fdr_bh")[1]
alc = N[N.Metric == "Alcohol consumption"].set_index("State")
# Accept legacy workbook labels without changing the underlying values.
alc = alc.rename(index={"Drinks per day": "Drinks per drinking occasion"})
alc["q_dist"] = multipletests(alc[P_DIST], method="fdr_bh")[1]
alc["q_icc"] = multipletests(alc[P_ICC], method="fdr_bh")[1]
ms = pd.read_excel(MASTER, sheet_name="Master")


def row(m, s):
    return brain[(brain.Metric == m) & (brain.State == f"State {s + 1}")].iloc[0]


def dist(r, z, scale=1.0):
    return tuple(r[c] / scale for c in (f"{z} mean distance", f"{z} mean distance 95% CI low",
                                        f"{z} mean distance 95% CI high"))


def h2row(r):
    return (r["Heritability h2 (Falconer)"], r["h2 95% CI low (bootstrap percentile)"],
            r["h2 95% CI high (bootstrap percentile)"])


def f3(x):
    return f"{x:.3f}".replace("-", "−")


# =============================================================================
# STYLE
# =============================================================================
COL = {"FO": "#A8417F", "DT": "#E59A12", "AR": "#E0605A",
       "Weekly alcohol consumption": "#5940C9", "Drinks per drinking occasion": "#D9306A"}
SIG_RED = "#D00000"
INK = "black"
FS_TICK, FS_LABEL, FS_TITLE, FS_SECTION, FS_LETTER, FS_ANNOT = 11.5, 12, 13, 15, 19, 10.5

mpl.rcParams.update({
    "font.family": ["Arial", "Liberation Sans"], "font.size": FS_TICK, "text.color": INK,
    "axes.labelcolor": INK, "axes.edgecolor": INK, "axes.linewidth": 1.0,
    "axes.labelsize": FS_LABEL, "axes.titlesize": FS_TITLE, "axes.titleweight": "bold",
    "axes.titlepad": 4, "axes.spines.top": False, "axes.spines.right": False,
    "xtick.color": INK, "ytick.color": INK, "xtick.labelsize": FS_TICK, "ytick.labelsize": FS_TICK,
    "xtick.major.width": 1.0, "ytick.major.width": 1.0, "xtick.major.size": 3.5, "ytick.major.size": 3.5,
    "xtick.major.pad": 2.5, "ytick.major.pad": 2.5, "legend.fontsize": FS_TICK,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
})

W, H = 11.4, 16.2           # inches
fig = plt.figure(figsize=(W, H))


def ax_at(left, top, width, height, **kw):
    """Axes placed in inches, measured from the top-left corner of the figure."""
    return fig.add_axes([left / W, 1 - (top + height) / H, width / W, height / H], **kw)


def fig_text(x, y_top, s, **kw):
    fig.text(x / W, 1 - y_top / H, s, **kw)


def letter(x, y_top, s):
    fig_text(x, y_top, s, fontsize=FS_LETTER, fontweight="bold", va="top", ha="left")


def section(y_top, title):
    fig_text(0.05, y_top, title, fontsize=FS_SECTION, fontweight="bold", va="center", ha="left")
    fig.add_artist(Line2D([0.05 / W, (W - 0.05) / W], [1 - (y_top + 0.17) / H] * 2, color=INK, lw=0.9))


def tint(color, amount):
    rgb = np.array(mpl.colors.to_rgb(color))
    return tuple(rgb + (1 - rgb) * amount)


def p_style(p):
    sig = p < 0.05
    return dict(color=SIG_RED if sig else INK, fontweight="bold" if sig else "normal")


class Bracket:
    """Legend handle drawn as a significance bracket."""


class HandlerBracket(HandlerBase):
    def create_artists(self, legend, orig, x0, y0, width, height, fontsize, trans):
        xs = [x0 + 0.15 * width, x0 + 0.15 * width, x0 + 0.85 * width, x0 + 0.85 * width]
        ys = [y0 + 0.15 * height, y0 + 0.85 * height, y0 + 0.85 * height, y0 + 0.15 * height]
        return [Line2D(xs, ys, color=INK, lw=1.0, transform=trans)]


def legend_row(y_top, x_center, handles, labels):
    fig.legend(handles, labels, loc="center", ncol=len(handles), handler_map={Bracket: HandlerBracket()},
               bbox_to_anchor=(x_center / W, 1 - y_top / H), frameon=False, handlelength=1.5,
               handletextpad=0.4, columnspacing=1.3, borderaxespad=0, fontsize=FS_TICK)


MZ_MARK = dict(marker="o", ms=7.5, mec=INK, mew=0.8)
DZ_MARK = dict(marker="o", ms=7.5, mew=1.4)
# MZ / DZ key in black only (the same symbols in every panel)
KEY_POINTS = [Line2D([], [], ls="none", mfc=INK, **MZ_MARK), Line2D([], [], ls="none", mfc="white", mec=INK, **DZ_MARK),
              Line2D([], [], color=INK, lw=1.0, marker="|", ms=9, mew=1.0)]
KEY_BARS = [Patch(fc=INK, ec=INK, lw=0.8), Patch(fc="white", ec=INK, lw=1.0), Bracket()]


def mz_dz_pair(ax, x, mz, dz, color, offset=0.14):
    (m1, lo1, hi1), (m2, lo2, hi2) = mz, dz
    ax.plot([x - offset, x + offset], [m1, m2], color=color, lw=1.5, zorder=2)
    for xx, (m, lo, hi) in ((x - offset, mz), (x + offset, dz)):
        ax.errorbar(xx, m, yerr=[[m - lo], [hi - m]], color=INK, lw=1.1, capsize=2.8, capthick=1.1, zorder=3,
                    fmt="none")
    ax.plot(x - offset, m1, mfc=color, zorder=4, **MZ_MARK)
    ax.plot(x + offset, m2, mfc=tint(color, 0.75), mec=color, zorder=4, **DZ_MARK)


def bracket(ax, x1, x2, y, h):
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], color=INK, lw=0.9, clip_on=False)


def pq_rows(ax, xs, ps, qs, y_p=1.13, y_q=1.03, fs=FS_ANNOT, lab_x=None):
    """Two aligned rows of numbers above an axis: p and q, one value per x position."""
    for x, p, q in zip(xs, ps, qs):
        ax.text(x, y_p, f"{p:.3f}", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=fs,
                **p_style(p))
        ax.text(x, y_q, f"{q:.3f}", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=fs,
                **p_style(q))
    if lab_x is not None:
        for y, t in ((y_p, "p"), (y_q, "q")):
            ax.text(lab_x, y, t, transform=ax.get_xaxis_transform(), ha="right", va="bottom", fontsize=fs,
                    fontstyle="italic", fontweight="bold")


def forest(ax, rows, colors, labels, sigs, xlim, xticks, fs_num, pq=None):
    """Falconer h2 with 95% CI; numbers column on the right (header 'h² [95% CI]')."""
    n = len(rows)
    y = np.arange(n)[::-1]
    ax.axvline(0, color=INK, lw=0.9, ls=(0, (3, 2)), zorder=1)
    for yi, (h2, lo, hi), c, sig in zip(y, rows, colors, sigs):
        ax.plot([lo, hi], [yi, yi], color=INK, lw=1.3, zorder=2, solid_capstyle="butt")
        for xx in (lo, hi):
            ax.plot([xx, xx], [yi - 0.2, yi + 0.2], color=INK, lw=1.3)
        ax.plot(h2, yi, marker="s", ms=8, mfc=c, mec=INK, mew=0.8, zorder=3)
        dy = 0.2 if pq is not None else 0
        ax.text(1.04, yi + dy, f"{f3(h2)} [{f3(lo)}, {f3(hi)}]", transform=ax.get_yaxis_transform(),
                va="center", ha="left", fontsize=fs_num, **p_style(0 if sig else 1))
        if pq is not None:                                   # p and q of h2 (r_MZ > r_DZ) under the numbers
            p, q = pq[n - 1 - yi]
            t = ax.text(1.04, yi - 0.2, f"p = {p:.3f}", transform=ax.get_yaxis_transform(), va="center",
                        ha="left", fontsize=fs_num, **p_style(p))
            ax.annotate(f",  q = {q:.3f}", xycoords=t, xy=(1, 0.5), va="center", ha="left", fontsize=fs_num,
                        annotation_clip=False, **p_style(q))
    ax.text(1.04, n - 0.35, "h² [95% CI]", transform=ax.get_yaxis_transform(), va="bottom", ha="left",
            fontsize=fs_num, fontweight="bold")
    ax.set_yticks(y, labels)
    ax.tick_params(axis="y", length=0, pad=3)
    ax.spines["left"].set_visible(False)
    ax.set_xlim(*xlim)
    ax.set_xticks(xticks, [f"{t:g}".replace("-", "−") for t in xticks])
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_xlabel("h² (95% CI)", labelpad=2)


# =============================================================================
# SECTION 1 - BRAIN-STATE DYNAMICS
# =============================================================================
section(0.2, "Brain-state dynamics")

# ---- Panel A ----------------------------------------------------------------
letter(0.02, 0.48, "A")
fig_text(W / 2, 0.62, "Mean within-pair difference, |twin 1 − twin 2| (95% CI)", fontsize=FS_TITLE,
         fontweight="bold", ha="center", va="center")
A_LEFT, A_W, A_GAP, A_TOP, A_H = 1.05, 2.8, 0.95, 1.55, 2.15
for i, m in enumerate(METRICS):
    ax = ax_at(A_LEFT + i * (A_W + A_GAP), A_TOP, A_W, A_H)
    top = 0
    for s in range(4):
        r = row(m, s)
        mz, dz = dist(r, "MZ"), dist(r, "DZ")
        mz_dz_pair(ax, s, mz, dz, COL[m])
        top = max(top, mz[2], dz[2])
    pq_rows(ax, range(4), [row(m, s)[P_DIST] for s in range(4)], [row(m, s).q_dist for s in range(4)],
            lab_x=-0.62)
    ax.set_xlim(-0.5, 3.5)
    ax.set_xticks(range(4), STATE_TICK, rotation=30, ha="right", rotation_mode="anchor", linespacing=0.95)
    ax.set_ylim(0, top * 1.08)
    if m == "FO":
        ax.set_yticks([0, 0.02, 0.04, 0.06])

    elif m == "DT":
        ax.set_yticks([0.00, 0.25, 0.50, 0.75, 1.00])

    elif m == "AR":
        ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
    ax.set_title(METRIC_TITLES[m], pad=40)
    ax.set_ylabel(A_YLABEL[m], labelpad=3, linespacing=0.95)

legend_row(4.95, W / 2, KEY_POINTS, ["MZ mean", "DZ mean", "95% CI"])
fig_text(W / 2, 5.2, "p: MZ vs DZ (Mann–Whitney U, two-sided); q: FDR across 12 tests; red, < 0.05",
         fontsize=FS_TICK, ha="center", va="center")

# ---- Panel B ----------------------------------------------------------------
B_TOP, B_H = 6.3, 2.6
letter(0.02, 5.6, "B")
B_LEFT, B_WIDTH = 0.75, 5.25
axB = ax_at(B_LEFT, B_TOP, B_WIDTH, B_H)
fig_text(B_LEFT + B_WIDTH / 2, 5.72, "Cross-twin correlation (ICC)", fontsize=FS_TITLE, fontweight="bold",
         ha="center", va="center")
bw, block_gap = 0.4, 0.7
centers, metric_mid = [], []
for mi, m in enumerate(METRICS):
    xs = np.arange(4) + mi * (4 + block_gap)
    metric_mid.append(xs.mean())
    for s, x in enumerate(xs):
        r = row(m, s)
        mz, dz = r["r_MZ (ICC)"], r["r_DZ (ICC)"]
        axB.bar(x - bw / 2, mz, bw, color=COL[m], edgecolor=INK, lw=0.7, zorder=2)
        axB.bar(x + bw / 2, dz, bw, color=tint(COL[m], 0.65), edgecolor=COL[m], lw=1.1, zorder=2)
        for xx, v in ((x - bw / 2, mz), (x + bw / 2, dz)):
            axB.text(xx, max(v, 0) + 0.012, f3(v), ha="center", va="bottom", fontsize=FS_ANNOT - 1.5,
                     rotation=90)
        bracket(axB, x - bw / 2, x + bw / 2, max(mz, dz) + 0.13, 0.025)
        centers.append(x)
    if mi < 2:
        axB.axvline(xs[-1] + (1 + block_gap) / 2, color=INK, lw=0.8, ls=(0, (2, 2)), zorder=0)
pq_rows(axB, centers, [row(m, s)[P_ICC] for m in METRICS for s in range(4)],
        [row(m, s).q_icc for m in METRICS for s in range(4)], y_p=1.09, y_q=1.01, fs=FS_ANNOT - 1.5,
        lab_x=-0.75)
axB.axhline(0, color=INK, lw=0.9, zorder=3)
axB.set_xlim(-0.6, centers[-1] + 0.6)
axB.set_ylim(-0.05, 0.85)
axB.set_yticks([0, 0.2, 0.4, 0.6, 0.8], ["0", "0.2", "0.4", "0.6", "0.8"])
axB.set_xticks(centers, STATE_TICK * 3, fontsize=FS_TICK - 1.5, rotation=55, ha="right", rotation_mode="anchor",
               linespacing=0.9)
axB.tick_params(axis="x", length=0, pad=3)
axB.spines["bottom"].set_visible(False)
axB.set_ylabel("Intraclass correlation (r)", labelpad=3)
for m, xm in zip(METRICS, metric_mid):
    axB.text(xm, -0.62, {"FO": "Fractional occupancy (FO)", "DT": "Dwell time (DT)",
                         "AR": "Appearance rate (AR)"}[m],
             transform=axB.get_xaxis_transform(), ha="center", va="top", fontsize=FS_LABEL, fontweight="bold")
legend_row(10.85, B_LEFT + B_WIDTH / 2, KEY_BARS, ["MZ ICC (r$_{MZ}$)", "DZ ICC (r$_{DZ}$)", "r$_{MZ}$ > r$_{DZ}$"])
fig_text(B_LEFT + B_WIDTH / 2, 11.1, "p: one-sided r$_{MZ}$ > r$_{DZ}$; q: FDR across 12 tests; red, < 0.05",
         fontsize=FS_TICK, ha="center", va="center")

# ---- Panel C ----------------------------------------------------------------
letter(6.2, 5.6, "C")
C_LEFT, C_WIDTH = 7.75, 1.6
axC = ax_at(C_LEFT, B_TOP - 0.15, C_WIDTH, B_H + 1.3)
fig_text(8.6, 5.72, "Falconer heritability (h²)", fontsize=FS_TITLE, fontweight="bold", ha="center", va="center")
rowsC = [h2row(row(m, s)) for m in METRICS for s in range(4)]
sigC = [row(m, s)[P_ICC] < 0.05 for m in METRICS for s in range(4)]
forest(axC, rowsC, [COL[m] for m in METRICS for s in range(4)],
       [f"{st} {m}" for m in METRICS for st in STATE_TICK], sigC, xlim=(-0.65, 1.7), xticks=[-0.5, 0, 0.5, 1, 1.5, 2.0],
       fs_num=FS_ANNOT - 1.5, pq=[(row(m, s)[P_ICC], row(m, s).q_icc) for m in METRICS for s in range(4)])
axC.tick_params(axis="y", labelsize=FS_TICK - 2.5)
for yb in (3.5, 7.5):
    axC.axhline(yb, color=INK, lw=0.8, ls=(0, (2, 2)), zorder=0, xmin=0, clip_on=False)
fig_text(8.9, 10.85, "p: h² > 0 (r$_{MZ}$ > r$_{DZ}$, as in B); red, < 0.05", fontsize=FS_TICK - 1,
         ha="center", va="center")

# =============================================================================
# SECTION 2 - ALCOHOL CONSUMPTION
# =============================================================================
section(11.55, "Alcohol consumption")
ROW_TOP, ROW_H = 12.75, 2.1

# ---- Panels D, E ------------------------------------------------------------
DE_LEFT, DE_W, DE_GAP = 1.05, 1.15, 1.15
for k, (name, letter_s) in enumerate(zip(ALC, "DE")):
    key, col = ALC_KEY[name]
    r, sd = alc.loc[key], ms[col].std()
    left = DE_LEFT + k * (DE_W + DE_GAP)
    ax = ax_at(left, ROW_TOP, DE_W, ROW_H)
    letter(0.02 if k == 0 else left - 0.95, 11.85, letter_s)
    fig_text(left + DE_W / 2, 12.0, "Weekly alcohol\nconsumption" if k == 0 else "Drinks per\ndrinking occasion",
             fontsize=FS_TITLE - 1, fontweight="bold", ha="center", va="center", linespacing=0.95)
    mz, dz = dist(r, "MZ", sd), dist(r, "DZ", sd)
    c = COL[name]
    ax.plot([0, 1], [mz[0], dz[0]], color=c, lw=1.6, zorder=2)
    for x, v in ((0, mz), (1, dz)):
        ax.errorbar(x, v[0], yerr=[[v[0] - v[1]], [v[2] - v[0]]], color=INK, lw=1.1, capsize=3, capthick=1.1,
                    fmt="none", zorder=3)
        ax.text(x + 0.17, v[0], f3(v[0]), va="center", ha="left", fontsize=FS_ANNOT, fontweight="bold")
    ax.plot(0, mz[0], mfc=c, zorder=4, **{**MZ_MARK, "ms": 9})
    ax.plot(1, dz[0], mfc=tint(c, 0.75), mec=c, zorder=4, **{**DZ_MARK, "ms": 9})
    pq_rows(ax, [0.5], [r[P_DIST]], [r.q_dist], y_p=1.13, y_q=1.03, lab_x=-0.15)
    ax.set_xlim(-0.4, 1.9)
    ax.set_xticks([0, 1], ["MZ", "DZ"], fontsize=FS_LABEL)
    ax.set_ylim(0, 2.0)
    ax.set_yticks([0, 0.5, 1.0, 1.5, 2.0], ["0", "0.5", "1.0", "1.5", "2.0"])
    ax.set_ylabel("Mean |twin 1 − twin 2|\n(z-score)", labelpad=3, linespacing=0.95)
legend_row(15.55, DE_LEFT + DE_W + DE_GAP / 2 - 0.2, KEY_POINTS, ["MZ mean", "DZ mean", "95% CI"])
fig_text(DE_LEFT + DE_W + DE_GAP / 2 - 0.2, 15.8, "p: Mann\u2013Whitney U (two-sided); q: FDR across 2", fontsize=FS_TICK - 1,
         ha="center", va="center")

# ---- Panel F ----------------------------------------------------------------
F_LEFT, F_W = 5.8, 1.45
letter(5.0, 11.85, "F")
axF = ax_at(F_LEFT, ROW_TOP, F_W, ROW_H)
fig_text(F_LEFT + F_W / 2, 12.0, "Cross-twin ICC", fontsize=FS_TITLE - 1, fontweight="bold", ha="center",
         va="center")
bwF = 0.38
for k, name in enumerate(ALC):
    r, c = alc.loc[ALC_KEY[name][0]], COL[name]
    mz, dz = r["r_MZ (ICC)"], r["r_DZ (ICC)"]
    axF.bar(k - bwF / 2, mz, bwF, color=c, edgecolor=INK, lw=0.7, zorder=2)
    axF.bar(k + bwF / 2, dz, bwF, color=tint(c, 0.65), edgecolor=c, lw=1.1, zorder=2)
    for xx, v in ((k - bwF / 2, mz), (k + bwF / 2, dz)):
        axF.text(xx, v + 0.015, f3(v), ha="center", va="bottom", fontsize=FS_ANNOT - 1, rotation=90)
    bracket(axF, k - bwF / 2, k + bwF / 2, max(mz, dz) + 0.17, 0.03)
pq_rows(axF, [0, 1], [alc.loc[ALC_KEY[n][0]][P_ICC] for n in ALC], [alc.loc[ALC_KEY[n][0]].q_icc for n in ALC],
        y_p=1.13, y_q=1.03, lab_x=-0.55)
axF.set_xlim(-0.6, 1.6)
axF.set_ylim(0, 0.85)
axF.set_yticks([0, 0.2, 0.4, 0.6, 0.8], ["0", "0.2", "0.4", "0.6", "0.8"])
axF.set_xticks([0, 1], ["Weekly\nalcohol\nconsumption", "Drinks per\ndrinking\noccasion"], fontsize=FS_TICK - 2.5,
               linespacing=0.95)
axF.tick_params(axis="x", length=0)
axF.set_ylabel("Intraclass correlation (r)", labelpad=3)
legend_row(15.55, F_LEFT + F_W / 2 - 0.15, KEY_BARS, ["MZ", "DZ", "r$_{MZ}$ > r$_{DZ}$"])
fig_text(F_LEFT + F_W / 2 - 0.15, 15.8, "p: one-sided; q: FDR across 2", fontsize=FS_TICK - 1, ha="center",
         va="center")

# ---- Panel G ----------------------------------------------------------------
G_LEFT, G_W = 8.1, 1.35
letter(7.35, 11.85, "G")
axG = ax_at(G_LEFT, ROW_TOP, G_W, ROW_H)
fig_text(G_LEFT + 0.65, 12.0, "Falconer h²", fontsize=FS_TITLE - 1, fontweight="bold", ha="center", va="center")
rowsG = [h2row(alc.loc[ALC_KEY[n][0]]) for n in ALC]
forest(axG, rowsG, [COL[n] for n in ALC], ["Weekly", "Per\noccasion"],
       [alc.loc[ALC_KEY[n][0]][P_ICC] < 0.05 for n in ALC], xlim=(-0.65, 1.15), xticks=[-0.5, 0, 0.5, 1, 1.5],
       fs_num=FS_ANNOT - 1.5)
for j, n in enumerate(ALC):               # p and q under each estimate (test as in F)
    r = alc.loc[ALC_KEY[n][0]]
    yj = 1 - j
    axG.text(1.04, yj - 0.28, f"p = {r[P_ICC]:.3f}", transform=axG.get_yaxis_transform(), va="center",
             ha="left", fontsize=FS_ANNOT - 1.5, **p_style(r[P_ICC]))
    axG.text(1.04, yj - 0.52, f"q = {r.q_icc:.3f}", transform=axG.get_yaxis_transform(), va="center",
             ha="left", fontsize=FS_ANNOT - 1.5, **p_style(r.q_icc))
axG.tick_params(axis="y", labelsize=FS_TICK)

# =============================================================================
OUT.mkdir(parents=True, exist_ok=True)
stem = OUT / "Figure2"
fig.savefig(f"{stem}.pdf")
fig.savefig(f"{stem}.png", dpi=480)          # 5472 x 7776 px (<= 8000 px)
print(f"Saved {stem}.pdf/.png")
