"""Build Supplementary Figure S2 (A-D) and Table Sx from a consensus-metrics .mat file.

Works with either
  * consensus_metrics_corrected_2026_11_12.mat  (MATLAB script output), or
  * the Python re-run export (same field names).
Usage:  python make_figS2_2026_11_12.py <metrics.mat> <out_dir> [null_metrics.mat]
"""
import sys, os
import numpy as np
import scipy.io as sio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import linkage, leaves_list, fcluster
from scipy.spatial.distance import squareform

INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
HL = {3: "#2a78d6", 4: "#eb6834", 5: "#1baf7a"}        # categorical slots 1-3 (fixed order)
NULL_C = "#8a8984"
SEQ = ["#ffffff", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
CMAP_CONS = LinearSegmentedColormap.from_list("cons", SEQ)
GRAYS = LinearSegmentedColormap.from_list("grays", ["#d6d5d0", "#3d3c39"])

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.edgecolor": INK2, "axes.labelcolor": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.8, "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 600,
})


def load(path):
    d = sio.loadmat(path, squeeze_me=True)
    r = {k: d[k] for k in d if not k.startswith("__")}
    r["K_range"] = np.atleast_1d(r["K_range"]).astype(int)
    for k in ["PAC", "CDF_area", "DeltaK", "meanC", "sdC"]:
        r[k] = np.atleast_1d(r[k]).astype(float)
    r["cdf_x"] = [np.atleast_1d(x).astype(float) for x in r["cdf_x"]]
    r["cdf_F"] = [np.atleast_1d(x).astype(float) for x in r["cdf_F"]]
    return r


def style_k_axis(ax, K):
    ax.set_xlim(K.min() - 0.5, K.max() + 0.5)
    ax.set_xticks(np.arange(2, K.max() + 1, 2))
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.set_xlabel("Number of clusters (K)")


def line_with_highlights(ax, K, y, label):
    ax.plot(K, y, color=INK, lw=1.4, zorder=2, label=label)
    ax.scatter(K, y, s=14, color=INK, zorder=3, edgecolor="white", linewidth=0.6)
    for k, c in HL.items():
        i = int(np.where(K == k)[0][0])
        ax.scatter([k], [y[i]], s=44, color=c, zorder=4, edgecolor="white", linewidth=1.0)


def panel_A(ax, r, null=None):
    K = r["K_range"]
    if null is not None:
        ax.plot(null["K_range"], null["PAC"], color=NULL_C, lw=1.2, ls="--", zorder=1,
                label="Null reference (ROI-permuted)")
    line_with_highlights(ax, K, r["PAC"], "Observed")
    style_k_axis(ax, K)
    ax.set_ylabel("PAC (0.1 < consensus < 0.9)")
    top = np.nanmax(r["PAC"]) if null is None else max(np.nanmax(r["PAC"]), np.nanmax(null["PAC"]))
    ax.set_ylim(0, np.ceil(top * 1.25 * 20) / 20)
    for k, dy in ((4, -11), (5, 7)):
        i = int(np.where(K == k)[0][0])
        ax.annotate(f"{r['PAC'][i]:.3f}", (k, r["PAC"][i]), xytext=(6, dy), textcoords="offset points",
                    fontsize=6.5, color=INK, ha="left", va="center")
    ax.set_title("A  Proportion of ambiguous clustering", loc="left", fontweight="bold")
    handles = [Line2D([], [], color=INK, lw=1.4, marker="o", ms=3.5, label="Observed")]
    if null is not None:
        handles.append(Line2D([], [], color=NULL_C, lw=1.2, ls="--", label="Null reference (ROI-permuted)"))
    handles += [Line2D([], [], ls="", marker="o", ms=6, color=c, label=f"K = {k}") for k, c in HL.items()]
    ax.legend(handles=handles, frameon=False, loc="upper right", fontsize=7)


def panel_B(ax, r, cax=None):
    K = r["K_range"]
    norm = plt.Normalize(K.min(), K.max())
    for i, k in enumerate(K):
        if k in HL:
            continue
        x = np.r_[0, r["cdf_x"][i]]; F = np.r_[0, r["cdf_F"][i]]
        ax.step(x, F, where="post", color=GRAYS(norm(k)), lw=0.7, zorder=1)
    for k, c in HL.items():
        i = int(np.where(K == k)[0][0])
        x = np.r_[0, r["cdf_x"][i]]; F = np.r_[0, r["cdf_F"][i]]
        ax.step(x, F, where="post", color=c, lw=1.8, zorder=3, label=f"K = {k}")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.01)
    ax.set_xlabel("Consensus index"); ax.set_ylabel("Cumulative fraction of pairs")
    ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True)
    ax.set_title("B  Empirical CDF of consensus, K = 2–20", loc="left", fontweight="bold")
    ax.legend(frameon=False, loc="lower right", fontsize=7)
    if cax is not None:
        sm = plt.cm.ScalarMappable(cmap=GRAYS, norm=norm)
        cb = plt.colorbar(sm, cax=cax)
        cb.set_label("K (grey lines)", fontsize=7); cb.set_ticks([2, 5, 10, 15, 20]); cb.outline.set_visible(False)


def panel_C(ax, r, null=None):
    K = r["K_range"]
    if null is not None:
        ax.plot(null["K_range"], null["DeltaK"], color=NULL_C, lw=1.2, ls="--", zorder=1)
    line_with_highlights(ax, K, r["DeltaK"], "ΔK")
    style_k_axis(ax, K)
    ax.set_ylabel("ΔK (relative change in CDF area)")
    ax.set_yscale("log")
    ax.set_title("C  ΔK profile", loc="left", fontweight="bold")
    i2 = int(np.where(K == 2)[0][0])
    ax.annotate("K = 2: ΔK = A(2)", (K[i2], r["DeltaK"][i2]), xytext=(8, -2), textcoords="offset points",
                fontsize=6.5, color=INK2, va="center")
    if null is not None:
        ax.text(0.04, 0.05, "Null-reference ΔK (dashed)\nlies under the observed curve", transform=ax.transAxes,
                fontsize=6.5, color=INK2, ha="left", va="bottom")
        ax.legend(handles=[Line2D([], [], color=INK, lw=1.4, marker="o", ms=3.5, label="Observed"),
                           Line2D([], [], color=NULL_C, lw=1.2, ls="--", label="Null reference")],
                  frameon=False, loc="upper right", fontsize=7)


def ordered_consensus(C, K):
    D = 1.0 - C
    np.fill_diagonal(D, 0.0)
    D = (D + D.T) / 2
    Z = linkage(squareform(D, checks=False), method="average")
    order = leaves_list(Z)
    cl = fcluster(Z, K, criterion="maxclust")[order]
    return C[np.ix_(order, order)], cl


def panel_D(axes, r, Ks=(3, 4, 5), cax=None):
    show = list(np.atleast_1d(r["K_show"]).astype(int))
    subC = r["subC"]
    for ax, k in zip(axes, Ks):
        s = show.index(k)
        C = np.asarray(subC[:, :, s], dtype=float)
        Co, cl = ordered_consensus(C, k)
        im = ax.imshow(Co, cmap=CMAP_CONS, vmin=0, vmax=1, interpolation="nearest", rasterized=True)
        b = np.where(np.diff(cl) != 0)[0] + 0.5
        for x in b:
            ax.axhline(x, color=HL.get(k, INK), lw=0.6); ax.axvline(x, color=HL.get(k, INK), lw=0.6)
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(True); sp.set_color(HL.get(k, INK)); sp.set_linewidth(1.5)
        ax.set_title(f"K = {k}", color=INK, fontsize=9)
    axes[0].text(0, 1.14, "D  Consensus matrices", transform=axes[0].transAxes, fontsize=9,
                 fontweight="bold", ha="left", va="bottom", color=INK)
    if cax is not None:
        cb = plt.colorbar(im, cax=cax); cb.set_label("Consensus index", fontsize=7); cb.outline.set_visible(False)


def main(path, out, null_path=None):
    os.makedirs(out, exist_ok=True)
    r = load(path)
    null = load(null_path) if null_path else None
    K = r["K_range"]
    n_sub = np.atleast_1d(r["subIdx"]).size

    # ---- Combined figure S2 (A-D) ----
    fig = plt.figure(figsize=(7.2, 6.3))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.0], hspace=0.42)
    top = gs[0].subgridspec(1, 3, width_ratios=[1, 1, 1], wspace=0.42)
    axA = fig.add_subplot(top[0]); axB = fig.add_subplot(top[1]); axC = fig.add_subplot(top[2])
    bot = gs[1].subgridspec(1, 4, width_ratios=[1, 1, 1, 0.05], wspace=0.08)
    axD = [fig.add_subplot(bot[i]) for i in range(3)]; caxD = fig.add_subplot(bot[3])
    panel_A(axA, r, null); panel_B(axB, r); panel_C(axC, r, null); panel_D(axD, r, cax=caxD)
    for a in (axA, axC):
        if a.get_legend():
            a.get_legend().remove()
    axA.set_title("A  PAC", loc="left", fontweight="bold"); axB.set_title("B  CDF, K = 2–20", loc="left", fontweight="bold")
    hl = [Line2D([], [], color=c, lw=2, marker="o", ms=4, label=f"K = {k}") for k, c in HL.items()]
    hl.append(Line2D([], [], color=GRAYS(0.5), lw=1, label="Other K (light→dark = 2→20)"))
    if null is not None:
        hl.append(Line2D([], [], color=NULL_C, lw=1.2, ls="--", label="Null reference (A, C)"))
    fig.legend(handles=hl, loc="upper center", ncol=len(hl), frameon=False, fontsize=7, bbox_to_anchor=(0.5, 0.955))
    fig.text(0.5, 0.06, f"D: random subset of {n_sub:,} prototypes, ordered by average-linkage clustering of 1 − consensus; lines mark the K consensus clusters.",
             ha="center", fontsize=6.5, color=INK2)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(out, f"2026-11-12-FigS2_combined.{ext}"), bbox_inches="tight")
    plt.close(fig)

    # ---- Individual panels ----
    f, ax = plt.subplots(figsize=(3.5, 2.8)); panel_A(ax, r, null)
    for ext in ("png", "pdf"): f.savefig(os.path.join(out, f"2026-11-12-FigS2A_PAC.{ext}"), bbox_inches="tight")
    plt.close(f)
    f = plt.figure(figsize=(3.8, 2.8)); g = f.add_gridspec(1, 2, width_ratios=[1, 0.04], wspace=0.08)
    ax = f.add_subplot(g[0]); cax = f.add_subplot(g[1]); panel_B(ax, r, cax=cax)
    for ext in ("png", "pdf"): f.savefig(os.path.join(out, f"2026-11-12-FigS2B_CDF.{ext}"), bbox_inches="tight")
    plt.close(f)
    f, ax = plt.subplots(figsize=(3.5, 2.8)); panel_C(ax, r, null)
    for ext in ("png", "pdf"): f.savefig(os.path.join(out, f"2026-11-12-FigS2C_DeltaK.{ext}"), bbox_inches="tight")
    plt.close(f)
    f = plt.figure(figsize=(7.2, 2.6)); g = f.add_gridspec(1, 4, width_ratios=[1, 1, 1, 0.05], wspace=0.08)
    axs = [f.add_subplot(g[i]) for i in range(3)]; panel_D(axs, r, cax=f.add_subplot(g[3]))
    for ext in ("png", "pdf"): f.savefig(os.path.join(out, f"2026-11-12-FigS2D_consensus_matrices.{ext}"), bbox_inches="tight")
    plt.close(f)

    # ---- Table Sx ----
    lines = ["K,PAC,CDF_area,DeltaK,Mean_consensus,SD_consensus" +
             (",PAC_null_reference,CDF_area_null_reference,DeltaK_null_reference" if null is not None else "")]
    for i, k in enumerate(K):
        row = f"{k},{r['PAC'][i]:.4f},{r['CDF_area'][i]:.4f},{r['DeltaK'][i]:.4f},{r['meanC'][i]:.4f},{r['sdC'][i]:.4f}"
        if null is not None:
            row += f",{null['PAC'][i]:.4f},{null['CDF_area'][i]:.4f},{null['DeltaK'][i]:.4f}"
        lines.append(row)
    open(os.path.join(out, "2026-11-12-TableSx_consensus_metrics.csv"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else None)
