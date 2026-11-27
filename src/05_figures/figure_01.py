#!/usr/bin/env python3

"""

Figure 1 - co-activation patterns (CAPs) of the twin sample.



  A    pipeline strip (whole-brain activation -> k-means, k = 4 -> CAP 1-4)

  B-E  one row per CAP: inflated surfaces (left/right, lateral/medial), sagittal and

       coronal sections, and next to them a 17-network radar plot of the MZ and DZ CAPs

  F-H  fractional occupancy, appearance rate and mean dwell time of each CAP in MZ and

       DZ twins (violins drawn from the master sheet; p = Welch's t-test, MZ vs DZ)



Inputs (all in inputs/):

  centroids/state_1..4.nii.gz            pooled (full-cohort) CAP centroids = CAP 1-4

  centroids/centroids-mono_4D.nii.gz     MZ k-means centroids (4 volumes, own label order)

  centroids/centroids-di_4D.nii.gz       DZ k-means centroids (4 volumes, own label order)

  2026-11-13-Figure1-Yeo17-correlations.xlsx   Pearson r of each MZ / DZ centroid with each

                                              Yeo-17 network mask (sheets MZ, DZ, Network names)

  2026-11-13-master-sheet.xlsx           FO, DT and AR of every participant (violins)

  figure1-pipeline-images/*.png          the two pictures of panel A

Steps:

  1. the MZ and DZ centroid volumes are matched to CAP 1-4 by spatial correlation with the

     pooled maps (k-means label numbers differ between runs)

  2. maps are shown from their native 2-mm grid (the 1-mm files repeat each value in 2x2x2

     blocks), thresholded at |value| > 2334.09 a.u., colour scale up to +/-7500 a.u.

  3. violins: kernel density of each group's values, drawn over the observed range,

     median marked; p-values: two-sided Welch's t-test, MZ (n = 160) vs DZ (n = 80)

Output: outputs/Figure1.pdf (vector text), outputs/Figure1.tif (1200 dpi), outputs/Figure1.png

Run   : python 2026-11-13-Figure1.py

"""

import glob

import inspect

import warnings

from pathlib import Path

import os

import matplotlib as mpl

import matplotlib.pyplot as plt

import nibabel as nib

import numpy as np

import pandas as pd

from matplotlib.backends.backend_agg import FigureCanvasAgg

from matplotlib.colors import LinearSegmentedColormap, Normalize

from matplotlib.figure import Figure

from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from nilearn import datasets, image, plotting, surface

from PIL import Image

from scipy import stats



warnings.filterwarnings("ignore", category=FutureWarning)

warnings.filterwarnings("ignore", category=DeprecationWarning)



# =============================================================================

# SETTINGS

# =============================================================================

INPUTS = Path("inputs")

STATE_DIR = Path(r"C:\Users\Tanya\Desktop\Figures and Tables\Excels\centroids")

state_files = sorted(glob.glob(os.path.join(STATE_DIR, "state_*.nii.gz")))



for file in state_files:

    print(file)



print(f"\nTotal files found: {len(state_files)}")

MZ_CENTROIDS = Path(r"C:\Users\Tanya\Desktop\Figures and Tables\Excels\centroids\centroids-mono_4D.nii.gz")

DZ_CENTROIDS = Path(r"C:\Users\Tanya\Desktop\Figures and Tables\Excels\centroids\centroids-di_4D.nii.gz")

YEO_SHEET = Path(r"C:\Users\Tanya\Desktop\Figures and Tables\Excels\Yeo17_network_correlations_figure1.xlsx")

MASTER = Path(r"C:\Users\Tanya\Desktop\Figures and Tables\Excels\Master-sheet.xlsx")

PIPE_WHOLEBRAIN = Path(r"C:\Users\Tanya\Desktop\Figures and Tables\Excels\figure1-pipeline-images\pipeline_wholebrain.png")

PIPE_KMEANS = Path(r"C:\Users\Tanya\Desktop\Figures and Tables\Excels\figure1-pipeline-images\pipeline_kmeans.png")

OUT_STEM = Path(r"C:\Users\Tanya\Desktop\Figures and Tables\Figure1")

DPI = 580                       # raster export (text and lines stay vector in the PDF)



CAP_NAMES = ["Visual Network-A", "DMN", "Visual Network-B", "Hippocampal–Temporoparietal"]

# violin panels: y-range and ticks (same as the published figure)

VIOLIN_AXES = {"FO": ((0.0267, 0.4265), [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4]),

               "AR": ((2.6047, 7.8198), [1, 2, 3, 4, 5, 6, 7, 8]),

               "DT": ((0.5728, 5.3608), [1, 2, 3, 4, 5])}

YLABELS = {"FO": "Fractional occupancy\n(proportion)", "AR": "Appearance rate\n(per min)",

           "DT": "Mean dwell time (s)"}

VIOLIN_POS, VIOLIN_WIDTH = (-0.17, 0.18), 0.31       # MZ / DZ offset from each CAP, max width



# Network names from the authors' table (sheet "Network names"), with corrected atlas labels.
# NAME_OVERRIDE changes labels only; correlation values stay with their IDs.

NAME_OVERRIDE = {10: "LIM-B", 11: "FPN-C", 12: "FPN-A", 13: "FPN-B", 14: "TempPar", 15: "DMN-C", 16: "DMN-A", 17: "DMN-B"}

# spoke order: grouped by family (IDs) - LIM-A/B = 9, 10; FPN-A/B/C = 11, 12, 13

SPOKE_ORDER = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]



THRESH, VMAX = 2334.09, 7500.0   # shared display cut-off / scale (as Fig. S4B)

SAG_X, COR_Y = -2, -54

R_MIN, R_MAX, R_RINGS = -0.5, 0.6, [-0.4, -0.2, 0.0, 0.2, 0.4]

MZ_COL, DZ_COL = "#2a78d6", "#eb6834"      # MZ blue, DZ orange - the same colours in Figures 1-3



FS_LETTER, FS_LABEL, FS_SMALL, FS_TINY = 21, 13.5, 12, 11.5   # extra large (journal version: 10, 7, 6.5, 6)

FS_SPOKE, FS_RING = 10, 10   # all 17 full network names fit at 8 pt

FS_PVAL = 10.5                  # p-values (APA style, one row above the violins)

VIOLIN_FILL = [MZ_COL, DZ_COL]          # MZ, DZ

INK, INK_SOFT, GREY = "#000000", "#000000", "#d9d9d9"   # GREY only = sub-threshold colour on the brain maps

MINUS, MM = "−", 1 / 25.4

plt.rcParams.update({

    "font.family": "sans-serif",

    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],

    "font.size": FS_LABEL, "font.weight": "bold", "axes.labelweight": "bold",

    "text.color": INK, "axes.edgecolor": INK, "axes.linewidth": 0.5,

    "xtick.labelsize": FS_TINY, "ytick.labelsize": FS_TINY,

    "axes.unicode_minus": True, "pdf.fonttype": 42, "ps.fonttype": 42,

    "figure.facecolor": "white", "savefig.facecolor": "white",

})

CROP = {"x": (-74, 74), "y": (-109, 75), "z": (-74, 84)}

LIMITS = {"x": (CROP["y"], CROP["z"]), "y": (CROP["x"], CROP["z"])}





# =============================================================================

# helpers

# =============================================================================

def fmt(v, nd=0):

    if nd and round(v, nd) == 0:

        v = 0.0

    return (f"{v:.{nd}f}" if nd else f"{v:g}").replace("-", MINUS)





def signed_cmap(t, vmax):

    pos, neg = ["#c0161b", "#f0601e", "#ffe14d"], ["#1c4fc9", "#3ea0f5", "#a8edff"]

    f, e = t / vmax, 1e-4

    lo, hi = 0.5 - f / 2, 0.5 + f / 2

    cm = LinearSegmentedColormap.from_list("signed", [

        (0, neg[2]), (lo / 2, neg[1]), (lo - e, neg[0]), (lo, GREY), (hi, GREY),

        (hi + e, pos[0]), ((1 + hi) / 2, pos[1]), (1, pos[2])], N=1024)

    cm.set_over(pos[2]); cm.set_under(neg[2])

    return cm





def render_surface(tex, hemi, view, fs, cmap):

    f = Figure(figsize=(4, 3), dpi=600)

    FigureCanvasAgg(f)

    ax = f.add_subplot(111, projection="3d")

    kw = dict(hemi=hemi, view=view, bg_map=fs[f"sulc_{hemi}"], threshold=THRESH, vmax=VMAX,

              cmap=cmap, colorbar=False, symmetric_cbar=True, axes=ax, figure=f)

    if "vmin" in inspect.signature(plotting.plot_surf_stat_map).parameters:

        kw["vmin"] = -VMAX

    plotting.plot_surf_stat_map(fs[f"infl_{hemi}"], tex, **kw)

    f.canvas.draw()

    rgba = np.asarray(f.canvas.buffer_rgba()).copy()

    ink = (rgba[..., :3] < 245).any(-1)

    r, c = np.where(ink.any(1))[0], np.where(ink.any(0))[0]

    rgba = rgba[max(r[0] - 6, 0):r[-1] + 6, max(c[0] - 6, 0):c[-1] + 6].copy()

    rgb = rgba[..., :3].astype(float)

    grey = (rgb.max(-1) - rgb.min(-1)) < 14

    rgb[grey] = 255 - 0.6 * (255 - rgb[grey])

    rgba[..., :3] = rgb.astype(np.uint8)

    return rgba





class MMFigure:

    def __init__(self, w, h):

        self.W, self.H = w, h

        self.fig = plt.figure(figsize=(w * MM, h * MM))



    def ax(self, x, y, w, h, **kw):

        return self.fig.add_axes([x / self.W, 1 - (y + h) / self.H, w / self.W, h / self.H], **kw)



    def text(self, x, y, s, **kw):

        kw.setdefault("color", INK)

        return self.fig.text(x / self.W, 1 - y / self.H, s, **kw)





def letter(mf, x, y, s):

    mf.text(x, y, s, ha="left", va="top", fontsize=FS_LETTER, fontweight="bold", color="black")





# =============================================================================

# main

# =============================================================================

def main():

    # ---- data ------------------------------------------------------------------

    files = [sorted(glob.glob(str(STATE_DIR / f"*state_{k}.nii*")))[0] for k in (1, 2, 3, 4)]

    ref = nib.load(files[0])

    maps = [nib.Nifti1Image(nib.load(f).get_fdata(dtype=np.float32), ref.affine) for f in files]

    # native 2-mm grid: each value is repeated in 2x2x2 blocks with phase (0, 1, 0)

    to2 = np.array([[2, 0, 0, 0.5], [0, 2, 0, 1.5], [0, 0, 2, 0.5], [0, 0, 0, 1.0]])

    maps2 = [nib.Nifti1Image(np.asarray(m.dataobj)[0::2, 1::2, 0::2].copy(), ref.affine @ to2)

             for m in maps]

    tmpl = datasets.load_mni152_template(resolution=1)



    def slab(img, mode, cut, interp, step=0.25):

        """Fine 2-D slice (3 planes thick) through world coordinate `cut`."""

        (h0, h1), (v0, v1) = LIMITS[mode]

        n_h, n_v = int((h1 - h0) / step) + 1, int((v1 - v0) / step) + 1

        aff = np.eye(4) * step; aff[3, 3] = 1

        if mode == "x":

            aff[:3, 3] = [cut - step, h0, v0]; shape = (3, n_h, n_v)

        else:

            aff[:3, 3] = [h0, cut - step, v0]; shape = (n_h, 3, n_v)

        return image.resample_img(img, target_affine=aff, target_shape=shape, interpolation=interp,

                                  force_resample=True, copy_header=True)

    Path(OUT_STEM).parent.mkdir(exist_ok=True)



    # step 1: match the MZ / DZ k-means centroid volumes to CAP 1-4 (spatial correlation)

    pooled = [np.asarray(m.dataobj)[0::2, 1::2, 0::2].astype(float) for m in maps]

    in_brain = np.zeros(pooled[0].shape, bool)

    for p_ in pooled:

        in_brain |= p_ != 0

    match_rows, cols = [], {}

    for grp, path in (("MZ", MZ_CENTROIDS), ("DZ", DZ_CENTROIDS)):

        img4 = nib.load(path)

        r = np.array([[np.corrcoef(np.asarray(img4.dataobj[..., v])[0::2, 1::2, 0::2].astype(float)[in_brain],

                                   pooled[c][in_brain])[0, 1] for c in range(4)] for v in range(4)])

        best = r.argmax(axis=0)                                   # centroid volume for each CAP

        assert len(set(best)) == 4, f"{grp}: centroids do not match the CAPs one-to-one"

        cols[grp] = [int(b) + 1 for b in best]

        for c in range(4):

            match_rows.append({"Group": grp, "CAP": c + 1, "Centroid volume": cols[grp][c],

                               "Spatial r": round(r[best[c], c], 4)})

    print(pd.DataFrame(match_rows).to_string(index=False))



    # Yeo-17 network correlations (authors' values) in CAP order

    mz = pd.read_excel(YEO_SHEET, sheet_name="MZ").set_index("Network_ID")

    dz = pd.read_excel(YEO_SHEET, sheet_name="DZ").set_index("Network_ID")

    names = pd.read_excel(YEO_SHEET, sheet_name="Network names").set_index("Network_ID").Network_Name

    names = {i: NAME_OVERRIDE.get(i, n) for i, n in names.items()}

    R_mz = np.array([[mz.loc[i, f"Centroid_{cols['MZ'][c]}"] for i in SPOKE_ORDER] for c in range(4)])

    R_dz = np.array([[dz.loc[i, f"Centroid_{cols['DZ'][c]}"] for i in SPOKE_ORDER] for c in range(4)])

    labels = [names[i] for i in SPOKE_ORDER]

    pd.DataFrame(np.vstack([R_mz, R_dz]).T, index=labels,

                 columns=[f"CAP{c + 1}_MZ" for c in range(4)] + [f"CAP{c + 1}_DZ" for c in range(4)]) \
        .to_csv(f"{OUT_STEM}_radar_values.csv", float_format="%.4f")



    # dynamic metrics from the master sheet

    master = pd.read_excel(MASTER, sheet_name="Master")

    master["Zygosity"] = master["Zygosity"].str.strip().str.upper()

    vstats = []

    for key in ("FO", "AR", "DT"):

        for c in range(4):

            a = master.loc[master.Zygosity == "MZ", f"{key}_State{c + 1}"].dropna()

            b = master.loc[master.Zygosity == "DZ", f"{key}_State{c + 1}"].dropna()

            vstats.append({"metric": key, "CAP": c + 1, "MZ n": len(a), "DZ n": len(b),

                           "MZ median": a.median(), "DZ median": b.median(),

                           "p (Welch t-test)": stats.ttest_ind(a, b, equal_var=False).pvalue})

    vstats = pd.DataFrame(vstats)

    vstats.to_csv(f"{OUT_STEM}_violin_stats.csv", index=False)



    cmap = signed_cmap(THRESH, VMAX)

    try:

        fs = datasets.fetch_surf_fsaverage("fsaverage")       # full-resolution mesh (downloaded once)

    except Exception as err:                                  # offline: nilearn's bundled fsaverage5 mesh

        print(f"WARNING: full fsaverage unavailable ({type(err).__name__}); using fsaverage5")

        fs = datasets.fetch_surf_fsaverage("fsaverage5")

    views = [("left", "lateral"), ("right", "lateral"), ("left", "medial"), ("right", "medial")]

    surf = {}

    for c, m in enumerate(maps2):

        for hemi in ("left", "right"):

            tex = np.nan_to_num(surface.vol_to_surf(m, fs[f"pial_{hemi}"]))

            for view in ("lateral", "medial"):

                surf[c, hemi, view] = render_surface(tex, hemi, view, fs, cmap)

    a_surf = np.median([im.shape[1] / im.shape[0] for im in surf.values()])



    # ---- layout (mm) -------------------------------------------------------------

    W = 227.0

    pipe_wb = np.asarray(Image.open(PIPE_WHOLEBRAIN).convert("RGB"))

    pipe_km = np.asarray(Image.open(PIPE_KMEANS).convert("RGB"))

    h_A = 31.0

    y_rows = h_A + 11.5

    row_h, row_gap = 50.0, 2.0

    lab_w = 33.0

    h_s = min((row_h - 1.0 - 2.0) / 2, 22.5)  # 2 x 2 surfaces

    w_s = h_s * a_surf

    x_surf = lab_w

    h_sl = 27.0

    w_sag = h_sl * np.ptp(LIMITS["x"][0]) / np.ptp(LIMITS["x"][1])

    w_cor = h_sl * np.ptp(LIMITS["y"][0]) / np.ptp(LIMITS["y"][1])

    x_sag = x_surf + 2 * w_s + 0.8 + 2.0

    x_cor = x_sag + w_sag + 0.8

    rad_d = 37.0

    x_radc = W - 36.0                         # full names need ~28 mm right, ~24 mm left

    x_rad0 = x_radc - 37.0

    y_leg = y_rows + 4 * row_h + 3 * row_gap + 2.5

    h_leg = 11.0                              # colour scale + MZ/DZ key only

    y_V = y_leg + h_leg + 1.5                 # top of the violin block

    v_left, v_gap, v_plot_h = 21.0, 3.0, 48.0  # y-label + tick-label margin, gap, plot height

    v_plot_w = (W - 3 * v_left - 2 * v_gap) / 3

    v_title = 8.0

    H = y_V + v_title + v_plot_h + 22.0          # room for the rotated state names

    mf = MMFigure(W, H)

    ov = mf.ax(0, 0, W, H)

    ov.set_xlim(0, W); ov.set_ylim(H, 0); ov.set_axis_off(); ov.patch.set_alpha(0)

    arrow = dict(arrowstyle="-|>", mutation_scale=11, lw=1.4, color=INK_SOFT)



    # ---- A: pipeline (pictures from the previous Figure 1) ------------------------

    letter(mf, 0, 0.3, "A")

    hb = 24.0; wb = hb * pipe_wb.shape[1] / pipe_wb.shape[0]

    xb, yb = 12.0, 1.5

    ax = mf.ax(xb, yb, wb, hb); ax.imshow(pipe_wb, interpolation="antialiased"); ax.set_axis_off()

    mf.text(xb + wb / 2, yb + hb + 1.2, "Whole-brain activation", ha="center", va="top",

            fontsize=FS_LABEL)

    hk = 15.0; wk = hk * pipe_km.shape[1] / pipe_km.shape[0]

    xk, yk = xb + wb + 24.0, yb + (hb - hk) / 2

    ov.add_patch(FancyArrowPatch((xb + wb + 2.0, yb + hb / 2), (xk - 2.0, yb + hb / 2), **arrow))

    ax = mf.ax(xk, yk, wk, hk); ax.imshow(pipe_km, interpolation="antialiased"); ax.set_axis_off()

    mf.text(xk + wk / 2, yb + hb + 1.2, "k-means clustering (k = 4)", ha="center", va="top",

            fontsize=FS_LABEL)

    xc = xk + wk + 14.0

    ov.add_patch(FancyArrowPatch((xk + wk + 2.0, yb + hb / 2), (xc - 2.0, yb + hb / 2), **arrow))

    mf.text(xc, yb + hb / 2 - 3.0, "CAP 1\u20134", ha="left", va="center", fontsize=FS_LABEL + 2)

    mf.text(xc, yb + hb / 2 + 3.4, "maps & network profiles (B\u2013E)\ndynamic metrics (F\u2013H)",

            ha="left", va="center", fontsize=FS_SMALL, color=INK_SOFT, linespacing=1.25)

    ov.plot([0, W], [h_A - 0.5, h_A - 0.5], color=INK, lw=0.5)



    # ---- column headers -----------------------------------------------------------

    for c, t in enumerate(["Left", "Right"]):

        mf.text(x_surf + c * (w_s + 0.8) + w_s / 2, y_rows - 5.4, t, ha="center", va="bottom",

                fontsize=FS_LABEL)

    mf.text(x_surf + w_s + 0.4, y_rows - 1.0, "top lateral \u00b7 bottom medial", ha="center",

            va="bottom", fontsize=FS_SMALL, color=INK_SOFT)

    for t, sub, x, w in (("Sagittal", f"x = {fmt(SAG_X)} mm", x_sag, w_sag),

                         ("Coronal", f"y = {fmt(COR_Y)} mm", x_cor, w_cor)):

        mf.text(x + w / 2, y_rows - 5.4, t, ha="center", va="bottom", fontsize=FS_LABEL)

        mf.text(x + w / 2, y_rows - 1.0, sub, ha="center", va="bottom", fontsize=FS_SMALL,

                color=INK_SOFT)

    mf.text(x_radc, y_rows - 5.4, "Yeo-17 network correlation", ha="center", va="bottom",

            fontsize=FS_LABEL)

    mf.text(x_radc, y_rows - 1.0, "MZ vs DZ", ha="center", va="bottom", fontsize=FS_SMALL,

            color=INK_SOFT)



    # ---- B-E: CAP rows -----------------------------------------------------------------

    stat_kw = dict(threshold=THRESH, cmap=cmap, vmin=-VMAX, vmax=VMAX, symmetric_cbar=True,

                   colorbar=False, annotate=False, draw_cross=False, black_bg=False, dim=-0.25)

    n = len(SPOKE_ORDER)

    theta = np.linspace(0, 2 * np.pi, n, endpoint=False)

    step = 2 * np.pi / n

    fam = [lab.split("-")[0] for lab in labels]

    ring_labels = [-0.4, 0.0, 0.4]

    for c in range(4):

        y0 = y_rows + c * (row_h + row_gap)

        letter(mf, 0, y0 + 0.2, "BCDE"[c])

        mf.text(0.3, y0 + row_h / 2 - 1.0, f"CAP {c + 1}", ha="left", va="center",

                fontsize=FS_LABEL + 1)

        mf.text(0.3, y0 + row_h / 2 + 4.6 + (2.2 if "–" in CAP_NAMES[c] else 0), CAP_NAMES[c].replace("–", "–\n"), ha="left", va="center",

                fontsize=FS_SMALL, color=INK_SOFT)

        if c < 3:

            sep = y0 + row_h + row_gap / 2

            ov.plot([0, W], [sep, sep], color=INK, lw=0.5)

        for i, (hemi, view) in enumerate(views):

            r, k = divmod(i, 2)

            ax = mf.ax(x_surf + k * (w_s + 0.8), y0 + 1.0 + r * (h_s + 1.0), w_s, h_s)

            ax.imshow(surf[c, hemi, view], interpolation="antialiased")

            ax.set_axis_off()

        for mode, cut, x, w in (("x", SAG_X, x_sag, w_sag), ("y", COR_Y, x_cor, w_cor)):

            ax = mf.ax(x, y0 + (row_h - h_sl) / 2, w, h_sl)

            kw = dict(stat_kw, bg_img=slab(tmpl, mode, cut, "continuous"))

            d = plotting.plot_stat_map(slab(maps2[c], mode, cut, "linear"), display_mode=mode,

                                       cut_coords=[cut], axes=ax, figure=mf.fig, **kw)

            cax = list(d.axes.values())[0].ax

            (h0, h1), (v0, v1) = LIMITS[mode]

            cax.set_xlim(h0, h1); cax.set_ylim(v0, v1)

            if c == 0 and mode == "y":

                for xx, t, ha in ((0.0, "L", "left"), (1.0, "R", "right")):

                    cax.text(xx, 0.0, t, transform=cax.transAxes, ha=ha, va="bottom",

                             fontsize=FS_SMALL, color=INK_SOFT)

        # radar

        rax = mf.ax(x_radc - rad_d / 2, y0 + (row_h - rad_d) / 2, rad_d, rad_d, projection="polar")

        rax.set_theta_offset(np.pi / 2); rax.set_theta_direction(-1)

        rax.set_ylim(R_MIN, R_MAX)

        shade = False

        for f in dict.fromkeys(fam):                               # family sectors

            idx = [i for i, ff in enumerate(fam) if ff == f]

            if shade:

                rax.bar(theta[idx].mean(), R_MAX - R_MIN, width=step * len(idx), bottom=R_MIN,

                        color="#ffffff", lw=0, zorder=0)          # no grey sector shading

            shade = not shade

        tt = np.linspace(0, 2 * np.pi, 361)

        for rr in R_RINGS + [R_MAX]:

            rax.plot(tt, np.full_like(tt, rr), color=INK, ls="-" if rr == 0 else (0, (1.5, 1.5)),

                     lw=1.3 if rr == 0 else 0.4, zorder=1)

        rax.set_yticks([]); rax.spines["polar"].set_visible(False)

        rax.set_xticks(theta); rax.set_xticklabels([])

        rax.xaxis.grid(True, color=INK, lw=0.25)

        for th, lab in zip(theta, labels):                       # full name on every spoke

            deg = np.degrees(th) % 360

            ha = "center" if (deg < 5 or deg > 355) else ("left" if deg < 180 else "right")

            va = "bottom" if (deg < 30 or deg > 330) else ("top" if 150 < deg < 210 else "center")

            rax.text(th, R_MAX + 0.05, lab, ha=ha, va=va, fontsize=FS_SPOKE, clip_on=False)

        for vals, col, ls in ((R_mz[c], MZ_COL, "-"), (R_dz[c], DZ_COL, (0, (2.4, 1.3)))):

            v = np.append(vals, vals[0]); t = np.append(theta, theta[0])

            rax.plot(t, v, color=col, lw=1.4, ls=ls, zorder=3)

            rax.plot(theta, vals, "o", ms=2.2, color=col, zorder=4)

        best, th_lab = np.inf, -step / 2                          # ring labels: emptiest gap

        for g in range(n):

            j = (g + 1) % n

            ang = (theta[g] + step / 2) % (2 * np.pi)

            deg = np.degrees(ang)

            if not (deg <= 55 or deg >= 305 or 125 <= deg <= 235):

                continue

            rg = [(R_mz[c][g] + R_mz[c][j]) / 2, (R_dz[c][g] + R_dz[c][j]) / 2]

            cost = sum(np.exp(-((r - ring) / 0.08) ** 2) for r in rg for ring in ring_labels)

            cost += 0.15 * min(ang, 2 * np.pi - ang) / np.pi

            if cost < best:

                best, th_lab = cost, theta[g] + step / 2

        for rr in ring_labels:

            rax.text(th_lab, rr, fmt(rr, 1) if rr else "0", fontsize=FS_RING, ha="center",

                     va="center", color=INK, zorder=5,

                     bbox=dict(boxstyle="square,pad=0.05", fc="white", ec="none", alpha=0.85))



    # ---- legend row: colour scale + MZ/DZ key (explanations go in the caption) ----

    cb_w, cb_h = 56.0, 2.8

    cb_x, cb_y = x_surf + 6.0, y_leg + 1.6

    cax = mf.ax(cb_x, cb_y, cb_w, cb_h)

    cb = mf.fig.colorbar(mpl.cm.ScalarMappable(norm=Normalize(-VMAX, VMAX), cmap=cmap), cax=cax,

                         orientation="horizontal", extend="both", extendfrac=0.04)

    ticks = [-VMAX, -THRESH, THRESH, VMAX]

    cb.set_ticks(ticks); cb.set_ticklabels([fmt(t / 1000, 1) for t in ticks])

    cb.outline.set_linewidth(0.6)

    cb.ax.tick_params(width=0.6, length=2.2, pad=1.2, labelsize=FS_SMALL)

    lx = x_rad0 - 8.0

    for i, (lab, col, ls) in enumerate((("MZ (monozygotic)", MZ_COL, "-"),

                                        ("DZ (dizygotic)", DZ_COL, (0, (2.4, 1.3))))):

        yy = y_leg + 2.2 + i * 5.0

        ov.plot([lx, lx + 7], [yy, yy], color=col, lw=1.5, ls=ls)

        ov.plot([lx + 3.5], [yy], "o", ms=2.6, color=col)

        ov.add_patch(mpl.patches.Rectangle((lx + 9.0, yy - 1.5), 4.0, 3.0, fc=VIOLIN_FILL[i],

                                           ec=INK, lw=0.4))

        mf.text(lx + 14.8, yy, lab, ha="left", va="center", fontsize=FS_SMALL)



    # ---- F-H: dynamic metrics, violins drawn from the master sheet ---------------------

    ov.plot([0, W], [y_V - 0.8, y_V - 0.8], color=INK, lw=0.5)

    titles = {"FO": "Fractional occupancy", "AR": "Appearance rate", "DT": "Mean dwell time"}

    for k, key in enumerate(("FO", "AR", "DT")):

        (y_lo, y_hi), yticks = VIOLIN_AXES[key]

        x0 = k * (v_left + v_plot_w + v_gap)

        letter(mf, x0, y_V + 0.3, "FGH"[k])

        mf.text(x0 + 7.0, y_V + 3.6, titles[key], ha="left", va="center", fontsize=FS_LABEL)

        ax = mf.ax(x0 + v_left, y_V + v_title, v_plot_w, v_plot_h)

        for c in range(4):

            for g_i, grp in enumerate(("MZ", "DZ")):

                vals = master.loc[master.Zygosity == grp, f"{key}_State{c + 1}"].dropna().to_numpy()

                pos = c + VIOLIN_POS[g_i]

                vp = ax.violinplot(vals, positions=[pos], widths=VIOLIN_WIDTH, showextrema=False,

                                   points=200)

                for body in vp["bodies"]:

                    body.set_facecolor(VIOLIN_FILL[g_i]); body.set_edgecolor(INK)

                    body.set_linewidth(0.5); body.set_alpha(1.0); body.set_zorder(1)

                ax.plot([pos - 0.07, pos + 0.07], [np.median(vals)] * 2, lw=1.0, solid_capstyle="butt",

                        color=INK, zorder=2)

        head = 0.13 * (y_hi - y_lo)              # headroom: one row of p-values

        ax.set_xlim(-0.5, 3.5); ax.set_ylim(y_lo, y_hi + head)

        ax.set_xticks(range(4))

        ax.set_xticklabels([f"CAP {c + 1}\n" + nm.replace("–", "–\n") for c, nm in enumerate(CAP_NAMES)], rotation=30, ha="right",

                           rotation_mode="anchor", linespacing=0.95, fontsize=10)

        ax.set_ylabel(YLABELS[key], fontsize=FS_SMALL, labelpad=2)

        ax.set_yticks(yticks)

        ax.set_yticklabels([f"{t:.2f}" if key == "FO" else f"{t:g}" for t in yticks])

        ax.tick_params(axis="y", labelsize=FS_SMALL, width=0.8, length=2.5, pad=1.5)

        ax.tick_params(axis="x", labelsize=10, width=0.8, length=2.5, pad=1.5)

        for side in ("top", "right"):

            ax.spines[side].set_visible(False)

        for side in ("left", "bottom"):

            ax.spines[side].set_linewidth(0.8)

        for c in range(4):                                 # MZ vs DZ p-value per CAP

            pv = vstats.loc[(vstats.metric == key) & (vstats.CAP == c + 1), "p (Welch t-test)"].iloc[0]

            ax.text(c, y_hi + 0.55 * head, f"{pv:.3f}", ha="center",

                    va="center", fontsize=FS_PVAL, zorder=3)

        ax.text(-0.62, y_hi + 0.55 * head, "p", ha="right", va="center", fontsize=FS_PVAL,

                fontstyle="italic", fontweight="bold", zorder=3)



    mf.fig.savefig(f"{OUT_STEM}.pdf", dpi=DPI)

    mf.fig.savefig(f"{OUT_STEM}.png", dpi=DPI)

    Image.MAX_IMAGE_PIXELS = None

    with Image.open(f"{OUT_STEM}.png") as im:

        im.convert("RGB").save(f"{OUT_STEM}.tif", dpi=(DPI, DPI), compression="tiff_lzw")

    print(f"Figure 1: {W:.0f} x {H:.1f} mm -> {OUT_STEM}.pdf/.tif/.png")





if __name__ == "__main__":

    main()
