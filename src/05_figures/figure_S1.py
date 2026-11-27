#!/usr/bin/env python3
"""
Supplementary Figure S1 - representative ICA components, brain-state maps and their
correspondence with the Yeo 17-network atlas.

  A  six representative group-ICA components (one per canonical network): sagittal,
     coronal and axial slices through the peak of the largest cluster + cortical surface
  B  the four brain-state centroid maps on the inflated surface and on sagittal, coronal
     and axial slices at the same MNI coordinates, one shared colour scale
  C  overlap of each state's strongest voxels with the Yeo-17 networks (% of voxels,
     positive and negative tails) and the threshold-free mean value per network

Inputs (all in inputs/):
  melodic_IC300_64_final.nii.gz            group-ICA spatial maps (64 components) - add this file
  centroids/state_1..4.nii.gz             brain-state centroid maps (State 1-4 = CAP 1-4)
  2026-11-13-SuppFigureS4-selected-ICs.xlsx  sheet "Selected ICs": the six ICs shown in panel A
                                           (sheet "Screening (all 64 ICs)" documents the choice)
  The Yeo-17 atlas (CBIG release) and the fsaverage surfaces are downloaded on the first run.
Steps:
  panel A  each IC thresholded at z > 10 for display (clusters < 100 mm3 hidden), shared colour
           scale z = 10-60, slices through the peak of the largest positive cluster
  panel B  one display cut-off for all states: |value| > 90th percentile of the pooled absolute
           values of the four maps (2334.09 a.u.); same cut-off defines the voxels in panel C
  panel C  % overlap = voxels(state tail AND network) / voxels(state tail in Yeo cortex) x 100
Output: outputs/SuppFigureS4.pdf (vector text), outputs/SuppFigureS4.tif (600 dpi), .png;
        panel files and the overlap tables (xlsx/csv) in outputs/
Run   : python 2026-11-13-SuppFigureS4.py
"""
import csv
import glob
import inspect
import json
import re
import sys
import urllib.request
import warnings
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import pandas as pd
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, Normalize
from matplotlib.figure import Figure
from nilearn import datasets, plotting, surface
from PIL import Image
from pypdf import PageObject, PdfReader, PdfWriter, Transformation
from scipy import ndimage

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

# =============================================================================
# 1. SETTINGS - input files and output folder
# =============================================================================
INPUTS = Path("inputs")
ICA_FILE = INPUTS / "melodic_IC300_64_final.nii.gz"
STATE_DIR = INPUTS / "centroids"                    # state_1..4.nii.gz
SELECTED_IC_SHEET = INPUTS / "2026-11-13-SuppFigureS4-selected-ICs.xlsx"
OUTPUT_DIR = Path("outputs")
PANEL_A_STEM = "SuppFigureS4_panelA"
PANELS_BC_STEM = "SuppFigureS4_panel"            # -> ..._panelB, ..._panelC
FINAL_STEM = "SuppFigureS4"

# ICs shown in panel A (1-based volume numbers in the 64-volume file, in display order)
SELECTED_ICS = (pd.read_excel(SELECTED_IC_SHEET, sheet_name="Selected ICs")["IC"]
                .astype(int).tolist())

# =============================================================================
# 2. PANEL A SETTINGS
# =============================================================================
# Display threshold on z (a visualisation choice, NOT a statistical test).
# State it in the legend, e.g. "maps thresholded at z > 10 for display".
THRESHOLD = 10.0

# Upper end of the shared colour scale. Values above it saturate (shown by the
# arrow end of the colour bar). None = automatic: the largest 99th percentile
# of the displayed values across the selected ICs, rounded up to a multiple of 5.
VMAX = 60.0

SHOW_NEGATIVE = False     # False = positive tail only (standard for RSN figures)

# Hide isolated clusters smaller than this (mm^3) - display only, state it in
# the legend ("clusters < 100 mm^3 not shown"). 0 = show everything.
MIN_CLUSTER_MM3 = 100

# Slice positions. By default every IC is cut through the peak voxel of its
# largest positive cluster. Override any IC with MNI mm coordinates (x, y, z):
CUT_COORDS_OVERRIDE = {
    # 13: (38, 10, 8),
}

# 4th column:
#   "surface" cortical surface (fsaverage, inflated) in the view that shows
#             most of the component: L/R x lateral/medial, labelled under it
#             <- recommended (clean, no outline mismatch)
#   "glass"   glass-brain projection (see GLASS_VIEW); maximum-intensity
#             projections pile up every deep cluster, so they look noisy
#   None      three slices only
FOURTH_COLUMN = "surface"
SURFACE_VIEW_OVERRIDE = {
    # 35: ("left", "medial"),
}
SURFACE_MESH = "fsaverage6"   # 41k vertices/hemisphere (downloaded once, ~5 MB);
                              # falls back to the bundled, coarser fsaverage5
SURFACE_BG_CONTRAST = 0.6     # 1 = nilearn default (dark sulci); lower = lighter
GLASS_VIEW = "x"   # "x" whole-brain sagittal; "l"/"r" ONE hemisphere only

# Layout (mm). 180 mm = full page width for most journals.
FIG_WIDTH_MM = 180
ICS_PER_ROW = 2           # 2 = compact (~80 mm tall); 1 = one IC per row
MAX_ROW_H_MM = 30         # cap on brain-image height (matters for ICS_PER_ROW = 1)
PANEL_LETTER = "A"
ROW_LABEL = "IC {ic}"     # e.g. "IC {ic}" -> "IC 2"

# Background anatomy: None = nilearn's MNI152 (2009c) template.
BG_IMG = None
BG_DIM = -0.25            # <0 lightens the anatomy so the colours stand out


SHOW_FIGURE = False

# =============================================================================
# 3. STYLE (shared by all panels)
# =============================================================================
FS_LETTER = 18  # panel letter          (journal version: 9)
FS_LABEL = 11   # titles, row labels    (journal version: 7)
FS_SMALL = 10   # coordinates, ticks    (journal version: 6)
FONT_WEIGHT = "bold"   # all text bold for viewing at a distance
INK = "#1a1a1a"
INK_SOFT = "#262626"   # near-black so secondary text reads from afar
SUBTHRESHOLD_GREY = "#d9d9d9"


# matplotlib style of each panel (the panels were designed with slightly different tick sizes)
RC_PANEL_A = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": FS_LABEL,
    "font.weight": FONT_WEIGHT,
    "axes.labelweight": FONT_WEIGHT,
    "axes.titleweight": FONT_WEIGHT,
    "text.color": INK,
    "axes.edgecolor": INK,
    "axes.linewidth": 0.5,
    "xtick.color": INK,
    "xtick.major.width": 0.5,
    "xtick.major.size": 2.0,
    "xtick.major.pad": 1.5,
    "xtick.labelsize": FS_SMALL,
    "axes.unicode_minus": True,
    "pdf.fonttype": 42,          # TrueType -> text stays editable
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
}
RC_PANELS_BC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": FS_LABEL, "font.weight": FONT_WEIGHT, "axes.labelweight": FONT_WEIGHT,
    "axes.titleweight": FONT_WEIGHT, "text.color": INK, "axes.edgecolor": INK,
    "axes.linewidth": 0.5, "xtick.color": INK, "ytick.color": INK,
    "xtick.major.width": 0.4, "ytick.major.width": 0.4,
    "xtick.major.size": 1.8, "ytick.major.size": 1.8,
    "xtick.labelsize": FS_SMALL, "ytick.labelsize": FS_SMALL,
    "axes.unicode_minus": True, "pdf.fonttype": 42, "ps.fonttype": 42,
    "svg.fonttype": "none", "figure.facecolor": "white", "savefig.facecolor": "white",
}

# Display bounds (MNI mm) - the MNI152 brain plus a 2 mm margin. Cropping to
# these removes the empty space nilearn leaves around every slice.
CROP = {"x": (-74, 74), "y": (-109, 75), "z": (-74, 84)}
SLICE_LIMITS = {   # (horizontal axis, vertical axis) of each slice view
    "x": (CROP["y"], CROP["z"]),   # sagittal: anterior to the right
    "y": (CROP["x"], CROP["z"]),   # coronal : neurological (L on left)
    "z": (CROP["x"], CROP["y"]),   # axial   : neurological, anterior up
}
MINUS = "\u2212"

MM = 1 / 25.4  # inches per mm



# =============================================================================
# 4. PANEL B-C SETTINGS
# =============================================================================
STATE_PATTERN = "*state_{k}.nii*"
N_STATES = 4
INCLUSION_PERCENTILE = 90        # shared |value| cut-off (pooled over states)
SENSITIVITY_PERCENTILES = (80, 95)
VMAX_PERCENTILE = 99             # colour-scale end (pooled), rounded up to 500 a.u.

SLICE_COORDS = (-2, -54, 28)     # x, y, z in MNI mm - identical for every state
                                 # (through the posterior cingulate / DMN hub)

# Yeo-17 (Schaefer/Kong network names, Yeo 2011 colours read from the CBIG LUT)
NETWORKS = ["VisCent", "VisPeri", "SomMotA", "SomMotB", "DorsAttnA", "DorsAttnB",
            "SalVentAttnA", "SalVentAttnB", "LimbicA", "LimbicB", "ContA", "ContB",
            "ContC", "DefaultA", "DefaultB", "DefaultC", "TempPar"]
FAMILIES = [("Visual", 2), ("SomMot", 2), ("DorsAttn", 2), ("SalVenAttn", 2),
            ("Limbic", 2), ("Control", 3), ("Default", 3), ("TP", 1)]
CBIG = ("https://raw.githubusercontent.com/ThomasYeoLab/CBIG/master/stable_projects/"
        "brain_parcellation/Yeo2011_fcMRI_clustering/1000subjects_reference/"
        "Yeo_JNeurophysiol11_SplitLabels/")
YEO_VOLUME = "Yeo2011_17Networks_N1000.split_components.FSL_MNI152_FreeSurferConformed_1mm.nii.gz"
YEO_LUT = "17Networks_ColorLUT_freeview.txt"

# =============================================================================
# 3. HELPERS
# =============================================================================

def fmt_num(v):
    """Integer-looking numbers without '.0', with a typographic minus."""
    s = f"{v:g}"
    return s.replace("-", MINUS)


def load_component(img, ic):
    """Read ONE volume lazily (the full 4D file never enters memory)."""
    data = np.array(img.dataobj[..., ic - 1], dtype=np.float32)  # writable copy
    data[~np.isfinite(data)] = 0
    return nib.Nifti1Image(data, img.affine)


def display_map(comp_img, threshold, show_negative, min_cluster_mm3):
    """Map exactly as displayed: optional positive-only, small clusters removed."""
    data = np.asarray(comp_img.dataobj).copy()
    if not show_negative:
        data[data < 0] = 0
    if min_cluster_mm3 > 0:
        vox_mm3 = abs(np.linalg.det(comp_img.affine[:3, :3]))
        min_vox = int(np.ceil(min_cluster_mm3 / vox_mm3))
        for supra in ([data > threshold, data < -threshold] if show_negative
                      else [data > threshold]):
            labels, _ = ndimage.label(supra, structure=np.ones((3, 3, 3)))
            sizes = np.bincount(labels.ravel())
            small = (sizes < min_vox)[labels] & supra
            data[small] = 0
    return nib.Nifti1Image(data, comp_img.affine)


def peak_coords(comp_img, threshold):
    """Peak voxel (MNI mm) of the largest positive cluster above threshold."""
    vol = np.asarray(comp_img.dataobj)
    supra = vol > threshold
    if not supra.any():                       # fall back to |z|
        vol = np.abs(vol)
        supra = vol > threshold
    if not supra.any():
        ijk = np.unravel_index(np.argmax(np.abs(vol)), vol.shape)
    else:
        labels, n = ndimage.label(supra, structure=np.ones((3, 3, 3)))
        sizes = np.bincount(labels.ravel())[1:]
        biggest = int(np.argmax(sizes)) + 1
        ijk = np.unravel_index(
            np.argmax(np.where(labels == biggest, vol, -np.inf)), vol.shape
        )
    xyz = nib.affines.apply_affine(comp_img.affine, ijk)
    return tuple(int(round(c)) for c in xyz)


def make_cmap(threshold, vmax, show_negative):
    """FSL-style map that is already bright AT the threshold.

    Positive: red -> orange -> yellow. Negative: blue -> sky -> pale cyan.
    Lightness rises with |z| on both sides, and the sub-threshold band
    (never drawn on the brain) is neutral grey on the colour bar.
    """
    pos = ["#c0161b", "#f0601e", "#ffe14d"]
    neg = ["#1c4fc9", "#3ea0f5", "#a8edff"]
    eps = 1e-4
    if not show_negative:
        f = threshold / vmax
        nodes = [(0, SUBTHRESHOLD_GREY), (f, SUBTHRESHOLD_GREY),
                 (f + eps, pos[0]), ((1 + f) / 2, pos[1]), (1, pos[2])]
    else:
        f = threshold / vmax
        lo, hi = 0.5 - f / 2, 0.5 + f / 2
        nodes = [(0, neg[2]), (lo / 2, neg[1]), (lo - eps, neg[0]),
                 (lo, SUBTHRESHOLD_GREY), (hi, SUBTHRESHOLD_GREY),
                 (hi + eps, pos[0]), ((1 + hi) / 2, pos[1]), (1, pos[2])]
    cmap = LinearSegmentedColormap.from_list("fsl_like", nodes, N=1024)
    cmap.set_over(pos[2])
    cmap.set_under(neg[2] if show_negative else SUBTHRESHOLD_GREY)
    return cmap


def first_cut_ax(display):
    return list(display.axes.values())[0].ax


def _mesh_coords(mesh):
    m = surface.load_surf_mesh(mesh)
    return np.asarray(m[0] if isinstance(m, (list, tuple)) else m.coordinates)


def surface_view(comp_img, fs, threshold, override=None):
    """Project to fsaverage5 and pick the hemisphere/view showing most signal."""
    tex = {h: np.nan_to_num(surface.vol_to_surf(comp_img, fs[f"pial_{h}"]))
           for h in ("left", "right")}
    if override:
        return tex, override
    best, best_score = ("left", "lateral"), -1.0
    for h in ("left", "right"):
        x = _mesh_coords(fs[f"infl_{h}"])[:, 0]
        lateral = x < np.median(x) if h == "left" else x > np.median(x)
        w = np.clip(np.abs(tex[h]) - threshold, 0, None)
        for view, sel in (("lateral", lateral), ("medial", ~lateral)):
            if w[sel].sum() > best_score:
                best, best_score = (h, view), w[sel].sum()
    return tex, best


def render_surface(tex, hemi, view, fs, threshold, vmin, vmax, cmap):
    """Render one inflated hemisphere off-screen; return a cropped RGBA image.

    Rasterising here keeps the PDF small (a 3D mesh as vector paths is
    ~20k polygons per brain and slow to open).
    """
    f = Figure(figsize=(4, 3), dpi=450)
    FigureCanvasAgg(f)
    ax = f.add_subplot(111, projection="3d")
    kw = dict(hemi=hemi, view=view, bg_map=fs[f"sulc_{hemi}"],
              threshold=threshold, vmax=vmax, cmap=cmap, colorbar=False,
              symmetric_cbar=SHOW_NEGATIVE, axes=ax, figure=f)
    if "vmin" in inspect.signature(plotting.plot_surf_stat_map).parameters:
        kw["vmin"] = vmin
    plotting.plot_surf_stat_map(fs[f"infl_{hemi}"], tex, **kw)
    f.canvas.draw()
    rgba = np.asarray(f.canvas.buffer_rgba()).copy()
    ink = (rgba[..., :3] < 245).any(-1) & (rgba[..., 3] > 0)
    rows, cols = np.where(ink.any(1))[0], np.where(ink.any(0))[0]
    pad = 6
    rgba = rgba[max(rows[0] - pad, 0):rows[-1] + pad, max(cols[0] - pad, 0):cols[-1] + pad]
    # Soften the sulcal shading (grey pixels only) so the overlay stands out;
    # white stays white. SURFACE_BG_CONTRAST = 1 keeps nilearn's default look.
    rgb = rgba[..., :3].astype(float)
    grey = (rgb.max(-1) - rgb.min(-1)) < 14
    rgb[grey] = 255 - SURFACE_BG_CONTRAST * (255 - rgb[grey])
    rgba[..., :3] = np.clip(rgb, 0, 255).astype(np.uint8)
    return rgba


def glass_aspect(sample_img, view):
    """Width/height of nilearn's glass-brain outline for this view."""
    f = plt.figure(figsize=(2, 2))
    d = plotting.plot_glass_brain(sample_img, display_mode=view, figure=f,
                                  colorbar=False, annotate=False,
                                  black_bg=False, threshold=1e9)
    ax = first_cut_ax(d)
    w, h = abs(np.diff(ax.get_xlim())[0]), abs(np.diff(ax.get_ylim())[0])
    plt.close(f)
    return w / h




# =============================================================================
# PANEL A
# =============================================================================
def make_panel_a(ica_file=ICA_FILE, output_dir=OUTPUT_DIR):
    """Panel A: representative ICs (slices through the peak + surface view)."""
    ica_file, output_dir = Path(ica_file), Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    img = nib.load(str(ica_file))
    print(f"ICA file : {ica_file.name}  shape {img.shape}")
    if len(img.shape) != 4:
        raise ValueError("Expected a 4D ICA NIfTI image.")
    if max(SELECTED_ICS) > img.shape[3] or min(SELECTED_ICS) < 1:
        raise ValueError(f"IC numbers must be 1..{img.shape[3]} (1-based).")

    # ---- load the selected components, find slice positions ---------------
    comps, coords, p99 = {}, {}, []
    for ic in SELECTED_ICS:
        raw = load_component(img, ic)
        n_brain = np.count_nonzero(np.asarray(raw.dataobj))
        if n_brain == 0:
            raise ValueError(f"IC {ic} is empty.")
        comp = display_map(raw, THRESHOLD, SHOW_NEGATIVE, MIN_CLUSTER_MM3)
        shown = np.asarray(comp.dataobj)
        shown = np.abs(shown[np.abs(shown) > THRESHOLD])
        p99.append(np.percentile(shown, 99) if shown.size else THRESHOLD)
        comps[ic] = comp
        coords[ic] = tuple(CUT_COORDS_OVERRIDE.get(ic, peak_coords(comp, THRESHOLD)))
        print(f"IC {ic:>2}: slices at x={coords[ic][0]:>4} y={coords[ic][1]:>4} "
              f"z={coords[ic][2]:>4} mm | {100 * shown.size / n_brain:4.1f}% of "
              f"brain voxels shown")
        del raw

    vmax = VMAX if VMAX is not None else 5 * np.ceil(max(p99) / 5)
    if vmax <= THRESHOLD:
        raise ValueError("VMAX must be larger than THRESHOLD.")
    cmap = make_cmap(THRESHOLD, vmax, SHOW_NEGATIVE)
    vmin = -vmax if SHOW_NEGATIVE else 0
    print(f"Colour scale: |z| {THRESHOLD:g} to {vmax:g} (shared by all ICs)")

    # ---- 4th column content ------------------------------------------------
    surf_imgs, surf_views = {}, {}
    if FOURTH_COLUMN == "surface":
        try:
            fs = datasets.fetch_surf_fsaverage(SURFACE_MESH)
        except Exception as err:                           # offline, etc.
            print(f"Could not get {SURFACE_MESH} ({err}); using fsaverage5")
            fs = datasets.fetch_surf_fsaverage("fsaverage5")  # ships with nilearn
        for ic in SELECTED_ICS:
            tex, (hemi, view) = surface_view(comps[ic], fs, THRESHOLD,
                                             SURFACE_VIEW_OVERRIDE.get(ic))
            surf_imgs[ic] = render_surface(tex[hemi], hemi, view, fs,
                                           THRESHOLD, vmin, vmax, cmap)
            surf_views[ic] = f"{hemi[0].upper()} {view}"
            print(f"IC {ic:>2}: surface view {surf_views[ic]}")

    # ---- layout in millimetres ---------------------------------------------
    aspects = [np.ptp(SLICE_LIMITS[v][0]) / np.ptp(SLICE_LIMITS[v][1])
               for v in ("x", "y", "z")]
    if FOURTH_COLUMN == "surface":
        aspects.append(float(np.median([im.shape[1] / im.shape[0]
                                        for im in surf_imgs.values()])))
    elif FOURTH_COLUMN == "glass":
        aspects.append(glass_aspect(comps[SELECTED_ICS[0]], GLASS_VIEW))

    label_w, view_gap, block_gap, row_gap = 11.5, 1.0, 5.0, 1.4
    top_h, header_h, coord_h, cbar_band = 1.0, 7.5, 4.4, 13.0
    n_rows = int(np.ceil(len(SELECTED_ICS) / ICS_PER_ROW))
    block_w = (FIG_WIDTH_MM - (ICS_PER_ROW - 1) * block_gap) / ICS_PER_ROW
    views_w = block_w - label_w - (len(aspects) - 1) * view_gap
    row_h = min(views_w / sum(aspects), MAX_ROW_H_MM)
    view_ws = [a * row_h for a in aspects]
    # centre each block's content if the height cap left spare width
    x_pad = (views_w - sum(view_ws)) / 2
    pitch = row_h + coord_h + row_gap
    grid_bottom = top_h + header_h + n_rows * pitch - row_gap
    fig_h = grid_bottom + cbar_band

    W, H = FIG_WIDTH_MM, fig_h
    fig = plt.figure(figsize=(W * MM, H * MM))

    def rect(x, y_top, w, h):     # mm from top-left -> figure fraction
        return [x / W, 1 - (y_top + h) / H, w / W, h / H]

    def text_mm(x, y_top, s, **kw):
        fig.text(x / W, 1 - y_top / H, s, **kw)

    stat_kw = dict(threshold=THRESHOLD, cmap=cmap, vmin=vmin, vmax=vmax,
                   symmetric_cbar=SHOW_NEGATIVE, colorbar=False, annotate=False,
                   draw_cross=False, black_bg=False, dim=BG_DIM)
    if BG_IMG is not None:
        stat_kw["bg_img"] = str(BG_IMG)

    # ---- panels ------------------------------------------------------------
    for k, ic in enumerate(SELECTED_ICS):
        r, b = divmod(k, ICS_PER_ROW)
        x0 = b * (block_w + block_gap) + x_pad
        y0 = top_h + header_h + r * pitch
        comp, (cx, cy, cz) = comps[ic], coords[ic]
        print(f"Plotting IC {ic} ...")

        # row label
        text_mm(x0 + label_w - 1.2, y0 + row_h / 2, ROW_LABEL.format(ic=ic),
                ha="right", va="center", fontsize=FS_LABEL, color=INK)

        x = x0 + label_w
        for v, (mode, cut) in enumerate(zip("xyz", (cx, cy, cz))):
            ax = fig.add_axes(rect(x, y0, view_ws[v], row_h))
            d = plotting.plot_stat_map(comp, display_mode=mode, cut_coords=[cut],
                                       axes=ax, figure=fig, **stat_kw)
            cax = first_cut_ax(d)
            (h0, h1), (v0, v1) = SLICE_LIMITS[mode]
            cax.set_xlim(h0, h1)
            cax.set_ylim(v0, v1)
            text_mm(x + view_ws[v] / 2, y0 + row_h + 0.35, f"{mode} = {fmt_num(cut)}",
                    ha="center", va="top", fontsize=FS_SMALL, color=INK_SOFT)
            if k == 0 and mode == "y":                  # orientation, once (bottom
                for xx, s, ha in ((0.0, "L", "left"), (1.0, "R", "right")):  # corners are empty)
                    cax.text(xx, 0.0, s, transform=cax.transAxes, ha=ha,
                             va="bottom", fontsize=FS_SMALL, color=INK_SOFT)
            x += view_ws[v] + view_gap

        if FOURTH_COLUMN == "surface":
            ax = fig.add_axes(rect(x, y0, view_ws[3], row_h))
            ax.imshow(surf_imgs[ic], interpolation="antialiased")
            ax.set_axis_off()
            text_mm(x + view_ws[3] / 2, y0 + row_h + 0.35, surf_views[ic],
                    ha="center", va="top", fontsize=FS_SMALL, color=INK_SOFT)
        elif FOURTH_COLUMN == "glass":
            ax = fig.add_axes(rect(x, y0, view_ws[3], row_h))
            plotting.plot_glass_brain(comp, display_mode=GLASS_VIEW, axes=ax,
                                      figure=fig, threshold=THRESHOLD, cmap=cmap,
                                      vmin=vmin, vmax=vmax,
                                      plot_abs=not SHOW_NEGATIVE,
                                      symmetric_cbar=SHOW_NEGATIVE, colorbar=False,
                                      annotate=False, black_bg=False, alpha=0.4)

        # column titles above the first row of each block
        if r == 0:
            xc = x0 + label_w
            titles = ["Sagittal", "Coronal", "Axial"] + {
                "surface": ["Surface"], "glass": ["Glass brain"]}.get(FOURTH_COLUMN, [])
            for v, title in enumerate(titles):
                text_mm(xc + view_ws[v] / 2, top_h + header_h - 1.0, title,
                        ha="center", va="bottom", fontsize=FS_LABEL, color=INK)
                xc += view_ws[v] + view_gap

    # ---- panel letter --------------------------------------------------------
    text_mm(0.0, 0.3, PANEL_LETTER, ha="left", va="top",
            fontsize=FS_LETTER, fontweight="bold", color="black")

    # ---- shared colour bar ---------------------------------------------------
    cb_w, cb_h = 50.0, 2.6
    cb_y = grid_bottom + 2.4
    cb_x = (W - cb_w) / 2
    cb_ax = fig.add_axes(rect(cb_x, cb_y, cb_w, cb_h))
    sm = mpl.cm.ScalarMappable(norm=Normalize(vmin, vmax), cmap=cmap)
    cb = fig.colorbar(sm, cax=cb_ax, orientation="horizontal",
                      extend="both" if SHOW_NEGATIVE else "max",
                      extendfrac=0.04)
    step = 10 if vmax <= 50 else 20
    mids = [t for t in np.arange(step, vmax, step) if t - THRESHOLD >= step / 2]
    ticks = ([-vmax, -THRESHOLD, THRESHOLD, vmax] if SHOW_NEGATIVE
             else [THRESHOLD, *mids, vmax])
    cb.set_ticks(ticks)
    cb.set_ticklabels([fmt_num(t) for t in ticks])
    cb.outline.set_linewidth(0.6)
    cb.ax.tick_params(width=0.6, length=2.2, pad=1.2, labelsize=FS_SMALL)
    text_mm(cb_x - 2.5, cb_y + cb_h / 2, "z-score", ha="right", va="center",
            fontsize=FS_LABEL, color=INK)

    # ---- export --------------------------------------------------------------
    stem = output_dir / PANEL_A_STEM
    fig.savefig(stem.with_suffix(".pdf"), dpi=600,
                metadata={"Title": f"{PANEL_A_STEM}", "Creator": "matplotlib/nilearn"})
    fig.savefig(stem.with_suffix(".png"), dpi=600)
    with Image.open(stem.with_suffix(".png")) as im:      # journal TIFF: RGB + LZW
        im.convert("RGB").save(stem.with_suffix(".tif"), dpi=(600, 600),
                               compression="tiff_lzw")

    with open(f"{stem}_slice_coordinates.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["IC", "x_mm", "y_mm", "z_mm", "source"])
        for ic in SELECTED_ICS:
            w.writerow([ic, *coords[ic],
                        "manual" if ic in CUT_COORDS_OVERRIDE else "cluster peak"])

    print(f"\nSaved ({W:.0f} x {H:.1f} mm):")
    for ext in (".pdf", ".tif", ".png"):
        print("  ", stem.with_suffix(ext).resolve())
    print("  ", Path(f"{stem}_slice_coordinates.csv").resolve())

    plt.close(fig)



# =============================================================================
# 3. DATA + ANALYSIS
# =============================================================================

def load_states(state_dir):
    files = []
    for k in range(1, N_STATES + 1):
        hits = sorted(glob.glob(str(Path(state_dir) / STATE_PATTERN.format(k=k))))
        if not hits:
            raise FileNotFoundError(f"No file for state {k} in {state_dir}")
        files.append(hits[0])
    imgs = [nib.load(f) for f in files]
    for f, im in zip(files, imgs):
        if im.shape != imgs[0].shape or not np.allclose(im.affine, imgs[0].affine):
            raise ValueError(f"{f} is not on the same grid as {files[0]}")
    data = np.stack([im.get_fdata(dtype=np.float32) for im in imgs], -1)
    mask = (data != 0).any(-1)
    print(f"States: {[Path(f).name for f in files]}\n  grid {imgs[0].shape}, "
          f"orientation {''.join(nib.aff2axcodes(imgs[0].affine))}, "
          f"voxel {imgs[0].header.get_zooms()[:3]} mm, {mask.sum()} in-mask voxels")
    return imgs[0], data, mask


def load_yeo17(atlas_dir, ref):
    atlas_dir = Path(atlas_dir)
    atlas_dir.mkdir(parents=True, exist_ok=True)
    for name in (YEO_VOLUME, YEO_LUT):
        if not (atlas_dir / name).exists():
            print(f"Downloading {name} (CBIG repository) ...")
            urllib.request.urlretrieve(CBIG + "MNI152/" + name, atlas_dir / name)
    img = nib.load(str(atlas_dir / YEO_VOLUME))
    lab = np.asarray(img.dataobj).squeeze().astype(np.int32)
    if lab.shape != ref.shape[:3] or not np.allclose(img.affine, ref.affine):
        raise ValueError("Atlas and state maps are on different grids: resample the "
                         "atlas with nearest-neighbour interpolation first.")
    merged, colours = np.zeros(lab.shape, np.int16), {}
    for line in open(atlas_dir / YEO_LUT):
        p = line.split()
        if len(p) < 5 or not p[0].isdigit() or int(p[0]) == 0:
            continue
        net = re.match(r"17Networks_[LR]H_([A-Za-z]+?)(?:_|$)", p[1]).group(1)
        merged[lab == int(p[0])] = NETWORKS.index(net) + 1
        colours.setdefault(net, tuple(int(c) / 255 for c in p[2:5]))
    return merged, [colours[n] for n in NETWORKS]


def overlap_tables(data, mask, yeo):
    pooled = np.abs(data[mask])
    cortex = (yeo > 0) & mask
    net_vox = np.array([(yeo == j + 1).sum() for j in range(17)])
    expected = 100 * net_vox / net_vox.sum()
    rows, thresholds = [], {}
    for pct in sorted({INCLUSION_PERCENTILE, *SENSITIVITY_PERCENTILES}):
        T = float(np.percentile(pooled, pct))
        thresholds[pct] = T
        for s in range(data.shape[-1]):
            v = data[..., s]
            for tail, sel in (("positive", v > T), ("negative", v < -T)):
                sel &= mask
                n_all, n_ctx = int(sel.sum()), int((sel & cortex).sum())
                counts = np.bincount(yeo[sel & cortex], minlength=18)[1:18]
                for j, net in enumerate(NETWORKS):
                    pc = 100 * counts[j] / n_ctx if n_ctx else np.nan
                    rows.append(dict(
                        inclusion_percentile=pct, threshold_abs_value=round(T, 2),
                        state=s + 1, tail=tail, network=net,
                        voxels_state_tail=n_all, voxels_state_tail_in_cortex=n_ctx,
                        voxels_overlap=int(counts[j]), pct_of_cortical_state_voxels=pc,
                        pct_of_all_state_voxels=100 * counts[j] / n_all if n_all else np.nan,
                        pct_state_voxels_outside_yeo_cortex=(100 * (n_all - n_ctx) / n_all
                                                             if n_all else np.nan),
                        network_share_of_cortex_pct=expected[j],
                        enrichment_obs_over_expected=pc / expected[j] if n_ctx else np.nan))
    long = pd.DataFrame(rows)
    signed = pd.DataFrame([[float(data[..., s][yeo == j + 1].mean()) for j in range(17)]
                           for s in range(data.shape[-1])],
                          index=[f"State {s + 1}" for s in range(data.shape[-1])],
                          columns=NETWORKS)
    return long, signed, thresholds


def write_tables(long, signed, params, stem):
    main = long[long.inclusion_percentile == INCLUSION_PERCENTILE]
    piv = lambda col: main.pivot_table(index=["state", "tail"], columns="network",
                                       values=col)[NETWORKS]
    piv("pct_of_cortical_state_voxels").to_csv(f"{stem}C_overlap_4x17_pct.csv",
                                               float_format="%.2f")
    signed.to_csv(f"{stem}C_signed_mean_4x17.csv", float_format="%.3f")
    long.to_csv(f"{stem}C_overlap_all_thresholds_long.csv", index=False,
                float_format="%.4f")
    with pd.ExcelWriter(f"{stem}C_overlap_results.xlsx") as xw:
        piv("pct_of_cortical_state_voxels").round(2).to_excel(xw, sheet_name="pct_overlap")
        main.groupby(["state", "tail"]).pct_state_voxels_outside_yeo_cortex.first() \
            .round(2).to_frame().to_excel(xw, sheet_name="pct_outside_cortex")
        piv("enrichment_obs_over_expected").round(3).to_excel(xw, sheet_name="enrichment")
        piv("pct_of_all_state_voxels").round(3).to_excel(xw, sheet_name="pct_of_whole_mask")
        signed.round(2).to_excel(xw, sheet_name="signed_mean")
        long.to_excel(xw, sheet_name="all_thresholds_long", index=False)
        pd.DataFrame([(k, str(v)) for k, v in params.items()],
                     columns=["parameter", "value"]).to_excel(xw, sheet_name="parameters",
                                                              index=False)
    with open(f"{stem}C_parameters.json", "w") as fh:
        json.dump(params, fh, indent=2, default=str)


# =============================================================================
# 4. DRAWING HELPERS
# =============================================================================

def fmt(v, nd=0):
    if nd and round(v, nd) == 0:
        v = 0.0
    s = f"{v:.{nd}f}" if nd else f"{v:g}"
    return s.replace("-", MINUS)


def signed_cmap(threshold, vmax):
    """FSL-style, bright at threshold; grey band = not shown (as panel A)."""
    pos = ["#c0161b", "#f0601e", "#ffe14d"]
    neg = ["#1c4fc9", "#3ea0f5", "#a8edff"]
    f, eps = threshold / vmax, 1e-4
    lo, hi = 0.5 - f / 2, 0.5 + f / 2
    nodes = [(0, neg[2]), (lo / 2, neg[1]), (lo - eps, neg[0]), (lo, SUBTHRESHOLD_GREY),
             (hi, SUBTHRESHOLD_GREY), (hi + eps, pos[0]), ((1 + hi) / 2, pos[1]), (1, pos[2])]
    cm = LinearSegmentedColormap.from_list("signed", nodes, N=1024)
    cm.set_over(pos[2]); cm.set_under(neg[2])
    return cm


def diverging_cmap():
    """Continuous blue - neutral - red for the threshold-free means."""
    return LinearSegmentedColormap.from_list(
        "div", ["#1c3f95", "#4f86c6", "#b9d3ea", "#f4f4f4", "#f3c1a7", "#d6604d", "#8e0c1e"],
        N=256)


def sequential_cmap():
    return LinearSegmentedColormap.from_list("seq", ["#f7f7f7", "#bdbdbd", "#636363",
                                                     "#1f1f1f"], N=256)


def crop_and_soften(rgba, soften=True):
    ink = (rgba[..., :3] < 245).any(-1) & (rgba[..., 3] > 0)
    rows, cols = np.where(ink.any(1))[0], np.where(ink.any(0))[0]
    p = 6
    rgba = rgba[max(rows[0] - p, 0):rows[-1] + p, max(cols[0] - p, 0):cols[-1] + p].copy()
    if soften:
        rgb = rgba[..., :3].astype(float)
        grey = (rgb.max(-1) - rgb.min(-1)) < 14
        rgb[grey] = 255 - SURFACE_BG_CONTRAST * (255 - rgb[grey])
        rgba[..., :3] = np.clip(rgb, 0, 255).astype(np.uint8)
    return rgba


def _offscreen():
    f = Figure(figsize=(4, 3), dpi=450)
    FigureCanvasAgg(f)
    return f, f.add_subplot(111, projection="3d")


def render_stat_surface(tex, hemi, view, fs, threshold, vmax, cmap):
    f, ax = _offscreen()
    kw = dict(hemi=hemi, view=view, bg_map=fs[f"sulc_{hemi}"], threshold=threshold,
              vmax=vmax, cmap=cmap, colorbar=False, symmetric_cbar=True, axes=ax, figure=f)
    if "vmin" in inspect.signature(plotting.plot_surf_stat_map).parameters:
        kw["vmin"] = -vmax
    plotting.plot_surf_stat_map(fs[f"infl_{hemi}"], tex, **kw)
    f.canvas.draw()
    return crop_and_soften(np.asarray(f.canvas.buffer_rgba()))


def render_roi_surface(labels, hemi, view, fs, colours):
    f, ax = _offscreen()
    kw = dict(hemi=hemi, view=view, bg_map=fs[f"sulc_{hemi}"],
              cmap=ListedColormap(colours), colorbar=False, axes=ax, figure=f)
    params = inspect.signature(plotting.plot_surf_roi).parameters
    if "vmin" in params:
        kw.update(vmin=0.5, vmax=len(colours) + 0.5)
    plotting.plot_surf_roi(fs[f"infl_{hemi}"], labels, **kw)
    f.canvas.draw()
    return crop_and_soften(np.asarray(f.canvas.buffer_rgba()))


def get_fsaverage():
    try:
        return datasets.fetch_surf_fsaverage(SURFACE_MESH)
    except Exception as err:
        print(f"Could not get {SURFACE_MESH} ({err}); using bundled fsaverage5")
        return datasets.fetch_surf_fsaverage("fsaverage5")


class MMFigure:
    """Figure laid out in millimetres from the top-left corner."""

    def __init__(self, w, h):
        self.W, self.H = w, h
        self.fig = plt.figure(figsize=(w * MM, h * MM))

    def ax(self, x, y, w, h, **kw):
        return self.fig.add_axes([x / self.W, 1 - (y + h) / self.H, w / self.W, h / self.H], **kw)

    def text(self, x, y, s, **kw):
        kw.setdefault("color", INK)
        return self.fig.text(x / self.W, 1 - y / self.H, s, **kw)

    def save(self, stem):
        self.fig.savefig(f"{stem}.pdf", dpi=600)
        self.fig.savefig(f"{stem}.png", dpi=600)
        with Image.open(f"{stem}.png") as im:
            im.convert("RGB").save(f"{stem}.tif", dpi=(600, 600), compression="tiff_lzw")
        plt.close(self.fig)


def colourbar(mf, x, y, w, h, cmap, vmin, vmax, ticks, labels, label,
              orientation="horizontal", extend="neither"):
    cax = mf.ax(x, y, w, h)
    cb = mf.fig.colorbar(mpl.cm.ScalarMappable(norm=Normalize(vmin, vmax), cmap=cmap),
                         cax=cax, orientation=orientation, extend=extend, extendfrac=0.04)
    cb.set_ticks(ticks)
    cb.set_ticklabels(labels)
    cb.outline.set_linewidth(0.6)
    cb.ax.tick_params(width=0.6, length=2.2, pad=1.2, labelsize=FS_SMALL)
    return cb


# =============================================================================
# 5. PANEL B
# =============================================================================

SURF_VIEWS = [("left", "lateral"), ("left", "medial"), ("right", "medial"), ("right", "lateral")]
SURF_TITLES = ["Left lateral", "Left medial", "Right medial", "Right lateral"]


def build_panel_b(ref, data, T, vmax, fs, stem, states=None, letter="B"):
    states = list(range(data.shape[-1])) if states is None else states
    cmap = signed_cmap(T, vmax)
    imgs = {s: nib.Nifti1Image(data[..., s], ref.affine) for s in states}

    surf = {}
    for s in states:
        for hemi in ("left", "right"):
            tex = np.nan_to_num(surface.vol_to_surf(imgs[s], fs[f"pial_{hemi}"]))
            for view in ("lateral", "medial"):
                surf[s, hemi, view] = render_stat_surface(tex, hemi, view, fs, T, vmax, cmap)
    a_lat = np.median([surf[k].shape[1] / surf[k].shape[0] for k in surf if k[2] == "lateral"])
    a_med = np.median([surf[k].shape[1] / surf[k].shape[0] for k in surf if k[2] == "medial"])
    aspects = [a_lat, a_med, a_med, a_lat] + [np.ptp(SLICE_LIMITS[v][0]) /
                                              np.ptp(SLICE_LIMITS[v][1]) for v in "xyz"]
    gaps = [0.6, 0.6, 0.6, 3.0, 0.8, 0.8]

    W, label_w, top, header_h, row_gap, cbar_band = 180.0, 14.5, 1.0, 10.0, 1.4, 13.0
    row_h = (W - label_w - sum(gaps)) / sum(aspects)
    ws = [a * row_h for a in aspects]
    n = len(states)
    grid_bottom = top + header_h + n * row_h + (n - 1) * row_gap
    mf = MMFigure(W, grid_bottom + cbar_band)

    xs, x = [], label_w
    for i, w in enumerate(ws):
        xs.append(x)
        x += w + (gaps[i] if i < len(gaps) else 0)

    stat_kw = dict(threshold=T, cmap=cmap, vmin=-vmax, vmax=vmax, symmetric_cbar=True,
                   colorbar=False, annotate=False, draw_cross=False, black_bg=False, dim=BG_DIM)
    for r, s in enumerate(states):
        y0 = top + header_h + r * (row_h + row_gap)
        mf.text(label_w - 1.5, y0 + row_h / 2, f"State {s + 1}", ha="right", va="center",
                fontsize=FS_LABEL)
        for c, (hemi, view) in enumerate(SURF_VIEWS):
            ax = mf.ax(xs[c], y0, ws[c], row_h)
            ax.imshow(surf[s, hemi, view], interpolation="antialiased")
            ax.set_axis_off()
        for c, (mode, cut) in enumerate(zip("xyz", SLICE_COORDS)):
            ax = mf.ax(xs[4 + c], y0, ws[4 + c], row_h)
            d = plotting.plot_stat_map(imgs[s], display_mode=mode, cut_coords=[cut],
                                       axes=ax, figure=mf.fig, **stat_kw)
            cax = list(d.axes.values())[0].ax
            (h0, h1), (v0, v1) = SLICE_LIMITS[mode]
            cax.set_xlim(h0, h1); cax.set_ylim(v0, v1)
            if r == 0 and mode == "y":                   # orientation once, in the
                for xx, t, ha in ((0.0, "L", "left"), (1.0, "R", "right")):  # empty corners
                    cax.text(xx, 0.0, t, transform=cax.transAxes, ha=ha, va="bottom",
                             fontsize=FS_SMALL, color=INK_SOFT)

    # column headers (slice coordinates once - they are identical for all states)
    for c, t in enumerate(SURF_TITLES):
        mf.text(xs[c] + ws[c] / 2, top + header_h - 0.8, t, ha="center", va="bottom",
                fontsize=FS_LABEL)
    for c, (t, mode, cut) in enumerate(zip(("Sagittal", "Coronal", "Axial"), "xyz",
                                           SLICE_COORDS)):
        xc = xs[4 + c] + ws[4 + c] / 2
        mf.text(xc, top + header_h - 4.9, t, ha="center", va="bottom", fontsize=FS_LABEL)
        mf.text(xc, top + header_h - 0.8, f"{mode} = {fmt(cut)}", ha="center", va="bottom",
                fontsize=FS_SMALL, color=INK_SOFT)
    mf.text(0.0, 0.3, letter, ha="left", va="top", fontsize=FS_LETTER, fontweight="bold",
            color="black")

    # shared colour bar, in thousands of a.u.
    k = 1000.0
    cbw, cbh = 54.0, 2.6
    cbx, cby = (W - cbw) / 2 + 10, grid_bottom + 2.6
    colourbar(mf, cbx, cby, cbw, cbh, cmap, -vmax, vmax,
              [-vmax, -T, T, vmax],
              [fmt(-vmax / k, 1), fmt(-T / k, 1), fmt(T / k, 1), fmt(vmax / k, 1)],
              None, extend="both")
    mf.text(cbx - 2.5, cby + cbh / 2, "Centroid value (a.u. \u00d7 10\u00b3)", ha="right",
            va="center", fontsize=FS_LABEL)
    mf.save(stem)
    return mf.H


def yeo_surface_labels(atlas_dir, fs):
    """Official Yeo-17 split-component annotation (fsaverage5, same CBIG release),
    merged to 17 networks and carried to the display mesh by nearest neighbour on
    the sphere. Medial wall and boundary vertices removed by the release = NaN."""
    from nibabel.freesurfer import read_annot
    from scipy.spatial import cKDTree
    atlas_dir = Path(atlas_dir)
    fs5 = datasets.fetch_surf_fsaverage("fsaverage5")

    def xyz(mesh):
        m = surface.load_surf_mesh(mesh)
        return np.asarray(m[0] if isinstance(m, (list, tuple)) else m.coordinates)

    out = {}
    for hemi, h in (("left", "lh"), ("right", "rh")):
        name = f"{h}.Yeo2011_17Networks_N1000.split_components.annot"
        if not (atlas_dir / name).exists():
            urllib.request.urlretrieve(CBIG + "fsaverage5/label/" + name, atlas_dir / name)
        idx, _, names = read_annot(str(atlas_dir / name))
        net = np.zeros(len(names))
        for i, n in enumerate(names):
            m = re.match(r"17Networks_[LR]H_([A-Za-z]+?)(?:_|$)", n.decode())
            net[i] = NETWORKS.index(m.group(1)) + 1 if m else np.nan
        lab5 = net[idx]
        lab5[lab5 == 0] = np.nan
        _, nn = cKDTree(xyz(fs5[f"sphere_{hemi}"])).query(xyz(fs[f"sphere_{hemi}"]))
        out[hemi] = lab5[nn]
    return out


# =============================================================================
# 6. PANEL C
# =============================================================================

def build_panel_c(ref, yeo, colours, long, signed, T, fs, stem, letter="C"):
    main = long[long.inclusion_percentile == INCLUSION_PERCENTILE]
    rows = [(s, t) for s in range(1, N_STATES + 1) for t in ("positive", "negative")]
    pct = np.array([[main[(main.state == s) & (main["tail"] == t) & (main.network == n)]
                     .pct_of_cortical_state_voxels.iloc[0] for n in NETWORKS] for s, t in rows])
    outside = np.array([main[(main.state == s) & (main["tail"] == t)]
                        .pct_state_voxels_outside_yeo_cortex.iloc[0] for s, t in rows])
    sm = signed.values / 1000.0

    # Yeo-17 brain map: the official surface annotation of the same release
    labels = yeo_surface_labels(Path(stem).parent / "yeo17_atlas", fs)
    roi = {(h, v): render_roi_surface(labels[h], h, v, fs, colours)
           for h in ("left", "right") for v in ("lateral", "medial")}

    W = 180.0
    brain_w, brain_gap = 20.0, 1.5
    a = np.median([im.shape[1] / im.shape[0] for im in roi.values()])
    brain_h = brain_w / a
    top, chip_h = 1.0, 2.6
    hm_x, label_col = 45.0, 20.5           # heatmap block starts at hm_x
    grid_x = hm_x + label_col
    out_gap, cb_gap, cb_w, cb_lab = 1.6, 2.4, 2.2, 10.6
    cell_w = (W - grid_x - out_gap - cb_gap - cb_w - cb_lab) / 18   # 17 networks + outside
    cell_h = 4.6
    grid_w = 17 * cell_w
    y_chip = top + 10.0
    y_h1 = y_chip + chip_h + 0.8
    h1 = len(rows) * cell_h
    y_h2 = y_h1 + h1 + 7.0
    h2 = N_STATES * cell_h
    lab_h = 21.0
    H = y_h2 + h2 + lab_h
    mf = MMFigure(W, H)
    fam_edges = np.cumsum([n for _, n in FAMILIES])[:-1]   # family boundaries (grouping)

    # -- Yeo-17 brain map (2 x 2) -----------------------------------------------
    mf.text(0.0, 0.3, letter, ha="left", va="top", fontsize=FS_LETTER, fontweight="bold",
            color="black")
    bx0 = 2.0
    mf.text(bx0 + brain_w + brain_gap / 2, y_chip - 1.0, "Yeo-17 networks",
            ha="center", va="bottom", fontsize=FS_LABEL)
    by0 = y_chip + 0.6
    row_step = brain_h + 5.0
    for i, (h, v) in enumerate([("left", "lateral"), ("right", "lateral"),
                                ("left", "medial"), ("right", "medial")]):
        r, c = divmod(i, 2)
        ax = mf.ax(bx0 + c * (brain_w + brain_gap), by0 + r * row_step, brain_w, brain_h)
        ax.imshow(roi[h, v], interpolation="antialiased")
        ax.set_axis_off()
        mf.text(bx0 + c * (brain_w + brain_gap) + brain_w / 2, by0 + r * row_step + brain_h + 0.4,
                f"{h[0].upper()} {v}", ha="center", va="top", fontsize=FS_SMALL,
                color=INK_SOFT)

    # -- colour chips (network key; families separated by white gaps below) ------
    for j in range(17):
        ax = mf.ax(grid_x + j * cell_w + 0.3, y_chip, cell_w - 0.6, chip_h)
        ax.set_facecolor(colours[j])
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_linewidth(0.4); sp.set_color("#7f7f7f")

    # -- network-family names over brackets; at this text size they alternate
    #    between two rows so wide names (e.g. SalVenAttn) never collide -----------
    ov = mf.ax(0, 0, W, H)                    # transparent overlay in mm
    ov.set_xlim(0, W); ov.set_ylim(H, 0); ov.set_axis_off(); ov.patch.set_alpha(0)
    yb = y_chip - 0.9                          # bracket line
    j0 = 0
    for k, (name, nfam) in enumerate(FAMILIES):
        xa, xb = grid_x + j0 * cell_w + 0.5, grid_x + (j0 + nfam) * cell_w - 0.5
        xm = (xa + xb) / 2
        ov.plot([xa, xa, xb, xb], [yb + 0.6, yb, yb, yb + 0.6], color=INK_SOFT, lw=0.8,
                solid_joinstyle="miter", clip_on=False)
        if k % 2 == 0:                         # upper row + leader down to the bracket
            y_lab = yb - 4.9
            ov.plot([xm, xm], [y_lab + 0.3, yb], color=INK_SOFT, lw=0.6, clip_on=False)
        else:                                  # lower row, right on the bracket
            y_lab = yb - 0.5
        mf.text(xm, y_lab, name, ha="center", va="bottom", fontsize=FS_SMALL, color=INK_SOFT)
        j0 += nfam

    # -- heatmap 1: % overlap ----------------------------------------------------
    seq = sequential_cmap()
    pmax = 10 * np.ceil(np.nanmax(pct) / 10)
    ax1 = mf.ax(grid_x, y_h1, grid_w, h1)
    ax1.imshow(pct, cmap=seq, vmin=0, vmax=pmax, aspect="auto", interpolation="nearest")
    axo = mf.ax(grid_x + grid_w + out_gap, y_h1, cell_w, h1)
    axo.imshow(outside[:, None], cmap=seq, vmin=0, vmax=pmax, aspect="auto",
               interpolation="nearest", alpha=0.55)
    for axx, M in ((ax1, pct), (axo, outside[:, None])):
        for (i, j), val in np.ndenumerate(M):
            dark = val / pmax > 0.55
            axx.text(j, i, f"{val:.0f}", ha="center", va="center", fontsize=FS_SMALL,
                     color="white" if (dark and axx is ax1) else INK)
        axx.set_xticks(np.arange(M.shape[1] + 1) - 0.5, minor=True)
        axx.set_yticks(np.arange(M.shape[0] + 1) - 0.5, minor=True)
        axx.grid(which="minor", color="white", lw=0.8)
        axx.tick_params(which="both", length=0)
        axx.set_xticks([]); axx.set_yticks([])
        for sp in axx.spines.values():
            sp.set_visible(False)
        for i in range(2, M.shape[0], 2):                  # separate states
            axx.axhline(i - 0.5, color="white", lw=2.6)
        if axx is ax1:
            for e in fam_edges:                            # separate network families
                axx.axvline(e - 0.5, color="white", lw=2.6)
    mf.text(grid_x + grid_w + out_gap, y_chip - 0.4, "Outside\ncortex",   # left-aligned:
            ha="left", va="bottom", fontsize=FS_SMALL, color=INK_SOFT,      # clears "TP"
            linespacing=0.95)
    for i, (s, t) in enumerate(rows):
        yc = y_h1 + (i + 0.5) * cell_h
        if t == "positive":
            mf.text(hm_x + 0.3, yc + cell_h / 2, f"State {s}", ha="left", va="center",
                    fontsize=FS_LABEL)
        sq = mf.ax(grid_x - 5.4, yc - 1.0, 2.0, 2.0)
        sq.set_facecolor("#c0161b" if t == "positive" else "#1c4fc9")
        sq.set_xticks([]); sq.set_yticks([])
        for sp in sq.spines.values():
            sp.set_visible(False)
        mf.text(grid_x - 0.7, yc, "+" if t == "positive" else MINUS, ha="right",
                va="center", fontsize=FS_SMALL)
    mf.text(hm_x + 0.3, y_chip - 0.6, "% overlap", ha="left", va="bottom",
            fontsize=FS_LABEL, fontweight="bold")
    colourbar(mf, grid_x + grid_w + out_gap + cell_w + cb_gap, y_h1, cb_w, h1, seq, 0, pmax,
              [0, pmax / 2, pmax], [fmt(0), fmt(pmax / 2), fmt(pmax)], None,
              orientation="vertical")
    mf.text(W - 1.9, y_h1 + h1 / 2,
            "% of tail voxels", rotation=90, ha="center", va="center", fontsize=FS_SMALL)

    # -- heatmap 2: threshold-free signed mean --------------------------------------
    div = diverging_cmap()
    smax = np.ceil(np.nanmax(np.abs(sm)))
    ax2 = mf.ax(grid_x, y_h2, grid_w, h2)
    ax2.imshow(sm, cmap=div, vmin=-smax, vmax=smax, aspect="auto", interpolation="nearest")
    # colour only at this text size; exact values are in the CSV / Excel output
    ax2.set_xticks(np.arange(18) - 0.5, minor=True)
    ax2.set_yticks(np.arange(N_STATES + 1) - 0.5, minor=True)
    ax2.grid(which="minor", color="white", lw=0.8)
    for e in fam_edges:
        ax2.axvline(e - 0.5, color="white", lw=2.6)
    ax2.tick_params(which="both", length=0)
    ax2.set_yticks([])
    for sp in ax2.spines.values():
        sp.set_visible(False)
    ax2.set_xticks(range(17))
    ax2.set_xticklabels(NETWORKS, rotation=55, ha="right", rotation_mode="anchor",
                        fontsize=FS_SMALL)
    ax2.tick_params(axis="x", pad=1.0)
    for i in range(N_STATES):
        mf.text(hm_x + 0.3, y_h2 + (i + 0.5) * cell_h, f"State {i + 1}", ha="left",
                va="center", fontsize=FS_LABEL)
    mf.text(hm_x + 0.3, y_h2 - 3.2, "Mean centroid value (threshold-free)", ha="left",
            va="center", fontsize=FS_LABEL, fontweight="bold")
    colourbar(mf, grid_x + grid_w + out_gap + cell_w + cb_gap, y_h2, cb_w, h2, div, -smax, smax,
              [-smax, 0, smax], [fmt(-smax), "0", fmt(smax)], None, orientation="vertical")
    mf.text(W - 1.9, y_h2 + h2 / 2, "a.u. \u00d7 10\u00b3",
            rotation=90, ha="center", va="center", fontsize=FS_SMALL)
    mf.save(stem)
    return mf.H




# =============================================================================
# PANELS B-C
# =============================================================================
def make_panels_bc(state_dir=STATE_DIR, output_dir=OUTPUT_DIR):
    """Panels B and C: state maps and their overlap with the Yeo-17 networks."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    stem = str(out / PANELS_BC_STEM)

    ref, data, mask = load_states(state_dir)
    yeo, colours = load_yeo17(out / "yeo17_atlas", ref)
    long, signed, thresholds = overlap_tables(data, mask, yeo)
    T = thresholds[INCLUSION_PERCENTILE]
    vmax = 500 * np.ceil(np.percentile(np.abs(data[mask]), VMAX_PERCENTILE) / 500)
    print(f"Shared cut-off |value| > {T:.1f} a.u. (pooled {INCLUSION_PERCENTILE}th pct); "
          f"colour scale \u00b1{vmax:.0f}")
    params = dict(
        state_files_grid=str(ref.shape), affine=ref.affine.round(3).tolist(),
        orientation="".join(nib.aff2axcodes(ref.affine)),
        reference_space="FSL MNI152 1 mm, FreeSurfer-conformed 256^3 (sform_code 2)",
        in_mask_voxels=int(mask.sum()),
        yeo_atlas="CBIG Yeo2011 17Networks N1000 split_components "
                  "FSL_MNI152_FreeSurferConformed_1mm, merged to 17 networks (no resampling)",
        inclusion_percentile=INCLUSION_PERCENTILE,
        inclusion_threshold_abs=round(T, 2),
        sensitivity_thresholds={str(k): round(v, 2) for k, v in thresholds.items()},
        display_vmax=float(vmax), slice_coords_mm=list(SLICE_COORDS),
        surface=f"{SURFACE_MESH} inflated; vol_to_surf on pial (nilearn defaults)",
        overlap_definition="voxels(state tail AND network) / voxels(state tail in Yeo cortex) x 100",
    )
    write_tables(long, signed, params, stem)

    fs = get_fsaverage()
    print("Building panel B ...")
    hb = build_panel_b(ref, data, T, vmax, fs, stem + "B")
    print("Building panel C ...")
    hc = build_panel_c(ref, yeo, colours, long, signed, T, fs, stem + "C")
    print(f"Saved to {out.resolve()}  (B: 180 x {hb:.1f} mm, C: 180 x {hc:.1f} mm)")




# =============================================================================
# ASSEMBLY OF THE THREE PANELS
# =============================================================================
GAP_MM = 3.0
PT_PER_MM = 72 / 25.4


def assemble(panel_pdfs, out_stem, gap_mm=GAP_MM):
    """Stack panels A, B and C (vector PDF + 600-dpi TIFF/PNG)."""
    pages = [PdfReader(str(p)).pages[0] for p in panel_pdfs]
    W = max(float(p.mediabox.width) for p in pages)
    heights = [float(p.mediabox.height) for p in pages]
    gap = gap_mm * PT_PER_MM
    H = sum(heights) + gap * (len(pages) - 1)
    canvas = PageObject.create_blank_page(width=W, height=H)
    y = H
    for p, h in zip(pages, heights):
        y -= h
        canvas.merge_transformed_page(p, Transformation().translate(0, y))
        y -= gap
    w = PdfWriter()
    w.add_page(canvas)
    w.add_metadata({"/Title": "Supplementary Figure S4"})
    with open(f"{out_stem}.pdf", "wb") as fh:
        w.write(fh)

    ims = [Image.open(Path(p).with_suffix(".png")).convert("RGB") for p in panel_pdfs]
    dpi = 600
    gap_px = round(gap_mm / 25.4 * dpi)
    Wpx = max(im.width for im in ims)
    out = Image.new("RGB", (Wpx, sum(im.height for im in ims) + gap_px * (len(ims) - 1)), "white")
    y = 0
    for im in ims:
        out.paste(im, (0, y))
        y += im.height + gap_px
    out.save(f"{out_stem}.tif", dpi=(dpi, dpi), compression="tiff_lzw")
    out.save(f"{out_stem}.png", dpi=(dpi, dpi))
    print(f"Figure S4: {W / PT_PER_MM:.1f} x {H / PT_PER_MM:.1f} mm -> {out_stem}.pdf/.tif/.png")




# =============================================================================
# MAIN: panel A, panels B-C, then one combined figure
# =============================================================================
def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with plt.rc_context(RC_PANEL_A):
        make_panel_a(ICA_FILE, OUTPUT_DIR)
    with plt.rc_context(RC_PANELS_BC):
        make_panels_bc(STATE_DIR, OUTPUT_DIR)
    assemble([OUTPUT_DIR / f"{PANEL_A_STEM}.pdf", OUTPUT_DIR / f"{PANELS_BC_STEM}B.pdf",
              OUTPUT_DIR / f"{PANELS_BC_STEM}C.pdf"], str(OUTPUT_DIR / FINAL_STEM))


if __name__ == "__main__":
    main()
