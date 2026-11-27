"""Supplementary Figure S3: IC300 CAPs vs IC100 centroids.
Usage: python 2026-11-14-make_figS3.py <state_1..4 dir (IC300)> <centroidsIC100_4D.nii.gz> <out_dir>
A: IC300 CAPs 1-4 (left lateral fsaverage), B: IC100 centroids 1-4 (same view/scale),
C: Pearson spatial correlation between them (whole-brain volumes on the 2 mm MNI grid).
Each map is z-scored over its non-zero voxels so all eight maps share one colour scale.
"""
import os, sys
import numpy as np
import nibabel as nib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from nilearn import datasets, surface
from matplotlib.collections import PolyCollection
from nilearn.image import resample_to_img, index_img

state_dir, c100_path, out = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(out, exist_ok=True)
tmp = os.path.join(out, "_panels"); os.makedirs(tmp, exist_ok=True)

INK, INK2 = "#0b0b0b", "#52514e"
# Figure 1 style: cyan -> dark blue (negative), red -> yellow (positive), no black midpoint
CAPMAP = LinearSegmentedColormap.from_list(
    "cap", [(0.0, "#00f0ff"), (0.2, "#1aa3ff"), (0.4, "#1f5cff"), (0.4999, "#2a4fe0"),
            (0.5001, "#e8230c"), (0.68, "#f2420f"), (0.86, "#ff9a00"), (1.0, "#ffe600")])
DIV = LinearSegmentedColormap.from_list("div", ["#104281", "#3987e5", "#f0efec", "#e66767", "#a31f1f"])
VMAX = 2.5
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "pdf.fonttype": 42})


def zscore_img(img):
    d = img.get_fdata(dtype=np.float32)
    m = d != 0
    z = np.zeros_like(d)
    z[m] = (d[m] - d[m].mean()) / d[m].std()
    return nib.Nifti1Image(z, img.affine), m


c100 = nib.load(c100_path)
ref = index_img(c100, 0)
maps300 = [zscore_img(nib.load(os.path.join(state_dir, f"state_{k}.nii.gz")))[0] for k in range(1, 5)]
maps100 = [zscore_img(index_img(c100, k))[0] for k in range(4)]

# ---- Panel C: spatial correlation on the 2 mm grid, voxels non-zero in all eight maps ----
v300 = [resample_to_img(m, ref, interpolation="continuous", force_resample=True, copy_header=True).get_fdata() for m in maps300]
v100 = [m.get_fdata() for m in maps100]
mask = np.all(np.stack([v != 0 for v in v300 + v100]), axis=0)
R = np.array([[np.corrcoef(a[mask], b[mask])[0, 1] for b in v100] for a in v300])
np.savetxt(os.path.join(out, "2026-11-14-FigS3_correlation_IC300_vs_IC100.csv"), R, delimiter=",", fmt="%.4f",
           header="rows: IC300 CAP1-4; columns: IC100 centroid 1-4; Pearson r over %d voxels (2 mm MNI)" % mask.sum())

# ---- Surface renders (left lateral, fsaverage pial, orthographic view with soft lighting) ----
fs = datasets.load_fsaverage("fsaverage")
pial, white = fs["pial"].parts["left"], fs["white_matter"].parts["left"]
PIAL_XYZ, PIAL_F = np.asarray(pial.coordinates, float), np.asarray(pial.faces)
def face_normals(coords, faces):
    v0, v1, v2 = coords[faces[:, 0]], coords[faces[:, 1]], coords[faces[:, 2]]
    n = np.cross(v1 - v0, v2 - v0)
    n /= np.linalg.norm(n, axis=1, keepdims=True) + 1e-12
    # make normals point outward
    cen = coords.mean(0)
    if np.mean(np.sum(n * ((v0 + v1 + v2) / 3 - cen), axis=1)) < 0:
        n = -n
    return n


def vertex_normals(coords, faces, fn):
    vn = np.zeros_like(coords)
    for k in range(3):
        np.add.at(vn, faces[:, k], fn)
    vn /= np.linalg.norm(vn, axis=1, keepdims=True) + 1e-12
    return vn


def render_lateral_left(coords, faces, values, cmap, vmin, vmax, path, width_in=4.0, dpi=300,
                        light=(-1.0, 0.35, 0.65), ambient=0.38, diffuse=0.62, spec=0.12):
    fn = face_normals(coords, faces)
    vn = vertex_normals(coords, faces, fn)
    L = np.asarray(light, float); L /= np.linalg.norm(L)
    V = np.array([-1.0, 0.0, 0.0])                    # viewer looks from -x (left lateral)
    H = (L + V) / np.linalg.norm(L + V)
    lam = np.clip(vn @ L, 0, 1)
    sp = np.clip(vn @ H, 0, 1) ** 30
    shade_v = ambient + diffuse * lam
    # per-face colour = mean of vertex colours, shaded
    norm = np.clip((values - vmin) / (vmax - vmin), 0, 1)
    vcol = cmap(norm)[:, :3] * shade_v[:, None] + spec * sp[:, None]
    fcol = np.clip(vcol[faces].mean(axis=1), 0, 1)
    visible = fn @ V > 0.0                             # faces turned toward the viewer
    f = faces[visible]; fc = fcol[visible]
    depth = coords[f].mean(axis=1)[:, 0]               # larger x = farther from viewer
    order = np.argsort(-depth)
    f, fc = f[order], fc[order]
    xy = np.stack([-coords[:, 1], coords[:, 2]], axis=1)   # anterior to the left, superior up
    polys = xy[f]
    xr = np.ptp(xy[:, 0]); yr = np.ptp(xy[:, 1])
    fig = plt.figure(figsize=(width_in, width_in * yr / xr))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_axis_off()
    pc = PolyCollection(polys, facecolors=fc, edgecolors=fc, linewidths=0.15, antialiaseds=False)
    ax.add_collection(pc)
    ax.set_xlim(xy[:, 0].min() - 1, xy[:, 0].max() + 1); ax.set_ylim(xy[:, 1].min() - 1, xy[:, 1].max() + 1)
    ax.set_aspect("equal")
    fig.savefig(path, dpi=dpi, facecolor="white")
    plt.close(fig)


def render(img, path):
    tex = np.nan_to_num(surface.vol_to_surf(img, pial, inner_mesh=white))
    render_lateral_left(PIAL_XYZ, PIAL_F, tex, CAPMAP, -VMAX, VMAX, path, width_in=4.0, dpi=300)


paths300, paths100 = [], []
for k in range(4):
    p = os.path.join(tmp, f"ic300_{k+1}.png"); render(maps300[k], p); paths300.append(p)
    p = os.path.join(tmp, f"ic100_{k+1}.png"); render(maps100[k], p); paths100.append(p)

# ---- Compose ----
from PIL import Image


def cropped(path):
    im = Image.open(path).convert("RGB"); a = np.asarray(im)
    nonwhite = np.where(a.min(axis=2) < 245)
    y0, y1, x0, x1 = nonwhite[0].min(), nonwhite[0].max(), nonwhite[1].min(), nonwhite[1].max()
    return np.asarray(im.crop((x0, y0, x1 + 1, y1 + 1)))


FW, FH = 7.2, 4.3
fig = plt.figure(figsize=(FW, FH))
colw, gap, x0 = 0.168, 0.006, 0.01
rows = [(paths300, 0.675, "A  IC300 (300-component group-ICA): CAPs in the twin cohort", "CAP"),
        (paths100, 0.27, "B  IC100 (100-component group-ICA): centroids", "Centroid")]
imh = 0.235
for paths, y, title, label in rows:
    fig.text(x0, y + imh + 0.075, title, fontsize=9, fontweight="bold", ha="left", va="bottom")
    for k in range(4):
        ax = fig.add_axes([x0 + k * (colw + gap), y, colw, imh])
        ax.imshow(cropped(paths[k])); ax.axis("off")
        ax.set_title(f"{label} {k+1}", fontsize=8.5, fontweight="bold", color=INK, pad=3)
# colour bar for maps
cax = fig.add_axes([0.13, 0.15, 0.44, 0.03])
cb = fig.colorbar(plt.cm.ScalarMappable(cmap=CAPMAP, norm=plt.Normalize(-VMAX, VMAX)), cax=cax, orientation="horizontal")
cb.outline.set_visible(False); cb.set_ticks([-VMAX, 0, VMAX]); cb.ax.tick_params(labelsize=7, length=2)
cb.set_label("Map value (z, standardized within each map)", fontsize=7, labelpad=2)
fig.text(0.125, 0.165, "Negatively activated", fontsize=6.8, color=INK2, ha="right", va="center")
fig.text(0.575, 0.165, "Positively activated", fontsize=6.8, color=INK2, ha="left", va="center")

# ---- Panel C heatmap ----
hw = 0.205; hh = hw * FW / FH
axc = fig.add_axes([0.775, 0.40, hw, hh])
im = axc.imshow(R, cmap=DIV, vmin=-1, vmax=1)
for i in range(4):
    for j in range(4):
        best = j == int(np.argmax(R[i]))
        axc.text(j, i, f"{R[i, j]:.2f}", ha="center", va="center", fontsize=7.5,
                 color="white" if abs(R[i, j]) > 0.55 else INK, fontweight="bold" if best else "normal")
axc.set_xticks(range(4)); axc.set_xticklabels([f"C{j+1}" for j in range(4)], fontsize=7)
axc.set_yticks(range(4)); axc.set_yticklabels([f"CAP {i+1}" for i in range(4)], fontsize=7)
axc.set_xlabel("IC100 centroid", fontsize=7.5, labelpad=2); axc.set_ylabel("IC300 CAP", fontsize=7.5)
axc.tick_params(length=0)
for sp in axc.spines.values(): sp.set_visible(False)
fig.text(0.70, 0.40 + hh + 0.07, "C  Spatial correlation", fontsize=9, fontweight="bold", ha="left", va="bottom")
cax2 = fig.add_axes([0.775, 0.30, hw, 0.025])
cb2 = fig.colorbar(im, cax=cax2, orientation="horizontal"); cb2.set_ticks([-1, 0, 1])
cb2.set_label("Pearson r (bold = best match per CAP)", fontsize=6.5, labelpad=2); cb2.outline.set_visible(False)
cb2.ax.tick_params(labelsize=7, length=2)

for ext in ("png", "pdf"):
    fig.savefig(os.path.join(out, f"2026-11-14-FigS3_IC300_vs_IC100_centroids.{ext}"), dpi=300, bbox_inches="tight", facecolor="white")
print(np.round(R, 2)); print("voxels:", int(mask.sum()))
