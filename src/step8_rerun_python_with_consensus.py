"""Corrected consensus clustering (Python re-run used for the 2026-11-12 figures).
items = 48,000 prototypes (240 subjects x 200), features = 64 ROIs; K = 2..20;
20 resamples of 80% of items; correlation k-means, k-means++ start, 5 replicates.

Usage (pick ONE run; each reproduces the values in the paper exactly):
  python run_corrected_consensus.py 2026-11-12-X_prototypes_extracted.mat out_dir            # main run (seed 4)
  python run_corrected_consensus.py 2026-11-12-X_prototypes_extracted.mat out_dir --repeat   # independent repeat (seed 99)
  python run_corrected_consensus.py 2026-11-12-X_prototypes_extracted.mat out_dir --null     # ROI-permuted null (seed 5)
Writes consensus_metrics_<real|repeat|null>_python.mat, readable by make_figS2_2026_11_12.py.
Runtime ~30 min per run on 2 CPU cores. Requires numpy, scipy, numba.
"""
import sys, os, time
import numpy as np
import scipy.io as sio
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from consensus_core import run_consensus_labels, pair_count_table, table_metrics, subset_consensus

src, out_dir = sys.argv[1], sys.argv[2]
null = "--null" in sys.argv
repeat = "--repeat" in sys.argv
if null and repeat:
    sys.exit("Choose either --null or --repeat, not both.")
tag, seed = ("null", 5) if null else (("repeat", 99) if repeat else ("real", 4))
os.makedirs(out_dir, exist_ok=True)
P = sio.loadmat(src)["P"]                       # 240 x 200 x 64
X = P.reshape(-1, P.shape[2])                   # 48000 x 64, subject-major
if null:
    rng = np.random.default_rng(2026)
    X = np.column_stack([rng.permutation(X[:, j]) for j in range(X.shape[1])])
K_range = list(range(2, 21))
t0 = time.time()
labels = run_consensus_labels(X, K_range, Pcc=0.8, N=20, seed=seed, replicates=5,
                              log=lambda s: print(s, f"[{time.time()-t0:.0f}s]", flush=True))
sub_idx = np.sort(np.random.default_rng(7).choice(X.shape[0], 2000, replace=False))
mt, tables = [], []
for K in K_range:
    T = pair_count_table(labels[K], 8); tables.append(T); mt.append(table_metrics(T))
    print(f"K={K:2d} PAC={mt[-1]['pac']:.4f} area={mt[-1]['area']:.4f}", flush=True)
area = np.array([m["area"] for m in mt])
cdf_x = np.empty(len(K_range), dtype=object); cdf_F = np.empty(len(K_range), dtype=object)
for i, m in enumerate(mt):
    cdf_x[i] = m["cdf_x"]; cdf_F[i] = m["cdf_F"]
out = dict(K_range=np.array(K_range), PAC=np.array([m["pac"] for m in mt]), CDF_area=area,
           DeltaK=np.r_[area[0], (area[1:] - area[:-1]) / area[:-1]],
           meanC=np.array([m["mean"] for m in mt]), sdC=np.array([m["sd"] for m in mt]),
           F_edges=np.column_stack([m["F_edges"] for m in mt]), cdf_x=cdf_x, cdf_F=cdf_F,
           T_all=np.stack(tables, axis=2).astype(float), subIdx=sub_idx + 1, N=20, Pcc=0.8, kmeans_replicates=5,
           seed=seed, run=tag)
if not null:
    K_show = [2, 3, 4, 5, 6]
    out["K_show"] = np.array(K_show)
    out["subC"] = np.stack([subset_consensus(labels[K], sub_idx).astype(np.float32) for K in K_show], axis=2)
    out["L_all"] = np.stack([labels[K] for K in K_range], axis=2)
sio.savemat(os.path.join(out_dir, f"consensus_metrics_{tag}_python.mat"), out, do_compression=True)
print("done in %.1f s" % (time.time() - t0))
