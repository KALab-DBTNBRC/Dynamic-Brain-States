"""Consensus clustering (Monti et al. 2003) mirroring TbCAPs CAP_ConsensusClustering
with 'items' resampling and MATLAB kmeans(..., 'Distance','correlation') defaults
('Start','plus', 1 replicate, batch updates, EmptyAction 'singleton').

Metrics are computed EXACTLY over all item pairs (i<j) from a table T[a, b] of
pair counts, where a = #resamples in which the pair was clustered together and
b = #resamples in which both items were drawn. Consensus index = a / b.
"""
import numpy as np
from numba import njit, prange


def corr_normalize(X):
    Xc = X - X.mean(axis=1, keepdims=True)
    nrm = np.linalg.norm(Xc, axis=1, keepdims=True)
    nrm[nrm == 0] = 1.0
    return Xc / nrm


def _unit_centroid(v):
    v = v - v.mean()
    n = np.linalg.norm(v)
    return v / n if n > 0 else v


def kmeans_corr(Xn, K, rng, max_iter=100, replicates=1):
    """Xn: rows already centered + unit norm. Returns labels 0..K-1 of the replicate
    with the lowest total within-cluster correlation distance (MATLAB 'Replicates')."""
    best, best_cost = None, np.inf
    for _ in range(replicates):
        lab, cost = _kmeans_corr_once(Xn, K, rng, max_iter)
        if cost < best_cost:
            best, best_cost = lab, cost
    return best


def _kmeans_corr_once(Xn, K, rng, max_iter=100):
    n, d = Xn.shape
    C = np.empty((K, d))
    C[0] = Xn[rng.integers(n)]
    dmin = np.maximum(1.0 - Xn @ C[0], 0.0)
    for k in range(1, K):
        s = dmin.sum()
        idx = rng.choice(n, p=dmin / s) if s > 0 else rng.integers(n)
        C[k] = Xn[idx]
        dmin = np.minimum(dmin, np.maximum(1.0 - Xn @ C[k], 0.0))
    labels = np.full(n, -1)
    for _ in range(max_iter):
        S = Xn @ C.T
        new = S.argmax(axis=1)
        if np.array_equal(new, labels):
            break
        labels = new
        dist = 1.0 - S[np.arange(n), labels]
        for k in range(K):
            m = labels == k
            if not m.any():  # EmptyAction 'singleton': farthest point seeds the empty cluster
                far = int(dist.argmax())
                labels[far] = k
                dist[far] = -np.inf
                m = labels == k
            C[k] = _unit_centroid(Xn[m].mean(axis=0))
    cost = float((1.0 - (Xn * C[labels]).sum(axis=1)).sum())
    return labels, cost


def run_consensus_labels(X, K_range, Pcc=0.8, N=20, seed=0, log=print, replicates=1):
    """Returns dict K -> int8 array (n_items x N); 0 = not drawn, else cluster id 1..K."""
    rng = np.random.default_rng(seed)
    Xn = corr_normalize(np.asarray(X, dtype=np.float64))
    n = Xn.shape[0]
    m = int(np.floor(Pcc * n))
    out = {}
    for K in K_range:
        L = np.zeros((n, N), dtype=np.int8)
        for h in range(N):
            ss = rng.choice(n, m, replace=False)
            L[ss, h] = kmeans_corr(Xn[ss], K, rng, replicates=replicates) + 1
        out[K] = L
        log(f"  K={K}: {N} resamples clustered")
    return out


@njit(parallel=True, cache=True)
def pair_count_table(L, nchunks):
    n, N = L.shape
    T = np.zeros((nchunks, N + 1, N + 1), dtype=np.int64)
    for c in prange(nchunks):
        li = np.empty(N, dtype=np.int8)
        for i in range(c, n, nchunks):
            for r in range(N):
                li[r] = L[i, r]
            for j in range(i + 1, n):
                a = 0
                b = 0
                for r in range(N):
                    x = li[r]
                    y = L[j, r]
                    if x > 0 and y > 0:
                        b += 1
                        if x == y:
                            a += 1
                T[c, a, b] += 1
    return T.sum(axis=0)


def table_metrics(T, lo=0.1, hi=0.9, breaks=100):
    """Exact metrics from pair-count table T[a, b]."""
    N = T.shape[0] - 1
    a_idx, b_idx = np.nonzero(T)
    cnt = T[a_idx, b_idx].astype(np.float64)
    valid = b_idx > 0
    never = cnt[~valid].sum()
    a_idx, b_idx, cnt = a_idx[valid], b_idx[valid], cnt[valid]
    v = a_idx / b_idx
    tot = cnt.sum()
    mean = (cnt * v).sum() / tot
    sd = np.sqrt((cnt * (v - mean) ** 2).sum() / tot)
    # PAC with strict inequalities lo < v < hi (as in the user's MATLAB code), exact rational test
    amb = (a_idx * 10 > b_idx * int(lo * 10)) & (a_idx * 10 < b_idx * int(hi * 10))
    pac = cnt[amb].sum() / tot
    # ConsensusClusterPlus CDF area: R hist() with breaks seq(0,1,1/breaks), right-closed bins,
    # area = sum_b cumulative_fraction(bin b) * bin width
    edges = np.arange(1, breaks + 1)  # right edges as integers e/breaks
    F_edges = np.array([cnt[a_idx * breaks <= e * b_idx].sum() / tot for e in edges])
    area_ccp = (F_edges * (1.0 / breaks)).sum()
    # exact step CDF over distinct values for plotting
    order = np.argsort(v, kind="stable")
    vs, cs = v[order], cnt[order]
    uv, inv = np.unique(vs, return_inverse=True)
    w = np.bincount(inv, weights=cs)
    F = np.cumsum(w) / tot
    return dict(mean=mean, sd=sd, pac=pac, area=area_ccp, n_pairs=tot, never_cosampled=never,
                cdf_x=uv, cdf_F=F, F_edges=F_edges)


def subset_consensus(L, idx):
    """Exact consensus matrix for a subset of items."""
    S = L[idx].astype(np.int16)
    k = len(idx)
    a = np.zeros((k, k), dtype=np.int32)
    b = np.zeros((k, k), dtype=np.int32)
    for r in range(S.shape[1]):
        s = S[:, r]
        drawn = s > 0
        b += np.outer(drawn, drawn)
        a += (s[:, None] == s[None, :]) & np.outer(drawn, drawn)
    with np.errstate(invalid="ignore", divide="ignore"):
        C = a / b
    np.fill_diagonal(C, 1.0)
    return C
