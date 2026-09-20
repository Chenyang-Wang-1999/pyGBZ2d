'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-08-27
Copyright © Department of Physics, Tsinghua University. All rights reserved

EXPERIMENTAL — research-grade code, no stability guarantee.

Band clustering of GBZ sweep point clouds by radius-graph connected
components (single-linkage; DBSCAN with ``min_samples=1``).

Idea
----
Bands whose gaps never close are disjoint sheets in the full
``(E, mu1, theta1, mu2, theta2)`` space, while points of one sheet are
chained at the scale of the E-grid step.  A radius between those two scales
recovers each sheet as one connected component.

Two geometry facts baked into the embedding:

* Angles enter as ``(cos, sin)`` so the Euclidean metric of ``cKDTree``
  respects the theta = 0 = 2pi seam (chord ~ arc at small distance).
* E enters SCALED DOWN by the grid step (factor ``ALPHA_E``): the E block
  must be small enough that adjacent slices of one sheet stay within eps,
  but big enough that two sheets crossing in the (mu, theta) projection at
  different energies stay separated (that separation is ``alpha_E * gap/dE``
  and the well-separated assumption promises gap >> dE).

No decimation: line points are used at the solver's native sampling
density.  A theta1-width-based stride was tried and REJECTED — it deletes
the body of near-vertical arcs (theta1 ~ const, theta2 sweeping) and breaks
chaining; see log/2026-08-27.  The performance cost is paid instead by the
scan cap below.

Scan cap: the eps stability ladder stops early at the first radius where
the cloud is fully merged — beyond that point larger radii only multiply
the query_pairs output (quadratically many pairs on dense 1-D E-families)
without adding information.  The absolute cap EPS_SCAN_MAX is a safety
net on top of the early stop.
'''

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from ..core import PointSubset, LineSubset, TWO_PI

# ---------------------------------------------------------------------------
# Module constants (single-consumer, home-module locality per doc/constants.md
# convention; assign to override process-wide, pass kwargs per call).
# ---------------------------------------------------------------------------

#: E-block downscale: one E-grid step contributes this many distance units.
#: 0.25 keeps the chaining floor well below the inter-band separations seen
#: in the playground validation datasets (0.44 - 12.5 in scaled units).
ALPHA_E: float = 0.25

#: eps ladder bounds / resolution for the stability scan.  The cap is a
#: safety net only — the early stop at the first fully-merged radius is the
#: real blowup guard for dense 1-D E-families, and legitimate band
#: separations can exceed the old 1.5 cap (y-SGBZ folded bands merge only
#: around 1.3-3.0).
EPS_SCAN_MIN: float = 0.02
EPS_SCAN_MAX: float = 3.0
N_SCAN_STEPS: int = 21

#: A plateau of the cluster-count curve counts as a clustering window only
#: if its largest cluster holds at least this fraction of all points —
#: rejects the trivial all-singletons plateaus at tiny eps.
MIN_LARGEST_FRAC: float = 0.1

#: How many of the biggest clusters enter the pairwise margin computation.
MAX_MARGIN_PAIRS: int = 12

#: kNN order for the chaining-floor diagnostic (k=2: nearest non-self).
KNN_K: int = 2

#: Clusters smaller than this are reported separately by summarize_clusters.
MIN_CLUSTER_SIZE: int = 5

#: A plateau needs at least this many clusters to count as a window.
MIN_CLUSTERS_PLATEAU: int = 2


# ---------------------------------------------------------------------------
# 1. Flattening: GBZResult list -> point cloud in (E, mu, theta) coordinates
# ---------------------------------------------------------------------------

@dataclass
class BandPoints:
    """Flat GBZ point cloud with provenance back into the source GBZResults.

    Attributes:
        E: (n,) complex reference energy of each point.
        mu1: (n,) ln|beta1|.
        theta1: (n,) arg(beta1) in [0, 2pi).
        mu2: (n,) ln|beta2|.
        theta2: (n,) arg(beta2) in [0, 2pi).
        slice_idx: (n,) index into the source results list.
        subset_idx: (n,) subset index within its result.
        is_line: (n,) True if the point came from a LineSubset sample.
        n_dropped: points dropped for non-finite mu/theta (0/infinity beta2
            padding leaking into subsets).  Recorded instead of silently
            vanishing.
    """
    E: np.ndarray
    mu1: np.ndarray
    theta1: np.ndarray
    mu2: np.ndarray
    theta2: np.ndarray
    slice_idx: np.ndarray
    subset_idx: np.ndarray
    is_line: np.ndarray
    n_dropped: int = 0


def flatten_results(results) -> BandPoints:
    """Flatten every PointSubset / LineSubset sample of a sweep into a cloud.

    Failed and empty GBZResults are skipped.  LineSubsets contribute ALL of
    their sampled points (see module docstring for why no decimation).
    """
    buf = {k: [] for k in
           ("E", "mu1", "theta1", "mu2", "theta2",
            "slice_idx", "subset_idx", "is_line")}
    n_dropped = 0

    for i, res in enumerate(results):
        if not (getattr(res, "success", False) and res.subsets):
            continue
        for j, sub in enumerate(res.subsets):
            if isinstance(sub, PointSubset):
                E = np.asarray([sub.E], dtype=complex)
                mu1 = np.asarray([np.log(abs(sub.beta1))])
                th1 = np.asarray([np.angle(sub.beta1)])
                mu2 = np.asarray([np.log(abs(sub.beta2))])
                th2 = np.asarray([np.angle(sub.beta2)])
                is_line = False
            elif isinstance(sub, LineSubset):
                E = np.full(len(sub.theta1_arr), complex(sub.E))
                mu1 = np.full(len(sub.theta1_arr), float(sub.mu1))
                th1 = np.asarray(sub.theta1_arr, dtype=float) % TWO_PI
                # log(0) = -inf on leaked padding roots is EXPECTED — the
                # non-finite filter below is the handler, so the warning is
                # noise; suppress it locally rather than pre-screening.
                with np.errstate(divide="ignore", invalid="ignore"):
                    mu2 = np.log(np.abs(sub.beta2_arr))
                th2 = np.angle(sub.beta2_arr)
                is_line = True
            else:
                raise TypeError(f"unknown subset type: {type(sub)}")

            finite = (np.isfinite(mu1) & np.isfinite(mu2)
                      & np.isfinite(th1) & np.isfinite(th2))
            n_dropped += int((~finite).sum())
            buf["E"].append(E[finite])
            buf["mu1"].append(mu1[finite])
            buf["theta1"].append(th1[finite] % TWO_PI)
            buf["mu2"].append(mu2[finite])
            buf["theta2"].append(th2[finite] % TWO_PI)
            buf["slice_idx"].append(np.full(int(finite.sum()), i, dtype=int))
            buf["subset_idx"].append(np.full(int(finite.sum()), j, dtype=int))
            buf["is_line"].append(np.full(int(finite.sum()), is_line))

    out = {k: np.concatenate(v) if v else np.array([]) for k, v in buf.items()}
    return BandPoints(n_dropped=n_dropped, **out)


# ---------------------------------------------------------------------------
# 2. Embedding: (E, mu, theta) -> R^8 with torus-aware angle blocks
# ---------------------------------------------------------------------------

def grid_steps(E: np.ndarray) -> tuple[float, float]:
    """Median grid step of the sweep in (Re E, Im E).

    Per-axis because playground sweeps are strongly anisotropic (e.g. 7.7 x
    1.0 extent over the same 201 points); a single scale would let one axis
    dominate.  A 1-D sweep (constant Im) falls back to the Re step.
    """
    re_u = np.unique(E.real)
    im_u = np.unique(E.imag)
    d_re = float(np.median(np.diff(re_u))) if len(re_u) > 1 else 1.0
    d_im = float(np.median(np.diff(im_u))) if len(im_u) > 1 else d_re
    return d_re, d_im


def embed(bp: BandPoints, *, alpha_E: Optional[float] = None,
          d_re: Optional[float] = None, d_im: Optional[float] = None,
          w_mu: Optional[float] = None) -> np.ndarray:
    """Map the cloud into R^8 for cKDTree.

    Blocks: ``[alpha_E*(ReE/dRe, ImE/dIm), w_mu*(mu1, mu2),
    (cos t1, sin t1), (cos t2, sin t2)]``.

    Parameters:
        bp: the point cloud.
        alpha_E: E-block downscale; ``None`` reads :data:`ALPHA_E`.
        d_re, d_im: grid steps; ``None`` infers from the cloud's energies.
        w_mu: mu-block scale; ``None`` auto-scales to 1/max(mu span, 1).
    """
    if alpha_E is None:
        alpha_E = ALPHA_E
    if d_re is None or d_im is None:
        d_re_, d_im_ = grid_steps(bp.E)
        d_re = d_re if d_re is not None else d_re_
        d_im = d_im if d_im is not None else d_im_
    if w_mu is None:
        span = max(np.ptp(bp.mu1), np.ptp(bp.mu2), 1.0) if len(bp.mu1) else 1.0
        w_mu = 1.0 / span

    X = np.empty((len(bp.E), 8))
    X[:, 0] = alpha_E * bp.E.real / d_re
    X[:, 1] = alpha_E * bp.E.imag / d_im
    X[:, 2] = w_mu * bp.mu1
    X[:, 3] = w_mu * bp.mu2
    X[:, 4] = np.cos(bp.theta1)
    X[:, 5] = np.sin(bp.theta1)
    X[:, 6] = np.cos(bp.theta2)
    X[:, 7] = np.sin(bp.theta2)
    return X


# ---------------------------------------------------------------------------
# 3. Clustering: radius graph + connected components
# ---------------------------------------------------------------------------

def radius_graph_labels(X: np.ndarray, eps: float) -> np.ndarray:
    """Single-linkage labels: connected components of the eps-ball graph."""
    pairs = cKDTree(X).query_pairs(eps, output_type="ndarray")
    graph = csr_matrix(
        (np.ones(len(pairs), dtype=np.int8), (pairs[:, 0], pairs[:, 1])),
        shape=(len(X), len(X)),
    )
    _, labels = connected_components(graph, directed=False)
    return labels


def knn_distance_stats(X: np.ndarray, k: Optional[int] = None) -> np.ndarray:
    """Distances to the k-th neighbour (default k=2: nearest non-self)."""
    if k is None:
        k = KNN_K
    d, _ = cKDTree(X).query(X, k=k, workers=-1)
    return d[:, -1]


def eps_stability_scan(X: np.ndarray, eps_list: np.ndarray) -> list[dict]:
    """Cluster count vs eps — the diagnostic that verifies the well-separated
    assumption (a usable window shows up as a wide plateau).

    Stops early at the first radius where the whole cloud is one component:
    beyond the full-merge radius nothing new appears, and on dense 1-D
    E-families the pair count explodes quadratically.
    """
    rows = []
    for eps in eps_list:
        labels = radius_graph_labels(X, float(eps))
        sizes = np.bincount(labels)
        rows.append({
            "eps": float(eps),
            "n_clusters": int(len(sizes)),
            "largest": int(sizes.max()),
            "n_singletons": int((sizes == 1).sum()),
            "n_pts": int(len(X)),
        })
        if len(sizes) == 1 and len(X) > 1:
            break
    return rows


def widest_plateau(rows: list[dict], min_clusters: Optional[int] = None,
                   min_largest_frac: Optional[float] = None):
    """Widest plateau (in log-eps) of the cluster-count curve.

    A plateau counts only if it has at least ``min_clusters`` (default
    :data:`MIN_CLUSTERS_PLATEAU`) clusters AND its largest cluster holds
    >= ``min_largest_frac`` (default :data:`MIN_LARGEST_FRAC`) of all
    points.  Returns ``(lo, hi, n_clusters)`` or ``None``.
    """
    if min_clusters is None:
        min_clusters = MIN_CLUSTERS_PLATEAU
    if min_largest_frac is None:
        min_largest_frac = MIN_LARGEST_FRAC
    n_pts = rows[0]["n_pts"]
    best = None
    i = 0
    while i < len(rows):
        j = i
        while j + 1 < len(rows) and rows[j + 1]["n_clusters"] == rows[i]["n_clusters"]:
            j += 1
        k = rows[i]["n_clusters"]
        if k >= min_clusters and rows[i]["largest"] >= min_largest_frac * n_pts:
            lo, hi = rows[i]["eps"], rows[j]["eps"]
            width = np.log(hi) - np.log(lo)
            if best is None or width > best[2]:
                best = (lo, hi, width, k)
        i = j + 1
    if best is None:
        return None
    return best[0], best[1], best[3]


def inter_cluster_margins(X: np.ndarray, labels: np.ndarray,
                          max_pairs: Optional[int] = None) -> list[tuple]:
    """Exact minimum distance between pairs of the biggest clusters.

    The ratio margin/eps is the safety factor of the clustering: how much
    further than the chaining radius the bands actually stay apart.
    """
    if max_pairs is None:
        max_pairs = MAX_MARGIN_PAIRS
    sizes = np.bincount(labels)
    big = np.argsort(-sizes)[:max_pairs]
    margins = []
    for a in range(len(big)):
        for b in range(a + 1, len(big)):
            Xa, Xb = X[labels == big[a]], X[labels == big[b]]
            d, _ = cKDTree(Xa).query(Xb, k=1, workers=-1)
            margins.append((int(big[a]), int(big[b]), float(d.min())))
    return sorted(margins, key=lambda t: t[2])


def summarize_clusters(bp: BandPoints, labels: np.ndarray,
                       min_size: Optional[int] = None) -> list[dict]:
    """Per-cluster statistics (sizes, slice counts, footprint ranges)."""
    if min_size is None:
        min_size = MIN_CLUSTER_SIZE
    stats = []
    for lab in range(int(labels.max()) + 1):
        m = labels == lab
        if not m.any():
            continue
        per_slice = np.bincount(bp.slice_idx[m])
        stats.append({
            "label": int(lab),
            "size": int(m.sum()),
            "n_slices": int(len(np.unique(bp.slice_idx[m]))),
            "pts_per_slice_med": float(np.median(per_slice[per_slice > 0])),
            "frac_line": float(bp.is_line[m].mean()),
            "E_re": (float(bp.E[m].real.min()), float(bp.E[m].real.max())),
            "E_im": (float(bp.E[m].imag.min()), float(bp.E[m].imag.max())),
            "mu1": (float(bp.mu1[m].min()), float(bp.mu1[m].max())),
            "mu2": (float(bp.mu2[m].min()), float(bp.mu2[m].max())),
        })
    stats.sort(key=lambda s: -s["size"])
    return [s for s in stats if s["size"] >= min_size] + \
           [s for s in stats if s["size"] < min_size]


# ---------------------------------------------------------------------------
# 4. Top-level entry point
# ---------------------------------------------------------------------------

@dataclass
class BandClustering:
    """Result of :func:`cluster_bands`.

    Attributes:
        points: the flattened cloud (BandPoints).
        X: the R^8 embedding the clustering ran on.
        labels: (n,) connected-component id per point.
        eps: the radius actually used.
        eps_window: (lo, hi) of the plateau ``eps`` was chosen from,
            or ``None`` when no plateau existed (fallback path).
        n_clusters: number of clusters.
        scan: the eps-stability ladder (may be empty when eps was explicit).
        margins: (cluster_a, cluster_b, min_dist) for the biggest cluster
            pairs, sorted by distance.
        knn_percentiles: nearest-neighbour distance percentiles of the
            embedding (the chaining-floor scale).
    """
    points: BandPoints
    X: np.ndarray
    labels: np.ndarray
    eps: float
    eps_window: Optional[tuple[float, float]]
    n_clusters: int
    scan: list[dict]
    margins: list[tuple]
    knn_percentiles: dict


def cluster_bands(results=None, *, points: Optional[BandPoints] = None,
                  eps: Optional[float] = None,
                  alpha_E: Optional[float] = None,
                  w_mu: Optional[float] = None,
                  max_scan_eps: Optional[float] = None,
                  n_scan_steps: Optional[int] = None,
                  d_re: Optional[float] = None,
                  d_im: Optional[float] = None) -> BandClustering:
    """Cluster an E-sweep of GBZ results into bands (radius graph).

    Provide exactly one of ``results`` (a list of :class:`GBZResult`) or
    ``points`` (a pre-flattened :class:`BandPoints`).

    Parameters:
        eps: explicit radius; ``None`` selects the geometric middle of the
            widest plateau of the eps-stability scan.
        alpha_E: E-block downscale (default :data:`ALPHA_E`).
        w_mu: mu-block scale (default auto: 1/max(mu span, 1)).
        max_scan_eps, n_scan_steps: eps ladder bounds
            (defaults :data:`EPS_SCAN_MAX` / :data:`N_SCAN_STEPS`).
        d_re, d_im: E-block grid steps.  Pass these whenever the cloud is
            NOT a plain rectangular sweep — e.g. boundary-enriched data,
            whose probe energies sit between grid nodes and would dilute
            the inferred median step, inflating the E block and
            fragmenting the clustering.  ``None`` infers from the cloud.

    Returns:
        BandClustering with labels and diagnostics.
    """
    if (results is None) == (points is None):
        raise ValueError("provide exactly one of 'results' and 'points'")
    bp = points if points is not None else flatten_results(results)
    if len(bp.E) == 0:
        raise ValueError("empty point cloud (no in-GBZ results?)")

    if max_scan_eps is None:
        max_scan_eps = EPS_SCAN_MAX
    if n_scan_steps is None:
        n_scan_steps = N_SCAN_STEPS

    X = embed(bp, alpha_E=alpha_E, w_mu=w_mu, d_re=d_re, d_im=d_im)
    knn = knn_distance_stats(X)
    knn_pct = {f"p{q}": float(np.percentile(knn, q))
               for q in (10, 25, 50, 75, 90, 99)}

    eps_window = None
    scan = []
    if eps is None:
        eps_list = np.geomspace(min(EPS_SCAN_MIN, max_scan_eps),
                                max_scan_eps, n_scan_steps)
        scan = eps_stability_scan(X, eps_list)
        plat = widest_plateau(scan)
        if plat is not None:
            eps_window = (plat[0], plat[1])
            eps = float(np.sqrt(plat[0] * plat[1]))  # geometric middle
        elif any(r["n_clusters"] == 1 for r in scan):
            # No >=2-cluster plateau: the cloud merges before separating.
            # Cluster at the first fully-merged radius — the honest answer
            # is "one component", not a fake split.
            eps = next(r["eps"] for r in scan if r["n_clusters"] == 1)
        else:
            eps = scan[-1]["eps"]

    labels = radius_graph_labels(X, float(eps))
    n_clusters = int(labels.max()) + 1
    margins = inter_cluster_margins(X, labels) if n_clusters >= 2 else []

    return BandClustering(points=bp, X=X, labels=labels, eps=float(eps),
                          eps_window=eps_window, n_clusters=n_clusters,
                          scan=scan, margins=margins, knn_percentiles=knn_pct)
