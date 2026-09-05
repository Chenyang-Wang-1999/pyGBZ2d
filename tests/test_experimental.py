"""Tests for pygbz2d.experimental.band_clustering.

All fixtures are SYNTHETIC point clouds — no solver runs, so the suite
stays fast.  Each fixture targets one geometric promise of the embedding:

* mu-separated sheets: identical (E, theta) trajectories, offset ONLY in
  mu -> the mu block must do the separating (control: w_mu=0 must fail);
* near-vertical arcs: theta1 ~ const, theta2 sweeping — the shape that a
  rejected theta1-width decimation used to destroy;
* seam crossing: theta1 trajectories winding through 0 = 2pi must not
  split (the (cos, sin) embedding is load-bearing).
"""

import numpy as np
import pytest

from pygbz2d.core import GBZResult, PointSubset, LineSubset, TWO_PI
from pygbz2d.experimental import (
    BandPoints,
    cluster_bands,
    embed,
    eps_stability_scan,
    flatten_results,
    grid_steps,
    radius_graph_labels,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _mu_separated_cloud(n_re=25, n_im=25, n_branch=3, d_mu=1.5, seed=0):
    """Two sheets with IDENTICAL (E, theta) trajectories, offset only in mu.

    Branch angles are shared between the sheets, so at w_mu = 0 same-branch
    points coincide exactly — the control test then MUST merge the sheets.
    """
    rng = np.random.default_rng(seed)
    re = np.linspace(-2.0, 2.0, n_re)
    im = np.linspace(-0.5, 0.5, n_im)
    RE, IM = np.meshgrid(re, im)
    E = (RE + 1j * IM).ravel()
    u, v = E.real / 4.0, E.imag
    n = len(E)

    chunks = []
    shared = []  # (th1, th2, warp1, warp2) per branch, shared by both sheets
    for b in range(n_branch):
        shared.append((
            (2 * np.pi * b / n_branch + 0.6 * np.sin(2 * u)
             + 0.5 * v + 0.02 * rng.standard_normal(n)) % TWO_PI,
            (2 * np.pi * (b + 1) / n_branch + 0.5 * np.cos(3 * u)
             + 0.4 * np.sin(4 * v)) % TWO_PI,
            0.15 * np.sin(3 * u + b) + 0.05 * v,
            0.12 * np.cos(2 * u + b) + 0.04 * np.sin(v),
        ))
    for off in (-d_mu / 2.0, +d_mu / 2.0):
        for th1, th2, w1, w2 in shared:
            chunks.append(dict(
                E=E, theta1=th1, theta2=th2,
                mu1=off + w1, mu2=0.8 * off + w2,
                slice_idx=np.arange(n),
                subset_idx=np.full(n, len(chunks)),
                is_line=np.zeros(n, dtype=bool),
            ))
    return BandPoints(**{k: np.concatenate([c[k] for c in chunks])
                         for k in chunks[0]})


def _vertical_arc_cloud(n_E=40, n_th2=100, separation=1.2):
    """Two bands of near-vertical arcs (theta1 ~ const, theta2 sweeping),
    offset by ``separation`` in theta1.  The Hermitian-Haldane failure
    shape: chaining must happen along the arc BODY, not just endpoints."""
    Es = np.linspace(-1.0, 1.0, n_E)
    th2 = np.linspace(0.0, TWO_PI, n_th2, endpoint=False)
    chunks = []
    for band, c0 in enumerate((0.3, 0.3 + separation)):
        for j, E in enumerate(Es):
            th1 = c0 + 0.02 * np.sin(3 * th2) + 0.004 * j
            chunks.append(dict(
                E=np.full(n_th2, E), theta1=th1 % TWO_PI, theta2=th2,
                mu1=np.zeros(n_th2), mu2=0.1 * np.sin(th2),
                slice_idx=np.full(n_th2, j + band * n_E),
                subset_idx=np.full(n_th2, j),
                is_line=np.ones(n_th2, dtype=bool),
            ))
    return BandPoints(**{k: np.concatenate([c[k] for c in chunks])
                         for k in chunks[0]})


def _seam_crossing_cloud(n_E=60, n_th2=80):
    """One band whose theta1 drifts DOWNWARD through 1.2 rad — it crosses
    the theta = 0 = 2pi seam near j = n_E/4.  Must stay a single cluster
    (a raw-radian embedding would tear it apart at the seam)."""
    Es = np.linspace(-1.0, 1.0, n_E)
    th2 = np.linspace(0.0, TWO_PI, n_th2, endpoint=False)
    chunks = []
    for j, E in enumerate(Es):
        th1 = np.full(n_th2, (0.3 - 1.2 * j / (n_E - 1)) % TWO_PI)
        chunks.append(dict(
            E=np.full(n_th2, E), theta1=th1, theta2=th2,
            mu1=np.zeros(n_th2), mu2=np.zeros(n_th2),
            slice_idx=np.full(n_th2, j),
            subset_idx=np.full(n_th2, j),
            is_line=np.ones(n_th2, dtype=bool),
        ))
    return BandPoints(**{k: np.concatenate([c[k] for c in chunks])
                         for k in chunks[0]})


def _fake_results():
    """One PointSubset + one clean + one dirty LineSubset + junk results."""
    line = LineSubset(E=1.0 + 0j, mu1=0.1,
                      theta1_arr=np.linspace(0.0, 1.0, 10),
                      beta2_arr=np.exp(0.2 + 1j * np.linspace(0.0, 1.0, 10)))
    dirty = LineSubset(E=2.0 + 0j, mu1=0.1,
                       theta1_arr=np.linspace(0.0, 1.0, 5),
                       beta2_arr=np.exp(0.2 + 1j * np.linspace(0.0, 1.0, 5)))
    dirty.beta2_arr[3] = 0.0  # ln|0| = -inf -> must be dropped, counted
    point = PointSubset(E=1.0 + 0j, beta1=np.exp(0.1 + 0.3j),
                        beta2=np.exp(0.2 + 0.4j))
    ok = GBZResult(E_ref=1.0 + 0j, subsets=[point, line], index=(1, 1))
    ok2 = GBZResult(E_ref=2.0 + 0j, subsets=[dirty], index=(0, 1))
    failed = GBZResult(E_ref=3.0 + 0j, success=False, error="x")
    empty = GBZResult(E_ref=4.0 + 0j)
    return [ok, ok2, failed, empty]


# ---------------------------------------------------------------------------
# flatten / grid_steps
# ---------------------------------------------------------------------------

class TestFlatten:
    def test_counts_and_provenance(self):
        bp = flatten_results(_fake_results())
        # 1 point + 10 line + (5-1 dropped) line = 15
        assert len(bp.E) == 15
        assert bp.n_dropped == 1
        assert not bp.is_line[0]
        assert bp.is_line[1:].all()
        # failed/empty results contribute no slices
        assert set(bp.slice_idx) == {0, 1}
        assert bp.subset_idx[0] == 0 and bp.subset_idx[1] == 1

    def test_coordinates(self):
        bp = flatten_results(_fake_results())
        assert bp.mu1[0] == pytest.approx(0.1)
        assert bp.theta1[0] == pytest.approx(0.3)
        assert bp.mu2[0] == pytest.approx(0.2)
        assert bp.theta2[0] == pytest.approx(0.4)
        assert np.all((bp.theta1 >= 0) & (bp.theta1 < TWO_PI))

    def test_grid_steps_1d_fallback(self):
        E = np.linspace(0, 1, 11) + 0j  # constant Im
        d_re, d_im = grid_steps(E)
        assert d_im == d_re == pytest.approx(0.1)


# ---------------------------------------------------------------------------
# Geometry of the embedding
# ---------------------------------------------------------------------------

class TestClustering:
    def test_mu_separated_sheets_split(self):
        """Default embedding: every cluster lives on one mu side."""
        cl = cluster_bands(points=_mu_separated_cloud())
        assert cl.n_clusters == 6  # 2 sheets x 3 non-braiding branches
        assert cl.eps_window is not None
        for c in range(cl.n_clusters):
            mu1 = cl.points.mu1[cl.labels == c]
            assert mu1.max() < -0.3 or mu1.min() > 0.3, (
                f"cluster {c} spans both mu sides: [{mu1.min()}, {mu1.max()}]"
            )

    def test_mu_block_control_merges(self):
        """w_mu = 0 strips the only separating block: sheets must merge."""
        bp = _mu_separated_cloud()
        X = embed(bp, w_mu=0.0)
        labels = radius_graph_labels(X, 0.5)
        assert labels.max() + 1 == 3  # branches only; sheets glued
        for c in range(3):
            mu1 = bp.mu1[labels == c]
            assert mu1.min() < -0.3 < 0.3 < mu1.max()

    def test_vertical_arcs_two_bands(self):
        """Near-vertical arcs chain along their bodies into exactly the
        two bands (the decimation-bug regression shape)."""
        cl = cluster_bands(points=_vertical_arc_cloud())
        assert cl.n_clusters == 2
        for c in range(2):
            th1 = cl.points.theta1[cl.labels == c]
            assert th1.max() - th1.min() < 0.5  # each band keeps its column

    def test_seam_crossing_single_cluster(self):
        cl = cluster_bands(points=_seam_crossing_cloud())
        assert cl.n_clusters == 1

    def test_scan_early_stops_at_full_merge(self):
        bp = _seam_crossing_cloud()
        X = embed(bp)
        rows = eps_stability_scan(X, np.geomspace(0.02, 1.5, 21))
        assert rows[-1]["n_clusters"] == 1
        assert len(rows) < 21  # ladder stopped early

    def test_explicit_eps_skips_scan(self):
        cl = cluster_bands(points=_mu_separated_cloud(), eps=0.6)
        assert cl.scan == []
        assert cl.n_clusters == 6

    def test_requires_exactly_one_input(self):
        with pytest.raises(ValueError):
            cluster_bands()
        with pytest.raises(ValueError):
            cluster_bands([], points=_seam_crossing_cloud())

    def test_empty_cloud_rejected(self):
        with pytest.raises(ValueError):
            cluster_bands(points=flatten_results([]))
