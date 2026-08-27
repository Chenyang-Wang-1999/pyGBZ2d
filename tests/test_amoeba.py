"""Tests for pygbz2d.amoeba module using 2D HN model analytic solution."""

import numpy as np
import pytest
from cmath import exp

import pygbz2d.amoeba as bfa
from pygbz2d.amoeba.bisect import _resolve_continuum
from pygbz2d.amoeba.zm_extract import AmoebaZeroManager
from pygbz2d.core import PointSubset, LineSubset, GBZResult, CharPoly, TWO_PI

from conftest import build_HN2D_polynomial


@pytest.fixture
def params_A():
    return {"J1": 1.0, "J2": 1.0, "gamma_1": 0.2, "gamma_2": 0.3,
            "delta_1": 0.0, "delta_2": 0.0}


@pytest.fixture
def poly_A(params_A):
    return build_HN2D_polynomial(**params_A)


class TestAmoeba:
    def test_returns_gbzresult(self, poly_A):
        coeffs, degs = poly_A
        gbz = bfa.collect_GBZ_subsets(coeffs, degs, 1.0 + 0j, 0.0)
        assert isinstance(gbz, GBZResult)

    def test_inside_spectrum(self, poly_A):
        coeffs, degs = poly_A
        gbz = bfa.collect_GBZ_subsets(coeffs, degs, 1.0 + 0j, 0.0)
        assert gbz.is_gbz
        assert not gbz.is_empty
        assert gbz.index != (0, 0)

    def test_outside_spectrum(self, poly_A):
        coeffs, degs = poly_A
        gbz = bfa.collect_GBZ_subsets(coeffs, degs, 5.0 + 0j, 0.0)
        # Amoeba check may return non-empty but is_amoeba check should fail
        # Result: empty GBZResult
        assert not gbz.is_gbz or gbz.is_empty
        # At minimum, outside E should have index (0,0)
        if not gbz.is_gbz:
            assert gbz.index == (0, 0)

    def test_beta_magnitudes_analytic(self, poly_A, params_A):
        """At E=0 the HN amoeba GBZ is a full-circle continuum: cos θ₁ ≥ -1
        holds for every θ₁, so every track sits at |β₂| = exp(γ₂) over the
        whole circle.  The result is LineSubsets only (no discrete points)."""
        coeffs, degs = poly_A
        gbz = bfa.collect_GBZ_subsets(coeffs, degs, 0.0 + 0j, 0.0)
        assert gbz.is_gbz, (
            f"E=0 full-circle continuum should be GBZ, got index={gbz.index}"
        )
        # Full-circle continuum → all subsets are LineSubsets (no PointSubsets).
        assert gbz.index[0] == 0, (
            f"E=0 should have no discrete points, got {gbz.index[0]}"
        )
        assert gbz.index[1] > 0, (
            f"E=0 should have continuum LineSubsets, got {gbz.index[1]}"
        )
        gamma_1, gamma_2 = params_A["gamma_1"], params_A["gamma_2"]
        for s in gbz.subsets:
            if isinstance(s, PointSubset):
                assert abs(s.beta1) == pytest.approx(exp(gamma_1), rel=1e-3)
                assert abs(s.beta2) == pytest.approx(exp(gamma_2), rel=1e-3)
            elif isinstance(s, LineSubset):
                assert s.mu1 == pytest.approx(gamma_1, rel=1e-3)

    def test_subsets_valid_types(self, poly_A):
        coeffs, degs = poly_A
        gbz = bfa.collect_GBZ_subsets(coeffs, degs, 1.0 + 0j, 0.0)
        for s in gbz.subsets:
            assert isinstance(s, (PointSubset, LineSubset))

    def test_index_consistent(self, poly_A):
        coeffs, degs = poly_A
        gbz = bfa.collect_GBZ_subsets(coeffs, degs, 1.0 + 0j, 0.0)
        n_0d = sum(1 for s in gbz.subsets if isinstance(s, PointSubset))
        n_1d = sum(1 for s in gbz.subsets if isinstance(s, LineSubset))
        assert gbz.index == (n_0d, n_1d)


def test_all_exports():
    expected = [
        "bisect_amoeba_ronkin_min",
        "AmoebaZeroManager",
        "extract_amoeba_subsets",
        "amoeba_windings",
        "collect_GBZ_subsets",
        "PointSubset", "LineSubset", "GBZResult", "ConnectedSubset",
    ]
    for name in expected:
        assert hasattr(bfa, name), f"Missing export: {name}"


# ---- continuum-through-MR joining tests ----
#
# A segment boundary is a multiple root, but only the roots in its
# ``cluster_indices`` are genuinely degenerate there — every other track is a
# regular root passing straight through.  A continuum LineSubset whose
# endpoint root is NOT in the cluster must not be truncated at the MR; it
# continues into the adjacent segment on the matched track.  These tests build
# a synthetic 2-segment ZeroManager topology (no polynomial solving) and verify
# the joining logic in ``_join_continuum_across_mrs`` directly.


def _build_two_segment_zm(mu2, cont_on_seg1=True, cluster_tracks=(0, 1),
                          boundary_mr=False):
    """Two segments [0, π] and [π, 2π] sharing an interior MR at θ₁=π.

    Track 2 is a continuum (|β₂| = exp(μ₂)) over both segments; tracks 0/1 form
    a degenerate cluster at the MR.  ``cont_on_seg1=False`` breaks the
    continuum on segment 1 (no matching continuation → the join must raise).
    ``cluster_tracks`` selects which tracks are listed as the MR cluster.
    ``boundary_mr=True`` also places a boundary MR at θ₁=0/2π (so the seam is a
    cluster terminator and the only joinable boundary is the interior MR).
    """
    from cmath import exp, pi
    from pygbz2d.continuation.zero_manager import SegmentData
    from pygbz2d.continuation.multiple_roots import MultipleRootInfo
    from types import SimpleNamespace

    K = 3
    mu1, E = 0.2, 0.5 + 0.1j
    th0 = np.linspace(0, pi, 11)
    tr0 = np.zeros((11, K), dtype=complex)
    for i, t in enumerate(th0):
        tr0[i, 0] = 0.5 * exp(1j * t)
        tr0[i, 1] = 0.5 * exp(-1j * t)
        tr0[i, 2] = exp(mu2 + 1j * t * 2)
    m = (tr0[-1, 0] + tr0[-1, 1]) / 2
    tr0[-1, 0] = tr0[-1, 1] = m
    seg0_left_mr = 0 if boundary_mr else -1
    seg0 = SegmentData(
        theta1_arr=th0, tracked_roots=tr0,
        abs_argsort=np.tile(np.arange(K), (11, 1)),
        left_mr=seg0_left_mr, right_mr=1 if boundary_mr else 0,
        tangents=np.zeros_like(tr0),
    )

    th1 = np.linspace(pi, TWO_PI, 11)
    tr1 = np.zeros((11, K), dtype=complex)
    mr_roots = np.array(sorted([m, m, tr0[-1, 2]], key=abs))
    tr1[0] = mr_roots
    cont_col = int(np.argmin(np.abs(mr_roots - tr0[-1, 2])))
    others = [c for c in range(K) if c != cont_col]
    for i in range(1, 11):
        if cont_on_seg1:
            tr1[i, cont_col] = exp(mu2 + 1j * th1[i] * 2)
        else:
            tr1[i, cont_col] = 0.9 * exp(1j * th1[i])  # off the continuum
        tr1[i, others[0]] = 0.5 * exp(1j * th1[i])
        tr1[i, others[1]] = 0.5 * exp(-1j * th1[i])
    seg1_right_mr = 0 if boundary_mr else -1
    seg1 = SegmentData(
        theta1_arr=th1, tracked_roots=tr1,
        abs_argsort=np.tile(np.arange(K), (11, 1)),
        left_mr=1 if boundary_mr else 0, right_mr=seg1_right_mr,
        tangents=np.zeros_like(tr1),
    )

    # Continuum mask: track whose |β₂| stays at exp(μ₂) for the whole segment.
    masks = []
    for seg in (seg0, seg1):
        la = np.log(np.abs(seg.tracked_roots))
        masks.append(np.mean(np.abs(la - mu2) < 1e-6, axis=0) > 0.9)

    cluster = [tuple(cluster_tracks)] if cluster_tracks else []
    mr_list = []
    if boundary_mr:
        # Boundary MR cluster includes track 2 so the seam is a terminator for
        # the continuum too — the only joinable boundary is then the interior MR.
        mr_list.append(MultipleRootInfo(
            theta1=0.0, cluster_indices=[(0, 1, 2)],
            roots=mr_roots, cluster_stds=(0.0,),
        ))
    mr_list.append(MultipleRootInfo(
        theta1=pi, cluster_indices=cluster,
        roots=mr_roots, cluster_stds=(0.0,),
    ))
    zm = SimpleNamespace(
        segments=[seg0, seg1],
        multiple_roots=mr_list,
        has_boundary_mr=boundary_mr,
        boundary_perm=np.arange(K),
        K=K,
    )
    return zm, masks, mu1, E


class TestContinuumThroughMR:
    def test_continuum_passing_through_mr_joins(self):
        """A continuum track that passes through the interior MR as a
        non-cluster root AND through the θ₁=0/2π seam forms a closed loop.
        Both passthroughs are non-cluster, so neither is a real endpoint;
        the only real terminator is the interior MR's cluster (tracks 0/1),
        but the continuum track (track 2) is NOT in it.  The loop is therefore
        opened at the interior MR θ₁=π: the merged array has both ends at π
        (same cluster value, since the track is regular-but-degenerate-valued
        there) and the seam lands inside, where β₂ is continuous through 0/2π.

        (A 2-segment full-circle continuum has no genuine endpoint anywhere, so
        it is a closed loop; we open it at the interior MR and keep the seam
        interior rather than fabricating endpoints at θ₁=0/2π.)"""
        from pygbz2d.amoeba.zm_extract import (
            _join_continuum_across_mrs, _LinePiece,
        )
        zm, masks, mu1, E = _build_two_segment_zm(mu2=0.3)
        pieces = []
        for s, seg in enumerate(zm.segments):
            for j in np.where(masks[s])[0]:
                pieces.append(_LinePiece(
                    E=E, mu1=mu1,
                    theta1_arr=seg.theta1_arr.copy(),
                    beta2_arr=seg.tracked_roots[:, j].copy(),
                    ml=s, mr=s,
                ))
        pieces = _join_continuum_across_mrs(zm, masks, pieces)

        assert len(pieces) == 1, f"expected 1 joined LineSubset, got {len(pieces)}"
        p = pieces[0]
        # Both ends fall on the interior MR (θ₁=π); the seam is interior.
        assert p.theta1_arr[0] == pytest.approx(np.pi)
        assert p.theta1_arr[-1] == pytest.approx(np.pi)
        # 11 + 11 rows, minus the dropped seam duplicate.
        assert len(p.theta1_arr) == 21
        # β₂ continuous across the seam: the step across the θ₁ 2π→0 wrap
        # must match a typical within-segment step (same track, same sampling),
        # not jump discontinuously.
        dtheta = np.diff(p.theta1_arr)
        wrap = int(np.argmin(dtheta))
        steps = np.abs(np.diff(p.beta2_arr))
        seam_step = steps[wrap]
        # Within-segment steps exclude the wrap; compare to their median.
        in_seg = np.delete(steps, wrap)
        assert seam_step < 2 * np.median(in_seg), (
            f"β₂ step at seam ({seam_step:.2e}) >> within-segment "
            f"median ({np.median(in_seg):.2e}) → wrong join direction"
        )

    def test_seam_join_direction(self):
        """When the interior MR is a cluster terminator for the continuum (so
        the interior boundary does NOT join) but the θ₁=0/2π seam is a
        non-cluster passthrough, the seam join alone connects the two segments.

        The shared seam point is segment 0's left end (θ=0) and the last
        segment's right end (θ=2π).  The join must align THESE — i.e. the
        merged array is ``[last_seg, segment0[1:]]`` so the seam lands inside
        the array (where β₂ is continuous through 0/2π) and both array ends
        fall on the real terminator at θ₁=π.  Concatenating in the other order
        would splice segment 0's right end (θ=π) onto the last segment's
        second point — a different physical point — and break β₂ continuity.
        """
        from pygbz2d.amoeba.zm_extract import (
            _join_continuum_across_mrs, _LinePiece,
        )
        # Track 2 in the interior-MR cluster → interior boundary does not join.
        # No boundary MR → the seam is the only joinable boundary.
        zm, masks, mu1, E = _build_two_segment_zm(
            mu2=0.3, cluster_tracks=(0, 1, 2), boundary_mr=False)
        pieces = []
        for s, seg in enumerate(zm.segments):
            for j in np.where(masks[s])[0]:
                pieces.append(_LinePiece(
                    E=E, mu1=mu1,
                    theta1_arr=seg.theta1_arr.copy(),
                    beta2_arr=seg.tracked_roots[:, j].copy(),
                    ml=s, mr=s,
                ))
        pieces = _join_continuum_across_mrs(zm, masks, pieces)

        assert len(pieces) == 1, f"seam join should yield 1 piece, got {len(pieces)}"
        p = pieces[0]
        # Both ends fall on the interior terminator θ₁=π (the real endpoint);
        # the seam (θ₁=0 ≡ 2π) is interior, where β₂ must be continuous.
        assert p.theta1_arr[0] == pytest.approx(np.pi)
        assert p.theta1_arr[-1] == pytest.approx(np.pi)
        # 11 + 11 rows, minus the dropped seam duplicate.
        assert len(p.theta1_arr) == 21
        # β₂ continuous across the seam: the step across the θ₁ 2π→0 wrap
        # must match a typical within-segment step, not jump discontinuously.
        steps = np.abs(np.diff(p.beta2_arr))
        wrap = int(np.argmin(np.diff(p.theta1_arr)))
        seam_step = steps[wrap]
        in_seg = np.delete(steps, wrap)
        assert seam_step < 2 * np.median(in_seg), (
            f"β₂ step at seam ({seam_step:.2e}) >> within-segment "
            f"median ({np.median(in_seg):.2e}) → wrong join direction"
        )

    def test_missing_continuation_raises(self):
        """If a non-cluster endpoint has no matching continuum in the adjacent
        segment, the join must raise (topology inconsistency), not silently
        truncate."""
        from pygbz2d.amoeba.zm_extract import (
            _join_continuum_across_mrs, _LinePiece,
        )
        zm, masks, mu1, E = _build_two_segment_zm(
            mu2=0.3, cont_on_seg1=False)
        pieces = []
        for s, seg in enumerate(zm.segments):
            for j in np.where(masks[s])[0]:
                pieces.append(_LinePiece(
                    E=E, mu1=mu1,
                    theta1_arr=seg.theta1_arr.copy(),
                    beta2_arr=seg.tracked_roots[:, j].copy(),
                    ml=s, mr=s,
                ))
        with pytest.raises(ValueError):
            _join_continuum_across_mrs(zm, masks, pieces)

    def test_cluster_terminator_does_not_join(self):
        """When the continuum track IS part of the MR cluster, it genuinely
        terminates at the MR and is not joined through it.

        With a boundary MR at θ₁=0/2π the seam is a cluster terminator, so the
        only joinable boundary is the interior MR at θ₁=π.  Track 2 outside
        the cluster → joins (1 piece); track 2 inside the cluster → does not
        (2 pieces)."""
        from pygbz2d.amoeba.zm_extract import (
            _join_continuum_across_mrs, _LinePiece,
        )

        def run(cluster_tracks):
            zm, masks, mu1, E = _build_two_segment_zm(
                mu2=0.3, cluster_tracks=cluster_tracks, boundary_mr=True)
            pieces = []
            for s, seg in enumerate(zm.segments):
                for j in np.where(masks[s])[0]:
                    pieces.append(_LinePiece(
                        E=E, mu1=mu1,
                        theta1_arr=seg.theta1_arr.copy(),
                        beta2_arr=seg.tracked_roots[:, j].copy(),
                        ml=s, mr=s,
                    ))
            return _join_continuum_across_mrs(zm, masks, pieces)

        # Track 2 NOT in the interior-MR cluster → the continuum passes
        # through θ₁=π and the two pieces join.
        assert len(run(cluster_tracks=(0, 1))) == 1
        # Track 2 IS in the cluster → it terminates at θ₁=π, no join.
        assert len(run(cluster_tracks=(0, 1, 2))) == 2



# ---- _resolve_continuum tests ----


def _build_char_poly(coeffs, degs):
    """Build a CharPoly from coefficient / degree arrays."""
    return CharPoly(coeffs, degs)


class TestResolveContinuum:
    """Tests for _resolve_continuum using 2D HN model.

    HN model: J1=J2=1, gamma_1=0.2, gamma_2=0.3.
    Analytic SGBZ / amoeba GBZ at mu1=0.2, mu2=0.3.
    """

    @pytest.fixture
    def hn_char_poly(self, poly_A):
        coeffs, degs = poly_A
        return _build_char_poly(coeffs, degs)

    def test_at_known_boundary(self, hn_char_poly):
        """At (mu1=0.2, mu2=0.3) with E=0, the amoeba boundary is a
        continuum.  Both w2 and w1 limits should straddle zero."""
        E_ref = 0.0 + 0j
        mu1 = 0.2
        mu2 = 0.3

        zm = AmoebaZeroManager(hn_char_poly, E_ref, mu1)
        zm.run()
        result = _resolve_continuum(
            hn_char_poly, E_ref, mu1, mu2, zm,
            continuum_perturb=1e-4,
        )

        # w2 limits should exist and straddle zero
        assert result["w2_left"] is not None
        assert result["w2_right"] is not None
        assert result["w2_opposite"] is True, (
            f"w2_left={result['w2_left']}, w2_right={result['w2_right']}"
        )

        # w1 should have been resolved and also straddle zero
        assert result["w1_resolved"] is True
        assert result["w1_left"] is not None
        assert result["w1_right"] is not None
        assert result["w1_opposite"] is True, (
            f"w1_left={result['w1_left']}, w1_right={result['w1_right']}"
        )

        assert result["is_boundary"] is True

    def test_outside_spectrum_no_crash(self, hn_char_poly):
        """E=5 is outside the HN spectrum.  The function should not crash
        and should return a well-formed dict regardless of the physics."""
        E_ref = 5.0 + 0j
        mu1 = 0.2
        mu2 = 0.3

        zm = AmoebaZeroManager(hn_char_poly, E_ref, mu1)
        zm.run()
        result = _resolve_continuum(
            hn_char_poly, E_ref, mu1, mu2, zm,
            continuum_perturb=1e-4,
        )

        # Verify well-formed return dict
        for key in ("w2_left", "w2_right", "w1_left", "w1_right"):
            assert key in result
        assert isinstance(result["w2_opposite"], bool)
        assert isinstance(result["is_boundary"], bool)
        # w1_opposite and w1_resolved may be None when w2 not opposite
        for key in ("w1_opposite", "w1_resolved"):
            assert key in result

        # w2 limits should be valid (not None) — even without continuum,
        # amoeba_windings returns normal winding values
        assert result["w2_left"] is not None
        assert result["w2_right"] is not None

    def test_w2_not_opposite_no_w1(self, hn_char_poly):
        """When w2 limits don't straddle 0, w1 should not be resolved.

        At mu2=0 (far below the Ronkin minimum at mu2=0.3), the w2
        winding is on the same side of zero at mu2 ± ε, so the limits
        don't straddle 0 and w1 resolution is skipped."""
        E_ref = 0.0 + 0j
        mu1 = 0.2
        mu2 = 0.0

        zm = AmoebaZeroManager(hn_char_poly, E_ref, mu1)
        zm.run()
        result = _resolve_continuum(
            hn_char_poly, E_ref, mu1, mu2, zm,
            continuum_perturb=1e-4,
        )

        assert result["w2_opposite"] is False
        assert result["w1_resolved"] is False
        assert result["w1_left"] is None
        assert result["w1_right"] is None
        assert result["w1_opposite"] is None
        assert result["is_boundary"] is False


# ---- plateau edge detection tests ----
# Reference case from diagnostics/plateau-check-commit-report.md:
# the next-nearest-coupling model at E = 1.648 + 0.0294j exposes a
# zero-plateau edge where the bisection finds 4 clustered zeros with
# canceling jumps (net_zero_count = 0).  The plateau pre-check must
# correctly trigger the probe and classify the point as outside the
# amoeba spectrum.


@pytest.fixture
def nnc_char_poly():
    """Build the next-nearest-coupling model from playground/imaginary-degeneracy-splitting.py."""
    # README lists BerryPy as OPTIONAL; skip (not error) when it is absent.
    tb = pytest.importorskip("BerryPy").TightBinding

    gamma = 0.2
    u1, v1, w = 1.0, 0.8, 0.5
    u2, v2 = 0.1, 0.1
    u3, v3 = 0.05, 0.05
    inter_cell = [
        [0, 0, u1 * np.exp(gamma), (1, 0)],
        [0, 0, u1 * np.exp(-gamma), (-1, 0)],
        [0, 0, u2, (2, 0)],
        [0, 0, u2, (-2, 0)],
        [0, 0, u3, (3, 0)],
        [0, 0, u3, (-3, 0)],
        [0, 0, v1, (0, 1)],
        [0, 0, v1, (0, -1)],
        [0, 0, v2, (0, 2)],
        [0, 0, v2, (0, -2)],
        [0, 0, v3, (0, 3)],
        [0, 0, v3, (0, -3)],
        [0, 0, w, (1, 1)],
        [0, 0, w, (1, -1)],
        [0, 0, w, (-1, 1)],
        [0, 0, w, (-1, -1)],
    ]
    model = tb.TightBindingModel(2, 1, np.eye(2), [], inter_cell, [(0, 0)])
    coeffs, degs = model.get_characteristic_polynomial_data()
    return _build_char_poly(coeffs, degs), coeffs, degs


class TestPlateauEdge:
    """Verify the plateau pre-check correctly handles a known zero-plateau edge."""

    E_PLATEAU_EDGE = 1.648 + 0.0294j

    def test_plateau_edge_outside_spectrum(self, nnc_char_poly):
        """At the plateau-edge E, should return is_gbz=False with empty subsets."""
        char_poly, coeffs, degs = nnc_char_poly
        gbz = bfa.collect_GBZ_subsets(
            coeffs, degs, self.E_PLATEAU_EDGE, 0.0, plateau_check=True,
        )
        assert not gbz.is_gbz, (
            f"Plateau edge should be outside spectrum, got is_gbz={gbz.is_gbz}"
        )
        assert gbz.index == (0, 0)

    def test_plateau_pre_check_triggers_probe(self, nnc_char_poly):
        """At the plateau edge the bisection lands at a genuine zero-w2
        continuum-touch: two tracks merely touch |β₂|=exp(μ₂) from opposite
        sides without crossing, so w2≡0 uniformly and there are no discrete
        crossings.

        With ZeroManager tracking this is visible directly (zeros empty,
        w1/w2 areas tiny), so the point is classified as a non-GBZ plateau
        without needing the clustering-based probe.  This replaces the old
        fixed-grid behaviour, which produced 4 canceling spurious zeros at
        the edge; the functional outcome (edge → non-GBZ) is unchanged."""
        char_poly, coeffs, degs = nnc_char_poly
        from pygbz2d.amoeba.bisect import bisect_amoeba_ronkin_min
        from pygbz2d.amoeba.ronkin_winding import _get_average_winding_from_zeros

        res = bisect_amoeba_ronkin_min(char_poly, self.E_PLATEAU_EDGE)
        zeros = res["zeros"]

        # Condition (a): non-zero winding area must be tiny — w2 is uniformly 0.
        w1_area = res["_w1_area"]
        _, w2_area = _get_average_winding_from_zeros(
            char_poly, self.E_PLATEAU_EDGE, res["mu1"], res["mu2"],
            zeros, direction=2,
        )
        plateau_area_threshold = 1e-2
        assert w1_area < plateau_area_threshold, f"w1_area={w1_area} should be tiny at plateau edge"
        assert w2_area < plateau_area_threshold, f"w2_area={w2_area} should be tiny at plateau edge"

        # Condition (b): with ZM tracking the touch is resolved as a uniform
        # zero-w2 (no crossings), not as canceling clustered zeros.
        assert len(zeros) == 0, (
            f"Expected no discrete crossings at plateau edge, got {len(zeros)}"
        )

    def test_plateau_edge_net_zero_count(self, nnc_char_poly):
        """At the plateau edge the zero-w2 touch yields no discrete
        crossings, so the net jump count is trivially 0 — the point is a
        uniform zero-w2 plateau, not a genuine GBZ with canceling jumps."""
        char_poly, coeffs, degs = nnc_char_poly
        from pygbz2d.amoeba.bisect import bisect_amoeba_ronkin_min

        res = bisect_amoeba_ronkin_min(char_poly, self.E_PLATEAU_EDGE)
        zeros = res["zeros"]

        # No crossings → net jump 0 (uniform zero-w2 plateau).
        jump_sum = sum(z[2] for z in zeros)
        assert jump_sum == 0, (
            f"Plateau-edge zeros should have canceling jumps, got sum={jump_sum}"
        )
