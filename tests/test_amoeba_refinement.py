"""Analytic root tracks exercise persistent refinement independently of ZM.run."""

import numpy as np
import pytest

from pygbz2d.core import CharPoly, TWO_PI
from pygbz2d.continuation.arclength import compute_tangent
from pygbz2d.continuation.zero_manager import SegmentData
from pygbz2d.amoeba.ronkin_winding import _find_exact_crossing
from pygbz2d.amoeba.zm_extract import (
    AmoebaZeroManager, _calculate_a2_winding_and_zeros, find_crossings,
)


def affine_tracks(*tracks):
    """Make f = product(β₂ - a_j - b_j β₁) on an intentionally coarse mesh."""
    terms = {(0, 0): 1.0}
    for a, b in tracks:
        new = {}
        for (p, q), c in terms.items():
            for dp, dq, factor in ((0, 1, 1), (0, 0, -a), (1, 0, -b)):
                key = p + dp, q + dq
                new[key] = new.get(key, 0) + c * factor
        terms = new
    poly = CharPoly(np.array(list(terms.values())),
                    np.array([(0, p, q) for p, q in terms]))
    zm = AmoebaZeroManager(poly, 0j, 0.)
    theta = np.array([0., np.pi, TWO_PI])
    roots = np.array([[a + b * np.exp(1j * t) for a, b in tracks] for t in theta])
    tangent = np.array([compute_tangent(poly, 0j, np.exp(1j*t), r)[0]
                        for t, r in zip(theta, roots)])
    zm.segments = [SegmentData(theta, roots, np.argsort(abs(roots), axis=1),
                               -1, -1, tangent)]
    zm.refresh_logabs()
    return zm


def expected_angles(a, b, mu2):
    cosine = (np.exp(2*mu2) - abs(a)**2 - abs(b)**2) / (2*abs(a*b))
    alpha = np.arccos(cosine)
    phase = np.angle(b / a)
    return sorted([(-phase-alpha) % TWO_PI, (-phase+alpha) % TWO_PI])


def assert_valid_mesh(zm):
    for s, seg in enumerate(zm.segments):
        assert np.all(np.diff(seg.theta1_arr) > 0)
        np.testing.assert_array_equal(zm.seg_logabs[s], np.log(abs(seg.tracked_roots)))
        assert seg.tracked_roots.shape == seg.tangents.shape == seg.abs_argsort.shape
        for t, roots in zip(seg.theta1_arr, seg.tracked_roots):
            for root in roots:
                assert abs(zm.poly.eval_val((0j, np.exp(1j*t), root))) < 1e-10


def test_refinement_retains_actual_roots_and_reuses_same_level(monkeypatch):
    zm = affine_tracks((1.2, .6))
    original = zm._solve
    calls = []

    def counted(theta):
        calls.append(theta)
        return original(theta)

    monkeypatch.setattr(zm, '_solve', counted)
    first = find_crossings(zm, 0., 0., return_refined=True)
    assert len(first) == 2
    assert calls
    assert len(calls) == len(set(calls))
    assert len(zm.segments[0].theta1_arr) == 3 + len(calls)
    angles = sorted(np.angle(b1) % TWO_PI for b1, _ in first)
    np.testing.assert_allclose(angles, expected_angles(1.2, .6, 0.), atol=2e-12, rtol=0)
    for b1, b2 in first:
        assert abs(zm.poly.eval_val((0j, b1, b2))) < 1e-10
    calls.clear()
    second = find_crossings(zm, 0., 0., return_refined=True)
    assert calls == []
    order = lambda crossing: np.angle(crossing[0]) % TWO_PI
    np.testing.assert_allclose(sorted(second, key=order), sorted(first, key=order), atol=2e-12)
    assert_valid_mesh(zm)


def test_next_mu2_retains_previous_rows_and_correct_midpoint_winding():
    zm = affine_tracks((1.2, .6))
    _calculate_a2_winding_and_zeros(zm, 0., 0., return_refined=True)
    old_theta = zm.segments[0].theta1_arr.copy()
    old_roots = zm.segments[0].tracked_roots.copy()
    mu2 = .03
    winding, zeros = _calculate_a2_winding_and_zeros(zm, 0., mu2, return_refined=True)
    expected = expected_angles(1.2, .6, mu2)
    assert len(zeros) == 2
    np.testing.assert_allclose(sorted(t1 for t1, _ in zeros), expected, atol=2e-12, rtol=0)
    assert winding == pytest.approx((expected[1] - expected[0]) / TWO_PI, abs=1e-12)
    indices = np.searchsorted(zm.segments[0].theta1_arr, old_theta)
    np.testing.assert_array_equal(zm.segments[0].theta1_arr[indices], old_theta)
    np.testing.assert_array_equal(zm.segments[0].tracked_roots[indices], old_roots)
    assert_valid_mesh(zm)


def test_new_samples_expose_crossings_in_previously_scanned_track():
    hidden = (1.1j, .2j * np.exp(-1j*(2-np.pi)))
    visible = (1.2, .6)
    zm = affine_tracks(hidden, visible)
    assert len(find_crossings(zm, 0., 0.)) == 2
    refined = find_crossings(zm, 0., 0., return_refined=True)
    assert len(refined) == 4
    expected = sorted(expected_angles(*hidden, 0.) + expected_angles(*visible, 0.))
    np.testing.assert_allclose(sorted(np.angle(b1) % TWO_PI for b1, _ in refined),
                               expected, atol=2e-12, rtol=0)
    assert_valid_mesh(zm)


def test_existing_seam_touch_does_not_insert_or_duplicate():
    zm = affine_tracks((1.2, .6))
    mu2 = float(zm.seg_logabs[0][0, 0])
    crossings = find_crossings(zm, 0., mu2, return_refined=True)
    assert len(crossings) == 1
    assert crossings[0][0] == 1 + 0j
    assert len(zm.segments[0].theta1_arr) == 3
    # The right seam endpoint must be read as 2π, not relocated to row zero.
    t1, _ = _find_exact_crossing(zm, mu2, 0, 0, TWO_PI, TWO_PI)
    assert t1 == TWO_PI


def test_coarse_and_avoided_paths_do_not_refine():
    zm = affine_tracks((1.2, .6))
    assert len(find_crossings(zm, 0., 0.)) == 2
    assert find_crossings(zm, 0., 0., avoided_segments=[(0, 0)], return_refined=True) == []
    assert len(zm.segments[0].theta1_arr) == 3


def test_invalid_bracket_reports_context_without_linear_fallback():
    zm = affine_tracks((1.2, .6))
    with pytest.raises(RuntimeError, match=r'segment=0, track=0, bracket='):
        _find_exact_crossing(zm, 2., 0, 0, 0., np.pi)


@pytest.mark.parametrize('wtol', [None, 1e-10])
def test_inner_bisection_returns_the_partition_used_for_convergence(monkeypatch, wtol):
    from pygbz2d import core
    from pygbz2d.amoeba import bisect
    if wtol is None:
        wtol = core.WINDING_ZERO_TOL
    zm = affine_tracks((1.2, .6))
    target = .1
    zeros = [(.3, .4), (.5, .6)]
    monkeypatch.setattr(bisect, 'calculate_a2_average_winding',
                        lambda zm, mu1, mu2, **kw: mu2-target)
    fine_calls = []

    def fine(zm, mu1, mu2, **kwargs):
        fine_calls.append(mu2)
        return mu2-target, zeros

    def forbidden(*args, **kwargs):
        pytest.fail('converged crossings must not be recomputed')

    monkeypatch.setattr(bisect, '_calculate_a2_winding_and_zeros', fine)
    monkeypatch.setattr(bisect, 'find_crossings', forbidden)
    result = bisect._bisect_mu2_discrete(
        zm.poly, 0j, 0., -1., 1., max_iter=bisect.BISECT_MAX_ITER,
        wtol=wtol, coarse_xtol=bisect.BISECT_COARSE_XTOL,
        max_range_expansions=bisect.MAX_RANGE_EXPANSIONS,
        range_expand_factor=bisect.RANGE_EXPAND_FACTOR, _zm=zm)
    assert result['zeros'] is zeros
    assert result['mu2'] == fine_calls[-1]
    assert abs(result['winding']) < wtol


@pytest.mark.parametrize('options', [{}, {'xtol': None, 'max_iter': None},
                                   {'xtol': 1e-13, 'max_iter': 200}])
def test_crossing_angular_tolerance_live_defaults_and_overrides(monkeypatch, options):
    from pygbz2d.amoeba import ronkin_winding
    original = ronkin_winding.brentq
    seen = []

    def solve(func, lo, hi, **kwargs):
        seen.append((kwargs['xtol'], kwargs['maxiter']))
        return original(func, lo, hi, **kwargs)

    monkeypatch.setattr(ronkin_winding, 'brentq', solve)
    for xtol, max_iter in ((1e-12, 500), (2e-12, 400)):
        monkeypatch.setattr(ronkin_winding, 'CROSSING_XTOL', xtol)
        monkeypatch.setattr(ronkin_winding, 'CROSSING_MAXITER', max_iter)
        zm = affine_tracks((1.2, .6))
        t1, t2 = _find_exact_crossing(zm, 0., 0, 0, 0., np.pi, **options)
        expected_xtol = options.get('xtol') or xtol
        assert seen[-1] == (expected_xtol, options.get('max_iter') or max_iter)
        assert t1 == pytest.approx(expected_angles(1.2, .6, 0.)[0], abs=2*expected_xtol, rel=0)
        assert abs(zm.poly.eval_val((0j, np.exp(1j*t1), np.exp(1j*t2)))) < 1e-10
