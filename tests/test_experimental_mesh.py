"""Geometry and solved-sample guarantees of the extracted mesh algorithms."""
from types import SimpleNamespace

import numpy as np
import pytest

from pygbz2d.experimental import torus_mesh as tm, mesh_refinement as mr


def test_periodic_delaunay_closes_seams_and_covers_one_torus():
    verts, _ = tm.regular_torus_mesh(8)
    verts = (verts + np.random.default_rng(7).uniform(-.02, .02, verts.shape)) % (2*np.pi)
    tri = tm.periodic_delaunay(verts)
    top = tm.mesh_topology(tri)
    assert (top['n_vertices'], top['n_triangles'], top['chi']) == (64, 128, 0)
    assert top['n_boundary_edges'] == 0 and top['n_components'] == 1
    assert len(np.unique(tri)) == len(verts)
    _, incidence = np.unique(np.sort(tm.edge_lengths(verts, tri)[0], axis=1), axis=0, return_counts=True)
    np.testing.assert_array_equal(incidence, 2)
    assert tm.triangle_areas(verts, tri).sum() == pytest.approx(4*np.pi**2)


def test_short_arc_geometry_and_orientation_across_seam():
    verts = np.array([[2*np.pi-.1, 0], [.1, 0], [2*np.pi-.1, .2]])
    tri, flips, label = tm.orient_triangles(verts, [[0, 2, 1]])
    assert flips == 1 and label == 'positive'
    np.testing.assert_allclose(tm.edge_lengths(verts, tri)[1], [.2, np.sqrt(.08), .2])
    np.testing.assert_allclose(tm.triangle_areas(verts, tri), [.02])
    np.testing.assert_allclose(tm.torus_midpoint(verts[0], verts[1]), [0, 0], atol=1e-14)


def test_band_selection_retains_first_sample_and_provenance():
    verts, _ = tm.regular_torus_mesh(5)
    verts = np.vstack([verts, verts[0]])
    points = SimpleNamespace(theta1=verts[:, 0], theta2=verts[:, 1],
                             E=np.arange(26, dtype=complex), mu1=np.zeros(26),
                             mu2=np.zeros(26), slice_idx=np.arange(26))
    v, tri, data, merged = tm.build_band_mesh(points)
    assert merged == 1 and len(v) == 25
    np.testing.assert_array_equal(data['E'], np.arange(25))
    np.testing.assert_array_equal(data['slice_idx'], np.arange(25))
    cl = SimpleNamespace(points=points, labels=np.r_[np.ones(25), 2])
    selected = tm.build_cluster_mesh(cl, 1)
    np.testing.assert_array_equal(selected[0], v)
    np.testing.assert_array_equal(selected[1], tri)
    assert selected[3] == 0


def test_dedup_default_is_live_and_explicit_tolerance_wins(monkeypatch):
    t1, t2 = np.array([0., 1e-4, 1.]), np.zeros(3)
    assert tm.dedup_vertices(t1, t2)[3] == 0
    monkeypatch.setattr(tm, 'DEDUP_TOL', 1e-3)
    assert tm.dedup_vertices(t1, t2)[3] == 1
    assert tm.dedup_vertices(t1, t2, tol=1e-5)[3] == 0


def test_energy_predictor_preserves_complex_vertex_values():
    verts, _ = tm.regular_torus_mesh(6)
    energies = np.cos(verts[:, 0]) + 1j*np.sin(verts[:, 1])
    interp = mr.build_E_interpolator(verts, energies)
    np.testing.assert_allclose(interp(verts), energies, atol=1e-12)
    np.testing.assert_allclose(interp([[0, 0], [2*np.pi, 0]]), [1, 1], atol=1e-12)


def test_refinement_inserts_solved_states_not_predicted_states(monkeypatch):
    verts, tri = tm.regular_torus_mesh(5)
    data = dict(E=np.cos(verts[:, 0]).astype(complex), mu1=np.zeros(25), mu2=np.zeros(25))
    original = {k: v.copy() for k, v in data.items()}
    accepted = {}
    def solved_batch(coeffs, degs, pos, predicted, aL, aR, anchor_E, anchor_idx,
                     method, match_tol, **kwargs):
        accepted['pos'] = (pos[:2] + [.003, .004]) % (2*np.pi)
        assert np.max(abs(predicted)) < 2
        return (accepted['pos'], np.array([77+1j, 78+1j]),
                np.array([[.7, .8], [.9, 1.]]), np.zeros((2, 2), int),
                dict(n_direct=2, n_boundary=0, n_rejected=len(pos)-2,
                     n_probes=0, match_p50=0., boundary_p50=float('nan')))
    monkeypatch.setattr(mr, 'solve_batch', solved_batch)
    v, _, refined = mr.refine_mesh(verts, tri, data, None, None, 'amoeba',
                                   np.zeros((25, 2), int), edge_thresh=1.3, max_iters=1)
    assert len(v) == 27
    np.testing.assert_allclose(v[-2:], accepted['pos'])
    np.testing.assert_array_equal(refined['E'][-2:], [77+1j, 78+1j])
    np.testing.assert_array_equal(refined['mu1'][-2:], [.7, .9])
    for key in original:
        np.testing.assert_array_equal(data[key], original[key])


def test_index_boundary_bisection_preserves_solved_anchor_side():
    def solver(coeffs, degs, energy):
        return SimpleNamespace(success=True, index=(2, 0) if energy.real < .3 else (0, 0))
    probes, energy, result = mr.bisect_to_boundary(0, (2, 0), 1, solver, None, None, 6)
    assert len(probes) == 6 and result.index == (2, 0)
    assert 0 <= .3-energy.real < 1/64


def test_already_resolved_mesh_needs_no_solver(monkeypatch):
    verts, tri = tm.regular_torus_mesh(5)
    def forbidden(*args, **kwargs):
        raise AssertionError('An already resolved mesh must not call the GBZ solver')
    monkeypatch.setattr(mr, 'solve_batch', forbidden)
    data = dict(E=np.ones(25), mu1=np.zeros(25), mu2=np.zeros(25))
    v, t, d = mr.refine_mesh(verts, tri, data, None, None, 'amoeba',
                            np.zeros((25, 2), int), edge_thresh=2.)
    np.testing.assert_array_equal(v, verts)
    np.testing.assert_array_equal(t, tri)
    np.testing.assert_array_equal(d['E'], data['E'])
