"""Continuum, gauge, and orientation checks for the experimental LR integrator."""
import numpy as np
import pytest

from pygbz2d.experimental.chern import (
    biorthogonal_edge_links,
    integrate_chern,
    left_right_eigenvectors_on_mesh,
    triangle_flux_biorthogonal,
    triangle_flux_biorthogonal_complex,
    triangle_flux_right,
)
from pygbz2d.experimental.torus_mesh import regular_torus_mesh


def analytic_states(h, variable_y):
    """R=(1,x), Lbra=(1-s*x,s): constant s or s=(1+i)*y."""
    x = np.array([0., h, 0.])
    y = np.array([0., 0., h])
    s = (1 + 1j) * (y if variable_y else np.ones(3))
    right = np.column_stack([np.ones(3), x]).astype(complex)
    left = np.column_stack([1 - s*x, s]).conj()
    return right, left


def test_zero_curvature_does_not_acquire_complex_metric_phase():
    # H=2|R><L|-I is smooth and has eigenvalues +/-1. Its connection
    # a=(1+i)dx is flat, although the old three-overlap rule is not.
    h = 1e-3
    right, left = analytic_states(h, variable_y=False)
    old = -np.angle(np.vdot(left[0], right[1]) * np.vdot(left[1], right[2])
                    * np.vdot(left[2], right[0]))
    assert old / (.5*h*h) == pytest.approx(4., abs=1e-8)
    corrected = triangle_flux_biorthogonal_complex(right, left, [[0, 1, 2]])
    np.testing.assert_allclose(corrected, 0, atol=1e-15)


def test_complex_curvature_converges_to_analytic_one_minus_i():
    # a=(1+i)y dx => B=i(da)_{xy}=1-i, including its imaginary part.
    errors = []
    for h in [.08, .04, .02, .01]:
        right, left = analytic_states(h, variable_y=True)
        flux = triangle_flux_biorthogonal_complex(right, left, [[0, 1, 2]])[0]
        errors.append(abs(flux/(.5*h*h) - (1-1j)))
    assert all(fine < .26*coarse for coarse, fine in zip(errors, errors[1:]))
    assert errors[-1] < 1.01e-4


def qwz_states(mass, gamma=.4, n=17):
    verts, tri = regular_torus_mesh(n)
    rng = np.random.default_rng(781)
    verts = (verts + rng.uniform(-.12, .12, verts.shape)*(2*np.pi/n)) % (2*np.pi)
    x, y = verts.T
    dx, dy, dz = np.sin(x)+1j*gamma, np.sin(y)+1j*gamma, mass-np.cos(x)-np.cos(y)
    matrices = np.empty((len(verts), 2, 2), complex)
    matrices[:, 0, 0], matrices[:, 1, 1] = dz, -dz
    matrices[:, 0, 1], matrices[:, 1, 0] = dx-1j*dy, dx+1j*dy
    energies, all_right = np.linalg.eig(matrices)
    band = np.argmin(energies.real, axis=1)
    ids = np.arange(len(verts))
    right = all_right[ids, :, band]
    # An independent eigendecomposition/inverse supplies exact dual bras.
    left = np.linalg.inv(all_right)[ids, band, :].conj()
    return verts, tri, energies[ids, band], right, left


@pytest.mark.parametrize("mass, expected", [(1.2, -1), (2.6, 0)])
def test_irregular_qwz_gauge_orientation_and_relabeling(mass, expected):
    verts, tri, _, right, left = qwz_states(mass)
    flux = triangle_flux_biorthogonal_complex(right, left, tri)
    assert flux.sum()/(2*np.pi) == pytest.approx(expected, abs=2e-13)
    reverse = triangle_flux_biorthogonal_complex(right, left, tri[:, ::-1])
    np.testing.assert_allclose(reverse, -flux, atol=2e-14)
    rng = np.random.default_rng(91)
    gauge = np.exp(rng.uniform(-12, 12, len(verts)) + 1j*rng.uniform(-np.pi, np.pi, len(verts)))
    changed = triangle_flux_biorthogonal_complex(
        right*gauge[:, None], left/gauge.conj()[:, None], tri)
    np.testing.assert_allclose(changed, flux, atol=2e-14)
    order = rng.permutation(len(verts))
    remap = np.argsort(order)
    relabeled = triangle_flux_biorthogonal_complex(right[order], left[order], remap[tri])
    np.testing.assert_allclose(relabeled, flux, atol=2e-14)


def test_hermitian_reduction_matches_right_flux_on_each_triangle():
    _, tri, _, right, left = qwz_states(1.2, gamma=0.)
    rr, bad = triangle_flux_right(right, tri)
    lr, lr_bad = triangle_flux_biorthogonal(right, left, tri)
    complex_flux = triangle_flux_biorthogonal_complex(right, left, tri)
    assert bad == lr_bad == 0
    np.testing.assert_allclose(lr, rr, atol=2e-15)
    np.testing.assert_allclose(complex_flux.imag, 0, atol=2e-15)


def test_reciprocal_edge_and_gauge_independent_square_root():
    right, left = analytic_states(.2, variable_y=True)
    edges = np.array([[0, 1], [1, 0]])
    links = biorthogonal_edge_links(right, left, edges)
    assert links.prod() == pytest.approx(1)
    # This gauge changes the sign of principal sqrt(Uij/Uji), but should
    # multiply the transport by exp(2i), with no extra sign.
    gauge = np.array([1, np.exp(2j), 1])
    changed = biorthogonal_edge_links(right*gauge[:, None], left/gauge.conj()[:, None], edges)
    np.testing.assert_allclose(changed, links*gauge[edges[:, 1]]/gauge[edges[:, 0]], atol=1e-15)


@pytest.mark.parametrize("kind", ["zero", "branch"])
def test_undefined_edge_is_rejected(kind):
    if kind == "zero":
        right = left = np.eye(2, dtype=complex)
    else:
        right = np.array([[1, 0], [1, 1]], complex)
        left = np.array([[1, 0], [-1, 2]], complex)
    with pytest.raises(ValueError, match="overlap product|branch cut"):
        biorthogonal_edge_links(right, left, [[0, 1]])


def test_exceptional_point_cannot_be_biorthogonally_normalized():
    h = lambda b1, b2: np.array([[0., 1.], [0., 0.]])
    with pytest.raises(ValueError, match="vertex 0.*self-overlap"):
        left_right_eigenvectors_on_mesh(h, np.array([0.]), np.ones(1), np.ones(1), 1e-6)


def test_integrate_chern_returns_complex_flux_and_rule_metadata():
    # The same analytic B=1-i example now exercises eigenvector extraction
    # and the returned output, rather than supplying states to the flux rule.
    h = .01
    verts = np.array([[0, 0], [h, 0], [0, h]])
    def hamiltonian(b1, b2):
        x, y = np.angle(b1), np.angle(b2)
        right = np.array([1, x])
        left_bra = np.array([1-(1+1j)*y*x, (1+1j)*y])
        return 2*np.outer(right, left_bra)-np.eye(2)
    result, saved = integrate_chern(verts, [[0, 1, 2]], np.ones(3), np.zeros(3),
                                    np.zeros(3), hamiltonian, mode="biorthogonal")
    assert result.link_rule == "biorthogonal-reciprocal-v1"
    assert result.total_flux_imag/(.5*h*h) == pytest.approx(-1, abs=1.1e-4)
    assert saved["flux_complex"][0]/(.5*h*h) == pytest.approx(1-1j, abs=1.1e-4)
    np.testing.assert_array_equal(saved["flux"], saved["flux_complex"].real)
    assert saved["link_rule"] == result.link_rule
