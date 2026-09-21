"""Benchmark numerical GBZ subsets against the 2D Hatano-Nelson solution.

Run with ``--compare-only`` for the numerical checks without plotting.
Without this option, the same checks also display the point and line subsets.
The comparison covers every stored sample, but does not establish that the
solver has found every connected component of an equal-energy set.
Use ``--coarse-sweep`` to inspect a 20 x 20 amoeba scan, or ``--full-sweep``
to follow it with a 100 x 100 scan of all four GBZs saved in application/data.

author:        Wang Chenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-09-15
Copyright © Department of Physics, Tsinghua University. All rights reserved
"""

from typing import Literal
from pathlib import Path

RUN_IN_SRC = True
# Prefer the checkout when running this example without an editable install.
if RUN_IN_SRC:
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from pygbz2d import core
from pygbz2d import amoeba
from pygbz2d import sgbz

GBZKind = Literal["amoeba", "x-strip", "y-strip", "11-strip"]
GBZ_KINDS = ("amoeba", "x-strip", "y-strip", "11-strip")
# Acceptance thresholds for this benchmark; solver tolerances remain unchanged.
LOG_RADIUS_ATOL = 1e-5
ENERGY_RTOL = 1e-7
SPECTRUM_ATOL = 1e-8


def get_HN_charpoly(
    Jx1: complex,
    Jx2: complex,
    Jy1: complex,
    Jy2: complex,
    which: Literal["x-y", "y-x", "11-y"],
) -> tuple[np.ndarray, np.ndarray]:
    """Return coefficient and degree arrays for ``f = E - h``.

    ``Jx1``/``Jy1`` hop in the positive coordinate directions, and
    ``Jx2``/``Jy2`` in the negative directions. Each row of ``degs`` gives
    the powers of ``(E, beta1, beta2)`` for the matching coefficient.
    ``which`` selects ``(beta_x, beta_y)``, ``(beta_y, beta_x)``, or
    ``(beta_x * beta_y, beta_y)``. The first Bloch factor determines the
    major axis for the strip solver. Invalid basis names raise ValueError.
    """

    if which == "x-y" or which == "y-x":
        coeffs = np.array([1, -Jx1, -Jx2, -Jy1, -Jy2])
        degs = np.array([
            # E, betax, betay
            [1, 0, 0],
            [0, -1, 0],
            [0, 1, 0],
            [0, 0, -1],
            [0, 0, 1],
        ])

        if which == "y-x":
            degs = degs[:, [0, 2, 1]]
    elif which == "11-y":
        coeffs = np.array([1, -Jx1, -Jx2, -Jy1, -Jy2])
        degs = np.array([
            # E, beta_([11]), beta_(y)
            [1, 0, 0],
            [0, -1, 1],
            [0, 1, -1],
            [0, 0, -1],
            [0, 0, 1]
        ])
    else:
        raise ValueError(f"Unknown which: {which}")

    return (coeffs, degs)


def factorize_coefficients(
    J_alpha1: complex, J_alpha2: complex,
) -> tuple[float, float, complex]:
    """Return ``(gamma, delta, J)`` for a pair of nonzero finite hoppings.

    The convention ``delta = Arg(J_alpha1 * J_alpha2) / 2`` fixes the
    sign ambiguity between ``J`` and ``delta``. It reconstructs the pair
    as ``exp(gamma + 1j*delta)*J`` and ``exp(-gamma + 1j*delta)*conj(J)``.
    A vanishing hopping has no finite logarithmic radius and is rejected.
    """
    if not all(np.isfinite(z) and abs(z) > 0 for z in (J_alpha1, J_alpha2)):
        raise ValueError("The closed-form benchmark requires nonzero finite hoppings.")
    gamma_alpha = 0.5 * np.log(abs(J_alpha1) / abs(J_alpha2))
    delta_alpha = 0.5 * np.angle(J_alpha1 * J_alpha2)
    J_alpha = J_alpha1 / np.exp(gamma_alpha + 1j * delta_alpha)
    return (gamma_alpha, delta_alpha, J_alpha)


def check_factorization():
    """Check reconstruction with reproducible complex hopping coefficients."""
    rng = np.random.default_rng(0)
    J_alpha1, J_alpha2 = rng.standard_normal(2) + 1j * rng.standard_normal(2)

    gamma_alpha, delta_alpha, J_alpha = factorize_coefficients(J_alpha1, J_alpha2)

    J_alpha1_rec = np.exp(gamma_alpha + 1j * delta_alpha) * J_alpha
    J_alpha2_rec = np.exp(-gamma_alpha + 1j * delta_alpha) * J_alpha.conjugate()

    np.testing.assert_allclose(
        [J_alpha1_rec, J_alpha2_rec], [J_alpha1, J_alpha2], rtol=1e-13, atol=1e-13,
    )


def calculate_subsets(
    E_ref: complex, Jx1: complex, Jx2: complex, Jy1: complex, Jy2: complex,
    which: GBZKind,
) -> core.GBZResult:
    """Solve in the requested basis, propagating solver errors for diagnosis.

    Amoeba uses Cartesian coordinates. The strip cases use their major
    axis as beta1. A successful result can be empty for an energy outside
    the spectrum; callers must distinguish this from a failed solve.
    """
    bases = {"amoeba": "x-y", "x-strip": "x-y", "y-strip": "y-x", "11-strip": "11-y"}
    if which not in bases:
        raise ValueError(f"Unknown GBZ kind: {which}")
    coeffs, degs = get_HN_charpoly(Jx1, Jx2, Jy1, Jy2, bases[which])
    solver = amoeba if which == "amoeba" else sgbz
    return solver.collect_GBZ_subsets(coeffs, degs, E_ref, debug_mode=True)


def subset_samples(subset) -> tuple[np.ndarray, np.ndarray]:
    """Read all stored samples in the solver's basis, without line decimation."""
    if isinstance(subset, core.PointSubset):
        return np.array([subset.beta1]), np.array([subset.beta2])
    if isinstance(subset, core.LineSubset):
        return np.exp(subset.mu1 + 1j * subset.theta1_arr), subset.beta2_arr
    raise TypeError(f"Unknown subset type: {type(subset)}")


def closed_form_log_radii(
    theta1, Jx1: complex, Jx2: complex, Jy1: complex, Jy2: complex,
    which: GBZKind,
) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate the analytic ``(ln|beta1|, ln|beta2|)`` at major-axis phases.

    These are GBZ constraints, not a solution of the fixed-energy equation.
    The [11] expression is Eq. (S3.30) of arXiv:2506.22743v3. Its finite
    radius chart excludes zeros of the numerator or denominator; those
    singular cases raise ValueError instead of clipping the radius.
    """
    theta1 = np.asarray(theta1, dtype=float)
    gamma_x, delta_x, Jx = factorize_coefficients(Jx1, Jx2)
    gamma_y, delta_y, Jy = factorize_coefficients(Jy1, Jy2)
    if which in ("amoeba", "x-strip"):
        return np.full_like(theta1, gamma_x), np.full_like(theta1, gamma_y)
    if which == "y-strip":
        return np.full_like(theta1, gamma_y), np.full_like(theta1, gamma_x)
    if which == "11-strip":
        delta_xy = delta_x - delta_y
        numerator = Jx.conjugate() * np.exp(1j * (delta_xy + theta1)) + Jy
        denominator = Jx * np.exp(1j * (delta_xy - theta1)) + Jy.conjugate()
        if np.any(np.abs(numerator) == 0) or np.any(np.abs(denominator) == 0):
            raise ValueError("Singular [11]-strip radius at a requested phase.")
        mu2 = gamma_y + 0.5 * (np.log(np.abs(numerator)) - np.log(np.abs(denominator)))
        return np.full_like(theta1, gamma_x + gamma_y), mu2
    raise ValueError(f"Unknown GBZ kind: {which}")


def _distance_to_convex_polygon(point: complex, vertices: np.ndarray) -> float:
    """Distance to a cyclically ordered convex polygon, including collapsed ones.

    Segment projections cover both the exterior and degenerate line spectra;
    the half-plane test handles the interior without an angular sampling mesh.
    """
    ends = np.roll(vertices, -1)
    edges = ends - vertices
    lengths_squared = np.abs(edges)**2
    projection = np.divide(
        np.real((point - vertices) * edges.conjugate()), lengths_squared,
        out=np.zeros(len(vertices)), where=lengths_squared > 0,
    )
    closest = vertices + np.clip(projection, 0, 1) * edges
    distance = float(np.min(np.abs(point - closest)))
    cross = np.imag(edges.conjugate() * (point - vertices))
    signed_area = np.sum(np.imag(vertices.conjugate() * ends))
    if signed_area != 0 and (np.all(cross >= 0) or np.all(cross <= 0)):
        return 0.0
    return distance


def closed_form_spectrum_membership(
    E_ref: complex, Jx1: complex, Jx2: complex, Jy1: complex, Jy2: complex,
    which: GBZKind,
    *,
    spectrum_atol: float = SPECTRUM_ATOL,
) -> dict:
    """Determine spectral membership independently of the numerical GBZ solver.

    Put u = 2*exp(i*delta_x)*|Jx| and v = 2*exp(i*delta_y)*|Jy|.
    The Cartesian spectra are {u*s + v*t: s,t in [-1,1]}. For the
    [11] strip, E**2 belongs to conv{0, (u+v)**2, (u-v)**2}.
    Both tests include the boundary and handle spectra collapsed to lines.

    ``distance`` is measured after scaling E by L = |u| + |v|; the [11]
    test is in the squared-energy plane, so its scale is L**2. Distances
    no larger than ``spectrum_atol`` are treated as inside. This tolerance
    is for the analytic comparison and never changes the numerical solver.
    """
    if not np.isfinite(E_ref):
        raise ValueError("The reference energy must be finite.")
    if not np.isfinite(spectrum_atol) or spectrum_atol < 0:
        raise ValueError("The spectrum tolerance must be finite and nonnegative.")
    _, delta_x, Jx = factorize_coefficients(Jx1, Jx2)
    _, delta_y, Jy = factorize_coefficients(Jy1, Jy2)
    u = 2 * np.exp(1j * delta_x) * abs(Jx)
    v = 2 * np.exp(1j * delta_y) * abs(Jy)
    scale = abs(u) + abs(v)
    u, v, energy = u / scale, v / scale, E_ref / scale
    if which in ("amoeba", "x-strip", "y-strip"):
        vertices = np.array([u + v, -u + v, -u - v, u - v])
        target, power = energy, 1
    elif which == "11-strip":
        vertices = np.array([0, (u + v)**2, (u - v)**2])
        target, power = energy**2, 2
    else:
        raise ValueError(f"Unknown GBZ kind: {which}")
    distance = _distance_to_convex_polygon(target, vertices)
    if not np.isfinite(distance):
        raise ValueError("The analytic spectrum distance is nonfinite.")
    return {"is_gbz": distance <= spectrum_atol, "distance": distance, "energy_power": power}


def compare_to_closed_form(
    result: core.GBZResult,
    Jx1: complex, Jx2: complex, Jy1: complex, Jy2: complex,
    which: GBZKind,
    *,
    log_radius_atol: float = LOG_RADIUS_ATOL,
    energy_rtol: float = ENERGY_RTOL,
    spectrum_atol: float = SPECTRUM_ATOL,
) -> dict:
    """Check spectral membership, then every returned sample against the GBZ.

    The two log-radius errors test GBZ membership. The energy residual is
    evaluated independently from the original Cartesian hopping terms,
    normalized by ``|E_ref| + sum(|hopping term|)`` at each sample.
    A successful empty result passes only when the analytic spectrum test
    also places E_ref outside. Its sample errors are None, not zero.
    A disagreement in membership, or inconsistent subset counts, raises
    AssertionError. Failed solves and nonfinite/zero factors cannot pass.
    An error above either acceptance threshold raises AssertionError with
    the measured errors. Passing does not establish branch completeness.
    """
    if not result.success:
        raise ValueError(f"Cannot compare a failed solve: {result.error}")
    if not all(np.isfinite(t) and t >= 0 for t in (log_radius_atol, energy_rtol)):
        raise ValueError("Comparison tolerances must be finite and nonnegative.")
    spectrum = closed_form_spectrum_membership(
        result.E_ref, Jx1, Jx2, Jy1, Jy2, which, spectrum_atol=spectrum_atol,
    )
    counts = (
        sum(isinstance(s, core.PointSubset) for s in result.subsets),
        sum(isinstance(s, core.LineSubset) for s in result.subsets),
    )
    if result.index != counts or sum(counts) != len(result.subsets):
        raise AssertionError("GBZResult.index is inconsistent with its subset objects.")
    report = {
        "which": which, "E_ref": result.E_ref, "index": result.index,
        "analytic_is_gbz": spectrum["is_gbz"], "numerical_is_gbz": result.is_gbz,
        "spectrum_distance": spectrum["distance"],
        "spectrum_energy_power": spectrum["energy_power"],
        "n_samples": 0, "max_mu1_error": None,
        "max_mu2_error": None, "max_energy_residual": None,
    }
    if report["analytic_is_gbz"] != report["numerical_is_gbz"]:
        raise AssertionError(f"Spectrum membership mismatch: {report}")
    if result.is_empty:
        return report

    errors = []
    for subset in result.subsets:
        beta1, beta2 = subset_samples(subset)
        if beta1.size == 0 or not all(
            np.all(np.isfinite(beta) & (np.abs(beta) > 0)) for beta in (beta1, beta2)
        ):
            raise ValueError("The comparison requires nonempty, finite, nonzero samples.")
        mu1_exact, mu2_exact = closed_form_log_radii(
            np.angle(beta1), Jx1, Jx2, Jy1, Jy2, which,
        )
        if which == "y-strip":
            beta_x, beta_y = beta2, beta1
        elif which == "11-strip":
            beta_x, beta_y = beta1 / beta2, beta2
        else:
            beta_x, beta_y = beta1, beta2
        terms = np.array([Jx1 / beta_x, Jx2 * beta_x, Jy1 / beta_y, Jy2 * beta_y])
        energy_scale = abs(result.E_ref) + np.sum(np.abs(terms), axis=0)
        errors.append(np.column_stack((
            np.abs(np.log(np.abs(beta1)) - mu1_exact),
            np.abs(np.log(np.abs(beta2)) - mu2_exact),
            np.abs(result.E_ref - np.sum(terms, axis=0)) / energy_scale,
        )))

    errors = np.concatenate(errors)
    if not np.all(np.isfinite(errors)):
        raise ValueError("The closed-form comparison produced nonfinite errors.")
    maxima = np.max(errors, axis=0)
    report.update({
        "n_samples": len(errors), "max_mu1_error": float(maxima[0]),
        "max_mu2_error": float(maxima[1]), "max_energy_residual": float(maxima[2]),
    })
    if maxima[0] > log_radius_atol or maxima[1] > log_radius_atol or maxima[2] > energy_rtol:
        raise AssertionError(f"Closed-form comparison failed: {report}")
    return report


def print_comparison(report: dict) -> None:
    """Print membership and errors; a verified empty result has no sample errors."""
    membership = "inside" if report["analytic_is_gbz"] else "outside"
    if report["n_samples"] == 0:
        print(
            f"{report['which']:8s}  index={report['index']}  N=    0  "
            f"spectrum={membership} (agrees)  distance={report['spectrum_distance']:.3e}  "
            "dmu1=dmu2=rE=N/A", flush=True,
        )
        return
    print(
        f"{report['which']:8s}  index={report['index']}  N={report['n_samples']:5d}  "
        f"spectrum={membership} (agrees)  "
        f"dmu1={report['max_mu1_error']:.3e}  dmu2={report['max_mu2_error']:.3e}  "
        f"rE={report['max_energy_residual']:.3e}",
        flush=True,
    )


def calculate_subset_and_show(
    E_ref: complex,
    Jx1: complex,
    Jx2: complex,
    Jy1: complex,
    Jy2: complex,
    which: GBZKind,
    *,
    show: bool = True,
) -> core.GBZResult:
    """Solve, compare with the analytic GBZ, and optionally display two plots.

    Returns the solver's GBZResult in its native coordinate basis. ``show``
    can be disabled for automated runs without importing Matplotlib.
    Empty results are checked against the analytic spectrum without plotting.
    Nonempty results must pass the radius and energy checks before plotting.
    """
    res = calculate_subsets(E_ref, Jx1, Jx2, Jy1, Jy2, which)
    if not res.success:
        raise RuntimeError(res.error)
    if res.is_empty or not show:
        return res

    import matplotlib.pyplot as plt

    fig, (ax_phase, ax_radius) = plt.subplots(1, 2, figsize=(10, 4), layout="constrained")
    mu2_samples = []
    for subset in res.subsets:
        beta1, beta2 = subset_samples(subset)
        theta1 = (np.angle(beta1) / np.pi) % 2
        mu2_samples.append(np.log(np.abs(beta2)))
        # Dots avoid drawing spurious connections across the angular seam.
        ax_phase.plot(theta1, (np.angle(beta2) / np.pi) % 2, ".")
        ax_radius.plot(theta1, mu2_samples[-1], ".")
    theta1 = np.linspace(0, core.TWO_PI, 501)
    _, mu2_exact = closed_form_log_radii(theta1, Jx1, Jx2, Jy1, Jy2, which)
    ax_radius.plot(theta1 / np.pi, mu2_exact, "k--", label="Closed-form radius")
    ax_radius.legend()
    # Keep a constant-radius GBZ from being magnified to the roundoff scale.
    mu2_all = np.concatenate([mu2_exact, *mu2_samples])
    ax_radius.set_ylim(np.min(mu2_all) - 0.1, np.max(mu2_all) + 0.1)
    ax_phase.set(xlabel=r"$\theta_1/\pi$", ylabel=r"$\theta_2/\pi$", xlim=(0, 2), ylim=(0, 2))
    ax_radius.set(xlabel=r"$\theta_1/\pi$", ylabel=r"$\mu_2 = \ln|\beta_2|$", xlim=(0, 2))
    fig.suptitle(f"{which} GBZ subsets, E = {E_ref}")
    plt.show()

    print_comparison(compare_to_closed_form(res, Jx1, Jx2, Jy1, Jy2, which))
    return res


def demo_point_subsets(*, show: bool = True):
    """Compare isolated subsets at E = 1+i using the paper's hopping parameters."""
    Jx1 = 1 + 1j
    Jx2 = 1.5 + 1.2j
    Jy1 = -1 + 1j
    Jy2 = -1.2 - 0.5j
    E_ref = 1 + 1j
    print(f"Point-subset benchmark: E = {E_ref}", flush=True)
    for which in GBZ_KINDS:
        result = calculate_subset_and_show(E_ref, Jx1, Jx2, Jy1, Jy2, which, show=show)
        if result.is_empty or any(not isinstance(s, core.PointSubset) for s in result.subsets):
            raise AssertionError(f"Expected nonempty point subsets for {which} at E = {E_ref}.")


def demo_line_subsets(*, show: bool = True):
    """Compare continuum samples at E = 1 for a model with a real spectrum."""
    Jx1 = 1
    Jx2 = 1.5
    Jy1 = -1
    Jy2 = -1.2
    E_ref = 1
    print(f"Line-subset benchmark: E = {E_ref}", flush=True)
    for which in GBZ_KINDS:
        result = calculate_subset_and_show(E_ref, Jx1, Jx2, Jy1, Jy2, which, show=show)
        if result.is_empty or any(not isinstance(s, core.LineSubset) for s in result.subsets):
            raise AssertionError(f"Expected nonempty line subsets for {which} at E = {E_ref}.")


def demo_spectrum_membership():
    """Check exterior energies and one energy whose membership depends on geometry."""
    cases = (
        ("complex hoppings, far outside", 6 + 6j,
         (1 + 1j, 1.5 + 1.2j, -1 + 1j, -1.2 - 0.5j), ()),
        ("real spectrum, beyond endpoints", 6, (1, 1.5, -1, -1.2), ()),
        ("real spectrum, transverse", 1j, (1, 1.5, -1, -1.2), ()),
        ("geometry-dependent membership", 2 + 2j,
         (1 + 1j, 1.5 + 1.2j, -1 + 1j, -1.2 - 0.5j), ("amoeba", "x-strip", "y-strip")),
    )
    for label, E_ref, hoppings, expected_inside in cases:
        print(f"Spectrum-membership benchmark ({label}): E = {E_ref}", flush=True)
        for which in GBZ_KINDS:
            if closed_form_spectrum_membership(E_ref, *hoppings, which)["is_gbz"] != (which in expected_inside):
                raise AssertionError(f"Unexpected analytic membership for {which} at E = {E_ref}.")
            calculate_subset_and_show(E_ref, *hoppings, which, show=False)


def benchmark_random_coeffs(
    *,
    show: bool = True,
    zero_Delta: bool = False
):
    """Compare random hoppings and a Cartesian spectral energy for all GBZs.

    ``zero_Delta`` sets delta_x = delta_y. In the general case the sampled
    energy need not belong to the [11]-strip spectrum; the analytic
    membership check decides whether its numerical result should be empty.
    """
    Jx, Jy = np.random.randn(2) + 1j * np.random.randn(2)
    deltay = np.random.rand() * 2 * np.pi
    gammax, gammay = np.random.randn(2)
    if zero_Delta:
        Deltaxy = 0
    else:
        Deltaxy = np.random.rand() * 2 * np.pi
    deltax = deltay + Deltaxy
    Jx1 = Jx * np.exp(gammax + 1j * deltax)
    Jy1 = Jy * np.exp(gammay + 1j * deltay)
    Jx2 = Jx.conjugate() * np.exp(-gammax + 1j * deltax)
    Jy2 = Jy.conjugate() * np.exp(-gammay + 1j * deltay)

    ##### random eigenvalues #####
    samp_thetax, samp_thetay = np.random.rand(2) * 2 * np.pi
    samp_betax = np.exp(gammax + 1j * samp_thetax)
    samp_betay = np.exp(gammay + 1j * samp_thetay)
    E_samp = Jx1 / samp_betax + Jx2 * samp_betax + Jy1 / samp_betay + Jy2 * samp_betay

    print(f"E_samp = {E_samp}")
    print(f"Jx1 = {Jx1}")
    print(f"Jx2 = {Jx2}")
    print(f"Jy1 = {Jy1}")
    print(f"Jy2 = {Jy2}")

    ##### Solve and check #####
    for which in GBZ_KINDS:
        calculate_subset_and_show(E_samp, Jx1, Jx2, Jy1, Jy2, which, show=show)


def compare_sweep_to_closed_form(scan: dict) -> dict[str, list[dict]]:
    """Validate a completed scan without rerunning solvers or changing its data.

    This HN-specific benchmark is separate from the reusable scan below.
    Fine-scan files have already been saved when the CLI calls this function,
    so an analytic mismatch leaves the numerical evidence available.
    """
    reports = {}
    for which, results in scan["results"].items():
        reports[which] = []
        for index, result in enumerate(results):
            try:
                report = compare_to_closed_form(result, which=which, **scan["hoppings"])
            except (ValueError, AssertionError) as exc:
                raise AssertionError(
                    f"Post-scan comparison failed for {which} at E[{index}] = {result.E_ref}: {exc}"
                ) from exc
            reports[which].append(report)
        print(f"{which}: closed-form comparison passed for {len(results)} scanned energies", flush=True)
    return reports


#### Full GBZ sweep ####
# This section can be copied with get_HN_charpoly into a separate script.
# Repeat its imports here so none of the earlier benchmark helpers is needed.
from pathlib import Path
from typing import Literal

import numpy as np
from pygbz2d import GBZResult, amoeba, sgbz

DEMO_HOPPINGS = (1 + 1j, 1.5 + 1.2j, -1 + 1j, -1.2 - 0.5j)
N_PROCESS = 1
SWEEP_BASES = {"amoeba": "x-y", "x-strip": "x-y", "y-strip": "y-x", "11-strip": "11-y"}

def _sweep_worker(task):
    """Call the library directly; keep workers importable under Windows spawn."""
    index, E_ref, coeffs, degs, which = task
    try:
        solver = amoeba.collect_GBZ_subsets if which == "amoeba" else sgbz.collect_GBZ_subsets
        result = solver(coeffs, degs, E_ref)
    except Exception as exc:
        raise RuntimeError(f"{which} sweep failed at E[{index}] = {E_ref}: {exc}") from exc
    return index, result


def sweep_GBZ(
    Jx1: complex,
    Jx2: complex,
    Jy1: complex,
    Jy2: complex,
    E_list: np.ndarray,
    which: str,
    n_process: int = N_PROCESS,
) -> list[GBZResult]:
    """Solve every energy with the library, preserving the input order.

    Argument order is (Jx1, Jx2, Jy1, Jy2), as in get_HN_charpoly. The
    polynomial is built once in the selected basis and passed to workers.
    The parent reports progress; worker exceptions immediately abort the pool.
    No closed-form solution is needed for scanning. With multiple
    processes, each worker uses one BLAS thread to avoid oversubscribing
    the machine. Call from an ``if __name__ == '__main__'`` guard on Windows.
    """
    import os
    import time
    import multiprocessing as mp

    if which not in SWEEP_BASES:
        raise ValueError(f"Unknown GBZ kind: {which}")
    if not isinstance(n_process, int) or isinstance(n_process, bool) or n_process < 1:
        raise ValueError("n_process must be a positive integer.")
    energies = np.asarray(E_list, dtype=complex)
    if energies.ndim != 1 or not np.all(np.isfinite(energies)):
        raise ValueError("E_list must be a one-dimensional array of finite energies.")
    if not len(energies):
        return []
    coeffs, degs = get_HN_charpoly(Jx1, Jx2, Jy1, Jy2, SWEEP_BASES[which])
    tasks = ((i, complex(E), coeffs, degs, which) for i, E in enumerate(energies))
    results = [None] * len(energies)
    started = time.monotonic()
    progress_step = max(1, len(energies) // 20)

    def collect(iterator):
        for completed, (index, result) in enumerate(iterator, 1):
            results[index] = result
            if completed == 1 or completed % progress_step == 0 or completed == len(energies):
                elapsed = time.monotonic() - started
                print(f"{which}: {completed}/{len(energies)} energies, {elapsed:.1f} s", flush=True)

    if n_process > 1:
        # Set these before spawn imports NumPy; an initializer would be too late.
        thread_vars = ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")
        previous = {name: os.environ.get(name) for name in thread_vars}
        try:
            for name in thread_vars:
                os.environ[name] = "1"
            with mp.get_context("spawn").Pool(min(n_process, len(energies))) as pool:
                collect(pool.imap_unordered(_sweep_worker, tasks, chunksize=1))
        finally:
            for name, value in previous.items():
                if value is None:
                    os.environ.pop(name, None)
                else:
                    os.environ[name] = value
    else:
        collect(map(_sweep_worker, tasks))
    return results


def _energy_grid(real_bounds, imag_bounds, grid_size):
    """Use rows for Im(E), columns for Re(E), and C-order flattening throughout."""
    if not isinstance(grid_size, int) or isinstance(grid_size, bool) or grid_size < 2:
        raise ValueError("grid_size must be an integer of at least two.")
    for bounds in (real_bounds, imag_bounds):
        if len(bounds) != 2 or not np.all(np.isfinite(bounds)) or bounds[0] >= bounds[1]:
            raise ValueError("Each energy interval must have finite, increasing endpoints.")
    real_axis = np.linspace(*real_bounds, grid_size)
    imag_axis = np.linspace(*imag_bounds, grid_size)
    energies = (real_axis[None, :] + 1j * imag_axis[:, None]).ravel()
    return real_axis, imag_axis, energies


def fine_bounds_from_coarse(real_axis, imag_axis, results, *, flatten_order="C"):
    """Expand the occupied coarse bounding box by one original spacing per side.

    An empty scan or an occupied outer row/column cannot establish a scan
    range. Such cases are reported rather than silently substituting an
    analytic bound or interpreting failed solves as exterior energies.
    """
    real_axis, imag_axis = np.asarray(real_axis), np.asarray(imag_axis)
    for axis in (real_axis, imag_axis):
        if axis.ndim != 1 or len(axis) < 2 or not np.all(np.isfinite(axis)):
            raise ValueError("Coarse axes must be finite one-dimensional grids.")
        if np.any(np.diff(axis) <= 0) or not np.allclose(np.diff(axis), axis[1] - axis[0]):
            raise ValueError("Coarse axes must be uniformly spaced and increasing.")
    if len(results) != len(real_axis) * len(imag_axis):
        raise ValueError("Coarse results do not match the grid shape.")
    if any(not result.success for result in results):
        raise RuntimeError("The coarse scan contains failed solves; its bounds are not usable.")
    occupied = np.array([result.is_gbz for result in results]).reshape(
        (len(imag_axis), len(real_axis)), order=flatten_order,
    )
    rows, columns = np.nonzero(occupied)
    if not len(rows):
        raise RuntimeError("No amoebic spectral points were found; inspect the coarse grid before refining.")
    if occupied[0, :].any() or occupied[-1, :].any() or occupied[:, 0].any() or occupied[:, -1].any():
        raise RuntimeError("The occupied coarse grid touches the scan boundary; enlarge the coarse region first.")
    dx, dy = real_axis[1] - real_axis[0], imag_axis[1] - imag_axis[0]
    return (
        (float(real_axis[columns.min()] - dx), float(real_axis[columns.max()] + dx)),
        (float(imag_axis[rows.min()] - dy), float(imag_axis[rows.max()] + dy)),
    )


def load_sweep(path) -> dict:
    """Load a scan saved in the same structure returned by the sweep functions."""
    import pickle

    with Path(path).open("rb") as stream:
        return pickle.load(stream)


def plot_sweep_spectra(scan: dict, *, show: bool = True):
    """Plot numerical membership; a coarse plot also outlines the fine region.

    Accepts either a returned scan or a loaded file without extra metadata.
    Energies and coarse bounds are derived from the axes and results. Figures
    are never saved by this function; the returned figure can be customized.
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    kinds = tuple(scan["results"])
    fig, axes = plt.subplots(1 if len(kinds) == 1 else 2, 1 if len(kinds) == 1 else 2,
                             figsize=(6, 5) if len(kinds) == 1 else (10, 8),
                             squeeze=False, layout="constrained")
    energies = (np.asarray(scan["real_axis"])[None, :]
                + 1j * np.asarray(scan["imag_axis"])[:, None]).ravel(order=scan["flatten_order"])
    for ax, which in zip(axes.flat, kinds):
        results = scan["results"][which]
        inside = np.array([result.is_gbz for result in results])
        success = np.array([result.success for result in results])
        outside = success & ~inside
        ax.scatter(energies[outside].real, energies[outside].imag, s=3, c="0.85", label="Outside")
        ax.scatter(energies[inside].real, energies[inside].imag, s=10, label="Inside")
        if not success.all():
            ax.scatter(energies[~success].real, energies[~success].imag,
                       s=20, c="tab:red", marker="x", label="Failed")
        if scan["stage"] == "coarse":
            xr, yr = fine_bounds_from_coarse(
                scan["real_axis"], scan["imag_axis"], scan["results"]["amoeba"],
                flatten_order=scan["flatten_order"],
            )
            ax.add_patch(Rectangle((xr[0], yr[0]), xr[1] - xr[0], yr[1] - yr[0],
                                   fill=False, color="tab:red", linestyle="--", label="Fine-scan region"))
            ax.autoscale_view()
        ax.set(xlabel="Re(E)", ylabel="Im(E)", title=f"{which}: {inside.sum()} spectral energies")
        ax.set_aspect("equal", adjustable="box")
        ax.legend(loc="best")
    fig.suptitle(f"{scan['stage'].capitalize()} GBZ scan ({len(scan['real_axis'])} x {len(scan['imag_axis'])})")
    if show:
        plt.show()
    return fig


def coarse_sweep(
    Jx1: complex = DEMO_HOPPINGS[0], Jx2: complex = DEMO_HOPPINGS[1],
    Jy1: complex = DEMO_HOPPINGS[2], Jy2: complex = DEMO_HOPPINGS[3],
    *, n_process: int = N_PROCESS, grid_size: int = 20,
    real_bounds=(-5.0, 5.0), imag_bounds=(-5.0, 5.0), show: bool = True,
) -> dict:
    """Scan only amoeba, show its spectrum, and return the unsaved coarse data.

    The 20 points on each axis include both endpoints, so the default
    spacing is 10/19. Fine bounds extend the occupied bounding box by that
    spacing along each axis. ``show=False`` is available for batch runs.
    """
    real_axis, imag_axis, energies = _energy_grid(real_bounds, imag_bounds, grid_size)
    print(f"Coarse scan: amoeba only, {grid_size} x {grid_size}, {n_process} processes", flush=True)
    hoppings = {"Jx1": Jx1, "Jx2": Jx2, "Jy1": Jy1, "Jy2": Jy2}
    results = sweep_GBZ(**hoppings, E_list=energies, which="amoeba", n_process=n_process)
    fine_bounds = fine_bounds_from_coarse(real_axis, imag_axis, results)
    scan = {
        "stage": "coarse", "hoppings": hoppings,
        "real_axis": real_axis, "imag_axis": imag_axis,
        "flatten_order": "C", "results": {"amoeba": results},
    }
    print(f"Fine-scan bounds: Re(E) in {fine_bounds[0]}, Im(E) in {fine_bounds[1]}", flush=True)
    if show:
        plot_sweep_spectra(scan)
    return scan


def fine_sweep(coarse: dict, *, n_process: int = N_PROCESS, grid_size: int = 100,
               output_dir=None, show: bool = True) -> dict:
    """Scan every fine-grid energy with all four methods and save full results.

    Returned scans and per-method files use the same six-field structure;
    files simply contain one entry in the results mapping. The energy grid
    is represented by its axes and flatten order. Run-specific filenames
    prevent overwriting previous scans. A solver exception stops the scan;
    previously completed files remain usable and no incomplete file is
    presented as a completed scan.
    """
    import pickle
    from datetime import datetime, timezone

    real_bounds, imag_bounds = fine_bounds_from_coarse(
        coarse["real_axis"], coarse["imag_axis"], coarse["results"]["amoeba"],
        flatten_order=coarse["flatten_order"],
    )
    real_axis, imag_axis, energies = _energy_grid(real_bounds, imag_bounds, grid_size)
    output_dir = Path(output_dir) if output_dir is not None else Path(__file__).resolve().parent / "data"
    output_dir.mkdir(parents=True, exist_ok=True)
    run_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    scan = {
        "stage": "fine", "hoppings": dict(coarse["hoppings"]),
        "real_axis": real_axis, "imag_axis": imag_axis,
        "flatten_order": "C", "results": {},
    }
    print(f"Fine scan: four GBZs, {grid_size} x {grid_size}, {n_process} processes", flush=True)
    for which in SWEEP_BASES:
        results = sweep_GBZ(**scan["hoppings"], E_list=energies, which=which, n_process=n_process)
        scan["results"][which] = results
        payload = {**scan, "results": {which: results}}
        path = output_dir / f"HN2D-{run_id}-{which}.pkl"
        temporary = path.with_suffix(".pkl.tmp")
        with temporary.open("wb") as stream:
            pickle.dump(payload, stream, protocol=pickle.HIGHEST_PROTOCOL)
        temporary.replace(path)
        print(f"Saved {len(results)} {which} results -> {path}", flush=True)
    if show:
        plot_sweep_spectra(scan)
    return scan


#### Command-line demos and optional benchmark validation ####
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mode",
        type=str,
        default="demo",
        help="Demo mode: \n\t"
            "'demo': show PointSubset and LineSubset demos. \n\t"
            "'benchmark': run random coefficient benchmark only. \n\t"
            "'coarse-sweep': show the unsaved 20 x 20 amoeba scan only. \n\t"
            "'full-sweep': coarse scan, then save all four 100 x 100 fine scans. \n\t"
            "'plot': plot the GBZ scan, need --data-fname"
        )
    parser.add_argument("--compare-only", action="store_true", help="Run checks without opening plots.")
    parser.add_argument("--n-process", type=int, default=N_PROCESS, help="Sweep worker count (default: 1).")
    parser.add_argument("--no-show", action="store_true", help="Disable figures for batch sweep runs.")
    parser.add_argument("--zero-Delta", action="store_true", help="Set delta_x = delta_y.")
    parser.add_argument("--data-fname", type=str, help="Data filename for plot mode.")
    args = parser.parse_args()
    if args.mode == "coarse-sweep" or args.mode == "full-sweep":
        if args.zero_Delta or args.compare_only:
            parser.error("Sweep modes use the complex demo hoppings; use --no-show to disable their figures.")
        coarse = coarse_sweep(n_process=args.n_process, show=not args.no_show)
        if args.mode == "full-sweep":
            fine = fine_sweep(coarse, n_process=args.n_process, show=not args.no_show)
        compare_sweep_to_closed_form(coarse)
        if args.mode == "full-sweep":
            compare_sweep_to_closed_form(fine)
    elif args.mode == "benchmark":
        benchmark_random_coeffs(show=not args.compare_only, zero_Delta=args.zero_Delta)
    elif args.mode == "plot":
        if args.data_fname is None:
            raise ValueError("--data-fname is required for plot mode.")
        scan = load_sweep(args.data_fname)
        plot_sweep_spectra(scan, show=not args.no_show)
    elif args.mode == "demo":
        demo_point_subsets(show=not args.compare_only)
        demo_line_subsets(show=not args.compare_only)
        demo_spectrum_membership()
    else:
        raise ValueError(f"Unknown mode: {args.mode}")
