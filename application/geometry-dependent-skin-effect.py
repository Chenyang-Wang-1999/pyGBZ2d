"""Gain-loss Haldane model: directional SGBZs and finite open geometries.

Run without arguments for a small end-to-end demonstration. Use --help for
independent sweep, OBC, plotting, validation and legacy-data extraction commands.
All paths default to this checkout, so the working directory is immaterial.
BerryPy and SymPy are needed for modeling; plotting additionally uses Matplotlib.
Only load pickle files from trusted research data.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import pickle
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DATA = ROOT / "application/data/geometry-dependent-skin-effect"
FIGURES = ROOT / "application/Figures/geometry-dependent-skin-effect"
PARAMS = (1.0, 0.5, np.pi / 3, 0.5j, 0.0)
DIRECTIONS = ("a1", "a2", "x", "y")
TRANSFORMS = {
    "a1": np.eye(2, dtype=int),
    "a2": np.array([[0, 1], [1, 0]]),
    "x": np.array([[1, 1], [1, -1]]),
    "y": np.array([[1, 1], [-1, 1]]),
}
COLORS = {"a1": "#305f9e", "a2": "#269c95", "x": "#d29229", "y": "#c24752"}


def build_model(params=PARAMS):
    """BerryPy convention: [source, destination, amplitude, cell shift]."""
    from BerryPy import TightBinding as tb

    t1, t2, phi, mass, gamma = params
    lattice = np.array([[-0.5, -0.5], [-np.sqrt(3) / 2, np.sqrt(3) / 2]])
    intracell = [[0, 0, mass], [1, 1, -mass], [1, 0, t1], [0, 1, t1]]
    intercell = [[1, 0, t1, (0, -1)], [1, 0, t1, (1, 0)],
                 [0, 1, t1, (0, 1)], [0, 1, t1, (-1, 0)]]
    for sublattice, sign in ((0, 1), (1, -1)):
        for shift in ((-1, 0), (0, -1), (1, 1)):
            intercell.append([sublattice, sublattice,
                              t2 * np.exp(1j * (gamma + sign * phi)), shift])
            intercell.append([sublattice, sublattice,
                              t2 * np.exp(1j * (gamma - sign * phi)),
                              tuple(-np.array(shift))])
    model = tb.TightBindingModel(2, 2, lattice, intracell, intercell)
    sites = np.array([[0, 1 / (2 * np.sqrt(3))], [0, -1 / (2 * np.sqrt(3))]])
    model.SiteCoord = model.cart2lattice(sites.T).T
    return model


def direction_model(direction, params=PARAMS):
    model = build_model(params)
    if direction == "a1":
        return model
    cells = [(0, 0)] if direction == "a2" else [(0, 0), (1, 0)]
    # BerryPy uses COLUMNS of this integer matrix as the new primitive vectors.
    return model.get_supercell(cells, TRANSFORMS[direction])


def polynomial_data(model):
    coefficients, degrees = model.get_characteristic_polynomial_data()
    return np.asarray(coefficients, complex), np.asarray(degrees, int)


def bz_spectrum(model, nk=101):
    phases = np.linspace(-np.pi, np.pi, nk, endpoint=False)
    # get_bulk_Hamiltonian takes fractional reciprocal coordinates, not radians.
    return np.concatenate([
        np.linalg.eigvals(model.get_bulk_Hamiltonian_complex(np.exp(1j * np.array([a, b]))).toarray())
        for a in phases for b in phases
    ])


def finite_hamiltonian(model, nx, ny):
    """Assemble BerryPy's hopping list in O(number of bonds), discarding exits.

    This avoids constructing a huge periodic supercell just to delete its wrap
    bonds. check_model() compares it with BerryPy's supercell OBC Hamiltonian.
    """
    from scipy.sparse import coo_matrix

    if nx < 1 or ny < 1:
        raise ValueError("Both OBC sizes must be positive")
    cells = [(i, j) for i in range(nx) for j in range(ny)]
    lookup = {cell: k for k, cell in enumerate(cells)}
    n = model.SiteNum
    rows, cols, values = [], [], []
    bonds = [(*b, (0, 0)) for b in model.InCell] + list(model.InterCell)
    for index, (i, j) in enumerate(cells):
        for source, dest, amplitude, shift in bonds:
            target = lookup.get((i + shift[0], j + shift[1]))
            if target is not None:
                rows.append(target * n + dest)
                cols.append(index * n + source)
                values.append(amplitude)
    coords = np.vstack([(np.asarray(model.SiteCoord) + cell) @ model.LatticeVec.T for cell in cells])
    return coo_matrix((values, (rows, cols)), shape=(len(coords), len(coords))).tocsr(), coords


def atomic_pickle(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("wb") as stream:
        pickle.dump(value, stream, protocol=4)
    os.replace(temporary, path)


class LegacyUnpickler(pickle.Unpickler):
    """Map old class locations without modifying global sys.modules."""
    def find_class(self, module, name):
        if module in ("gbz_types", "bfgbz2d", "bfgbz2d.gbz_types", "bfgbz2d.core"):
            if name in ("GBZResult", "PointSubset", "LineSubset"):
                from pygbz2d import core
                return getattr(core, name)
        return super().find_class(module, name)


def load_sweep(path):
    with Path(path).open("rb") as stream:
        data = LegacyUnpickler(stream).load()
    expected = len(data["E_real"]) * len(data["E_imag"])
    if len(data["results"]) != expected:
        raise ValueError(f"{path}: result count does not match energy grid")
    energies = energy_grid(data["E_real"], data["E_imag"])
    for e, result in zip(energies, data["results"]):
        if result is not None and not np.isclose(e, result.E_ref, rtol=0, atol=1e-12):
            raise ValueError(f"{path}: result order does not match the energy grid")
    return data


def energy_grid(real, imag):
    er, ei = np.meshgrid(real, imag, indexing="xy")
    return (er + 1j * ei).ravel()


_POLYNOMIAL = None


def _init_worker(coefficients, degrees):
    global _POLYNOMIAL
    _POLYNOMIAL = coefficients, degrees


def _solve_energy(task):
    from pygbz2d.sgbz import collect_GBZ_subsets
    index, energy = task
    return index, collect_GBZ_subsets(*_POLYNOMIAL, complex(energy))


def sweep(direction, real, imag, output, workers=1, checkpoint_every=20, retry_failed=False):
    """Resume missing entries; retry failures only when explicitly requested."""
    coefficients, degrees = polynomial_data(direction_model(direction))
    output = Path(output)
    if output.exists():
        data = load_sweep(output)
        compatible = (data.get("direction") == direction
                      and np.array_equal(data["E_real"], real)
                      and np.array_equal(data["E_imag"], imag)
                      and np.array_equal(data["degs"], degrees)
                      and np.shape(data["coeffs"]) == coefficients.shape
                      and np.allclose(data["coeffs"], coefficients, rtol=1e-12, atol=1e-12)
                      and tuple(data["params"]) == PARAMS)
        if not compatible:
            raise ValueError(f"{output}: incompatible checkpoint; choose a new output path")
    else:
        data = dict(E_real=np.asarray(real), E_imag=np.asarray(imag),
                    results=[None] * (len(real) * len(imag)), params=PARAMS,
                    coeffs=coefficients, degs=degrees, direction=direction,
                    created_utc=datetime.now(timezone.utc).isoformat(), schema_version=1)
    tasks = [(i, e) for i, (e, r) in enumerate(zip(energy_grid(real, imag), data["results"]))
             if r is None or (retry_failed and not r.success)]
    print(f"{direction}: {len(tasks)} pending / {len(data['results'])} energies", flush=True)
    atomic_pickle(output, data)
    executor = None
    try:
        if workers == 1:
            _init_worker(coefficients, degrees)
            answers = map(_solve_energy, tasks)
        else:
            executor = ProcessPoolExecutor(max_workers=workers, initializer=_init_worker,
                                           initargs=(coefficients, degrees))
            answers = executor.map(_solve_energy, tasks, chunksize=1)
        for completed, (index, result) in enumerate(answers, 1):
            data["results"][index] = result
            if not result.success:
                print(f"FAILED {direction} E={result.E_ref}: {result.error}", flush=True)
            if completed % checkpoint_every == 0 or completed == len(tasks):
                atomic_pickle(output, data)
                print(f"{direction}: saved {completed}/{len(tasks)} -> {output}", flush=True)
    finally:
        atomic_pickle(output, data)
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
    return data


def subset_samples(data):
    """Keep every PointSubset and every stored row of each LineSubset."""
    from pygbz2d.core import PointSubset, LineSubset
    samples = []
    for result in data["results"]:
        if result is None or not result.success:
            continue
        for subset in result.subsets:
            if isinstance(subset, PointSubset):
                samples.append((result.E_ref, subset.beta1, subset.beta2))
            elif isinstance(subset, LineSubset):
                beta1 = np.exp(subset.mu1 + 1j * subset.theta1_arr)
                samples.extend(zip(np.full(len(beta1), result.E_ref), beta1, subset.beta2_arr))
            else:
                raise TypeError(f"Unsupported subset type: {type(subset)}")
    return np.asarray(samples, complex).reshape(-1, 3)


def sweep_report(data):
    samples = subset_samples(data)
    valid = np.all(np.isfinite(samples), axis=1) & np.all(np.abs(samples[:, 1:]) > 0, axis=1)
    mu = np.log(np.abs(samples[valid, 1:]))
    results = data["results"]
    report = dict(total=len(results), pending=sum(r is None for r in results),
                  failed=sum(r is not None and not r.success for r in results),
                  in_spectrum=sum(r is not None and r.is_gbz for r in results),
                  samples=len(samples), invalid_samples=int(np.sum(~valid)),
                  max_abs_mu=np.max(np.abs(mu), axis=0).tolist() if len(mu) else None,
                  p99_abs_mu=np.quantile(np.abs(mu), .99, axis=0).tolist() if len(mu) else None)
    report["failures"] = [dict(energy=[r.E_ref.real, r.E_ref.imag], error=r.error)
                          for r in results if r is not None and not r.success]
    return report


def summarize_eigensystem(eigenvalues, vectors, coords, window=(-2.5, -1.5), target=-2 + .2j):
    """Right-state densities, normalized separately for every eigenvector.

    Column batches avoid allocating another full N x N array for legacy runs.
    The energy window is a reproducible bulk-band selection, not a topological
    classifier. All-state density is also retained for comparison.
    """
    eigenvalues, coords = np.asarray(eigenvalues), np.asarray(coords)
    n = len(eigenvalues)
    if vectors.shape != (n, n) or coords.shape != (n, 2):
        raise ValueError("OBC data must contain N eigenvalues, N x N right vectors and N x 2 coordinates")
    mask = (eigenvalues.real >= window[0]) & (eigenvalues.real <= window[1])
    candidates = np.flatnonzero(mask)
    if not len(candidates):
        raise ValueError(f"No states in the selected Re(E) window {window}")
    selected = candidates[np.argmin(np.abs(eigenvalues[candidates] - target))]
    total = np.zeros(n)
    bulk = np.zeros(n)
    ipr = np.zeros(n)
    selected_density = None
    for start in range(0, n, 128):
        stop = min(start + 128, n)
        probabilities = np.abs(vectors[:, start:stop]) ** 2
        norms = probabilities.sum(axis=0)
        if np.any(~np.isfinite(norms)) or np.any(norms <= 0):
            raise ValueError("Zero or nonfinite eigenvector norm")
        probabilities /= norms
        total += probabilities.sum(axis=1)
        bulk += probabilities[:, mask[start:stop]].sum(axis=1)
        ipr[start:stop] = np.sum(probabilities ** 2, axis=0)
        if start <= selected < stop:
            selected_density = probabilities[:, selected - start].copy()
    return dict(eigenvalues=eigenvalues, coords=coords, density_all=total / n,
                density_window=bulk / len(candidates), ipr=ipr,
                selected_index=np.array(selected), selected_density=selected_density,
                window=np.asarray(window), window_count=np.array(len(candidates)),
                target=np.array(target))


def save_npz(path, data, metadata):
    payload = dict(data, metadata=np.array(json.dumps(metadata)))
    if str(path) == "-":
        np.savez_compressed(sys.stdout.buffer, **payload)
    else:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("wb") as stream:
            np.savez_compressed(stream, **payload)


def calculate_obc(geometry, nx, ny, output, window=(-2.5, -1.5), save_vectors=False):
    from scipy.linalg import eig
    direction = "a1" if geometry == "rhombus" else "x"
    model = direction_model(direction)
    matrix, coords = finite_hamiltonian(model, nx, ny)
    n = matrix.shape[0]
    print(f"OBC {geometry}: {n} sites; one dense complex matrix uses {16*n*n/2**30:.3f} GiB; "
          "eigensolver workspace needs additional memory", flush=True)
    started = time.monotonic()
    eigenvalues, vectors = eig(matrix.toarray(), check_finite=True)
    summary = summarize_eigensystem(eigenvalues, vectors, coords, window)
    residual = np.linalg.norm(matrix @ vectors[:, int(summary["selected_index"])]
                              - eigenvalues[int(summary["selected_index"])]
                              * vectors[:, int(summary["selected_index"])])
    metadata = dict(geometry=geometry, nx=nx, ny=ny, params=[str(v) for v in PARAMS],
                    source="computed by geometry-dependent-skin-effect.py",
                    created_utc=datetime.now(timezone.utc).isoformat(),
                    elapsed_seconds=time.monotonic()-started, selected_residual=float(residual))
    save_npz(output, summary, metadata)
    if save_vectors:
        atomic_pickle(Path(output).with_suffix(".pkl"), (eigenvalues, vectors, coords))
    print(f"saved {output}; selected eigenpair residual {residual:.3e}", flush=True)
    return summary


def summarize_legacy(path, output, window=(-2.5, -1.5)):
    path = Path(path).expanduser()
    with path.open("rb") as stream:
        eigenvalues, vectors, coords = pickle.load(stream)
    data = summarize_eigensystem(eigenvalues, vectors, coords, window)
    save_npz(output, data, dict(source=str(path.resolve()), source_bytes=path.stat().st_size,
                              source_mtime_utc=datetime.fromtimestamp(path.stat().st_mtime,
                                                                    timezone.utc).isoformat(),
                              parameter_provenance="Legacy companion script ALL_PARAMS; pickle has no parameter metadata",
                              extracted_utc=datetime.now(timezone.utc).isoformat()))


def _plot_setup():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                         "savefig.dpi": 180, "svg.fonttype": "none"})
    return plt


def _save_figure(figure, directory, name):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    for extension in ("svg", "png"):
        figure.savefig(directory / f"{name}.{extension}", bbox_inches="tight")


def plot_sgbz(sgbz_dir, output_dir, directions=DIRECTIONS, nk=81):
    plt = _plot_setup()
    from matplotlib.colors import Normalize

    # Four columns must remain legible after fitting the figure to an A4 page.
    plt.rcParams.update({"font.size": 12})
    reference = bz_spectrum(build_model(), nk)
    fig, axes = plt.subplots(2, len(directions), figsize=(2.5*len(directions)+.4, 4.8),
                             squeeze=False, constrained_layout=True)
    reports = {}
    color_normalization = Normalize(0, .6)
    for column, direction in enumerate(directions):
        path = Path(sgbz_dir) / f"Haldane-gain-loss-{direction}-SGBZ.pkl"
        data = load_sweep(path)
        if not np.allclose(np.asarray(data["params"], complex), np.asarray(PARAMS, complex)):
            raise ValueError(f"{path}: parameters differ from the plotted BZ reference")
        reports[direction] = dict(source=str(path.resolve()), **sweep_report(data))
        if reports[direction]["invalid_samples"]:
            raise ValueError(f"{path}: cannot plot zero or nonfinite Bloch factors")
        samples = subset_samples(data)
        energies = np.array([r.E_ref for r in data["results"] if r is not None and r.is_gbz])
        failed = np.array([r.E_ref for r in data["results"] if r is not None and not r.success], complex)
        ax = axes[0, column]
        ax.scatter(reference.real, reference.imag, s=1, c="#d5d9dc", rasterized=True, label="BZ")
        ax.scatter(energies.real, energies.imag, s=2, c=COLORS[direction], rasterized=True, label="SGBZ")
        if len(failed):
            ax.scatter(failed.real, failed.imag, marker="x", s=10, c="black", label="Failed")
        ax.set(xlabel=r"Re $E/t_1$", ylabel=r"Im $E/t_1$" if column == 0 else "", xlim=(-3.15, 4.65), ylim=(-.54, .54),
               title=f"{direction}: {'armchair' if direction == 'y' else 'zigzag'}"
               + (f" ({reports[direction]['pending']} pending)" if reports[direction]['pending'] else ""))
        ax.legend(loc="upper right", markerscale=3, fontsize=8, frameon=False)
        ax = axes[1, column]
        if len(samples):
            mu = np.log(np.abs(samples[:, 1:]))
            departure = np.max(np.abs(mu), axis=1)
            ax.scatter(np.angle(samples[:, 1]), np.angle(samples[:, 2]), c=departure,
                       cmap="magma", norm=color_normalization, s=1.3, rasterized=True)
        ax.set(xlabel=r"$\theta_1$", ylabel=r"$\theta_2$" if column == 0 else "", xlim=(-np.pi, np.pi), ylim=(-np.pi, np.pi))
        ax.set_xticks([-np.pi, 0, np.pi], [r"$-\pi$", "0", r"$\pi$"])
        ax.set_yticks([-np.pi, 0, np.pi], [r"$-\pi$", "0", r"$\pi$"])
    scalar = plt.cm.ScalarMappable(norm=color_normalization, cmap="magma")
    fig.colorbar(scalar, ax=list(axes[1]), shrink=.85, label=r"$\max_j |\ln|\beta_j||$")
    _save_figure(fig, output_dir, "sgbz-comparison")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.3), constrained_layout=True)
    for direction in directions:
        samples = subset_samples(load_sweep(Path(sgbz_dir) / f"Haldane-gain-loss-{direction}-SGBZ.pkl"))
        if len(samples):
            mu = np.log(np.abs(samples[:, 1:]))
            for j, ax in enumerate(axes):
                ax.scatter(samples[:, 0].real, mu[:, j], s=1.2, c=COLORS[direction],
                           label=direction, alpha=.65, rasterized=True)
    for j, ax in enumerate(axes):
        ax.axhline(0, color="black", lw=.5)
        ax.set(xlabel=r"Re $E/t_1$", ylabel=rf"$\mu_{j+1}=\ln|\beta_{j+1}|$")
        if ax.get_legend_handles_labels()[0]:
            ax.legend(frameon=False, markerscale=4)
        if j == 0:
            ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
    _save_figure(fig, output_dir, "sgbz-radii")
    plt.close(fig)
    Path(output_dir, "sgbz-summary.json").write_text(json.dumps(reports, indent=2), encoding="utf-8")
    return reports


def plot_geometries(output_dir):
    plt = _plot_setup()
    from matplotlib.collections import LineCollection
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.5), constrained_layout=True)
    for ax, geometry, direction, nx, ny in zip(axes, ("Rhombus", "Rectangle"), ("a1", "x"), (8, 10), (8, 5)):
        model = direction_model(direction)
        matrix, coords = finite_hamiltonian(model, nx, ny)
        coo = matrix.tocoo()
        bonds = [(coords[i], coords[j]) for i, j, v in zip(coo.row, coo.col, coo.data)
                 if i < j and np.isclose(abs(v), PARAMS[0])]
        ax.add_collection(LineCollection(bonds, colors="#9ba5b0", linewidths=.7))
        ax.scatter(coords[::2, 0], coords[::2, 1], s=12, color="#c24752", label="A: gain")
        ax.scatter(coords[1::2, 0], coords[1::2, 1], s=12, color="#305f9e", label="B: loss")
        for j, label in enumerate((r"$a_1$", r"$a_2$") if direction == "a1" else (r"$x$", r"$y$")):
            origin = coords.mean(axis=0)
            end = origin + 2.0 * model.LatticeVec[:, j] / np.linalg.norm(model.LatticeVec[:, j])
            ax.annotate("", xy=end, xytext=origin, arrowprops=dict(arrowstyle="->", color="black", lw=1.5))
            ax.annotate(label, xy=end, xytext=(-9, 6), textcoords="offset points", fontsize=12,
                        bbox=dict(facecolor="white", edgecolor="none", alpha=.75, pad=.5))
        ax.set_aspect("equal")
        ax.set(xlabel="Cartesian x", ylabel="Cartesian y", title=geometry +
               (": zigzag sides" if direction == "a1" else ": zigzag + armchair sides"))
        ax.legend(loc="upper right", fontsize=8, frameon=True, facecolor="white", framealpha=.95)
    _save_figure(fig, output_dir, "geometries")
    plt.close(fig)


def plot_obc(paths, output_dir, nk=81):
    plt = _plot_setup()
    from matplotlib.colors import LogNorm
    reference = bz_spectrum(build_model(), nk)
    summaries = []
    for path in paths:
        with np.load(path, allow_pickle=False) as archive:
            summaries.append({k: archive[k] for k in archive.files})
    fig, axes = plt.subplots(2, len(paths), figsize=(5*len(paths), 6.5), squeeze=False, constrained_layout=True)
    largest_density = max(float(np.max(len(d["coords"]) * d["density_window"])) for d in summaries)
    norm = LogNorm(vmin=.02, vmax=max(2, largest_density))
    for column, (path, data) in enumerate(zip(paths, summaries)):
        eigenvalues, coords = data["eigenvalues"], data["coords"]
        metadata = json.loads(str(data["metadata"]))
        label = metadata.get("geometry", "rectangle" if "square" in metadata["source"] else "rhombus")
        ax = axes[0, column]
        ax.scatter(reference.real, reference.imag, s=1, c="#d5d9dc", rasterized=True, label="BZ")
        ax.scatter(eigenvalues.real, eigenvalues.imag, s=2, c="#305f9e", rasterized=True, label="OBC")
        ax.axvspan(*data["window"], color="#d29229", alpha=.14, label="Density window")
        ax.set(title=f"{label.capitalize()}: {len(eigenvalues):,} sites", xlabel=r"Re $E/t_1$", ylabel=r"Im $E/t_1$",
               xlim=(-3.15, 4.65), ylim=(-.54, .54))
        ax.legend(frameon=False, fontsize=8, markerscale=3)
        ax = axes[1, column]
        artist = ax.scatter(coords[:, 0], coords[:, 1], c=len(coords) * data["density_window"],
                            cmap="inferno", norm=norm, s=5, rasterized=True)
        ax.set_aspect("equal")
        ax.set(xlabel="Cartesian x", ylabel="Cartesian y",
               title=f"Mean right-state density: {int(data['window_count']):,} states")
    fig.colorbar(artist, ax=list(axes[1]), label=r"$N\overline{|\psi_R(r)|^2}$", shrink=.85)
    _save_figure(fig, output_dir, "obc-comparison")
    plt.close(fig)


def check_model():
    """Independent determinant and finite-boundary checks, without an SGBZ sweep."""
    rng = np.random.default_rng(20260925)
    diagnostics = {}
    for direction in DIRECTIONS:
        model = direction_model(direction)
        coefficients, degrees = polynomial_data(model)
        errors = []
        for _ in range(8):
            variables = np.r_[rng.normal() + 1j*rng.normal(), np.exp(rng.normal(0, .15, 2) + 1j*rng.uniform(-np.pi, np.pi, 2))]
            actual = np.sum(coefficients * np.prod(variables ** degrees, axis=1))
            expected = np.linalg.det(variables[0] * np.eye(model.SiteNum)
                                     - model.get_bulk_Hamiltonian_complex(variables[1:]).toarray())
            errors.append(abs(actual-expected) / max(1, abs(expected)))
        assert max(errors) < 1e-10, (direction, errors)
        direct, coords = finite_hamiltonian(model, 3, 2)
        supercell = model.get_supercell([(0, j) for j in range(2)], np.diag([1, 2]))
        supercell = supercell.get_supercell([(i, 0) for i in range(3)], np.diag([3, 1]))
        berry = supercell.get_bulk_Hamiltonian_complex((None, None)).toarray()
        matrix_error = float(np.max(np.abs(direct.toarray()-berry)))
        coord_error = float(np.max(np.abs(coords-supercell.lattice2cart(supercell.SiteCoord.T).T)))
        assert matrix_error < 1e-12 and coord_error < 1e-12
        diagnostics[direction] = dict(sites_per_cell=model.SiteNum, terms=len(coefficients),
                                     determinant_relative_error=max(errors), obc_matrix_error=matrix_error,
                                     obc_coordinate_error=coord_error)
    print(json.dumps(diagnostics, indent=2))
    return diagnostics


def positive_int(value):
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError("must be positive")
    return number


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command")
    demo = commands.add_parser("demo", help="Small end-to-end run (default); does not overwrite legacy data")
    demo.add_argument("--output-dir", type=Path, default=DATA / "demo")
    demo.add_argument("--workers", type=positive_int, default=1)
    commands.add_parser("check", help="Compare BerryPy determinants and finite OBC matrices")
    sweep_parser = commands.add_parser("sweep", help="Compute/resume SGBZ energy grids")
    sweep_parser.add_argument("--directions", nargs="+", choices=DIRECTIONS, default=list(DIRECTIONS))
    sweep_parser.add_argument("--n-real", type=positive_int, default=7)
    sweep_parser.add_argument("--n-imag", type=positive_int, default=5)
    sweep_parser.add_argument("--real-range", type=float, nargs=2, default=(-3.1, 4.6), metavar=("MIN", "MAX"))
    sweep_parser.add_argument("--imag-range", type=float, nargs=2, default=(-.51, .51), metavar=("MIN", "MAX"))
    sweep_parser.add_argument("--workers", type=positive_int, default=1)
    sweep_parser.add_argument("--checkpoint-every", type=positive_int, default=20)
    sweep_parser.add_argument("--retry-failed", action="store_true")
    sweep_parser.add_argument("--output-dir", type=Path, default=DATA / "sweep")
    obc = commands.add_parser("obc", help="Compute finite OBC eigensystem and save a compact summary")
    obc.add_argument("--geometry", choices=("rhombus", "rectangle"), required=True)
    obc.add_argument("--nx", type=positive_int, default=8)
    obc.add_argument("--ny", type=positive_int, default=8)
    obc.add_argument("--output", type=Path, required=True)
    obc.add_argument("--save-vectors", action="store_true", help="Also retain full eigenvectors as a trusted pickle")
    obc.add_argument("--window", type=float, nargs=2, default=(-2.5, -1.5))
    extract = commands.add_parser("summarize-obc", help="Extract a legacy (eigenvalues, right vectors, coordinates) pickle")
    extract.add_argument("input", type=Path)
    extract.add_argument("--output", required=True, help="NPZ path; '-' writes binary NPZ to stdout for SSH streaming")
    extract.add_argument("--window", type=float, nargs=2, default=(-2.5, -1.5))
    plot = commands.add_parser("plot", help="Plot saved SGBZ and optional OBC summaries; no GBZ solving")
    plot.add_argument("--sgbz-dir", type=Path, default=ROOT / "data")
    plot.add_argument("--directions", nargs="+", choices=DIRECTIONS, default=list(DIRECTIONS))
    plot.add_argument("--obc", type=Path, nargs="*", default=[])
    plot.add_argument("--output-dir", type=Path, default=FIGURES)
    plot.add_argument("--nk", type=positive_int, default=81)
    args = parser.parse_args(argv)
    if args.command is None:
        args = parser.parse_args(["demo"])
    if args.command == "check":
        check_model()
    elif args.command == "sweep":
        real = np.linspace(*args.real_range, args.n_real)
        imag = np.linspace(*args.imag_range, args.n_imag)
        for direction in args.directions:
            result = sweep(direction, real, imag, args.output_dir / f"Haldane-gain-loss-{direction}-SGBZ.pkl",
                           args.workers, args.checkpoint_every, args.retry_failed)
            print(json.dumps(sweep_report(result), indent=2))
    elif args.command == "obc":
        calculate_obc(args.geometry, args.nx, args.ny, args.output, args.window, args.save_vectors)
    elif args.command == "summarize-obc":
        summarize_legacy(args.input, args.output, args.window)
    elif args.command == "plot":
        plot_geometries(args.output_dir)
        print(json.dumps(plot_sgbz(args.sgbz_dir, args.output_dir, args.directions, args.nk), indent=2))
        if args.obc:
            plot_obc(args.obc, args.output_dir, args.nk)
    elif args.command == "demo":
        check_model()
        for direction in DIRECTIONS:
            # Two interior energies and one exterior energy keep the demo small;
            # the real-energy point also exercises continuum LineSubsets for y.
            data = sweep(direction, np.array([-2., 1.5, 5.]), np.array([0., .3]),
                         args.output_dir / f"Haldane-gain-loss-{direction}-SGBZ.pkl", args.workers, 1)
            report = sweep_report(data)
            if (report["failed"] or report["pending"] or report["invalid_samples"]
                    or not 0 < report["in_spectrum"] < report["total"]):
                raise RuntimeError(f"Demo did not cover interior/exterior successfully: {direction}: {report}")
            if direction == "y" and not any(r.index[1] for r in data["results"]):
                raise RuntimeError("The y demo did not exercise a continuum LineSubset")
        paths = []
        for geometry, nx, ny in (("rhombus", 8, 8), ("rectangle", 8, 4)):
            path = args.output_dir / f"obc-{geometry}.npz"
            calculate_obc(geometry, nx, ny, path)
            paths.append(path)
        output_dir = args.output_dir / "figures"
        plot_geometries(output_dir)
        plot_sgbz(args.output_dir, output_dir, nk=31)
        plot_obc(paths, output_dir, nk=31)
        print(f"Demo complete: {output_dir}")


if __name__ == "__main__":
    main()
