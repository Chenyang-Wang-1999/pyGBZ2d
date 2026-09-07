"""Compute actual pyGBZ2d data for the companion plotting script.

Use --model hn2d for the documented real-coupling benchmark, or --polynomial
model.npz with arrays coeffs (complex, n) and degs (integer, n, 3).
No saved pickles are loaded. Failures and numerical warnings remain in JSON.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import sys
import time
import warnings

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
import numpy as np
import scipy
from pygbz2d import CharPoly, PointSubset
from pygbz2d.continuation import ZeroManager
from pygbz2d.sgbz import collect_GBZ_subsets as sgbz
from pygbz2d.amoeba import collect_GBZ_subsets as amoeba


def zpair(z):
    return [float(z.real), float(z.imag)]


def safe_json(obj):
    if isinstance(obj, dict):
        return {k: safe_json(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [safe_json(v) for v in obj]
    if isinstance(obj, float) and not np.isfinite(obj):
        return None
    return obj


def samples(sub):
    if isinstance(sub, PointSubset):
        return [(sub.beta1, sub.beta2)]
    return list(zip(np.exp(sub.mu1 + 1j * sub.theta1_arr), sub.beta2_arr))


def normalized_residual(coeffs, degs, energy, b1, b2):
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        terms = coeffs * np.prod(np.array([energy, b1, b2], complex) ** degs, axis=1)
        denominator = float(np.abs(terms).sum())
        return float(abs(terms.sum()) / denominator) if denominator > 0 else float("nan")


def axis(spec):
    low, high, n = spec
    if not np.isfinite([low, high, n]).all() or n != int(n) or n < 1 or high < low:
        raise ValueError("grid requires finite low <= high and a positive integer count")
    return np.linspace(low, high, int(n))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    model = p.add_mutually_exclusive_group(required=True)
    model.add_argument("--model", choices=["hn2d"])
    model.add_argument("--polynomial", type=Path)
    p.add_argument("--J1", type=float, default=1.)
    p.add_argument("--J2", type=float, default=1.)
    p.add_argument("--gamma1", type=float, default=.2)
    p.add_argument("--gamma2", type=float, default=.3)
    p.add_argument("--energy", type=complex, required=True, help="fixed energy for root tracks, e.g. 1+0j")
    p.add_argument("--mu1", type=float, required=True, help="fixed log radius for root tracks")
    p.add_argument("--sweep-re", nargs=3, type=float, metavar=("LOW", "HIGH", "COUNT"))
    p.add_argument("--sweep-im", nargs=3, type=float, default=(0., 0., 1.), metavar=("LOW", "HIGH", "COUNT"))
    p.add_argument("--methods", nargs="+", choices=["sgbz", "amoeba"], default=["sgbz", "amoeba"])
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--overwrite", action="store_true")
    args = p.parse_args()
    if args.output.exists() and not args.overwrite:
        p.error("output exists; use a new path or --overwrite")
    if not np.isfinite([args.energy.real, args.energy.imag, args.mu1]).all():
        p.error("energy and mu1 must be finite")
    if args.polynomial:
        with np.load(args.polynomial, allow_pickle=False) as data:
            coeffs = np.asarray(data["coeffs"], complex)
            raw_degs = np.asarray(data["degs"])
        if not np.isfinite(raw_degs).all() or not np.equal(raw_degs, np.round(raw_degs)).all():
            p.error("degs must contain finite integer exponents")
        degs = raw_degs.astype(int)
        model_name = str(args.polynomial.resolve())
    else:
        pars = np.array([args.J1, args.J2, args.gamma1, args.gamma2])
        if not np.isfinite(pars).all() or min(args.J1, args.J2) <= 0:
            p.error("HN benchmark requires finite parameters and positive J1, J2")
        coeffs = np.array([1, -args.J1*np.exp(args.gamma1), -args.J1*np.exp(-args.gamma1),
                           -args.J2*np.exp(args.gamma2), -args.J2*np.exp(-args.gamma2)], complex)
        degs = np.array([[1,0,0],[0,-1,0],[0,1,0],[0,0,-1],[0,0,1]], int)
        model_name = "HN2D, real J, zero delta, axis-aligned basis (10)"
    if coeffs.ndim != 1 or degs.shape != (len(coeffs), 3) or not np.isfinite(coeffs).all():
        p.error("expected finite coeffs (n,) and degs (n,3)")
    poly = CharPoly(coeffs, degs, backend="numpy")
    if min(poly.get_minor_degrees()) <= 0:
        p.error("this experiment requires M > 0 and N > 0 in beta2")
    try:
        energies = ([complex(r, i) for i in axis(args.sweep_im) for r in axis(args.sweep_re)]
                    if args.sweep_re else [args.energy])
    except ValueError as exc:
        p.error(str(exc))
    # Both public solver factories obey this environment setting.
    import os
    os.environ["POLY_BACKEND"] = "numpy"
    revision = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO,
                              capture_output=True, text=True)
    result = {"schema_version": 1, "model": model_name,
              "created_utc": datetime.now(timezone.utc).isoformat(),
              "revision": revision.stdout.strip() if revision.returncode == 0 else None,
              "numpy": np.__version__, "scipy": scipy.__version__,
              "package_source": str(REPO / "src" / "pygbz2d"),
              "backend": "numpy", "solver_options": "repository defaults; plateau_check=True",
              "command": sys.argv, "coeffs": [zpair(z) for z in coeffs], "degs": degs.tolist(),
              "tracks": {"energy": zpair(args.energy), "mu1": args.mu1, "segments": []},
              "sweep": []}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            zm = ZeroManager(poly, args.energy, args.mu1)
            zm.run()
            result["tracks"]["segments"] = [
                {"theta1": s.theta1_arr.tolist(),
                 "beta2": [[zpair(z) for z in row] for row in s.tracked_roots]}
                for s in zm.segments]
            result["tracks"]["multiple_root_angles"] = [float(m.theta1) for m in zm.multiple_roots]
        except Exception as exc:
            result["tracks"]["error"] = f"{type(exc).__name__}: {exc}"
        result["tracks"]["warnings"] = [str(w.message) for w in caught]
    for energy in energies:
        for name in dict.fromkeys(args.methods):
            started = time.perf_counter()
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                gbz = {"sgbz": sgbz, "amoeba": amoeba}[name](coeffs, degs, energy)
            entry = {"energy": zpair(energy), "method": name, "success": gbz.success,
                     "error": gbz.error, "index": list(gbz.index), "is_gbz": gbz.is_gbz,
                     "elapsed_s": time.perf_counter()-started,
                     "warnings": [str(w.message) for w in caught], "subsets": []}
            residuals = []
            for sub in gbz.subsets:
                values = samples(sub)
                entry["subsets"].append({"kind": "point" if isinstance(sub, PointSubset) else "line",
                    "beta1": [zpair(b1) for b1, _ in values], "beta2": [zpair(b2) for _, b2 in values]})
                residuals.extend(normalized_residual(coeffs, degs, energy, b1, b2) for b1, b2 in values)
            entry["max_normalized_residual"] = max(residuals) if residuals else None
            result["sweep"].append(entry)
            print(f"{name}: E={energy}, success={gbz.success}, index={gbz.index}", flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(safe_json(result), indent=2, allow_nan=False), encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()
