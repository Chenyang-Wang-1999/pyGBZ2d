'''
Diagnostic experiment for demos/imaginary-degeneracy-splitting.py.

This script intentionally does not modify the amoeba implementation.  It
measures the a1 winding curve obtained after solving a2(mu1, mu2) = 0, then
checks the two criteria discussed in the debugging note:

1. a1(mu1) should not decrease.
2. If a1(mu1) has a finite zero plateau, the reference energy is outside the
   amoeba spectrum; returning "is_amoeba=True" at the plateau edge indicates a
   numerical/search bug.
'''

from __future__ import annotations

import argparse
import csv
import importlib.util
import sys
from dataclasses import dataclass
from math import pi
from pathlib import Path
from typing import Iterable

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np
import poly_tools as pt

import brute_force_amoeba.amoeba as amo


DEMO_PATH = ROOT / "demos" / "imaginary-degeneracy-splitting.py"
DEFAULT_E_REF = 1.648 + 0.0294j


@dataclass
class CurvePoint:
    mu1: float
    mu2: float | None
    a1: float | None
    a2: float | tuple[float, float] | None
    zero_count: int
    net_zero_count: int
    is_continuum: bool
    error: str | None = None


def load_demo_module():
    spec = importlib.util.spec_from_file_location("imaginary_degeneracy_splitting", DEMO_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import demo from {DEMO_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def build_char_poly(demo_module):
    model = demo_module.get_model(**demo_module.DEFAULT_PARAMS)
    coeffs, degs = model.get_characteristic_polynomial_data()
    char_poly = pt.CLaurent(3)
    char_poly.set_Laurent_by_terms(
        pt.CScalarVec(coeffs),
        pt.CLaurentIndexVec(degs.flatten()),
    )
    return char_poly, coeffs, degs


def angular_distance(a: float, b: float) -> float:
    return abs(((a - b + pi) % (2 * pi)) - pi)


def net_zero_count(zeros: Iterable[tuple], tol: float = 1e-7) -> int:
    """Count crossing clusters whose jump directions do not cancel."""
    clusters: list[dict] = []
    for zero in zeros:
        t1, t2 = float(zero[0]), float(zero[1])
        jump = int(zero[2]) if len(zero) >= 3 else 1
        for cluster in clusters:
            if (
                angular_distance(t1, cluster["theta1"]) < tol
                and angular_distance(t2, cluster["theta2"]) < tol
            ):
                cluster["jump_sum"] += jump
                cluster["size"] += 1
                break
        else:
            clusters.append({"theta1": t1, "theta2": t2, "jump_sum": jump, "size": 1})
    return sum(1 for cluster in clusters if cluster["jump_sum"] != 0)


def solve_curve_point(
    char_poly,
    e_ref: complex,
    mu1: float,
    n_points: int,
    xtol: float,
    continuum_tol: float,
    min_continuum_pts: int,
) -> CurvePoint:
    try:
        inner = amo._find_mu2_for_a2_zero(
            char_poly,
            e_ref,
            mu1,
            -1.0,
            1.0,
            target_winding=0.0,
            N_points=n_points,
            continuum_tol=continuum_tol,
            min_continuum_pts=min_continuum_pts,
            xtol=xtol,
        )
        a1 = amo._get_average_winding_from_zeros(
            char_poly,
            e_ref,
            mu1,
            inner["mu2"],
            inner["zeros"],
            direction=1,
        )
        return CurvePoint(
            mu1=mu1,
            mu2=float(inner["mu2"]),
            a1=float(a1),
            a2=inner["winding"],
            zero_count=len(inner["zeros"]),
            net_zero_count=net_zero_count(inner["zeros"]),
            is_continuum=bool(inner["is_continuum"]),
        )
    except Exception as exc:
        return CurvePoint(
            mu1=mu1,
            mu2=None,
            a1=None,
            a2=None,
            zero_count=0,
            net_zero_count=0,
            is_continuum=False,
            error=repr(exc),
        )


def sweep_curve(
    char_poly,
    e_ref: complex,
    mu1_low: float,
    mu1_high: float,
    samples: int,
    n_points: int,
    xtol: float,
    continuum_tol: float,
    min_continuum_pts: int,
) -> list[CurvePoint]:
    return [
        solve_curve_point(
            char_poly,
            e_ref,
            float(mu1),
            n_points,
            xtol,
            continuum_tol,
            min_continuum_pts,
        )
        for mu1 in np.linspace(mu1_low, mu1_high, samples)
    ]


def is_zero_plateau_point(point: CurvePoint, winding_tol: float) -> bool:
    return (
        point.error is None
        and point.a1 is not None
        and abs(point.a1) <= winding_tol
        and point.zero_count == 0
        and not point.is_continuum
    )


def find_plateau_intervals(points: list[CurvePoint], winding_tol: float):
    intervals = []
    start = None
    for idx, point in enumerate(points):
        if is_zero_plateau_point(point, winding_tol):
            if start is None:
                start = idx
        elif start is not None:
            intervals.append((start, idx - 1))
            start = None
    if start is not None:
        intervals.append((start, len(points) - 1))
    return intervals


def find_decreases(points: list[CurvePoint], monotonic_tol: float):
    decreases = []
    valid = [point for point in points if point.error is None and point.a1 is not None]
    for left, right in zip(valid, valid[1:]):
        delta = right.a1 - left.a1
        if delta < -monotonic_tol:
            decreases.append((left.mu1, right.mu1, left.a1, right.a1, delta))
    return decreases


def refine_plateau_boundary(
    char_poly,
    e_ref: complex,
    left_mu1: float,
    right_mu1: float,
    left_is_plateau: bool,
    n_points: int,
    xtol: float,
    continuum_tol: float,
    min_continuum_pts: int,
    winding_tol: float,
    iterations: int = 28,
) -> float:
    left = left_mu1
    right = right_mu1
    for _ in range(iterations):
        mid = 0.5 * (left + right)
        point = solve_curve_point(
            char_poly,
            e_ref,
            mid,
            n_points,
            xtol,
            continuum_tol,
            min_continuum_pts,
        )
        mid_is_plateau = is_zero_plateau_point(point, winding_tol)
        if mid_is_plateau == left_is_plateau:
            left = mid
        else:
            right = mid
    return 0.5 * (left + right)


def write_csv(path: Path, points: list[CurvePoint]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fp:
        writer = csv.writer(fp)
        writer.writerow([
            "mu1",
            "mu2",
            "a1",
            "a2",
            "zero_count",
            "net_zero_count",
            "is_continuum",
            "error",
        ])
        for point in points:
            writer.writerow([
                f"{point.mu1:.16g}",
                "" if point.mu2 is None else f"{point.mu2:.16g}",
                "" if point.a1 is None else f"{point.a1:.16g}",
                repr(point.a2),
                point.zero_count,
                point.net_zero_count,
                int(point.is_continuum),
                "" if point.error is None else point.error,
            ])


def format_check_result(result: dict) -> str:
    zeros = result.get("zeros", [])
    fields = [
        f"success={result.get('success')}",
        f"is_amoeba={result.get('is_amoeba')}",
        f"classification={result.get('_classification_reason')}",
        f"plateau_check={result.get('_plateau_check')}",
        f"mu1={result.get('mu1')}",
        f"mu2={result.get('mu2')}",
        f"is_continuum={result.get('is_continuum')}",
        f"zero_count={len(zeros)}",
        f"net_zero_count={net_zero_count(zeros)}",
    ]
    return ", ".join(fields)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--energy-real", type=float, default=DEFAULT_E_REF.real)
    parser.add_argument("--energy-imag", type=float, default=DEFAULT_E_REF.imag)
    parser.add_argument("--mu1-low", type=float, default=0.09)
    parser.add_argument("--mu1-high", type=float, default=0.13)
    parser.add_argument("--samples", type=int, default=81)
    parser.add_argument("--n-points", type=int, default=901)
    parser.add_argument("--xtol", type=float, default=1e-10)
    parser.add_argument("--winding-tol", type=float, default=1e-8)
    parser.add_argument("--monotonic-tol", type=float, default=1e-6)
    parser.add_argument("--continuum-tol", type=float, default=1e-8)
    parser.add_argument("--min-continuum-pts", type=int, default=3)
    parser.add_argument(
        "--resolution-n-points",
        type=int,
        nargs="+",
        default=[301, 401, 601, 901, 1201],
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=ROOT / "diagnostics" / "imaginary-degeneracy-splitting",
    )
    args = parser.parse_args()

    e_ref = complex(args.energy_real, args.energy_imag)
    demo = load_demo_module()
    char_poly, coeffs, degs = build_char_poly(demo)

    check_result = amo.check_amoeba(coeffs, degs, e_ref, 0.0, True)
    returned_mu1 = float(check_result["mu1"]) if check_result.get("success") else None

    resolution_rows = []
    if returned_mu1 is not None:
        for n_points in args.resolution_n_points:
            point = solve_curve_point(
                char_poly,
                e_ref,
                returned_mu1,
                n_points,
                args.xtol,
                args.continuum_tol,
                args.min_continuum_pts,
            )
            resolution_rows.append((n_points, point))

    points = sweep_curve(
        char_poly,
        e_ref,
        args.mu1_low,
        args.mu1_high,
        args.samples,
        args.n_points,
        args.xtol,
        args.continuum_tol,
        args.min_continuum_pts,
    )

    decreases = find_decreases(points, args.monotonic_tol)
    plateau_intervals = find_plateau_intervals(points, args.winding_tol)

    refined_intervals = []
    for start, end in plateau_intervals:
        left_boundary = points[start].mu1
        right_boundary = points[end].mu1
        if start > 0:
            left_boundary = refine_plateau_boundary(
                char_poly,
                e_ref,
                points[start - 1].mu1,
                points[start].mu1,
                left_is_plateau=False,
                n_points=args.n_points,
                xtol=args.xtol,
                continuum_tol=args.continuum_tol,
                min_continuum_pts=args.min_continuum_pts,
                winding_tol=args.winding_tol,
            )
        if end + 1 < len(points):
            right_boundary = refine_plateau_boundary(
                char_poly,
                e_ref,
                points[end].mu1,
                points[end + 1].mu1,
                left_is_plateau=True,
                n_points=args.n_points,
                xtol=args.xtol,
                continuum_tol=args.continuum_tol,
                min_continuum_pts=args.min_continuum_pts,
                winding_tol=args.winding_tol,
            )
        refined_intervals.append((left_boundary, right_boundary))

    args.out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.out_dir / "a1_vs_mu1.csv"
    report_path = args.out_dir / "report.md"
    write_csv(csv_path, points)

    valid_deltas = [
        right.a1 - left.a1
        for left, right in zip(points, points[1:])
        if left.a1 is not None and right.a1 is not None
    ]
    min_delta = min(valid_deltas) if valid_deltas else None

    lines = [
        "# imaginary-degeneracy-splitting diagnostic report",
        "",
        f"Reference energy: `{e_ref}`",
        f"Curve definition: solve `a2(mu1, mu2)=0`, then evaluate `a1(mu1, mu2)`.",
        "",
        "## check_amoeba result",
        "",
        format_check_result(check_result),
        "",
        "## Resolution sensitivity at returned mu1",
        "",
        "| N_points | mu2 | a1 | zero_count | net_zero_count | is_continuum |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for n_points, point in resolution_rows:
        lines.append(
            "| "
            f"{n_points} | "
            f"{'' if point.mu2 is None else f'{point.mu2:.12g}'} | "
            f"{'' if point.a1 is None else f'{point.a1:.12g}'} | "
            f"{point.zero_count} | {point.net_zero_count} | {int(point.is_continuum)} |"
        )

    lines.extend([
        "",
        "## a1(mu1) sweep",
        "",
        f"mu1 range: `[{args.mu1_low}, {args.mu1_high}]`, samples: `{args.samples}`, N_points: `{args.n_points}`",
        f"Minimum adjacent delta: `{min_delta}`",
        f"Monotonic decreases below tolerance {args.monotonic_tol}: `{len(decreases)}`",
    ])

    if decreases:
        lines.extend([
            "",
            "First decreasing intervals:",
            "",
            "| mu1_left | mu1_right | a1_left | a1_right | delta |",
            "|---:|---:|---:|---:|---:|",
        ])
        for row in decreases[:10]:
            lines.append("| " + " | ".join(f"{value:.12g}" for value in row) + " |")

    lines.extend([
        "",
        "## Zero plateau test",
        "",
        f"Plateau predicate: `abs(a1) <= {args.winding_tol}` and `zero_count == 0`.",
    ])

    if plateau_intervals:
        lines.append("")
        lines.append("| sampled_start | sampled_end | refined_start | refined_end | width |")
        lines.append("|---:|---:|---:|---:|---:|")
        for (start, end), (left, right) in zip(plateau_intervals, refined_intervals):
            lines.append(
                f"| {points[start].mu1:.12g} | {points[end].mu1:.12g} | "
                f"{left:.12g} | {right:.12g} | {right - left:.12g} |"
            )
    else:
        lines.append("")
        lines.append("No zero plateau was found by this sampled test.")

    lines.extend([
        "",
        "## Diagnostic conclusion",
        "",
    ])
    if decreases:
        lines.append(
            "The measured a1(mu1) curve has decreasing intervals, which points to a winding/root-tracking bug."
        )
    elif plateau_intervals and check_result.get("is_amoeba"):
        lines.append(
            "No monotonicity violation was found, but a finite a1=0 plateau with zero crossings absent was found. "
            "Because check_amoeba still returned is_amoeba=True, this supports the plateau-edge/search-classification bug hypothesis."
        )
    elif plateau_intervals:
        lines.append(
            "No monotonicity violation was found, and the zero plateau criterion indicates the energy is outside the amoeba spectrum."
        )
    else:
        lines.append(
            "No monotonicity violation and no zero plateau were found by this test, so the measurement supports the counter-intuitive in-spectrum conclusion."
        )

    lines.extend([
        "",
        "## Files",
        "",
        f"CSV data: `{csv_path.relative_to(ROOT)}`",
    ])

    report_path.write_text("\n".join(lines) + "\n")

    print(f"Wrote {report_path}")
    print(f"Wrote {csv_path}")
    print(format_check_result(check_result))
    print(f"decreases={len(decreases)}, plateau_intervals={len(plateau_intervals)}")


if __name__ == "__main__":
    main()
