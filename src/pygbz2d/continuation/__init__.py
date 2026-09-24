"""
continuation — Adaptive-step root tracking (pseudo-arclength parameterization).

Provides adaptive-step root tracking along theta1 using analytic derivatives
from CharPoly and Hungarian matching.  The step size h is chosen by an
RK45-style *error* controller; the pseudo-arclength parameterization then maps
it to dtheta1 = h / ||V|| with V = [1, d ln beta2_j/dtheta1].

The two are independent.  The parameterization self-refines only where
||V|| grows (beta2 -> 0 or infinity, or the tangent diverges at an exact
branch point); a close approach elsewhere in the beta2 plane is resolved by
the error controller alone.  See doc/continuation.md §2.1.1.

Public API
----------
Low-level step:
    compute_tangent, predict_roots, estimate_error, arclength_step,
    StepResult, StepControl

Polynomial interpolation:
    cubic_hermite_poly, hermite_interp_poly

Multiple-root detection & refinement:
    multiple_root_point_trigger, MultipleRootIntervalTrigger,
    MRTriggerRecord, detect_cluster, solve_multiple_roots_in_interval,
    solve_multiple_roots_iterative, MultipleRootInfo

Segment integration (stops at end or MR):
    integrate_segment, SegmentResult, StopReason

ZeroManager (top-level orchestrator):
    ZeroManager, SegmentData
"""

from .arclength import (
    compute_tangent,
    predict_roots,
    estimate_error,
    arclength_step,
    StepResult,
    StepControl,
)

from .interpolation import (
    cubic_hermite_poly,
    hermite_interp_poly,
)

from .zero_manager import (
    ZeroManager,
    SegmentData,
    StopReason,
    integrate_segment,
    SegmentResult,
)

from .multiple_roots import (
    MultipleRootInfo,
    MRTriggerRecord,
    multiple_root_point_trigger,
    MultipleRootIntervalTrigger,
    detect_cluster,
    solve_multiple_roots_in_interval,
    solve_multiple_roots_iterative,
)

__all__ = [
    "compute_tangent",
    "predict_roots",
    "estimate_error",
    "arclength_step",
    "cubic_hermite_poly",
    "hermite_interp_poly",
    "multiple_root_point_trigger",
    "MultipleRootIntervalTrigger",
    "MRTriggerRecord",
    "detect_cluster",
    "solve_multiple_roots_in_interval",
    "solve_multiple_roots_iterative",
    "integrate_segment",
    "StepResult",
    "StepControl",
    "MultipleRootInfo",
    "SegmentResult",
    "StopReason",
    "ZeroManager",
    "SegmentData",
]
