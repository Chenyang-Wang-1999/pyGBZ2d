"""
continuation — Pseudo-arclength continuation for root tracking.

Provides adaptive-step root tracking along theta1 using analytic derivatives
from CharPoly, Hungarian matching, and pseudo-arclength step-size control
modelled after scipy's RK45 integrator.

Public API
----------
Low-level step:
    compute_tangent, predict_roots, estimate_error, arclength_step, StepResult

Multiple-root detection & refinement:
    multiple_root_point_trigger, MultipleRootIntervalTrigger,
    detect_cluster, solve_multiple_roots_in_interval,
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
    SAFETY,
    MIN_FACTOR,
    MAX_FACTOR,
    ERROR_EXPONENT,
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
    "multiple_root_point_trigger",
    "MultipleRootIntervalTrigger",
    "detect_cluster",
    "solve_multiple_roots_in_interval",
    "solve_multiple_roots_iterative",
    "integrate_segment",
    "StepResult",
    "MultipleRootInfo",
    "SegmentResult",
    "StopReason",
    "ZeroManager",
    "SegmentData",
    "SAFETY",
    "MIN_FACTOR",
    "MAX_FACTOR",
    "ERROR_EXPONENT",
]
