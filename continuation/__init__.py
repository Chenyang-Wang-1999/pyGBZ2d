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
    detect_cluster, refine_multiple_root_theta, MultipleRootInfo

Segment integration (stops at end or MR):
    integrate_segment, SegmentResult

ZeroManager (top-level orchestrator):
    ZeroManager, SegmentData
"""

from continuation.arclength import (
    compute_tangent,
    predict_roots,
    estimate_error,
    arclength_step,
    detect_possible_multiple_root,
    detect_cluster,
    refine_multiple_root_theta,
    integrate_segment,
    StepResult,
    MultipleRootInfo,
    SegmentResult,
    SAFETY,
    MIN_FACTOR,
    MAX_FACTOR,
    ERROR_EXPONENT,
)

from continuation.zero_manager import ZeroManager, SegmentData

__all__ = [
    "compute_tangent",
    "predict_roots",
    "estimate_error",
    "arclength_step",
    "detect_possible_multiple_root",
    "detect_cluster",
    "refine_multiple_root_theta",
    "integrate_segment",
    "StepResult",
    "MultipleRootInfo",
    "SegmentResult",
    "ZeroManager",
    "SegmentData",
    "SAFETY",
    "MIN_FACTOR",
    "MAX_FACTOR",
    "ERROR_EXPONENT",
]
