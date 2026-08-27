'''
Central registry of every numerical constant in bfGBZ2d.

Single definition point + live lookup.  Modules read values as
``config.X`` attribute accesses at *use* time, so a global override
takes effect on the next call without re-importing anything::

    import bfgbz2d as bz
    bz.config.CONTINUUM_TOL = 1e-8            # global default
    with bz.config.override(CROSSING_TOL=1e-12):   # temporary
        result = bz.sgbz.collect_GBZ_subsets(poly, E)

Per-call keyword arguments (where available) always take precedence over
the values here.

Grouping encodes how safe a constant is to tune:

  1. Model/algorithm scale — the knobs users actually want (tolerances
     expressed in physical units: μ-band widths, ln|β₂| distances, θ).
  2. Step-size & budget — convergence-speed knobs; wrong values trade
     accuracy against cost but keep results correct.
  3. Machine-precision-anchored — guard tolerances tied to float64
     resolution.  Do NOT retune these casually; they exist to make
     floating-point coincidence tests deterministic, and raising them
     can silently change discrete decisions.
'''

from contextlib import contextmanager
import functools
import inspect
from typing import Optional

# ---------------------------------------------------------------------------
# 1. Model / algorithm scale (user-facing)
# ---------------------------------------------------------------------------

#: Width of the continuum (degenerate-band) tolerance in μ-space.  SGBZ
#: tie detection and amoeba continuum detection share this single value
#: (both modules previously defined their own copy).
CONTINUUM_TOL: float = 1e-6

#: Fraction of in-band mesh rows above which an ItemView counts as a
#: continuum cluster.
CONTINUUM_FRAC: float = 0.9

#: brentq xtol when refining a pairwise ln|β₂| crossing.
CROSSING_TOL: float = 1e-10

#: |a₁| below which the SGBZ zero-plateau predicate fires.
ZERO_TOL: float = 1e-10

#: μ-perturbation used to escape a continuum band when probing the two
#: winding limits (SGBZ plateau probe perturbs μ₁, amoeba bisection μ₂).
#: The escape ladders (×1,2,4,8) handle wider bands.
CONTINUUM_PERTURB: float = 1e-4

#: Amoeba seam-snap radius (θ₁ scale).
SNAP_TOL: float = 1e-3

#: Amoeba root-deduplication radius (β₂-value scale).
ROOT_TOL: float = 1e-9

#: Default cluster tolerance for ZeroManager.run (β₂ chordal scale).
CLUSTER_TOL: float = 1e-4

#: Default cluster tolerance for detect_cluster in multiple-root search.
MR_CLUSTER_TOL: float = 1e-6

#: Amoeba plateau-probe cluster tolerance.
PLATEAU_CLUSTER_TOL: float = 1e-2

#: Amoeba plateau-probe area threshold (non-zero winding area).
PLATEAU_AREA_THRESHOLD: float = 1e-2

#: Minimum |Re(V_a) − Re(V_b)| below which a pair crossing is treated as
#: a tangent touch (crossing direction not trustworthy).
MIN_DIRECTION_DERIV: float = 1e-12

#: Clamp band for ln|β₂| in the μ₂_mid path pieces; |β₂| = e^±14 ≈ 1.2e6.
LOGABS_CLAMP: float = 14.0

#: Gauge tolerance for MR cluster consistency checks.
MR_GAUGE_TOL: float = 0.1

#: scipy.integrate.quad settings for the loop-winding integral (the
#: result is rounded to an integer, so these stay loose).
WINDING_QUAD_EPSABS: float = 1e-3
WINDING_QUAD_EPSREL: float = 1e-3
WINDING_QUAD_LIMIT: int = 200

#: Floor for the amoeba winding tolerance (never tighter than this).
WINDING_TOL_FLOOR: float = 1e-10

# ---------------------------------------------------------------------------
# 2. Step-size & budget (advanced)
# ---------------------------------------------------------------------------

# RK45-style predictive-controller knobs for the pseudo-arclength stepper.
SAFETY: float = 0.9
MIN_FACTOR: float = 0.2
MAX_FACTOR: float = 10.0
#: -1/(p+1) for the first-order tangent predictor (p=1).
ERROR_EXPONENT: float = -0.5

#: |β₂| below which a root is treated as a singular 0-padding root.
ZERO_THRESHOLD: float = 1e-6
#: |β₂| above which a root is treated as a singular ∞-padding root.
INF_THRESHOLD: float = 1e6

#: Step acceptance tolerances (StepControl).
STEP_ATOL: float = 1e-12
STEP_RTOL: float = 1e-3
STEP_MAX_ITER: int = 20

#: Arclength step bounds (StepControl).
MAX_STEP: float = 0.5
MIN_STEP: float = 1e-12

#: Default initial θ₁ step and floor for ZeroManager.run.
H0: float = 0.1
MIN_DTHETA: float = 1e-10

# Multi-crossing mesh refinement (pairwise pre-scan).
REFINE_MAX_ROUNDS: int = 3
REFINE_SAFETY_FACTOR: float = 4.0
REFINE_MAX_SUBINTERVALS: int = 64
REFINE_MAX_TOTAL_INSERTS: int = 2000

#: brentq iteration budget (the bracket is guaranteed by the sign scan).
BRENTQ_MAXITER: int = 100

#: Cap on μ₁ bracket-expansion steps in the SGBZ bisection.
MAX_BRACKET_EXPANSIONS: int = 10

# Amoeba bisection budget.
AMOEBA_MAX_ITER: int = 60
AMOEBA_XTOL: float = 1e-10
AMOEBA_MAX_RANGE_EXPANSIONS: int = 10
AMOEBA_RANGE_EXPAND_FACTOR: float = 2.0

# fsolve refinement of amoeba (θ₁, θ₂) crossings.
FSOLVE_XTOL: float = 1e-12
FSOLVE_MAXFEV: int = 500
#: Residual gate that accepts a refined crossing.
FSOLVE_RESIDUAL_TOL: float = 1e-10

# MR restart geometry: restart distance past a multiple root, with the
# branch-point-safe floors (a closer restart re-detects the same MR).
MR_JUMP: float = 1e-6
MR_RESTART_FACTOR_H0: float = 10.0
MR_RESTART_FACTOR_ABS: float = 100.0

# ---------------------------------------------------------------------------
# 3. Machine-precision-anchored guards (do not retune casually)
# ---------------------------------------------------------------------------

#: Exact-endpoint float comparison when reading mesh rows in refinement.
THETA_EQ_TOL: float = 1e-15

#: Real-root / duplicate-θ filter, in units of max(1, interval length).
REFINE_REL_TOL: float = 1e-12

#: Forward-progress guard for MR refinement (θ must advance past this).
MR_STUCK_TOL: float = 1e-12

#: |θ − 2π| below which a boundary MR is pinned to exactly 2π.
BOUNDARY_THETA_TOL: float = 1e-6

#: Hard cap on the number of tracked segments (runaway protection).
MAX_SEGMENTS: int = 10000


def live_defaults(**param_to_key):
    """Decorator: resolve ``None`` parameters from config at CALL time.

    Apply to public entry points whose keyword defaults would otherwise
    freeze the config value at import time (Python evaluates ``def``
    defaults once).  Usage::

        @live_defaults(continuum_tol="CONTINUUM_TOL")
        def detect_continuum_simple(zm, continuum_tol: Optional[float] = None):
            ...   # continuum_tol is never None inside

    Resolution order stays: per-call kwarg > config value (live) — the
    import-time built-in default no longer exists as a third tier, since
    config IS the single definition point.  Intended for top-level APIs
    (called ~once per solve), not per-step hot loops; those carry their
    knobs explicitly via ``StepControl``.
    """
    def decorator(fn):
        sig = inspect.signature(fn)

        @functools.wraps(fn)
        def wrapper(*args, **kwargs):
            bound = sig.bind(*args, **kwargs)
            bound.apply_defaults()
            for pname, ckey in param_to_key.items():
                if bound.arguments[pname] is None:
                    bound.arguments[pname] = globals()[ckey]
            return fn(*bound.args, **bound.kwargs)
        return wrapper
    return decorator


@contextmanager
def override(**kwargs):
    """Temporarily override config values, restored on exit.

    >>> import bfgbz2d as bz
    >>> with bz.config.override(CONTINUUM_TOL=1e-8):
    ...     ...   # computations here see the tighter tolerance
    """
    unknown = [k for k in kwargs
               if k.startswith("_")
               or k not in globals()
               or not isinstance(globals()[k], (int, float))]
    if unknown:
        raise AttributeError(
            f"unknown config key(s): {unknown}; see bfgbz2d.config for "
            f"the defined names"
        )
    saved = {k: globals()[k] for k in kwargs}
    globals().update(kwargs)
    try:
        yield
    finally:
        globals().update(saved)
