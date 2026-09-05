"""Constant-management contract tests (module-constant architecture).

The customization model has exactly two layers:
  1. per-call kwargs (highest precedence, resolved from home-module
     constants via ``live_defaults`` when left as None);
  2. direct assignment of a home-module constant — takes effect
     process-wide on the NEXT read (module attribute lookup).

There is no third channel: no central config module, no options dict,
no override() context manager.
"""

import ast
import importlib
import inspect
import pathlib

import pytest

import pygbz2d as bz
from pygbz2d.core import live_defaults

# Home-module map: constant -> module path (single source of truth for
# the AST lint rules below and for doc/constants.md).
HOME = {
    "CONTINUUM_TOL": "core", "CONTINUUM_FRAC": "core",
    "CONTINUUM_PERTURB": "core", "WINDING_ZERO_TOL": "core",
    "PLATEAU_CLUSTER_TOL": "core", "ESCAPE_LADDER": "core",
    "PROBE_XTOL": "core",
    "SAFETY": "continuation.arclength", "MIN_FACTOR": "continuation.arclength",
    "MAX_FACTOR": "continuation.arclength", "ERROR_EXPONENT": "continuation.arclength",
    "STEP_ATOL": "continuation.arclength", "STEP_RTOL": "continuation.arclength",
    "STEP_MAX_ITER": "continuation.arclength", "MAX_STEP": "continuation.arclength",
    "MIN_STEP": "continuation.arclength", "ZERO_THRESHOLD": "continuation.arclength",
    "INF_THRESHOLD": "continuation.arclength",
    "PREDICT_MAX_ABS_ARG": "continuation.arclength",
    "CLUSTER_TOL": "continuation.multiple_roots",
    "MIN_DIST_THRESHOLD": "continuation.multiple_roots",
    "H0": "continuation.zero_manager", "MIN_DTHETA": "continuation.zero_manager",
    "MR_JUMP": "continuation.zero_manager",
    "MR_RESTART_FACTOR_H0": "continuation.zero_manager",
    "MR_RESTART_FACTOR_ABS": "continuation.zero_manager",
    "MR_STUCK_TOL": "continuation.zero_manager",
    "BOUNDARY_THETA_TOL": "continuation.zero_manager",
    "MR_GAUGE_TOL": "continuation.zero_manager",
    "MAX_SEGMENTS": "continuation.zero_manager",
    "CROSSING_TOL": "sgbz.pairwise", "MIN_DIRECTION_DERIV": "sgbz.pairwise",
    "REFINE_MAX_ROUNDS": "sgbz.pairwise", "REFINE_SAFETY_FACTOR": "sgbz.pairwise",
    "REFINE_MAX_SUBINTERVALS": "sgbz.pairwise",
    "REFINE_MAX_TOTAL_INSERTS": "sgbz.pairwise", "REFINE_REL_TOL": "sgbz.pairwise",
    "THETA_EQ_TOL": "sgbz.pairwise", "BRENTQ_MAXITER": "sgbz.pairwise",
    "LOGABS_CLAMP": "sgbz.mu2mid",
    "WINDING_QUAD_EPSABS": "sgbz.winding", "WINDING_QUAD_EPSREL": "sgbz.winding",
    "WINDING_QUAD_LIMIT": "sgbz.winding", "SEED_N_PER_INTERVAL": "sgbz.winding",
    "MU1_MAX_ITER": "sgbz.sgbz_solver",
    "MAX_BRACKET_EXPANSIONS": "sgbz.sgbz_solver",
    "SNAP_TOL": "amoeba.amoeba",
    "BISECT_MAX_ITER": "amoeba.bisect", "BISECT_XTOL": "amoeba.bisect",
    "BISECT_COARSE_XTOL": "amoeba.bisect",
    "EXTREMUM_INSERT_REL_TOL": "amoeba.bisect",
    "MAX_RANGE_EXPANSIONS": "amoeba.bisect", "RANGE_EXPAND_FACTOR": "amoeba.bisect",
    "FSOLVE_XTOL": "amoeba.ronkin_winding", "FSOLVE_MAXFEV": "amoeba.ronkin_winding",
    "CROSSING_RESIDUAL_TOL": "amoeba.ronkin_winding",
    "PLATEAU_AREA_THRESHOLD": "amoeba.amoeba",
}

SRC = pathlib.Path(bz.__file__).parent


# ---------------------------------------------------------------------------
# live_defaults decorator mechanics
# ---------------------------------------------------------------------------

class TestLiveDefaults:
    def test_none_resolved_at_call_time(self):
        @live_defaults(a="core:PROBE_XTOL")
        def f(a=None):
            return a
        assert f() == bz.core.PROBE_XTOL
        assert f(a=5.0) == 5.0
        assert f(a=0.0) == 0.0          # 0.0 is a real value, not "missing"

    def test_module_assignment_is_live(self):
        from pygbz2d.sgbz import pairwise
        old = pairwise.CROSSING_TOL
        try:
            @live_defaults(a="sgbz.pairwise:CROSSING_TOL")
            def f(a=None):
                return a
            assert f() == old
            pairwise.CROSSING_TOL = 1e-13
            assert f() == 1e-13          # next call sees the assignment
        finally:
            pairwise.CROSSING_TOL = old

    def test_method_decoration(self):
        class C:
            @live_defaults(tol="core:CONTINUUM_TOL")
            def m(self, tol=None):
                return tol
        assert C().m() == bz.core.CONTINUUM_TOL
        assert C().m(tol=1e-9) == 1e-9

    def test_bad_key_format_raises(self):
        @live_defaults(a="no-colon-here")
        def f(a=None):
            return a
        with pytest.raises(ValueError):
            f()


# ---------------------------------------------------------------------------
# StepControl resolves from its home module at construction
# ---------------------------------------------------------------------------

class TestStepControl:
    def test_default_construction(self):
        from pygbz2d.continuation import arclength
        c = arclength.StepControl()
        assert c.safety == arclength.SAFETY
        assert c.max_iter == arclength.STEP_MAX_ITER

    def test_assignment_reaches_new_controllers(self):
        from pygbz2d.continuation import arclength
        old = arclength.SAFETY
        try:
            arclength.SAFETY = 0.95
            assert arclength.StepControl().safety == 0.95
            assert arclength.StepControl(safety=0.5).safety == 0.5
        finally:
            arclength.SAFETY = old

    def test_run_defaults_are_none_sentinels(self):
        from pygbz2d.continuation.zero_manager import ZeroManager
        sig = inspect.signature(ZeroManager.run)
        for name in ("h0", "min_dtheta", "cluster_tol", "mr_jump"):
            assert sig.parameters[name].default is None


# ---------------------------------------------------------------------------
# End-to-end: assignment propagates through a public entry point
# ---------------------------------------------------------------------------

class TestEndToEnd:
    def test_continuum_tol_assignment_reaches_analyze(self):
        from pygbz2d.sgbz import mu2mid
        from pygbz2d.sgbz.continuum_lines import detect_continuum_simple
        from conftest import build_HN2D_polynomial

        seen = {}
        orig = mu2mid.Mu2MidZM.analyze

        def spy(self, *args, **kwargs):
            seen["tie_tol"] = kwargs.get("tie_tol")
            return orig(self, *args, **kwargs)

        coeffs, degs = build_HN2D_polynomial(
            1.0, 1.0, 0.2, 0.3, 0.0, 0.0, basis="10")
        poly = bz.CharPoly(coeffs, degs)
        zm = mu2mid.Mu2MidZM(poly, 1.0 + 0j, 0.1)
        mu2mid.Mu2MidZM.analyze = spy
        old = bz.core.CONTINUUM_TOL
        try:
            bz.core.CONTINUUM_TOL = 1e-7
            detect_continuum_simple(zm, poly)
            assert seen["tie_tol"] == 1e-7
        finally:
            bz.core.CONTINUUM_TOL = old
            mu2mid.Mu2MidZM.analyze = orig


# ---------------------------------------------------------------------------
# AST lint rule 1: no numeric signature defaults outside the whitelist
# ---------------------------------------------------------------------------

WHITELIST = {"direction", "mu1_guess", "mu2_low", "mu2_high", "perc"}


class TestNoMagicDefaults:
    def test_signatures_have_no_numeric_defaults(self):
        def const_num(d):
            if isinstance(d, ast.Constant) and isinstance(d.value, (int, float)):
                return [d.value]
            if isinstance(d, (ast.Tuple, ast.List)):
                return [e.value for e in d.elts
                        if isinstance(e, ast.Constant)
                        and isinstance(e.value, (int, float))]
            return []
        offenders = []
        for p in SRC.rglob("*.py"):
            for node in ast.walk(ast.parse(p.read_text(encoding="utf-8"))):
                if not isinstance(node, ast.FunctionDef):
                    continue
                pairs = []
                if node.args.defaults:
                    pairs += list(zip(node.args.args[-len(node.args.defaults):],
                                      node.args.defaults))
                pairs += [(a, d) for a, d in
                          zip(node.args.kwonlyargs, node.args.kw_defaults)
                          if d is not None]
                for a, d in pairs:
                    nums = const_num(d)
                    if (nums and any(n not in (0, 1) for n in nums)
                            and a.arg not in WHITELIST):
                        offenders.append(f"{p.name}:{node.lineno}:"
                                         f"{node.name}({a.arg}={nums})")
        assert offenders == [], (
            "numeric defaults must live as module constants (or join the "
            f"whitelist): {offenders}")


# ---------------------------------------------------------------------------
# AST lint rule 2: constants are never imported BY VALUE across modules
# (a from-import copies the value at import time and silently freezes it)
# ---------------------------------------------------------------------------

class TestNoConstantValueImports:
    def test_no_from_import_of_constants(self):
        offenders = []
        for p in SRC.rglob("*.py"):
            if p.name == "__init__.py":
                continue  # package re-exports are fine
            for node in ast.walk(ast.parse(p.read_text(encoding="utf-8"))):
                if isinstance(node, ast.ImportFrom):
                    for alias in node.names:
                        if alias.name in HOME:
                            offenders.append(
                                f"{p.name}:{node.lineno}: "
                                f"'from {node.module} import {alias.name}' "
                                f"— use module attribute access instead")
        assert offenders == [], offenders


# ---------------------------------------------------------------------------
# Structural: every live_defaults param exists and its key resolves
# ---------------------------------------------------------------------------

class TestDecoratorStructure:
    def test_params_exist_and_keys_resolve(self):
        bad = []
        for p in SRC.rglob("*.py"):
            for node in ast.walk(ast.parse(p.read_text(encoding="utf-8"))):
                if not isinstance(node, ast.FunctionDef):
                    continue
                for d in node.decorator_list:
                    if not (isinstance(d, ast.Call)
                            and getattr(d.func, "id", "") == "live_defaults"):
                        continue
                    fn_params = {a.arg for a in
                                 node.args.args + node.args.kwonlyargs}
                    for kw in d.keywords:
                        if kw.arg not in fn_params:
                            bad.append(f"{p.name}:{node.lineno}:"
                                       f"{node.name}:{kw.arg} not a param")
                        key = ast.literal_eval(kw.value)
                        try:
                            mod_name, const_name = key.split(":")
                            module = importlib.import_module(
                                "pygbz2d." + mod_name)
                            getattr(module, const_name)
                        except Exception as e:
                            bad.append(f"{p.name}:{node.lineno}:"
                                       f"{node.name}:{kw.arg} key {key!r}: {e}")
        assert bad == [], bad

    def test_home_map_is_accurate(self):
        for const, mod_path in HOME.items():
            module = importlib.import_module("pygbz2d." + mod_path)
            assert hasattr(module, const), (mod_path, const)
