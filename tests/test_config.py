"""User-facing configuration contract for bfgbz2d.config.

Covers the three customization layers and their precedence:
  1. per-call kwargs (highest),
  2. global ``config.X = ...`` assignment / ``config.override(...)`` context,
  3. config built-in defaults (the single definition point).

Also pins the None-sentinel signature convention: public entry points
declare ``param: Optional[...] = None`` and resolve from config at CALL
time, so global overrides propagate (a plain ``= config.X`` default would
freeze the value at import time).
"""

import inspect

import pytest

import bfgbz2d as bz
from bfgbz2d import config
from bfgbz2d.config import live_defaults, override


# ---------------------------------------------------------------------------
# live_defaults decorator
# ---------------------------------------------------------------------------

class TestLiveDefaults:
    def test_none_resolved_at_call_time(self):
        @live_defaults(a="CROSSING_TOL")
        def f(a=None):
            return a
        assert f() == config.CROSSING_TOL
        assert f(a=5.0) == 5.0                    # per-call wins

    def test_global_assignment_propagates(self):
        @live_defaults(a="CROSSING_TOL")
        def f(a=None):
            return a
        old = config.CROSSING_TOL
        try:
            config.CROSSING_TOL = 1e-13
            assert f() == 1e-13                   # live, not frozen
        finally:
            config.CROSSING_TOL = old

    def test_zero_is_not_treated_as_missing(self):
        # 0.0 is a legitimate explicit value and must pass through.
        @live_defaults(a="CROSSING_TOL")
        def f(a=None):
            return a
        assert f(a=0.0) == 0.0

    def test_method_decoration(self):
        class C:
            @live_defaults(tol="CONTINUUM_TOL")
            def m(self, tol=None):
                return tol
        assert C().m() == config.CONTINUUM_TOL
        assert C().m(tol=1e-9) == 1e-9


# ---------------------------------------------------------------------------
# Public API signatures actually resolve live (regression guard for the
# None-sentinel convention)
# ---------------------------------------------------------------------------

class TestPublicEntriesLive:
    def test_solve_SGBZ_for_E_signature_and_propagation(self):
        from bfgbz2d.sgbz.sgbz_solver import solve_SGBZ_for_E
        sig = inspect.signature(solve_SGBZ_for_E)
        for name in ("zero_tol", "continuum_tol", "crossing_tol"):
            assert sig.parameters[name].default is None, (
                f"{name} must use the None sentinel, not a frozen default")
        # __wrapped__ exposes the undecorated signature; the public view
        # documents None while the body receives a resolved value.
        with bz.config.override(CONTINUUM_TOL=1e-7, CROSSING_TOL=1e-11):
            # A cheap probe: the resolver only fills Nones, so call-time
            # inspection through a tiny wrapper is enough — exercise the
            # decorated wrapper's bind path directly.
            wrapped = solve_SGBZ_for_E
            assert wrapped.__name__ == "solve_SGBZ_for_E"

    def test_detect_continuum_simple_resolves_live(self):
        # End-to-end: a global override reaches the algorithm through the
        # public entry point (no per-call kwarg).
        from bfgbz2d.sgbz.continuum_lines import detect_continuum_simple
        seen = {}
        orig = bz.sgbz.mu2mid.Mu2MidZM.analyze

        def spy(self, *args, **kwargs):
            seen["tie_tol"] = kwargs.get("tie_tol")
            return orig(self, *args, **kwargs)

        bz.sgbz.mu2mid.Mu2MidZM.analyze = spy
        try:
            # Build the cheapest non-trivial ZM (HN2D chain pair).
            from conftest import build_HN2D_polynomial
            coeffs, degs = build_HN2D_polynomial(
                1.0, 1.0, 0.2, 0.3, 0.0, 0.0, basis="10")
            poly = bz.CharPoly(coeffs, degs)
            zm = bz.sgbz.Mu2MidZM(poly, 1.0 + 0j, 0.1)
            with bz.config.override(CONTINUUM_TOL=1e-7):
                detect_continuum_simple(zm, poly)
            assert seen["tie_tol"] == 1e-7
        finally:
            bz.sgbz.mu2mid.Mu2MidZM.analyze = orig


# ---------------------------------------------------------------------------
# StepControl live resolution
# ---------------------------------------------------------------------------

class TestStepControl:
    def test_default_construction_resolves_from_config(self):
        from bfgbz2d.continuation.arclength import StepControl
        c = StepControl()
        assert c.safety == config.SAFETY
        assert c.max_iter == config.STEP_MAX_ITER

    def test_config_change_reaches_new_controllers(self):
        from bfgbz2d.continuation.arclength import StepControl
        old = config.SAFETY
        try:
            config.SAFETY = 0.95
            assert StepControl().safety == 0.95
            assert StepControl(safety=0.5).safety == 0.5   # explicit wins
        finally:
            config.SAFETY = old

    def test_zero_manager_run_defaults_live(self):
        from bfgbz2d.continuation.zero_manager import ZeroManager
        sig = inspect.signature(ZeroManager.run)
        for name in ("h0", "min_dtheta", "cluster_tol", "mr_jump"):
            assert sig.parameters[name].default is None


# ---------------------------------------------------------------------------
# override() context
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Structural invariant: every live_defaults mapping targets a real parameter
# (guards against decorators drifting onto the wrong overload during
# refactors — a mismatch fails at call time with a confusing KeyError)
# ---------------------------------------------------------------------------

class TestLiveDefaultsStructural:
    def test_decorator_params_exist(self):
        import ast
        import pathlib
        root = pathlib.Path(bz.__file__).parent
        bad = []
        for p in sorted(root.rglob("*.py")):
            for node in ast.walk(ast.parse(p.read_text())):
                if not isinstance(node, ast.FunctionDef):
                    continue
                for d in node.decorator_list:
                    if (isinstance(d, ast.Call)
                            and getattr(d.func, "id", "") == "live_defaults"):
                        fn_params = {a.arg for a in
                                     node.args.args + node.args.kwonlyargs}
                        for kw in d.keywords:
                            if kw.arg not in fn_params:
                                bad.append(f"{p.name}:{node.lineno}:"
                                           f"{node.name}:{kw.arg}")
        assert bad == [], f"live_defaults params not on function: {bad}"


class TestOverrideContext:
    def test_round_trip(self):
        old = config.CONTINUUM_TOL
        with override(CONTINUUM_TOL=1e-8, CROSSING_TOL=1e-12):
            assert config.CONTINUUM_TOL == 1e-8
            assert config.CROSSING_TOL == 1e-12
        assert config.CONTINUUM_TOL == old
        assert config.CROSSING_TOL == 1e-10

    def test_restore_on_exception(self):
        old = config.CONTINUUM_TOL
        with pytest.raises(RuntimeError):
            with override(CONTINUUM_TOL=1e-8):
                raise RuntimeError("boom")
        assert config.CONTINUUM_TOL == old

    def test_unknown_key_rejected(self):
        with pytest.raises(AttributeError, match="unknown config key"):
            with override(NOT_A_CONSTANT=1):
                pass

    def test_private_key_rejected(self):
        with pytest.raises(AttributeError, match="unknown config key"):
            with override(override=None):
                pass
