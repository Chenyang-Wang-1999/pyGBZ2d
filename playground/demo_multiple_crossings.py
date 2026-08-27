"""
Minimal replication of the two failed SGBZ cases.

Both cases fail with "Multiple crossings detected near θ₁=0.000000e+00"
at energies E = -1.56 and E = -1.406.
"""

import numpy as np
from pathlib import Path
import sys
import importlib.util
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from brute_force_SGBZ import collect_GBZ_subsets
from gbz_types import CharPoly

# Import from Haldane model (filename has hyphens, so use importlib)
import importlib.util
spec = importlib.util.spec_from_file_location(
    "haldane",
    Path(__file__).parent / "Haldane-model-gainloss.py"
)
haldane = importlib.util.module_from_spec(spec)
spec.loader.exec_module(haldane)
Haldane_non_Hermitian_phase = haldane.Haldane_non_Hermitian_phase
ALL_PARAMS = haldane.ALL_PARAMS


def test_failed_energy(E_ref: complex, direction: str = "a1"):
    """Test a single energy and report results."""
    print(f"\n{'='*70}")
    print(f"Testing E = {E_ref:.6f} (direction = {direction})")
    print(f"{'='*70}")

    # Build model
    model = Haldane_non_Hermitian_phase(*ALL_PARAMS)
    if direction == "a2":
        model = model.get_supercell(
            [(0,0)],
            np.array([[0, 1], [1, 0]], dtype=int)
        )
    elif direction == "x":
        model = model.get_supercell(
            [(0, 0), (1, 0)],
            np.array([[1, 1], [1, -1]], dtype=int)
        )

    coeffs, degs = model.get_characteristic_polynomial_data()
    poly = CharPoly(coeffs, degs)

    # Test with collect_GBZ_subsets
    print("\n1. collect_GBZ_subsets result:")
    result = collect_GBZ_subsets(coeffs, degs, E_ref)
    print(f"   success: {result.success}")
    print(f"   is_gbz: {result.is_gbz}")
    print(f"   is_empty: {result.is_empty}")
    if hasattr(result, 'error') and result.error:
        print(f"   error: {result.error}")

    # Try direct solve to see more details
    print("\n2. Direct solve_SGBZ_for_E:")
    try:
        from brute_force_SGBZ.sgbz_solver import solve_SGBZ_for_E
        solution = solve_SGBZ_for_E(poly, E_ref)
        print(f"   SUCCESS (should not reach here)")
        print(f"   mu1: {solution.get('mu1')}")
        print(f"   winding: {solution.get('winding')}")
    except Exception as e:
        print(f"   Exception: {type(e).__name__}: {e}")

        # If it's the multiple crossings error, investigate further
        if "Multiple crossings" in str(e):
            print("\n   --> This is the multiple crossings issue")
            print("   --> See traceback above for exact location in crossings.py")

            # Show where the error originates
            print(f"\n   Error originates from:")
            print(f"     crossings.py:412 in _check_no_duplicates()")
            print(f"     This function checks for duplicate crossings")
            print(f"     at the same theta1 value")

    return result


def main():
    print("="*70)
    print("MINIMAL REPLICATION OF FAILED SGBZ CASES")
    print("="*70)
    print("\nTwo energies fail with 'Multiple crossings detected near θ₁=0':")
    print("  E = -1.56")
    print("  E = -1.406")

    # Test both failed energies for all three directions
    failed_energies = [-1.56, -1.406]
    directions = ["a1", "a2", "x"]

    for E in failed_energies:
        for direction in directions:
            test_failed_energy(E, direction)

    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print("\nAll 6 cases (2 energies × 3 directions) should fail with")
    print("the same 'Multiple crossings' error at theta1=0.")
    print("\nNext steps: Investigate why these energies have degenerate")
    print("crossings at theta1=0 (θ₁=0 corresponds to beta1=1 on the unit circle).")


if __name__ == "__main__":
    main()