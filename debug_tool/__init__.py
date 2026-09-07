'''
debug_tool — fixed-(E_ref, mu1) GBZ / winding debugging utilities.

Entry points (see ``gbz_debug`` module docstring for the full picture):

  * :func:`collect_debug_subsets` — GBZ subsets of BOTH methods (amoeba +
    SGBZ) at a frozen (E_ref, mu1), with the SGBZ PointSubset charges.
  * :func:`compute_loop_windings` — loop winding numbers at given theta2
    values + SGBZ charge-consistency gap checks.
  * :func:`plot_winding_debug` — subsets, charges and winding loops on the
    theta1-theta2 torus.
  * :func:`check_mesh_orientation` — torus triangle-mesh orientation
    diagnostics (signed-area census + directed-edge balance).

Run ``python debug_tool/demo_debug_tool.py`` for a ready-made session.
'''

# These repository utilities must use the same source as the checkout tests,
# including when invoked with `python -m debug_tool.demo_debug_tool`.
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from .gbz_debug import (
    MethodDebug,
    GBZDebugReport,
    LoopWinding,
    GapCheck,
    MethodLoops,
    collect_debug_subsets,
    auto_theta2_grid,
    compute_loop_windings,
    plot_winding_debug,
)
from .mesh_orientation import MeshOrientationReport, check_mesh_orientation

__all__ = [
    "MethodDebug",
    "GBZDebugReport",
    "LoopWinding",
    "GapCheck",
    "MethodLoops",
    "collect_debug_subsets",
    "auto_theta2_grid",
    "compute_loop_windings",
    "plot_winding_debug",
    "MeshOrientationReport",
    "check_mesh_orientation",
]
