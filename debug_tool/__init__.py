'''
debug_tool — fixed-(E_ref, mu1) GBZ / winding debugging utilities.

Entry points (see ``gbz_debug`` module docstring for the full picture):

  * :func:`collect_debug_subsets` — GBZ subsets of BOTH methods (amoeba +
    SGBZ) at a frozen (E_ref, mu1), with the SGBZ PointSubset charges.
  * :func:`compute_loop_windings` — loop winding numbers at given theta2
    values + SGBZ charge-consistency gap checks.
  * :func:`plot_winding_debug` — subsets, charges and winding loops on the
    theta1-theta2 torus.

Run ``python debug_tool/demo_debug_tool.py`` for a ready-made session.
'''

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
]
