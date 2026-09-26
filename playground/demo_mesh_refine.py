# Copyright © Department of Physics, Tsinghua University. All rights reserved.
# Numerical routines extracted from the playground prototypes by wangchenyang.

"""Compatibility imports for the former mesh-refinement prototype.

The numerical implementation now lives in pygbz2d.experimental; application
scripts should import the package directly.
"""
from pygbz2d.experimental.mesh_refinement import (
    MATCH_TOL, BISECT_ITERS, build_E_interpolator, chord4,
    match_result_to_target, dyadic_geometry_refine, bisect_to_boundary,
    solve_batch,
)
from pygbz2d.experimental.torus_mesh import (
    torus_midpoint, torus_pair_dist, edge_lengths, triangle_areas,
    edge_len_percentiles,
)
