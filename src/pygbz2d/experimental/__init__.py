'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-08-27
Copyright © Department of Physics, Tsinghua University. All rights reserved

EXPERIMENTAL band clustering, torus meshes, refinement and Chern integration.
Research-grade: APIs may change without
notice, correctness is not guaranteed, and nothing here is used by the
solver packages.  Import explicitly:

    from pygbz2d.experimental import cluster_bands
'''

from .band_clustering import (
    ALPHA_E,
    EPS_SCAN_MAX,
    EPS_SCAN_MIN,
    N_SCAN_STEPS,
    KNN_K,
    MIN_CLUSTER_SIZE,
    MIN_CLUSTERS_PLATEAU,
    MAX_MARGIN_PAIRS,
    BandPoints,
    BandClustering,
    flatten_results,
    grid_steps,
    embed,
    radius_graph_labels,
    knn_distance_stats,
    eps_stability_scan,
    widest_plateau,
    inter_cluster_margins,
    summarize_clusters,
    cluster_bands,
)
from .torus_mesh import (
    build_band_mesh,
    build_cluster_mesh,
    periodic_delaunay,
    mesh_topology,
    edge_lengths,
    triangle_areas,
    regular_torus_mesh,
)
from .mesh_refinement import refine_mesh
from .chern import ChernResult, VertexDiagnostics, integrate_chern

__all__ = [
    "build_band_mesh",
    "build_cluster_mesh",
    "periodic_delaunay",
    "mesh_topology",
    "edge_lengths",
    "triangle_areas",
    "regular_torus_mesh",
    "refine_mesh",
    "ChernResult",
    "VertexDiagnostics",
    "integrate_chern",
    "ALPHA_E",
    "EPS_SCAN_MAX",
    "EPS_SCAN_MIN",
    "N_SCAN_STEPS",
    "KNN_K",
    "MIN_CLUSTER_SIZE",
    "MIN_CLUSTERS_PLATEAU",
    "MAX_MARGIN_PAIRS",
    "BandPoints",
    "BandClustering",
    "flatten_results",
    "grid_steps",
    "embed",
    "radius_graph_labels",
    "knn_distance_stats",
    "eps_stability_scan",
    "widest_plateau",
    "inter_cluster_margins",
    "summarize_clusters",
    "cluster_bands",
]
