'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-08-27
Copyright © Department of Physics, Tsinghua University. All rights reserved

EXPERIMENTAL band clustering and diagnostics for GBZ sweep data.
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

__all__ = [
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
