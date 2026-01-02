"""
Clustering module for AMR Thesis Project.

Exports:
- Clustering parameters (LINKAGE_METHOD, DISTANCE_METRIC_PRIMARY, etc.)
- run_clustering_pipeline: Main clustering pipeline
- get_cluster_summary: Get cluster summary statistics
- create_cluster_summary_table: Create mandatory cluster summary table
- create_environmental_distribution_table: Create environmental distribution table
- perform_robustness_check: Robustness analysis with alternative distance metrics
- save_clustering_artifacts: Save clustering objects for reproducibility
"""

from .hierarchical_clustering import (
    LINKAGE_METHOD,
    DISTANCE_METRIC_PRIMARY,
    DISTANCE_METRIC_ROBUSTNESS,
    DEFAULT_N_CLUSTERS,
    CLUSTER_CUT_CRITERION,
    run_clustering_pipeline,
    get_cluster_summary,
    get_cluster_profiles,
    create_cluster_summary_table,
    create_environmental_distribution_table,
    create_regional_distribution_table,
    perform_robustness_check,
    save_clustering_artifacts,
    assign_clusters,
    determine_optimal_clusters,
)
