"""
Visualization module for AMR Thesis Project.

Exports:
- create_resistance_heatmap: Basic resistance heatmap
- create_dendrogram: Hierarchical clustering dendrogram
- create_clustered_heatmap_with_dendrogram: Dendrogram-linked heatmap with annotations
- create_cluster_profile_heatmap: Mean resistance profile per cluster
- create_mdr_distribution_plot: MDR distribution visualization
- create_mar_distribution_plot: MAR index distribution visualization
- generate_all_visualizations: Generate all required visualizations
"""

from .visualization import (
    create_resistance_heatmap,
    create_dendrogram,
    create_clustered_heatmap_with_dendrogram,
    create_cluster_profile_heatmap,
    create_mdr_distribution_plot,
    create_mar_distribution_plot,
    generate_all_visualizations,
    RESISTANCE_COLORS,
    RESISTANCE_CMAP,
    CLUSTER_COLORS,
    MDR_COLORS,
)
