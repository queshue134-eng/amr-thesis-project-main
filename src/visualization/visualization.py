"""
Visualization Module for AMR Thesis Project
Phase 3.2 - Heatmaps, Dendrograms, and Resistance Pattern Visualization

VISUALIZATION PRINCIPLES:
-------------------------
1. Dendrogram-linked heatmap: Uses the SAME dendrogram ordering for all visualizations
   to ensure visual patterns are consistent with hierarchical structure.

2. Cluster annotation: Color bars and vertical separators clearly identify cluster
   membership with consistent labeling (C1, C2, ...).

3. MDR signal enhancement: Side bars and annotations highlight MDR proportion per
   isolate and cluster-level MDR percentages.

4. Clusters as phenotypes: Labels describe clusters as "resistance phenotypes"
   NOT as taxonomic or causal groups.
"""

import pandas as pd
import numpy as np
from typing import List, Tuple, Dict, Optional
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import dendrogram
import os


# Color scheme for resistance values
RESISTANCE_COLORS = {
    0: '#4CAF50',  # Susceptible - Green
    1: '#FFC107',  # Intermediate - Yellow/Amber
    2: '#F44336',  # Resistant - Red
}

# Colormap for heatmap
RESISTANCE_CMAP = mcolors.ListedColormap(['#4CAF50', '#FFC107', '#F44336'])

# Color palette for clusters (distinct colors for up to 10 clusters)
CLUSTER_COLORS = [
    '#1f77b4',  # Blue
    '#ff7f0e',  # Orange
    '#2ca02c',  # Green
    '#d62728',  # Red
    '#9467bd',  # Purple
    '#8c564b',  # Brown
    '#e377c2',  # Pink
    '#7f7f7f',  # Gray
    '#bcbd22',  # Yellow-green
    '#17becf',  # Cyan
]

# MDR color scheme
MDR_COLORS = {
    0: '#4CAF50',  # Non-MDR - Green
    1: '#F44336',  # MDR - Red
}


def create_resistance_heatmap(df: pd.DataFrame,
                             feature_cols: List[str],
                             cluster_col: str = 'CLUSTER',
                             figsize: Tuple[int, int] = (14, 10),
                             save_path: str = None) -> plt.Figure:
    """
    Create a heatmap of resistance patterns.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe with encoded resistance values and cluster labels
    feature_cols : list
        List of feature column names
    cluster_col : str
        Column name for cluster labels
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save the figure
    
    Returns:
    --------
    matplotlib.Figure
        Heatmap figure
    """
    # Prepare data
    existing_cols = [c for c in feature_cols if c in df.columns]
    
    # Sort by cluster
    if cluster_col in df.columns:
        df_sorted = df.sort_values(cluster_col)
    else:
        df_sorted = df.copy()
    
    # Extract feature matrix
    data_matrix = df_sorted[existing_cols].values
    
    # Clean column names for display
    display_cols = [c.replace('_encoded', '') for c in existing_cols]
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create heatmap
    im = ax.imshow(data_matrix, aspect='auto', cmap=RESISTANCE_CMAP,
                   vmin=0, vmax=2)
    
    # Set labels
    ax.set_xticks(np.arange(len(display_cols)))
    ax.set_xticklabels(display_cols, rotation=45, ha='right')
    ax.set_xlabel('Antibiotics')
    ax.set_ylabel('Isolates')
    ax.set_title('Resistance Profile Heatmap')
    
    # Add cluster boundaries
    if cluster_col in df.columns:
        cluster_labels = df_sorted[cluster_col].values
        cluster_boundaries = np.where(np.diff(cluster_labels) != 0)[0]
        for boundary in cluster_boundaries:
            ax.axhline(y=boundary + 0.5, color='black', linewidth=1.5)
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, ticks=[0, 1, 2])
    cbar.set_label('Resistance Status')
    cbar.ax.set_yticklabels(['S (Susceptible)', 'I (Intermediate)', 'R (Resistant)'])
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Heatmap saved to: {save_path}")
    
    return fig


def create_dendrogram(linkage_matrix: np.ndarray,
                     labels: List[str] = None,
                     figsize: Tuple[int, int] = (12, 8),
                     color_threshold: float = None,
                     save_path: str = None) -> plt.Figure:
    """
    Create a dendrogram from hierarchical clustering.
    
    Parameters:
    -----------
    linkage_matrix : np.ndarray
        Linkage matrix from hierarchical clustering
    labels : list, optional
        Labels for leaf nodes
    figsize : tuple
        Figure size
    color_threshold : float, optional
        Threshold for coloring clusters
    save_path : str, optional
        Path to save the figure
    
    Returns:
    --------
    matplotlib.Figure
        Dendrogram figure
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create dendrogram
    dendro = dendrogram(
        linkage_matrix,
        ax=ax,
        labels=labels,
        color_threshold=color_threshold,
        leaf_rotation=90,
        leaf_font_size=8
    )
    
    ax.set_xlabel('Isolates')
    ax.set_ylabel('Distance')
    ax.set_title('Hierarchical Clustering Dendrogram')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Dendrogram saved to: {save_path}")
    
    return fig


def create_clustered_heatmap_with_dendrogram(df: pd.DataFrame,
                                              feature_cols: List[str],
                                              linkage_matrix: np.ndarray,
                                              figsize: Tuple[int, int] = (18, 14),
                                              save_path: str = None,
                                              show_mdr: bool = True,
                                              show_cluster_bars: bool = True) -> Tuple[plt.Figure, List[int]]:
    """
    Create a DENDROGRAM-LINKED HEATMAP with cluster annotations and MDR signal.
    
    CRITICAL REQUIREMENTS IMPLEMENTED:
    1. Uses dendrogram ordering for heatmap rows (no independent reordering)
    2. Explicit cluster annotation with color bars
    3. Vertical separators between clusters
    4. Consistent cluster labeling (C1, C2, ...)
    5. MDR proportion sidebar per isolate
    6. Cluster-level MDR % annotation
    
    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe with encoded resistance values and cluster labels
    feature_cols : list
        List of feature column names
    linkage_matrix : np.ndarray
        Linkage matrix from hierarchical clustering
    figsize : tuple
        Figure size (default: 18x14 for high-resolution)
    save_path : str, optional
        Path to save the figure
    show_mdr : bool
        Whether to show MDR signal enhancement (default: True)
    show_cluster_bars : bool
        Whether to show cluster color bars (default: True)
    
    Returns:
    --------
    tuple
        (Figure object, Dendrogram order list)
    """
    # Prepare data
    existing_cols = [c for c in feature_cols if c in df.columns]
    data_matrix = df[existing_cols].values
    display_cols = [c.replace('_encoded', '') for c in existing_cols]
    
    # Create figure with gridspec
    fig = plt.figure(figsize=figsize)
    
    # Define width ratios based on components
    # Dendrogram | Cluster bar | MDR bar | Heatmap
    width_ratios = [1.5]  # Dendrogram
    if show_cluster_bars and 'CLUSTER' in df.columns:
        width_ratios.append(0.3)  # Cluster bar
    if show_mdr and 'MDR_FLAG' in df.columns:
        width_ratios.append(0.3)  # MDR bar
    width_ratios.append(6)  # Heatmap
    
    gs = fig.add_gridspec(1, len(width_ratios), width_ratios=width_ratios, wspace=0.02)
    
    col_idx = 0
    
    # 1. Dendrogram subplot (leftmost)
    ax_dendro = fig.add_subplot(gs[col_idx])
    col_idx += 1
    
    # Create dendrogram (rotated) - this determines row order
    dendro = dendrogram(
        linkage_matrix,
        ax=ax_dendro,
        orientation='left',
        no_labels=True,
        color_threshold=0  # All same color for cleaner look
    )
    ax_dendro.set_xticks([])
    ax_dendro.set_yticks([])
    ax_dendro.invert_yaxis()
    ax_dendro.set_title('Dendrogram', fontsize=10)
    
    # Get dendrogram order - CRITICAL: this order is used for ALL visualizations
    dendro_order = dendro['leaves']
    
    # Reorder all data according to dendrogram
    data_ordered = data_matrix[dendro_order, :]
    
    # Get cluster labels in dendrogram order
    if 'CLUSTER' in df.columns:
        cluster_labels_ordered = df['CLUSTER'].values[dendro_order]
    
    # Get MDR flags in dendrogram order
    if 'MDR_FLAG' in df.columns:
        mdr_flags_ordered = df['MDR_FLAG'].values[dendro_order]
    
    # 2. Cluster color bar (optional)
    if show_cluster_bars and 'CLUSTER' in df.columns:
        ax_cluster = fig.add_subplot(gs[col_idx])
        col_idx += 1
        
        # Create cluster color array
        unique_clusters = sorted(df['CLUSTER'].unique())
        cluster_cmap = {c: CLUSTER_COLORS[i % len(CLUSTER_COLORS)] 
                       for i, c in enumerate(unique_clusters)}
        
        cluster_colors_array = np.array([
            mcolors.to_rgb(cluster_cmap[c]) for c in cluster_labels_ordered
        ])
        
        # Display as vertical bar
        ax_cluster.imshow(cluster_colors_array[:, np.newaxis, :], 
                         aspect='auto', interpolation='nearest')
        ax_cluster.set_xticks([])
        ax_cluster.set_yticks([])
        ax_cluster.set_title('Cluster', fontsize=10, rotation=0)
        
        # Add cluster separators (horizontal lines at cluster boundaries)
        prev_cluster = cluster_labels_ordered[0]
        for i, c in enumerate(cluster_labels_ordered):
            if c != prev_cluster:
                ax_cluster.axhline(y=i-0.5, color='white', linewidth=2)
                prev_cluster = c
    
    # 3. MDR color bar (optional)
    if show_mdr and 'MDR_FLAG' in df.columns:
        ax_mdr = fig.add_subplot(gs[col_idx])
        col_idx += 1
        
        # Create MDR color array
        mdr_colors_array = np.array([
            mcolors.to_rgb(MDR_COLORS.get(int(m), '#808080')) 
            for m in mdr_flags_ordered
        ])
        
        # Display as vertical bar
        ax_mdr.imshow(mdr_colors_array[:, np.newaxis, :], 
                     aspect='auto', interpolation='nearest')
        ax_mdr.set_xticks([])
        ax_mdr.set_yticks([])
        ax_mdr.set_title('MDR', fontsize=10, rotation=0)
    
    # 4. Heatmap subplot (rightmost)
    ax_heatmap = fig.add_subplot(gs[col_idx])
    
    im = ax_heatmap.imshow(data_ordered, aspect='auto', cmap=RESISTANCE_CMAP,
                           vmin=0, vmax=2)
    
    ax_heatmap.set_xticks(np.arange(len(display_cols)))
    ax_heatmap.set_xticklabels(display_cols, rotation=45, ha='right', fontsize=9)
    ax_heatmap.set_yticks([])
    ax_heatmap.set_xlabel('Antibiotics', fontsize=11)
    ax_heatmap.set_title('Resistance Profile Heatmap\n(Row order matches dendrogram)', fontsize=12)
    
    # Add cluster boundary lines on heatmap
    if 'CLUSTER' in df.columns:
        prev_cluster = cluster_labels_ordered[0]
        for i, c in enumerate(cluster_labels_ordered):
            if c != prev_cluster:
                ax_heatmap.axhline(y=i-0.5, color='black', linewidth=1.5, linestyle='-')
                prev_cluster = c
    
    # Add colorbar for resistance
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.3])
    cbar = plt.colorbar(im, cax=cbar_ax, ticks=[0, 1, 2])
    cbar.ax.set_yticklabels(['S (Susceptible)', 'I (Intermediate)', 'R (Resistant)'])
    cbar.set_label('Resistance Status', fontsize=10)
    
    # Create legend for clusters and MDR
    legend_elements = []
    
    if show_cluster_bars and 'CLUSTER' in df.columns:
        for c in unique_clusters:
            legend_elements.append(
                Patch(facecolor=cluster_cmap[c], edgecolor='black',
                      label=f'C{c}')
            )
    
    if show_mdr and 'MDR_FLAG' in df.columns:
        legend_elements.append(Patch(facecolor=MDR_COLORS[0], edgecolor='black', label='Non-MDR'))
        legend_elements.append(Patch(facecolor=MDR_COLORS[1], edgecolor='black', label='MDR'))
    
    if legend_elements:
        ax_heatmap.legend(handles=legend_elements, loc='upper left', 
                         bbox_to_anchor=(1.1, 1.0), fontsize=9)
    
    # Add cluster-level MDR annotation
    if 'CLUSTER' in df.columns and 'MDR_FLAG' in df.columns:
        # Calculate MDR% per cluster
        cluster_mdr = df.groupby('CLUSTER')['MDR_FLAG'].mean() * 100
        
        # Add text annotation at right side
        annotation_text = "Cluster MDR%:\n"
        for c in sorted(cluster_mdr.index):
            annotation_text += f"C{c}: {cluster_mdr[c]:.1f}%\n"
        
        fig.text(0.92, 0.55, annotation_text, fontsize=9, 
                verticalalignment='bottom', family='monospace',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.suptitle('Dendrogram-Linked Resistance Heatmap with Cluster & MDR Annotations',
                fontsize=14, fontweight='bold', y=0.98)
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Dendrogram-linked heatmap saved to: {save_path}")
    
    return fig, dendro_order


def create_cluster_profile_heatmap(cluster_profiles: pd.DataFrame,
                                   figsize: Tuple[int, int] = (12, 6),
                                   save_path: str = None) -> plt.Figure:
    """
    Create a heatmap showing mean resistance profile per cluster.
    
    Uses consistent cluster labeling (C1, C2, ...) as "resistance phenotypes".
    
    Parameters:
    -----------
    cluster_profiles : pd.DataFrame
        Mean resistance values per cluster
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save the figure
    
    Returns:
    --------
    matplotlib.Figure
        Cluster profile heatmap figure
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # Clean column names
    display_cols = [c.replace('_encoded', '') for c in cluster_profiles.columns]
    # Use consistent C1, C2, ... labeling (resistance phenotypes)
    cluster_labels = [f'C{i} (Resistance Phenotype)' for i in cluster_profiles.index]
    
    im = ax.imshow(cluster_profiles.values, aspect='auto', cmap='RdYlGn_r',
                   vmin=0, vmax=2)
    
    ax.set_xticks(np.arange(len(display_cols)))
    ax.set_xticklabels(display_cols, rotation=45, ha='right')
    ax.set_yticks(np.arange(len(cluster_labels)))
    ax.set_yticklabels(cluster_labels)
    ax.set_xlabel('Antibiotics')
    ax.set_ylabel('Cluster (Resistance Phenotype)')
    ax.set_title('Mean Resistance Profile by Cluster\n(Clusters represent resistance phenotypes, NOT taxonomic groups)')
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Mean Resistance Level\n(0=S, 1=I, 2=R)')
    
    # Add text annotations
    for i in range(len(cluster_labels)):
        for j in range(len(display_cols)):
            value = cluster_profiles.values[i, j]
            color = 'white' if value > 1.5 else 'black'
            ax.text(j, i, f'{value:.2f}', ha='center', va='center', color=color, fontsize=8)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Cluster profile heatmap saved to: {save_path}")
    
    return fig


def create_mdr_distribution_plot(df: pd.DataFrame,
                                 group_col: str = 'CLUSTER',
                                 figsize: Tuple[int, int] = (10, 6),
                                 save_path: str = None) -> plt.Figure:
    """
    Create a bar plot showing MDR distribution by group.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe with MDR_FLAG and group column
    group_col : str
        Column to group by
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save the figure
    
    Returns:
    --------
    matplotlib.Figure
        MDR distribution figure
    """
    if 'MDR_FLAG' not in df.columns:
        print("MDR_FLAG column not found")
        return None
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # Calculate MDR proportions
    mdr_by_group = df.groupby(group_col)['MDR_FLAG'].agg(['sum', 'count'])
    mdr_by_group['proportion'] = mdr_by_group['sum'] / mdr_by_group['count']
    
    groups = [f'{group_col} {i}' for i in mdr_by_group.index]
    proportions = mdr_by_group['proportion'].values * 100
    
    bars = ax.bar(groups, proportions, color='#F44336', alpha=0.7, edgecolor='black')
    
    ax.set_xlabel(group_col)
    ax.set_ylabel('MDR Proportion (%)')
    ax.set_title(f'MDR Distribution by {group_col}')
    ax.set_ylim(0, 100)
    
    # Add value labels
    for bar, prop in zip(bars, proportions):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 2,
                f'{prop:.1f}%', ha='center', va='bottom')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"MDR distribution plot saved to: {save_path}")
    
    return fig


def create_mar_distribution_plot(df: pd.DataFrame,
                                group_col: str = 'CLUSTER',
                                figsize: Tuple[int, int] = (10, 6),
                                save_path: str = None) -> plt.Figure:
    """
    Create a box plot showing MAR index distribution by group.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe with MAR_INDEX_COMPUTED and group column
    group_col : str
        Column to group by
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save the figure
    
    Returns:
    --------
    matplotlib.Figure
        MAR distribution figure
    """
    if 'MAR_INDEX_COMPUTED' not in df.columns:
        print("MAR_INDEX_COMPUTED column not found")
        return None
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # Prepare data for boxplot
    groups = sorted(df[group_col].unique())
    data_by_group = [df[df[group_col] == g]['MAR_INDEX_COMPUTED'].dropna().values 
                     for g in groups]
    
    bp = ax.boxplot(data_by_group, labels=[f'{group_col} {g}' for g in groups],
                    patch_artist=True)
    
    # Color the boxes
    colors = plt.cm.Set3(np.linspace(0, 1, len(groups)))
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
    
    ax.set_xlabel(group_col)
    ax.set_ylabel('MAR Index')
    ax.set_title(f'MAR Index Distribution by {group_col}')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"MAR distribution plot saved to: {save_path}")
    
    return fig


def generate_all_visualizations(df: pd.DataFrame,
                               feature_cols: List[str],
                               linkage_matrix: np.ndarray,
                               output_dir: str,
                               clustering_info: Dict = None) -> Dict[str, str]:
    """
    Generate all visualizations and save to output directory.
    
    DELIVERABLES CHECKLIST:
    - [OK] High-resolution dendrogram figure
    - [OK] Dendrogram-linked heatmap
    - [OK] MDR-enhanced visualization
    - [OK] Cluster profile heatmap
    
    Parameters:
    -----------
    df : pd.DataFrame
        Clustered dataframe
    feature_cols : list
        List of feature column names
    linkage_matrix : np.ndarray
        Linkage matrix
    output_dir : str
        Directory to save figures
    clustering_info : dict, optional
        Clustering info for dendrogram order consistency
    
    Returns:
    --------
    dict
        Dictionary with paths to saved figures
    """
    os.makedirs(output_dir, exist_ok=True)
    
    saved_figures = {}
    
    print("=" * 70)
    print("PHASE 3.2: Generating Visualizations (Dendrogram-Linked)")
    print("=" * 70)
    
    # 1. Resistance heatmap (basic)
    print("\n1. Creating basic resistance heatmap...")
    heatmap_path = os.path.join(output_dir, 'resistance_heatmap.png')
    create_resistance_heatmap(df, feature_cols, save_path=heatmap_path)
    saved_figures['heatmap'] = heatmap_path
    
    # 2. High-resolution dendrogram (PRIMARY RESULT)
    print("\n2. Creating HIGH-RESOLUTION dendrogram...")
    dendro_path = os.path.join(output_dir, 'dendrogram_highres.png')
    fig = create_dendrogram(linkage_matrix, figsize=(16, 10), save_path=dendro_path)
    saved_figures['dendrogram'] = dendro_path
    plt.close(fig)
    
    # 3. Dendrogram-linked heatmap with cluster annotations and MDR signal (CRITICAL)
    print("\n3. Creating DENDROGRAM-LINKED heatmap with cluster & MDR annotations...")
    clustered_path = os.path.join(output_dir, 'dendrogram_linked_heatmap.png')
    fig, dendro_order = create_clustered_heatmap_with_dendrogram(
        df, feature_cols, linkage_matrix, 
        figsize=(18, 14),
        save_path=clustered_path,
        show_mdr=True,
        show_cluster_bars=True
    )
    saved_figures['dendrogram_linked_heatmap'] = clustered_path
    plt.close(fig)
    
    # 4. Cluster profiles with consistent labeling
    # Try relative import first, fall back to absolute import for different contexts
    try:
        from ..clustering.hierarchical_clustering import get_cluster_profiles, create_cluster_summary_table, create_environmental_distribution_table
    except (ImportError, ValueError):
        # Fall back for running as standalone or from different entry points
        try:
            from clustering.hierarchical_clustering import get_cluster_profiles, create_cluster_summary_table, create_environmental_distribution_table
        except ImportError:
            import sys
            import os as os_module
            src_path = os_module.path.dirname(os_module.path.dirname(os_module.path.abspath(__file__)))
            if src_path not in sys.path:
                sys.path.insert(0, src_path)
            from clustering.hierarchical_clustering import get_cluster_profiles, create_cluster_summary_table, create_environmental_distribution_table
    
    print("\n4. Creating cluster profile heatmap (resistance phenotypes)...")
    cluster_profiles = get_cluster_profiles(df, feature_cols)
    profile_path = os.path.join(output_dir, 'cluster_profiles.png')
    fig = create_cluster_profile_heatmap(cluster_profiles, save_path=profile_path)
    saved_figures['cluster_profiles'] = profile_path
    plt.close(fig)
    
    # 5. MDR distribution by cluster
    print("\n5. Creating MDR distribution plot...")
    mdr_path = os.path.join(output_dir, 'mdr_distribution.png')
    fig = create_mdr_distribution_plot(df, save_path=mdr_path)
    if fig:
        saved_figures['mdr_distribution'] = mdr_path
        plt.close(fig)
    
    # 6. MAR distribution by cluster
    print("\n6. Creating MAR distribution plot...")
    mar_path = os.path.join(output_dir, 'mar_distribution.png')
    fig = create_mar_distribution_plot(df, save_path=mar_path)
    if fig:
        saved_figures['mar_distribution'] = mar_path
        plt.close(fig)
    
    # 7. Generate cluster summary table visualization
    print("\n7. Creating cluster summary table visualization...")
    summary_table = create_cluster_summary_table(df, feature_cols)
    
    if not summary_table.empty:
        # Save as CSV
        table_csv_path = os.path.join(output_dir, 'cluster_summary_table.csv')
        summary_table.to_csv(table_csv_path, index=False)
        saved_figures['cluster_summary_table_csv'] = table_csv_path
        print(f"   Saved cluster summary table to: {table_csv_path}")
        
        # Create figure of the table
        fig, ax = plt.subplots(figsize=(14, 4 + len(summary_table) * 0.5))
        ax.axis('tight')
        ax.axis('off')
        
        table = ax.table(
            cellText=summary_table.values,
            colLabels=summary_table.columns,
            cellLoc='center',
            loc='center'
        )
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1.2, 1.5)
        
        # Style header
        for (row, col), cell in table.get_celld().items():
            if row == 0:
                cell.set_text_props(fontweight='bold')
                cell.set_facecolor('#4CAF50')
                cell.set_text_props(color='white')
        
        ax.set_title('MANDATORY CLUSTER SUMMARY TABLE\n(Clusters = Resistance Phenotypes)',
                    fontsize=12, fontweight='bold', pad=20)
        
        table_fig_path = os.path.join(output_dir, 'cluster_summary_table.png')
        plt.savefig(table_fig_path, dpi=300, bbox_inches='tight')
        saved_figures['cluster_summary_table_fig'] = table_fig_path
        print(f"   Saved cluster summary table figure to: {table_fig_path}")
        plt.close(fig)
    
    # 8. Generate environmental distribution table
    print("\n8. Creating environmental distribution table...")
    env_table = create_environmental_distribution_table(df)
    
    if not env_table.empty:
        env_csv_path = os.path.join(output_dir, 'environmental_distribution_table.csv')
        env_table.to_csv(env_csv_path)
        saved_figures['environmental_distribution_csv'] = env_csv_path
        print(f"   Saved environmental distribution table to: {env_csv_path}")
    
    plt.close('all')
    
    print(f"\n{'='*70}")
    print(f"All visualizations saved to: {output_dir}")
    print(f"{'='*70}")
    
    return saved_figures


if __name__ == "__main__":
    from pathlib import Path
    import sys
    sys.path.insert(0, str(Path(__file__).parent.parent))
    from clustering.hierarchical_clustering import run_clustering_pipeline
    
    project_root = Path(__file__).parent.parent.parent
    analysis_path = project_root / "data" / "processed" / "analysis_ready_dataset.csv"
    
    if analysis_path.exists():
        df = pd.read_csv(analysis_path)
        feature_cols = [c for c in df.columns if c.endswith('_encoded')]
        
        # Run clustering first
        df_clustered, linkage_matrix, _ = run_clustering_pipeline(df, feature_cols)
        
        # Generate visualizations
        output_dir = project_root / "data" / "processed" / "figures"
        generate_all_visualizations(df_clustered, feature_cols, linkage_matrix, str(output_dir))
    else:
        print(f"Analysis-ready dataset not found at {analysis_path}")
