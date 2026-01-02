"""
Visualization Module for AMR Pattern Recognition Dashboard
Phase 6 - Mandatory Visual Components

This module provides visualization functions for:
- Heatmap with dendrogram (read-only)
- PCA scatter plots with cluster/region/environment toggle
- Cluster summary tables
- Feature importance bar plots
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from typing import List, Tuple, Dict, Optional


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


def create_resistance_heatmap(df: pd.DataFrame,
                             feature_cols: List[str],
                             cluster_col: str = 'CLUSTER',
                             figsize: Tuple[int, int] = (14, 10)) -> plt.Figure:
    """
    Create a read-only resistance heatmap.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe with encoded resistance values
    feature_cols : list
        List of feature column names
    cluster_col : str
        Column name for cluster labels
    figsize : tuple
        Figure size
    
    Returns:
    --------
    matplotlib.Figure
        Heatmap figure (read-only)
    """
    existing_cols = [c for c in feature_cols if c in df.columns]
    
    # Sort by cluster if available
    if cluster_col in df.columns:
        df_sorted = df.sort_values(cluster_col)
    else:
        df_sorted = df.copy()
    
    data_matrix = df_sorted[existing_cols].values
    display_cols = [c.replace('_encoded', '') for c in existing_cols]
    
    fig, ax = plt.subplots(figsize=figsize)
    
    im = ax.imshow(data_matrix, aspect='auto', cmap=RESISTANCE_CMAP, vmin=0, vmax=2)
    
    ax.set_xticks(np.arange(len(display_cols)))
    ax.set_xticklabels(display_cols, rotation=45, ha='right', fontsize=9)
    ax.set_xlabel('Antibiotics', fontsize=11)
    ax.set_ylabel('Isolates', fontsize=11)
    ax.set_title('Resistance Profile Heatmap (Read-Only)', fontsize=12, fontweight='bold')
    
    # Add cluster boundaries
    if cluster_col in df.columns:
        cluster_labels = df_sorted[cluster_col].values
        cluster_boundaries = np.where(np.diff(cluster_labels) != 0)[0]
        for boundary in cluster_boundaries:
            ax.axhline(y=boundary + 0.5, color='black', linewidth=1.5)
    
    cbar = plt.colorbar(im, ax=ax, ticks=[0, 1, 2])
    cbar.set_label('Resistance Status')
    cbar.ax.set_yticklabels(['S (Susceptible)', 'I (Intermediate)', 'R (Resistant)'])
    
    plt.tight_layout()
    return fig


def create_heatmap_with_dendrogram(df: pd.DataFrame,
                                   feature_cols: List[str],
                                   figsize: Tuple[int, int] = (16, 12)) -> plt.Figure:
    """
    Create a read-only heatmap with dendrogram.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe with encoded resistance values
    feature_cols : list
        List of feature column names
    figsize : tuple
        Figure size
    
    Returns:
    --------
    matplotlib.Figure
        Heatmap with dendrogram figure (read-only)
    """
    existing_cols = [c for c in feature_cols if c in df.columns]
    data_matrix = df[existing_cols].values
    display_cols = [c.replace('_encoded', '') for c in existing_cols]
    
    # Handle missing values
    imputer = SimpleImputer(strategy='median')
    data_imputed = imputer.fit_transform(data_matrix)
    
    # Compute linkage
    Z = linkage(data_imputed, method='ward', metric='euclidean')
    
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 4], wspace=0.02)
    
    # Dendrogram
    ax_dendro = fig.add_subplot(gs[0])
    dendro = dendrogram(Z, ax=ax_dendro, orientation='left', no_labels=True, color_threshold=0)
    ax_dendro.set_xticks([])
    ax_dendro.set_yticks([])
    ax_dendro.invert_yaxis()
    ax_dendro.set_title('Dendrogram', fontsize=10)
    
    # Get dendrogram order
    dendro_order = dendro['leaves']
    data_ordered = data_imputed[dendro_order, :]
    
    # Heatmap
    ax_heatmap = fig.add_subplot(gs[1])
    im = ax_heatmap.imshow(data_ordered, aspect='auto', cmap=RESISTANCE_CMAP, vmin=0, vmax=2)
    ax_heatmap.set_xticks(np.arange(len(display_cols)))
    ax_heatmap.set_xticklabels(display_cols, rotation=45, ha='right', fontsize=9)
    ax_heatmap.set_yticks([])
    ax_heatmap.set_xlabel('Antibiotics', fontsize=11)
    ax_heatmap.set_title('Resistance Profile Heatmap (Read-Only)\nOrdered by Hierarchical Clustering', 
                         fontsize=12, fontweight='bold')
    
    # Colorbar
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    cbar = plt.colorbar(im, cax=cbar_ax, ticks=[0, 1, 2])
    cbar.ax.set_yticklabels(['S', 'I', 'R'])
    cbar.set_label('Resistance')
    
    plt.suptitle('Dendrogram-Linked Resistance Heatmap', fontsize=14, fontweight='bold', y=0.98)
    
    return fig


def perform_pca_analysis(df: pd.DataFrame, 
                        feature_cols: List[str]) -> Tuple[np.ndarray, PCA, List[str]]:
    """
    Perform PCA on resistance data.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe with encoded resistance values
    feature_cols : list
        List of feature column names
    
    Returns:
    --------
    tuple
        (PCA-transformed data, fitted PCA object, feature names used)
    """
    existing_cols = [c for c in feature_cols if c in df.columns]
    X = df[existing_cols].copy()
    
    imputer = SimpleImputer(strategy='median')
    X_imputed = imputer.fit_transform(X)
    
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_imputed)
    
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_scaled)
    
    return X_pca, pca, existing_cols


def create_pca_plot(X_pca: np.ndarray, 
                   df: pd.DataFrame, 
                   color_col: str, 
                   pca: PCA,
                   figsize: Tuple[int, int] = (10, 8)) -> plt.Figure:
    """
    Create a PCA scatter plot colored by a specified column.
    
    Parameters:
    -----------
    X_pca : np.ndarray
        PCA-transformed data
    df : pd.DataFrame
        Original dataframe
    color_col : str
        Column to use for coloring points
    pca : PCA
        Fitted PCA object
    figsize : tuple
        Figure size
    
    Returns:
    --------
    matplotlib.Figure
        PCA scatter plot
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    if color_col in df.columns:
        categories = df[color_col].dropna().unique()
        n_categories = len(categories)
        
        if n_categories <= 10:
            colors = CLUSTER_COLORS[:n_categories]
        else:
            colors = plt.cm.tab20(np.linspace(0, 1, n_categories))
        
        for i, category in enumerate(sorted(categories)):
            mask = df[color_col] == category
            ax.scatter(X_pca[mask, 0], X_pca[mask, 1],
                      c=[colors[i] if n_categories <= 10 else colors[i]], 
                      label=str(category), alpha=0.7, s=50, edgecolors='white', linewidth=0.5)
        
        ax.legend(title=color_col, bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)
    else:
        ax.scatter(X_pca[:, 0], X_pca[:, 1], alpha=0.7, s=50, c='steelblue', 
                  edgecolors='white', linewidth=0.5)
    
    ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)', fontsize=11)
    ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)', fontsize=11)
    ax.set_title('PCA of Resistance Profiles', fontsize=12, fontweight='bold')
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    return fig


def create_cluster_summary_table(df: pd.DataFrame, 
                                feature_cols: List[str]) -> pd.DataFrame:
    """
    Create a cluster summary table.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe with cluster labels
    feature_cols : list
        List of feature column names
    
    Returns:
    --------
    pd.DataFrame
        Cluster summary table
    """
    if 'CLUSTER' not in df.columns:
        return pd.DataFrame()
    
    existing_cols = [c for c in feature_cols if c in df.columns]
    summary_rows = []
    
    for cluster in sorted(df['CLUSTER'].unique()):
        cluster_df = df[df['CLUSTER'] == cluster]
        n_isolates = len(cluster_df)
        
        row = {
            'Cluster': f'C{int(cluster)}',
            'N Isolates': n_isolates,
            'Proportion (%)': round(n_isolates / len(df) * 100, 1)
        }
        
        # Add MDR proportion if available
        if 'MDR_FLAG' in cluster_df.columns:
            row['MDR (%)'] = round(cluster_df['MDR_FLAG'].mean() * 100, 1)
        
        # Add mean MAR index if available
        if 'MAR_INDEX_COMPUTED' in cluster_df.columns:
            row['Mean MAR'] = round(cluster_df['MAR_INDEX_COMPUTED'].mean(), 3)
        
        # Add dominant species if available
        if 'ISOLATE_ID' in cluster_df.columns:
            species_counts = cluster_df['ISOLATE_ID'].value_counts()
            row['Dominant Species'] = species_counts.index[0]
            row['Species (%)'] = round(species_counts.iloc[0] / n_isolates * 100, 1)
        
        # Add top resistant antibiotics
        if existing_cols:
            mean_resistance = cluster_df[existing_cols].mean()
            top_resistant = mean_resistance.sort_values(ascending=False).head(3)
            top_abs = [c.replace('_encoded', '') for c in top_resistant.index]
            row['Top Resistant'] = ', '.join(top_abs)
        
        summary_rows.append(row)
    
    return pd.DataFrame(summary_rows)


def create_feature_importance_plot(importance_scores: Dict[str, float],
                                   title: str = 'Feature Importance',
                                   top_n: int = 15,
                                   figsize: Tuple[int, int] = (10, 6)) -> plt.Figure:
    """
    Create a feature importance bar plot.
    
    Parameters:
    -----------
    importance_scores : dict
        Dictionary mapping feature names to importance scores
    title : str
        Plot title
    top_n : int
        Number of top features to show
    figsize : tuple
        Figure size
    
    Returns:
    --------
    matplotlib.Figure
        Feature importance bar plot
    """
    # Sort by importance and take top N
    sorted_items = sorted(importance_scores.items(), key=lambda x: x[1], reverse=True)[:top_n]
    
    if not sorted_items:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, 'No feature importance data available', 
                ha='center', va='center', fontsize=12)
        ax.axis('off')
        return fig
    
    features = [item[0].replace('_encoded', '') for item in sorted_items]
    scores = [item[1] for item in sorted_items]
    
    fig, ax = plt.subplots(figsize=figsize)
    
    bars = ax.barh(range(len(features)), scores, color='steelblue', edgecolor='black', alpha=0.8)
    ax.set_yticks(range(len(features)))
    ax.set_yticklabels(features, fontsize=10)
    ax.invert_yaxis()
    ax.set_xlabel('Importance Score', fontsize=11)
    ax.set_title(f'{title}\n(Top {top_n} Features)', fontsize=12, fontweight='bold')
    
    # Add value labels
    for bar, score in zip(bars, scores):
        ax.text(bar.get_width() + 0.001, bar.get_y() + bar.get_height()/2,
                f'{score:.4f}', va='center', fontsize=9)
    
    plt.tight_layout()
    return fig


def create_mdr_distribution_plot(df: pd.DataFrame, 
                                group_col: str = 'CLUSTER',
                                figsize: Tuple[int, int] = (10, 6)) -> Optional[plt.Figure]:
    """
    Create MDR distribution bar plot.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe with MDR_FLAG and group column
    group_col : str
        Column to group by
    figsize : tuple
        Figure size
    
    Returns:
    --------
    matplotlib.Figure or None
        MDR distribution plot
    """
    if 'MDR_FLAG' not in df.columns or group_col not in df.columns:
        return None
    
    fig, ax = plt.subplots(figsize=figsize)
    
    mdr_by_group = df.groupby(group_col)['MDR_FLAG'].mean() * 100
    
    groups = [f'C{int(g)}' if group_col == 'CLUSTER' else str(g) for g in mdr_by_group.index]
    proportions = mdr_by_group.values
    
    bars = ax.bar(groups, proportions, color='#F44336', alpha=0.8, edgecolor='black')
    
    ax.set_xlabel(group_col, fontsize=11)
    ax.set_ylabel('MDR Proportion (%)', fontsize=11)
    ax.set_title(f'MDR Distribution by {group_col}', fontsize=12, fontweight='bold')
    ax.set_ylim(0, 100)
    
    # Add value labels
    for bar, prop in zip(bars, proportions):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 2,
                f'{prop:.1f}%', ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    return fig


def create_regional_distribution_plot(df: pd.DataFrame,
                                      figsize: Tuple[int, int] = (10, 6)) -> Optional[plt.Figure]:
    """
    Create regional distribution bar plot.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe with REGION column
    figsize : tuple
        Figure size
    
    Returns:
    --------
    matplotlib.Figure or None
        Regional distribution plot
    """
    if 'REGION' not in df.columns:
        return None
    
    fig, ax = plt.subplots(figsize=figsize)
    
    region_counts = df['REGION'].value_counts()
    
    bars = ax.bar(region_counts.index, region_counts.values, 
                 color='steelblue', alpha=0.8, edgecolor='black')
    
    ax.set_xlabel('Region', fontsize=11)
    ax.set_ylabel('Number of Isolates', fontsize=11)
    ax.set_title('Isolates by Region', fontsize=12, fontweight='bold')
    
    plt.xticks(rotation=45, ha='right')
    
    # Add value labels
    for bar in bars:
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                f'{int(bar.get_height())}', ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    return fig
