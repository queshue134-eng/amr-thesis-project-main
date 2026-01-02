"""
Clustering Module for AMR Pattern Recognition Dashboard
Phase 6 - Clustering utilities for visualization

This module provides clustering utilities for the dashboard.
Note: Actual clustering is done by the main pipeline - this module
provides read-only visualization and interpretation support.
"""

import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from sklearn.impute import SimpleImputer
from typing import List, Tuple, Dict, Optional


def get_cluster_summary(df: pd.DataFrame,
                       feature_cols: List[str] = None,
                       cluster_col: str = 'CLUSTER') -> Dict:
    """
    Get comprehensive summary statistics for each cluster.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe with cluster labels
    feature_cols : list, optional
        Feature columns for resistance analysis
    cluster_col : str
        Column name for cluster labels
    
    Returns:
    --------
    dict
        Summary statistics per cluster
    """
    if cluster_col not in df.columns:
        return {}
    
    if feature_cols is None:
        feature_cols = [c for c in df.columns if c.endswith('_encoded')]
    
    summary = {}
    total_isolates = len(df)
    
    for cluster_id in sorted(df[cluster_col].unique()):
        cluster_df = df[df[cluster_col] == cluster_id]
        n_isolates = len(cluster_df)
        
        cluster_summary = {
            'cluster_label': f'C{int(cluster_id)}',
            'n_isolates': n_isolates,
            'percentage': round((n_isolates / total_isolates) * 100, 1)
        }
        
        # MDR proportion
        if 'MDR_FLAG' in cluster_df.columns:
            mdr_count = cluster_df['MDR_FLAG'].sum()
            cluster_summary['mdr_count'] = int(mdr_count)
            cluster_summary['mdr_proportion'] = round((mdr_count / n_isolates) * 100, 1)
        
        # Mean MAR index
        if 'MAR_INDEX_COMPUTED' in cluster_df.columns:
            cluster_summary['mean_mar_index'] = round(cluster_df['MAR_INDEX_COMPUTED'].mean(), 4)
        
        # Species composition
        if 'ISOLATE_ID' in cluster_df.columns:
            species_counts = cluster_df['ISOLATE_ID'].value_counts()
            if len(species_counts) > 0:
                cluster_summary['dominant_species'] = species_counts.index[0]
                cluster_summary['dominant_species_pct'] = round(
                    (species_counts.iloc[0] / n_isolates) * 100, 1
                )
        
        # Top resistant antibiotics
        existing_cols = [c for c in feature_cols if c in cluster_df.columns]
        if existing_cols:
            mean_resistance = cluster_df[existing_cols].mean()
            top_resistant = mean_resistance.sort_values(ascending=False).head(5)
            cluster_summary['top_resistant_antibiotics'] = [
                c.replace('_encoded', '') for c in top_resistant.index.tolist()
            ]
        
        # Regional distribution
        if 'REGION' in cluster_df.columns:
            region_counts = cluster_df['REGION'].value_counts()
            if len(region_counts) > 0:
                cluster_summary['major_region'] = region_counts.index[0]
                cluster_summary['major_region_pct'] = round(
                    (region_counts.iloc[0] / n_isolates) * 100, 1
                )
        
        summary[int(cluster_id)] = cluster_summary
    
    return summary


def get_cluster_profiles(df: pd.DataFrame,
                        feature_cols: List[str],
                        cluster_col: str = 'CLUSTER') -> pd.DataFrame:
    """
    Calculate mean resistance profile for each cluster.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe with cluster labels
    feature_cols : list
        List of feature column names
    cluster_col : str
        Column name for cluster labels
    
    Returns:
    --------
    pd.DataFrame
        Mean resistance profile per cluster
    """
    if cluster_col not in df.columns:
        return pd.DataFrame()
    
    existing_cols = [c for c in feature_cols if c in df.columns]
    if not existing_cols:
        return pd.DataFrame()
    
    cluster_profiles = df.groupby(cluster_col)[existing_cols].mean()
    
    return cluster_profiles


def compute_linkage_for_visualization(df: pd.DataFrame,
                                      feature_cols: List[str],
                                      method: str = 'ward',
                                      metric: str = 'euclidean') -> np.ndarray:
    """
    Compute linkage matrix for dendrogram visualization.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe with encoded resistance values
    feature_cols : list
        List of feature column names
    method : str
        Linkage method (ward, complete, average, single)
    metric : str
        Distance metric
    
    Returns:
    --------
    np.ndarray
        Linkage matrix
    """
    existing_cols = [c for c in feature_cols if c in df.columns]
    data_matrix = df[existing_cols].values
    
    # Handle missing values
    imputer = SimpleImputer(strategy='median')
    data_imputed = imputer.fit_transform(data_matrix)
    
    # Compute linkage (ward requires euclidean)
    if method == 'ward':
        metric = 'euclidean'
    
    Z = linkage(data_imputed, method=method, metric=metric)
    
    return Z


def get_dendrogram_order(linkage_matrix: np.ndarray) -> List[int]:
    """
    Get the leaf order from a dendrogram.
    
    Parameters:
    -----------
    linkage_matrix : np.ndarray
        Linkage matrix from hierarchical clustering
    
    Returns:
    --------
    list
        Order of leaves in dendrogram
    """
    import matplotlib.pyplot as plt
    
    fig, ax = plt.subplots()
    dendro = dendrogram(linkage_matrix, no_plot=True)
    plt.close(fig)
    
    return dendro['leaves']


def describe_archetype(cluster_summary: Dict,
                      cluster_id: int) -> str:
    """
    Generate human-readable archetype description.
    
    Parameters:
    -----------
    cluster_summary : dict
        Cluster summary from get_cluster_summary()
    cluster_id : int
        Cluster ID to describe
    
    Returns:
    --------
    str
        Archetype description
    """
    if cluster_id not in cluster_summary:
        return f"No data available for Cluster {cluster_id}"
    
    info = cluster_summary[cluster_id]
    
    description = f"**{info.get('cluster_label', f'C{cluster_id}')}** "
    description += f"({info.get('n_isolates', 0)} isolates, {info.get('percentage', 0):.1f}%)"
    
    parts = []
    
    if 'mdr_proportion' in info:
        mdr_pct = info['mdr_proportion']
        if mdr_pct >= 70:
            parts.append(f"High MDR ({mdr_pct:.0f}%)")
        elif mdr_pct >= 30:
            parts.append(f"Mixed MDR ({mdr_pct:.0f}%)")
        else:
            parts.append(f"Low MDR ({mdr_pct:.0f}%)")
    
    if 'dominant_species' in info:
        parts.append(f"Dominant: {info['dominant_species']} ({info.get('dominant_species_pct', 0):.0f}%)")
    
    if 'top_resistant_antibiotics' in info:
        top_abs = info['top_resistant_antibiotics'][:3]
        parts.append(f"Resistant to: {', '.join(top_abs)}")
    
    if parts:
        description += "\n- " + "\n- ".join(parts)
    
    return description
