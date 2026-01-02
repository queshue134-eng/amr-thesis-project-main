"""
Clustering Module for AMR Thesis Project
Phase 3.1 - Hierarchical Agglomerative Clustering for Structure Identification

This module performs unsupervised structure identification to discover natural
groupings in resistance patterns without pre-defined categories.

METHODOLOGY NOTES:
------------------
- Clustering is performed ONLY on resistance data (X_resistance matrix).
- Metadata (species, region, environment, MDR status) is EXCLUDED from clustering
  to prevent data leakage and circular analysis.
- Clusters should be interpreted as "resistance phenotypes" - characteristic
  patterns of antibiotic susceptibility/resistance, NOT taxonomic or causal groups.
- Post-hoc analysis may reveal associations between clusters and metadata,
  but these are correlational findings, not causal relationships.

CLUSTERING PARAMETERS (explicitly defined for reproducibility):
---------------------------------------------------------------
- LINKAGE_METHOD: "ward" - Ward's minimum variance method minimizes the total
  within-cluster variance. Suitable for numerical resistance vectors where we
  want compact, spherical clusters. Ward linkage requires Euclidean distance.
  
- DISTANCE_METRIC: "euclidean" (primary) - The standard choice for numerical
  data and required for Ward linkage. Manhattan distance can be used as a
  robustness check with other linkage methods.
  
- CLUSTER_CUT_RULE: "maxclust" (default) - Clusters are determined by specifying
  a fixed number of clusters based on visual elbow analysis of the dendrogram.
  Alternative: "distance" threshold or "inconsistency" coefficient can be used.

References:
-----------
- Ward, J. H. (1963). Hierarchical Grouping to Optimize an Objective Function.
  Journal of the American Statistical Association, 58(301), 236-244.
- Magiorakos et al. (2012). Multidrug-resistant, extensively drug-resistant and
  pandrug-resistant bacteria. Clinical Microbiology and Infection, 18(3), 268-281.
"""

import pandas as pd
import numpy as np
from typing import List, Tuple, Dict, Optional
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram, inconsistent
from scipy.spatial.distance import pdist, squareform
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
import warnings
import pickle
import os

# Import from centralized configuration
try:
    from config import CLUSTERING_CONFIG, RANDOM_STATE
except ImportError:
    # Fallback for standalone execution
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).parent.parent))
    from config import CLUSTERING_CONFIG, RANDOM_STATE


# =============================================================================
# CLUSTERING PARAMETERS - Now imported from centralized config.py
# =============================================================================

# Primary linkage method: Ward's minimum variance method
# Justification: Ward minimizes within-cluster variance, producing compact clusters
# that are appropriate for identifying distinct resistance phenotypes.
LINKAGE_METHOD = CLUSTERING_CONFIG.get('linkage_method', 'ward')

# Primary distance metric: Euclidean distance
# Justification: Required for Ward linkage; standard for numerical resistance data
DISTANCE_METRIC_PRIMARY = CLUSTERING_CONFIG.get('distance_metric', 'euclidean')

# Alternative distance metric for robustness checking
# Justification: Manhattan (cityblock) distance is more robust to outliers and can serve
# as a sensitivity analysis to verify cluster stability
# Note: scipy uses 'cityblock' for Manhattan distance
DISTANCE_METRIC_ROBUSTNESS = "cityblock"

# Default number of clusters (should be overridden by elbow/silhouette analysis)
DEFAULT_N_CLUSTERS = CLUSTERING_CONFIG.get('default_n_clusters', 4)

# Cluster cutting criterion options:
# - "maxclust": Fixed number of clusters (based on elbow/domain knowledge)
# - "distance": Distance threshold for cluster formation
# - "inconsistent": Inconsistency coefficient threshold
CLUSTER_CUT_CRITERION = "maxclust"

# Inconsistency threshold for 'inconsistent' criterion (if used)
# Values > 1.0 indicate significant discontinuity in cluster structure
INCONSISTENCY_THRESHOLD = 1.5


def prepare_clustering_data(df: pd.DataFrame,
                           feature_cols: List[str],
                           impute_strategy: str = 'median') -> Tuple[np.ndarray, pd.DataFrame]:
    """
    Prepare data for clustering by handling missing values.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe with encoded resistance values
    feature_cols : list
        List of feature column names
    impute_strategy : str
        Strategy for imputing missing values ('mean', 'median', 'most_frequent')
    
    Returns:
    --------
    tuple
        (Imputed feature matrix as numpy array, Original dataframe with valid rows)
    """
    # Extract feature matrix
    existing_cols = [c for c in feature_cols if c in df.columns]
    feature_matrix = df[existing_cols].copy()
    
    # Track which rows have data
    valid_mask = feature_matrix.notna().any(axis=1)
    feature_matrix_valid = feature_matrix[valid_mask].copy()
    df_valid = df[valid_mask].copy()
    
    # Impute missing values
    imputer = SimpleImputer(strategy=impute_strategy)
    imputed_data = imputer.fit_transform(feature_matrix_valid)
    
    print(f"Prepared {imputed_data.shape[0]} isolates with {imputed_data.shape[1]} features")
    
    return imputed_data, df_valid


def compute_distance_matrix(data: np.ndarray, metric: str = 'euclidean') -> np.ndarray:
    """
    Compute pairwise distance matrix.
    
    Parameters:
    -----------
    data : np.ndarray
        Feature matrix
    metric : str
        Distance metric ('euclidean', 'manhattan', 'cosine')
    
    Returns:
    --------
    np.ndarray
        Condensed distance matrix
    """
    return pdist(data, metric=metric)


def perform_hierarchical_clustering(data: np.ndarray,
                                   method: str = 'ward',
                                   metric: str = 'euclidean') -> np.ndarray:
    """
    Perform hierarchical agglomerative clustering.
    
    Parameters:
    -----------
    data : np.ndarray
        Feature matrix (samples x features)
    method : str
        Linkage method ('ward', 'complete', 'average', 'single')
    metric : str
        Distance metric ('euclidean', 'manhattan')
    
    Returns:
    --------
    np.ndarray
        Linkage matrix for dendrogram
    """
    if method == 'ward':
        # Ward's method requires Euclidean distance
        metric = 'euclidean'
    
    # Compute linkage
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        Z = linkage(data, method=method, metric=metric)
    
    return Z


def assign_clusters(linkage_matrix: np.ndarray,
                   n_clusters: int = None,
                   distance_threshold: float = None,
                   inconsistency_threshold: float = None,
                   criterion: str = None) -> Tuple[np.ndarray, Dict]:
    """
    Assign cluster labels based on linkage matrix.
    
    CLUSTER CUT RULE DOCUMENTATION:
    -------------------------------
    This function supports three methods for determining clusters:
    
    1. MAXCLUST (default): Fixed number of clusters
       - Uses n_clusters parameter
       - Appropriate when number of groups is known/estimated from elbow analysis
       - Most reproducible method
    
    2. DISTANCE: Distance threshold cut
       - Uses distance_threshold parameter
       - Cuts dendrogram at specified height
       - More interpretable but may vary with data scale
    
    3. INCONSISTENCY: Inconsistency coefficient
       - Uses inconsistency_threshold parameter
       - Identifies natural breaks in cluster structure
       - Adaptive but can be sensitive to noise
    
    Parameters:
    -----------
    linkage_matrix : np.ndarray
        Linkage matrix from hierarchical clustering
    n_clusters : int, optional
        Number of clusters to form (for 'maxclust' criterion)
    distance_threshold : float, optional
        Distance threshold for forming clusters (for 'distance' criterion)
    inconsistency_threshold : float, optional
        Inconsistency coefficient threshold (for 'inconsistent' criterion)
    criterion : str, optional
        Override automatic criterion selection ('maxclust', 'distance', 'inconsistent')
    
    Returns:
    --------
    tuple
        (Cluster labels array, Cut rule info dict)
    """
    cut_info = {
        'criterion': None,
        'threshold': None,
        'n_clusters_formed': None,
        'justification': None
    }
    
    # Determine criterion based on provided parameters
    if criterion is not None:
        selected_criterion = criterion
    elif n_clusters is not None:
        selected_criterion = 'maxclust'
    elif distance_threshold is not None:
        selected_criterion = 'distance'
    elif inconsistency_threshold is not None:
        selected_criterion = 'inconsistent'
    else:
        # Default: use maxclust with default number of clusters
        selected_criterion = 'maxclust'
        n_clusters = DEFAULT_N_CLUSTERS
    
    # Apply the selected criterion
    if selected_criterion == 'maxclust':
        if n_clusters is None:
            n_clusters = DEFAULT_N_CLUSTERS
        labels = fcluster(linkage_matrix, n_clusters, criterion='maxclust')
        cut_info['criterion'] = 'maxclust'
        cut_info['threshold'] = n_clusters
        cut_info['justification'] = (
            f"Fixed number of clusters (k={n_clusters}) based on visual elbow "
            "analysis of dendrogram heights and cluster coherence."
        )
    
    elif selected_criterion == 'distance':
        if distance_threshold is None:
            # Use median linkage distance as default
            distance_threshold = np.median(linkage_matrix[:, 2])
        labels = fcluster(linkage_matrix, distance_threshold, criterion='distance')
        cut_info['criterion'] = 'distance'
        cut_info['threshold'] = distance_threshold
        cut_info['justification'] = (
            f"Distance threshold cut at d={distance_threshold:.4f}. "
            "All clusters formed have within-cluster distances below this threshold."
        )
    
    elif selected_criterion == 'inconsistent':
        if inconsistency_threshold is None:
            inconsistency_threshold = INCONSISTENCY_THRESHOLD
        labels = fcluster(linkage_matrix, inconsistency_threshold, criterion='inconsistent')
        cut_info['criterion'] = 'inconsistent'
        cut_info['threshold'] = inconsistency_threshold
        cut_info['justification'] = (
            f"Inconsistency coefficient threshold of {inconsistency_threshold:.2f}. "
            "Clusters are formed at linkages where inconsistency exceeds this value, "
            "indicating natural discontinuities in the data structure."
        )
    
    else:
        raise ValueError(f"Unknown criterion: {selected_criterion}")
    
    cut_info['n_clusters_formed'] = len(np.unique(labels))
    
    return labels, cut_info


def determine_optimal_clusters(linkage_matrix: np.ndarray,
                              max_clusters: int = 10) -> Dict:
    """
    Analyze cluster quality for different numbers of clusters.
    
    Parameters:
    -----------
    linkage_matrix : np.ndarray
        Linkage matrix from hierarchical clustering
    max_clusters : int
        Maximum number of clusters to evaluate
    
    Returns:
    --------
    dict
        Dictionary with cluster quality metrics
    """
    from scipy.cluster.hierarchy import inconsistent
    
    # Compute inconsistency coefficients
    inconsist = inconsistent(linkage_matrix)
    
    # Analyze cluster sizes for different cuts
    cluster_analysis = {}
    
    for k in range(2, max_clusters + 1):
        labels = fcluster(linkage_matrix, k, criterion='maxclust')
        unique_labels, counts = np.unique(labels, return_counts=True)
        
        cluster_analysis[k] = {
            'n_clusters': k,
            'cluster_sizes': dict(zip(unique_labels.astype(int), counts.astype(int))),
            'min_size': int(counts.min()),
            'max_size': int(counts.max()),
            'size_std': float(counts.std())
        }
    
    return cluster_analysis


def perform_robustness_check(data: np.ndarray,
                            n_clusters: int,
                            primary_method: str = LINKAGE_METHOD,
                            primary_metric: str = DISTANCE_METRIC_PRIMARY,
                            robustness_metric: str = DISTANCE_METRIC_ROBUSTNESS) -> Dict:
    """
    Perform robustness check by comparing clustering with different distance metrics.
    
    This function re-runs clustering with Manhattan distance (robustness metric)
    and compares cluster stability and composition with the primary clustering.
    
    Parameters:
    -----------
    data : np.ndarray
        Feature matrix (samples x features)
    n_clusters : int
        Number of clusters to form
    primary_method : str
        Primary linkage method (default: ward)
    primary_metric : str
        Primary distance metric (default: euclidean)
    robustness_metric : str
        Alternative distance metric for robustness check (default: manhattan)
    
    Returns:
    --------
    dict
        Robustness analysis results including cluster stability metrics
    """
    results = {
        'primary_metric': primary_metric,
        'robustness_metric': robustness_metric,
        'cluster_stability': {},
        'agreement_score': None,
        'interpretation': []
    }
    
    # Primary clustering
    primary_linkage = linkage(data, method=primary_method, metric=primary_metric)
    primary_labels = fcluster(primary_linkage, n_clusters, criterion='maxclust')
    
    # For robustness with Manhattan, we need to use a different linkage method
    # since Ward requires Euclidean. Use average linkage for comparison.
    robustness_method = 'average' if primary_method == 'ward' else primary_method
    
    robustness_linkage = linkage(data, method=robustness_method, metric=robustness_metric)
    robustness_labels = fcluster(robustness_linkage, n_clusters, criterion='maxclust')
    
    # Create contingency table using pandas crosstab
    contingency_table = pd.crosstab(
        pd.Series(primary_labels, name='Primary'),
        pd.Series(robustness_labels, name='Robustness')
    )
    results['contingency_table'] = contingency_table.to_dict()
    
    # Calculate Adjusted Rand Index for cluster agreement
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
    
    ari = adjusted_rand_score(primary_labels, robustness_labels)
    nmi = normalized_mutual_info_score(primary_labels, robustness_labels)
    
    results['agreement_score'] = {
        'adjusted_rand_index': float(ari),
        'normalized_mutual_info': float(nmi)
    }
    
    # Interpretation
    if ari > 0.8:
        stability = "High"
        results['interpretation'].append(
            f"High cluster stability (ARI={ari:.3f}). Clusters are robust to "
            "distance metric choice, supporting the validity of the groupings."
        )
    elif ari > 0.5:
        stability = "Moderate"
        results['interpretation'].append(
            f"Moderate cluster stability (ARI={ari:.3f}). Most clusters are stable, "
            "but some samples may be borderline cases."
        )
    else:
        stability = "Low"
        results['interpretation'].append(
            f"Low cluster stability (ARI={ari:.3f}). Results may be sensitive to "
            "methodological choices. Consider examining cluster boundaries carefully."
        )
    
    results['stability_level'] = stability
    
    # Calculate per-cluster stability
    for cluster_id in sorted(np.unique(primary_labels)):
        primary_mask = primary_labels == cluster_id
        n_primary = np.sum(primary_mask)
        
        # Find best matching robustness cluster
        robustness_assignments = robustness_labels[primary_mask]
        most_common = np.bincount(robustness_assignments).argmax()
        match_count = np.sum(robustness_assignments == most_common)
        match_rate = match_count / n_primary
        
        results['cluster_stability'][int(cluster_id)] = {
            'n_samples': int(n_primary),
            'best_match_cluster': int(most_common),
            'match_rate': float(match_rate),
            'stable': match_rate > 0.7
        }
    
    return results


def save_clustering_artifacts(linkage_matrix: np.ndarray,
                            clustering_info: Dict,
                            dendrogram_order: List[int],
                            output_dir: str,
                            robustness_results: Dict = None) -> Dict[str, str]:
    """
    Save clustering artifacts for reproducibility.
    
    Saves:
    - Linkage matrix (pickle format)
    - Dendrogram order (pickle format)
    - Clustering parameters and info (pickle format)
    - Robustness analysis results (if provided)
    
    Parameters:
    -----------
    linkage_matrix : np.ndarray
        Linkage matrix from hierarchical clustering
    clustering_info : dict
        Dictionary with clustering parameters and results
    dendrogram_order : list
        Order of leaves in the dendrogram
    output_dir : str
        Directory to save artifacts
    robustness_results : dict, optional
        Results from robustness check
    
    Returns:
    --------
    dict
        Dictionary with paths to saved artifacts
    """
    os.makedirs(output_dir, exist_ok=True)
    saved_paths = {}
    
    # Save linkage matrix
    linkage_path = os.path.join(output_dir, 'linkage_matrix.pkl')
    with open(linkage_path, 'wb') as f:
        pickle.dump(linkage_matrix, f)
    saved_paths['linkage_matrix'] = linkage_path
    print(f"   Saved linkage matrix to: {linkage_path}")
    
    # Save dendrogram order
    dendro_order_path = os.path.join(output_dir, 'dendrogram_order.pkl')
    with open(dendro_order_path, 'wb') as f:
        pickle.dump(dendrogram_order, f)
    saved_paths['dendrogram_order'] = dendro_order_path
    print(f"   Saved dendrogram order to: {dendro_order_path}")
    
    # Save clustering info
    info_path = os.path.join(output_dir, 'clustering_info.pkl')
    with open(info_path, 'wb') as f:
        pickle.dump(clustering_info, f)
    saved_paths['clustering_info'] = info_path
    print(f"   Saved clustering info to: {info_path}")
    
    # Save robustness results if provided
    if robustness_results is not None:
        robustness_path = os.path.join(output_dir, 'robustness_analysis.pkl')
        with open(robustness_path, 'wb') as f:
            pickle.dump(robustness_results, f)
        saved_paths['robustness_analysis'] = robustness_path
        print(f"   Saved robustness analysis to: {robustness_path}")
    
    return saved_paths


def run_clustering_pipeline(df: pd.DataFrame,
                           feature_cols: List[str],
                           n_clusters: int = DEFAULT_N_CLUSTERS,
                           linkage_method: str = LINKAGE_METHOD,
                           distance_metric: str = DISTANCE_METRIC_PRIMARY,
                           perform_robustness: bool = True,
                           output_dir: str = None) -> Tuple[pd.DataFrame, np.ndarray, Dict]:
    """
    Main clustering pipeline with explicit parameter control and reproducibility features.
    
    METHODOLOGY NOTES:
    -----------------
    1. Clustering is performed ONLY on resistance features (feature_cols).
       Metadata is explicitly excluded to prevent data leakage.
    
    2. Ward linkage is used by default because it minimizes within-cluster
       variance, producing compact clusters suitable for resistance phenotypes.
    
    3. The dendrogram structure is preserved for visualization consistency.
    
    4. Optional robustness check with Manhattan distance validates cluster stability.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe with encoded resistance values
    feature_cols : list
        List of feature column names (resistance data ONLY - no metadata)
    n_clusters : int
        Number of clusters to form (default: 5, based on elbow analysis)
    linkage_method : str
        Linkage method for hierarchical clustering (default: 'ward')
        - 'ward': Minimizes within-cluster variance (recommended)
        - 'complete': Maximum distance between clusters
        - 'average': Mean distance between clusters
        - 'single': Minimum distance between clusters
    distance_metric : str
        Distance metric (default: 'euclidean')
        Note: Ward linkage requires Euclidean distance
    perform_robustness : bool
        Whether to perform robustness check with Manhattan distance (default: True)
    output_dir : str, optional
        Directory to save clustering artifacts for reproducibility
    
    Returns:
    --------
    tuple
        (Dataframe with cluster labels, Linkage matrix, Clustering info dict)
    """
    print("=" * 70)
    print("PHASE 3.1: Structure Identification via Hierarchical Clustering")
    print("=" * 70)
    
    # Explicit parameter documentation
    print("\n>>> CLUSTERING PARAMETERS (for reproducibility):")
    print(f"    Linkage method: {linkage_method}")
    print(f"    Distance metric: {distance_metric}")
    print(f"    Number of clusters: {n_clusters}")
    print(f"    Cluster cut criterion: {CLUSTER_CUT_CRITERION}")
    
    # No-leakage assertion
    print("\n>>> DATA LEAKAGE PREVENTION:")
    print(f"    Using {len(feature_cols)} resistance features for clustering")
    metadata_cols = ['CODE', 'ISOLATE_ID', 'REGION', 'SITE', 'SAMPLE_SOURCE', 
                    'MDR_FLAG', 'MDR_CATEGORY', 'SPECIES']
    excluded = [c for c in metadata_cols if c in df.columns]
    print(f"    Metadata EXCLUDED from clustering: {excluded}")
    print("    ASSERTION: Clusters represent resistance phenotypes, NOT taxonomic groups")
    
    # Prepare data
    print("\n1. Preparing resistance data for clustering...")
    imputed_data, df_valid = prepare_clustering_data(df, feature_cols)
    
    # Perform clustering
    print(f"\n2. Performing hierarchical clustering...")
    print(f"   Method: {linkage_method} (minimizes within-cluster variance)")
    print(f"   Metric: {distance_metric}")
    linkage_matrix = perform_hierarchical_clustering(
        imputed_data, method=linkage_method, metric=distance_metric
    )
    
    # Generate dendrogram to get leaf order (for visualization consistency)
    import matplotlib.pyplot as plt
    fig_temp, ax_temp = plt.subplots()
    dendro_result = dendrogram(linkage_matrix, no_plot=True)
    dendrogram_order = dendro_result['leaves']
    plt.close(fig_temp)
    
    # Assign clusters with documented cut rule
    print(f"\n3. Assigning clusters (criterion: {CLUSTER_CUT_CRITERION})...")
    cluster_labels, cut_info = assign_clusters(linkage_matrix, n_clusters=n_clusters)
    print(f"   {cut_info['justification']}")
    
    # Add cluster labels to dataframe
    df_clustered = df_valid.copy()
    df_clustered['CLUSTER'] = cluster_labels
    
    # Analyze cluster distribution
    print("\n4. Cluster distribution:")
    cluster_dist = df_clustered['CLUSTER'].value_counts().sort_index()
    
    for cluster_id, count in cluster_dist.items():
        pct = (count / len(df_clustered)) * 100
        print(f"   C{cluster_id}: {count} isolates ({pct:.1f}%)")
    
    # Determine optimal clusters analysis
    print("\n5. Analyzing cluster quality metrics...")
    optimal_analysis = determine_optimal_clusters(linkage_matrix)
    
    # Perform robustness check
    robustness_results = None
    if perform_robustness:
        print(f"\n6. Performing robustness check with {DISTANCE_METRIC_ROBUSTNESS} distance...")
        robustness_results = perform_robustness_check(
            imputed_data, n_clusters, linkage_method, distance_metric
        )
        print(f"   Cluster stability: {robustness_results['stability_level']}")
        print(f"   Adjusted Rand Index: {robustness_results['agreement_score']['adjusted_rand_index']:.3f}")
        for interp in robustness_results['interpretation']:
            print(f"   {interp}")
    
    # Clustering info
    clustering_info = {
        'method': linkage_method,
        'method_justification': (
            "Ward linkage minimizes within-cluster variance, producing compact, "
            "spherical clusters suitable for identifying distinct resistance phenotypes."
        ),
        'metric': distance_metric,
        'metric_justification': (
            "Euclidean distance is the standard for numerical data and is required "
            "for Ward linkage. It measures straight-line distance in resistance space."
        ),
        'n_clusters': n_clusters,
        'cluster_cut_rule': cut_info,
        'total_isolates': len(df_clustered),
        'cluster_distribution': cluster_dist.to_dict(),
        'linkage_matrix_shape': linkage_matrix.shape,
        'dendrogram_order': dendrogram_order,
        'optimal_cluster_analysis': optimal_analysis,
        'robustness_check': robustness_results,
        'no_leakage_statement': (
            "Clustering was performed using ONLY resistance phenotype data. "
            "Metadata (species, region, environment, MDR status) was excluded "
            "from the clustering algorithm to prevent circular analysis. "
            "Any post-hoc associations with metadata are correlational findings."
        )
    }
    
    # Save artifacts if output directory provided
    if output_dir:
        print(f"\n7. Saving clustering artifacts for reproducibility...")
        saved_paths = save_clustering_artifacts(
            linkage_matrix, clustering_info, dendrogram_order, 
            output_dir, robustness_results
        )
        clustering_info['artifact_paths'] = saved_paths
    
    print(f"\n{'='*70}")
    print(f"CLUSTERING COMPLETE: {len(df_clustered)} isolates -> {n_clusters} clusters")
    print(f"{'='*70}")
    
    return df_clustered, linkage_matrix, clustering_info


def get_cluster_profiles(df_clustered: pd.DataFrame,
                        feature_cols: List[str]) -> pd.DataFrame:
    """
    Calculate mean resistance profile for each cluster.
    
    Parameters:
    -----------
    df_clustered : pd.DataFrame
        Dataframe with cluster labels
    feature_cols : list
        List of feature column names
    
    Returns:
    --------
    pd.DataFrame
        Mean resistance profile per cluster
    """
    existing_cols = [c for c in feature_cols if c in df_clustered.columns]
    
    cluster_profiles = df_clustered.groupby('CLUSTER')[existing_cols].mean()
    
    return cluster_profiles


def get_cluster_summary(df_clustered: pd.DataFrame,
                       feature_cols: List[str] = None,
                       metadata_cols: List[str] = None) -> Dict:
    """
    Get comprehensive summary statistics for each cluster.
    
    Creates the MANDATORY CLUSTER SUMMARY TABLE with:
    - Cluster ID (C1, C2, ...)
    - N isolates
    - Dominant species (%)
    - MDR %
    - Top resistant antibiotics
    - Major region
    - Major environment
    
    Parameters:
    -----------
    df_clustered : pd.DataFrame
        Dataframe with cluster labels
    feature_cols : list, optional
        Feature columns for identifying top resistant antibiotics
    metadata_cols : list, optional
        Metadata columns for cross-tabulation
    
    Returns:
    --------
    dict
        Comprehensive summary statistics per cluster
    """
    if feature_cols is None:
        feature_cols = [c for c in df_clustered.columns if c.endswith('_encoded')]
    
    summary = {}
    total_isolates = len(df_clustered)
    
    for cluster_id in sorted(df_clustered['CLUSTER'].unique()):
        cluster_df = df_clustered[df_clustered['CLUSTER'] == cluster_id]
        n_isolates = len(cluster_df)
        
        cluster_summary = {
            'cluster_label': f'C{cluster_id}',
            'n_isolates': n_isolates,
            'percentage': (n_isolates / total_isolates) * 100
        }
        
        # Species composition - find dominant species and top 3
        if 'ISOLATE_ID' in cluster_df.columns:
            species_counts = cluster_df['ISOLATE_ID'].value_counts()
            if len(species_counts) > 0:
                dominant_species = species_counts.index[0]
                dominant_count = species_counts.iloc[0]
                dominant_pct = (dominant_count / n_isolates) * 100
                cluster_summary['dominant_species'] = dominant_species
                cluster_summary['dominant_species_pct'] = dominant_pct
                cluster_summary['species_distribution'] = species_counts.to_dict()
                # Enhanced: Top 3 species composition
                top_3_species = []
                for i, (sp, count) in enumerate(species_counts.head(3).items()):
                    pct = (count / n_isolates) * 100
                    top_3_species.append(f"{sp} ({pct:.1f}%)")
                cluster_summary['top_3_species'] = top_3_species
        
        # MDR proportion
        if 'MDR_FLAG' in cluster_df.columns:
            mdr_count = cluster_df['MDR_FLAG'].sum()
            cluster_summary['mdr_count'] = int(mdr_count)
            cluster_summary['mdr_proportion'] = (mdr_count / n_isolates) * 100 if n_isolates > 0 else 0
        
        # Top resistant antibiotics
        if feature_cols:
            existing_cols = [c for c in feature_cols if c in cluster_df.columns]
            if existing_cols:
                # Calculate mean resistance per antibiotic (2 = Resistant)
                mean_resistance = cluster_df[existing_cols].mean()
                # Sort by resistance level (descending)
                top_resistant = mean_resistance.sort_values(ascending=False).head(5)
                cluster_summary['top_resistant_antibiotics'] = [
                    c.replace('_encoded', '') for c in top_resistant.index.tolist()
                ]
                cluster_summary['top_resistance_values'] = {
                    c.replace('_encoded', ''): float(v) 
                    for c, v in top_resistant.items()
                }
        
        # Regional distribution - find major region
        if 'REGION' in cluster_df.columns:
            region_counts = cluster_df['REGION'].value_counts()
            if len(region_counts) > 0:
                major_region = region_counts.index[0]
                region_pct = (region_counts.iloc[0] / n_isolates) * 100
                cluster_summary['major_region'] = major_region
                cluster_summary['major_region_pct'] = region_pct
                cluster_summary['regional_distribution'] = region_counts.to_dict()
        
        # Environmental distribution - use ENVIRONMENT column (FIX: was looking for SAMPLE_SOURCE)
        if 'ENVIRONMENT' in cluster_df.columns:
            env_counts = cluster_df['ENVIRONMENT'].value_counts()
            if len(env_counts) > 0:
                major_env = env_counts.index[0]
                env_pct = (env_counts.iloc[0] / n_isolates) * 100
                cluster_summary['major_environment'] = major_env
                cluster_summary['major_environment_pct'] = env_pct
                cluster_summary['environmental_distribution'] = env_counts.to_dict()
        
        # Enhanced: Barangay/Local Site distribution
        if 'LOCAL_SITE' in cluster_df.columns:
            local_counts = cluster_df['LOCAL_SITE'].value_counts()
            if len(local_counts) > 0:
                major_local = local_counts.index[0]
                local_pct = (local_counts.iloc[0] / n_isolates) * 100
                cluster_summary['major_barangay'] = major_local
                cluster_summary['major_barangay_pct'] = local_pct
                cluster_summary['barangay_distribution'] = local_counts.to_dict()
        
        # Enhanced: Detailed Sampling Source distribution
        if 'SAMPLING_SOURCE' in cluster_df.columns:
            source_counts = cluster_df['SAMPLING_SOURCE'].value_counts()
            if len(source_counts) > 0:
                major_source = source_counts.index[0]
                source_pct = (source_counts.iloc[0] / n_isolates) * 100
                cluster_summary['major_source'] = major_source
                cluster_summary['major_source_pct'] = source_pct
                cluster_summary['source_distribution'] = source_counts.to_dict()
        
        # Enhanced: Per-region barangay breakdown (to clarify mixed-region clusters)
        if 'REGION' in cluster_df.columns and 'LOCAL_SITE' in cluster_df.columns:
            per_region_barangay = {}
            for region in cluster_df['REGION'].unique():
                region_subset = cluster_df[cluster_df['REGION'] == region]
                if len(region_subset) > 0:
                    local_counts = region_subset['LOCAL_SITE'].value_counts()
                    if len(local_counts) > 0:
                        major_local = local_counts.index[0]
                        local_pct = (local_counts.iloc[0] / len(region_subset)) * 100
                        # Abbreviated region name for display
                        if 'BARMM' in str(region):
                            region_abbrev = 'BARMM'
                        elif 'Eastern' in str(region) or 'VIII' in str(region):
                            region_abbrev = 'EVS'
                        elif 'Central' in str(region) or 'III' in str(region):
                            region_abbrev = 'CLZ'
                        else:
                            region_abbrev = str(region)[:10]
                        per_region_barangay[region_abbrev] = f"{major_local} ({local_pct:.0f}%)"
            cluster_summary['barangay_by_region'] = per_region_barangay
        
        # Mean MAR index
        if 'MAR_INDEX_COMPUTED' in cluster_df.columns:
            cluster_summary['mean_mar_index'] = float(cluster_df['MAR_INDEX_COMPUTED'].mean())
        
        summary[int(cluster_id)] = cluster_summary
    
    return summary


def create_cluster_summary_table(df_clustered: pd.DataFrame,
                                feature_cols: List[str] = None) -> pd.DataFrame:
    """
    Create the ENHANCED MANDATORY CLUSTER SUMMARY TABLE as a DataFrame.
    
    Table format (enhanced with 9 columns):
    Cluster | N isolates | Species Composition | MDR % | Top Resistant Antibiotics | 
    Major Region | Major Barangay | Major Environment | Major Source
    
    Parameters:
    -----------
    df_clustered : pd.DataFrame
        Dataframe with cluster labels
    feature_cols : list, optional
        Feature columns for identifying top resistant antibiotics
    
    Returns:
    --------
    pd.DataFrame
        Enhanced summary table ready for display or export
    """
    summary = get_cluster_summary(df_clustered, feature_cols)
    
    table_data = []
    for cluster_id, info in sorted(summary.items()):
        # Enhanced species composition showing top 3
        if 'top_3_species' in info:
            species_comp = '; '.join(info['top_3_species'])
        elif 'dominant_species' in info:
            species_comp = f"{info['dominant_species']} ({info.get('dominant_species_pct', 0):.1f}%)"
        else:
            species_comp = 'N/A'
        
        row = {
            'Cluster': info.get('cluster_label', f'C{cluster_id}'),
            'N Isolates': info.get('n_isolates', 0),
            'Species Composition': species_comp,
            'MDR %': f"{info.get('mdr_proportion', 0):.1f}%" if 'mdr_proportion' in info else 'N/A',
            'Top Resistant Antibiotics': ', '.join(info.get('top_resistant_antibiotics', [])[:3]) if 'top_resistant_antibiotics' in info else 'N/A',
            'Major Region': f"{info.get('major_region', 'N/A')} ({info.get('major_region_pct', 0):.1f}%)" if 'major_region' in info else 'N/A',
            'Major Barangay': f"{info.get('major_barangay', 'N/A')} ({info.get('major_barangay_pct', 0):.1f}%)" if 'major_barangay' in info else 'N/A',
            'Major Environment': f"{info.get('major_environment', 'N/A')} ({info.get('major_environment_pct', 0):.1f}%)" if 'major_environment' in info else 'N/A',
            'Major Source': f"{info.get('major_source', 'N/A')} ({info.get('major_source_pct', 0):.1f}%)" if 'major_source' in info else 'N/A',
            'Barangay by Region': '; '.join([f"{k}: {v}" for k, v in info.get('barangay_by_region', {}).items()]) if 'barangay_by_region' in info else 'N/A'
        }
        table_data.append(row)
    
    return pd.DataFrame(table_data)


def create_environmental_distribution_table(df_clustered: pd.DataFrame,
                                           environment_col: str = 'ENVIRONMENT',
                                           normalize: bool = True) -> pd.DataFrame:
    """
    Create ENVIRONMENTAL DISTRIBUTION TABLE (Cluster × Environment contingency table).
    
    Reports proportions (not just counts) for clearer interpretation.
    
    Parameters:
    -----------
    df_clustered : pd.DataFrame
        Dataframe with cluster labels
    environment_col : str
        Column name for environmental/sample source (default: 'SAMPLE_SOURCE')
    normalize : bool
        If True, report proportions; if False, report counts (default: True)
    
    Returns:
    --------
    pd.DataFrame
        Contingency table with cluster × environment distribution
    """
    if environment_col not in df_clustered.columns:
        print(f"Warning: {environment_col} column not found in dataframe")
        return pd.DataFrame()
    
    if 'CLUSTER' not in df_clustered.columns:
        print("Warning: CLUSTER column not found in dataframe")
        return pd.DataFrame()
    
    # Create contingency table
    contingency = pd.crosstab(
        df_clustered['CLUSTER'], 
        df_clustered[environment_col],
        margins=True,
        margins_name='Total'
    )
    
    if normalize:
        # Calculate proportions (row-wise normalization)
        # Each row sums to 1.0 (100%)
        contingency_proportions = contingency.div(contingency['Total'], axis=0) * 100
        contingency_proportions = contingency_proportions.round(1)
        
        # Format as percentages
        for col in contingency_proportions.columns:
            if col != 'Total':
                contingency_proportions[col] = contingency_proportions[col].apply(lambda x: f"{x:.1f}%")
        
        # Add count in parentheses
        contingency_proportions['N (total)'] = contingency['Total']
        
        return contingency_proportions
    
    return contingency


def create_regional_distribution_table(df_clustered: pd.DataFrame,
                                      region_col: str = 'REGION',
                                      normalize: bool = True) -> pd.DataFrame:
    """
    Create REGIONAL DISTRIBUTION TABLE (Cluster × Region contingency table).
    
    Reports proportions (not just counts) for clearer interpretation.
    
    Parameters:
    -----------
    df_clustered : pd.DataFrame
        Dataframe with cluster labels
    region_col : str
        Column name for region (default: 'REGION')
    normalize : bool
        If True, report proportions; if False, report counts (default: True)
    
    Returns:
    --------
    pd.DataFrame
        Contingency table with cluster × region distribution
    """
    if region_col not in df_clustered.columns:
        print(f"Warning: {region_col} column not found in dataframe")
        return pd.DataFrame()
    
    if 'CLUSTER' not in df_clustered.columns:
        print("Warning: CLUSTER column not found in dataframe")
        return pd.DataFrame()
    
    # Create contingency table
    contingency = pd.crosstab(
        df_clustered['CLUSTER'], 
        df_clustered[region_col],
        margins=True,
        margins_name='Total'
    )
    
    if normalize:
        # Calculate proportions (row-wise normalization)
        contingency_proportions = contingency.div(contingency['Total'], axis=0) * 100
        contingency_proportions = contingency_proportions.round(1)
        
        # Format as percentages
        for col in contingency_proportions.columns:
            if col != 'Total':
                contingency_proportions[col] = contingency_proportions[col].apply(lambda x: f"{x:.1f}%")
        
        # Add count in parentheses
        contingency_proportions['N (total)'] = contingency['Total']
        
        return contingency_proportions
    
    return contingency


if __name__ == "__main__":
    from pathlib import Path
    
    project_root = Path(__file__).parent.parent.parent
    analysis_path = project_root / "data" / "processed" / "analysis_ready_dataset.csv"
    
    if analysis_path.exists():
        df = pd.read_csv(analysis_path)
        
        # Get encoded columns
        feature_cols = [c for c in df.columns if c.endswith('_encoded')]
        
        # Run clustering
        df_clustered, linkage_matrix, info = run_clustering_pipeline(
            df, feature_cols, n_clusters=4
        )
        
        # Get cluster summary
        summary = get_cluster_summary(df_clustered)
        
        print("\n" + "=" * 50)
        print("CLUSTER SUMMARY")
        print("=" * 50)
        
        for cluster_id, cluster_info in summary.items():
            print(f"\nCluster {cluster_id}:")
            print(f"  Isolates: {cluster_info['n_isolates']} ({cluster_info['percentage']:.1f}%)")
            if 'mdr_proportion' in cluster_info:
                print(f"  MDR proportion: {cluster_info['mdr_proportion']:.2%}")
            if 'mean_mar_index' in cluster_info:
                print(f"  Mean MAR index: {cluster_info['mean_mar_index']:.4f}")
        
        # Save clustered data
        clustered_path = project_root / "data" / "processed" / "clustered_dataset.csv"
        df_clustered.to_csv(clustered_path, index=False)
        print(f"\nClustered dataset saved to: {clustered_path}")
    else:
        print(f"Analysis-ready dataset not found at {analysis_path}")
        print("Run feature_engineering.py first.")
