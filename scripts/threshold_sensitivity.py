"""
Threshold Sensitivity Analysis Script
=====================================

This script performs sensitivity analysis on data cleaning thresholds to validate
the robustness of clustering results. It addresses the methodological requirement
to justify the choice of min_antibiotic_coverage (70%) and max_isolate_missing (30%)
thresholds.

METHODOLOGY:
-----------
1. Test multiple threshold combinations: (50%, 60%, 70%, 80%) for antibiotic coverage
2. For each combination:
   - Apply cleaning with specified thresholds
   - Perform hierarchical clustering (Ward's linkage, k=5)
   - Record cluster assignments
3. Compare cluster stability across thresholds using Adjusted Rand Index (ARI)
4. Generate comprehensive report with recommendations

INTERPRETATION:
--------------
- ARI > 0.90: High stability - results robust to threshold choice
- ARI 0.70-0.90: Moderate stability - results reasonably robust
- ARI < 0.70: Low stability - results are threshold-dependent (concerning)

Reference:
---------
Hubert, L., & Arabie, P. (1985). Comparing partitions.
Journal of Classification, 2(1), 193-218.

Author: AMR Thesis Project
Date: December 2025
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from itertools import combinations
import warnings
import sys

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.metrics import adjusted_rand_score, silhouette_score
from sklearn.impute import SimpleImputer

try:
    from config import (
        RANDOM_STATE, CLUSTERING_CONFIG, 
        MIN_ANTIBIOTIC_COVERAGE, MAX_ISOLATE_MISSING
    )
    from preprocessing.data_cleaning import (
        filter_antibiotics_by_coverage,
        remove_isolates_with_excessive_missing,
        compute_antibiotic_test_coverage
    )
except ImportError:
    # Fallback defaults
    RANDOM_STATE = 42
    MIN_ANTIBIOTIC_COVERAGE = 70.0
    MAX_ISOLATE_MISSING = 30.0


# =============================================================================
# THRESHOLD COMBINATIONS TO TEST
# =============================================================================

# Primary threshold pairs: (antibiotic_coverage%, isolate_missing%)
THRESHOLD_PAIRS = [
    (50.0, 40.0),  # Lenient
    (60.0, 35.0),  # Moderate-lenient
    (70.0, 30.0),  # Default (current implementation)
    (80.0, 25.0),  # Strict
    (90.0, 20.0),  # Very strict
]

# Baseline threshold (current implementation)
BASELINE_THRESHOLD = (70.0, 30.0)


def compute_cleaning_impact(
    df: pd.DataFrame,
    antibiotic_cols: List[str],
    ab_coverage_threshold: float,
    isolate_missing_threshold: float
) -> Tuple[pd.DataFrame, Dict]:
    """
    Apply cleaning thresholds and return cleaned data with impact statistics.
    
    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe with resistance data
    antibiotic_cols : list
        List of antibiotic column names
    ab_coverage_threshold : float
        Minimum coverage percentage for antibiotics (0-100)
    isolate_missing_threshold : float
        Maximum missing percentage for isolates (0-100)
    
    Returns
    -------
    tuple
        (Cleaned DataFrame, Impact statistics dict)
    """
    initial_isolates = len(df)
    initial_antibiotics = len(antibiotic_cols)
    
    # Step 1: Filter antibiotics by coverage
    retained_abs, excluded_abs, coverage_stats = filter_antibiotics_by_coverage(
        df, antibiotic_cols, min_coverage=ab_coverage_threshold
    )
    
    # Step 2: Remove isolates with excessive missing data
    df_clean, removed_count, removed_info = remove_isolates_with_excessive_missing(
        df, retained_abs, max_missing_pct=isolate_missing_threshold
    )
    
    impact = {
        'ab_coverage_threshold': ab_coverage_threshold,
        'isolate_missing_threshold': isolate_missing_threshold,
        'initial_isolates': initial_isolates,
        'final_isolates': len(df_clean),
        'isolates_removed': initial_isolates - len(df_clean),
        'isolate_retention_pct': (len(df_clean) / initial_isolates) * 100,
        'initial_antibiotics': initial_antibiotics,
        'final_antibiotics': len(retained_abs),
        'antibiotics_removed': initial_antibiotics - len(retained_abs),
        'antibiotic_retention_pct': (len(retained_abs) / initial_antibiotics) * 100,
        'retained_antibiotics': retained_abs,
        'excluded_antibiotics': excluded_abs
    }
    
    return df_clean, retained_abs, impact


def perform_clustering(
    df: pd.DataFrame,
    feature_cols: List[str],
    n_clusters: int = 5,
    linkage_method: str = 'ward',
    distance_metric: str = 'euclidean'
) -> Tuple[np.ndarray, float]:
    """
    Perform hierarchical clustering and return cluster labels with silhouette score.
    
    Parameters
    ----------
    df : pd.DataFrame
        Cleaned dataframe
    feature_cols : list
        List of feature column names
    n_clusters : int
        Number of clusters to form
    linkage_method : str
        Linkage method for hierarchical clustering
    distance_metric : str
        Distance metric
    
    Returns
    -------
    tuple
        (Cluster labels array, Silhouette score)
    """
    # Extract features and handle missing values
    existing_cols = [c for c in feature_cols if c in df.columns]
    X = df[existing_cols].values
    
    # Impute missing values with median
    imputer = SimpleImputer(strategy='median')
    X_imputed = imputer.fit_transform(X)
    
    # Perform hierarchical clustering
    Z = linkage(X_imputed, method=linkage_method, metric=distance_metric)
    labels = fcluster(Z, n_clusters, criterion='maxclust')
    
    # Calculate silhouette score
    if len(np.unique(labels)) > 1:
        sil_score = silhouette_score(X_imputed, labels)
    else:
        sil_score = -1.0  # Invalid clustering
    
    return labels, sil_score


def run_sensitivity_analysis(
    df: pd.DataFrame,
    antibiotic_cols: List[str],
    threshold_pairs: List[Tuple[float, float]] = None,
    n_clusters: int = 5,
    baseline_threshold: Tuple[float, float] = None
) -> Dict:
    """
    Run comprehensive threshold sensitivity analysis.
    
    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe with resistance data
    antibiotic_cols : list
        List of antibiotic column names
    threshold_pairs : list of tuples
        List of (ab_coverage, isolate_missing) threshold pairs to test
    n_clusters : int
        Number of clusters for comparison
    baseline_threshold : tuple
        Baseline threshold pair for ARI comparison
    
    Returns
    -------
    dict
        Comprehensive sensitivity analysis results
    """
    if threshold_pairs is None:
        threshold_pairs = THRESHOLD_PAIRS
    
    if baseline_threshold is None:
        baseline_threshold = BASELINE_THRESHOLD
    
    print("=" * 70)
    print("THRESHOLD SENSITIVITY ANALYSIS")
    print("=" * 70)
    print(f"\nTesting {len(threshold_pairs)} threshold combinations...")
    print(f"Baseline threshold: {baseline_threshold}")
    print(f"Number of clusters: k={n_clusters}")
    print("-" * 70)
    
    results = {
        'threshold_pairs_tested': threshold_pairs,
        'baseline_threshold': baseline_threshold,
        'n_clusters': n_clusters,
        'per_threshold_results': [],
        'pairwise_ari': {},
        'ari_vs_baseline': {},
        'stability_assessment': None,
        'recommendation': None
    }
    
    # Store cluster labels for each threshold
    cluster_labels_dict = {}
    valid_indices_dict = {}
    
    # Run clustering for each threshold pair
    for ab_thresh, iso_thresh in threshold_pairs:
        print(f"\n>>> Testing threshold: AB coverage ≥{ab_thresh}%, Isolate missing ≤{iso_thresh}%")
        
        # Apply cleaning
        df_clean, retained_abs, impact = compute_cleaning_impact(
            df, antibiotic_cols, ab_thresh, iso_thresh
        )
        
        print(f"    Retained: {impact['final_isolates']}/{impact['initial_isolates']} isolates "
              f"({impact['isolate_retention_pct']:.1f}%)")
        print(f"    Retained: {impact['final_antibiotics']}/{impact['initial_antibiotics']} antibiotics "
              f"({impact['antibiotic_retention_pct']:.1f}%)")
        
        if len(df_clean) < n_clusters * 2:
            print(f"    ⚠️ Insufficient samples for clustering (need at least {n_clusters * 2})")
            results['per_threshold_results'].append({
                'threshold': (ab_thresh, iso_thresh),
                'impact': impact,
                'clustering_valid': False,
                'error': 'Insufficient samples'
            })
            continue
        
        # Get encoded columns for clustering
        encoded_cols = [c + '_encoded' if not c.endswith('_encoded') else c 
                       for c in retained_abs]
        existing_encoded = [c for c in encoded_cols if c in df_clean.columns]
        
        if len(existing_encoded) < 3:
            print(f"    ⚠️ Insufficient features for clustering (need at least 3)")
            results['per_threshold_results'].append({
                'threshold': (ab_thresh, iso_thresh),
                'impact': impact,
                'clustering_valid': False,
                'error': 'Insufficient features'
            })
            continue
        
        # Perform clustering
        try:
            labels, sil_score = perform_clustering(
                df_clean, existing_encoded, n_clusters
            )
            
            print(f"    Silhouette score: {sil_score:.4f}")
            
            # Store results
            threshold_key = (ab_thresh, iso_thresh)
            cluster_labels_dict[threshold_key] = labels
            valid_indices_dict[threshold_key] = df_clean.index.tolist()
            
            # Cluster distribution
            unique, counts = np.unique(labels, return_counts=True)
            cluster_dist = dict(zip(unique.astype(int), counts.astype(int)))
            
            results['per_threshold_results'].append({
                'threshold': threshold_key,
                'impact': impact,
                'clustering_valid': True,
                'silhouette_score': sil_score,
                'cluster_distribution': cluster_dist,
                'n_samples': len(labels)
            })
            
        except Exception as e:
            print(f"    ❌ Clustering failed: {str(e)}")
            results['per_threshold_results'].append({
                'threshold': (ab_thresh, iso_thresh),
                'impact': impact,
                'clustering_valid': False,
                'error': str(e)
            })
    
    # Compute pairwise ARI between all valid threshold pairs
    print("\n" + "=" * 70)
    print("PAIRWISE ADJUSTED RAND INDEX (ARI) COMPARISON")
    print("=" * 70)
    
    valid_thresholds = list(cluster_labels_dict.keys())
    
    for thresh1, thresh2 in combinations(valid_thresholds, 2):
        # Find common indices
        indices1 = set(valid_indices_dict[thresh1])
        indices2 = set(valid_indices_dict[thresh2])
        common_indices = list(indices1.intersection(indices2))
        
        if len(common_indices) < n_clusters:
            print(f"\n{thresh1} vs {thresh2}: Insufficient common samples")
            continue
        
        # Get labels for common indices
        idx_map1 = {idx: i for i, idx in enumerate(valid_indices_dict[thresh1])}
        idx_map2 = {idx: i for i, idx in enumerate(valid_indices_dict[thresh2])}
        
        labels1 = [cluster_labels_dict[thresh1][idx_map1[idx]] for idx in common_indices]
        labels2 = [cluster_labels_dict[thresh2][idx_map2[idx]] for idx in common_indices]
        
        ari = adjusted_rand_score(labels1, labels2)
        
        results['pairwise_ari'][(thresh1, thresh2)] = {
            'ari': ari,
            'n_common_samples': len(common_indices)
        }
        
        print(f"\n{thresh1} vs {thresh2}:")
        print(f"  ARI = {ari:.4f} (n={len(common_indices)} common samples)")
    
    # Compute ARI vs baseline for all thresholds
    print("\n" + "=" * 70)
    print(f"ARI VS BASELINE ({baseline_threshold})")
    print("=" * 70)
    
    if baseline_threshold in cluster_labels_dict:
        baseline_labels = cluster_labels_dict[baseline_threshold]
        baseline_indices = set(valid_indices_dict[baseline_threshold])
        baseline_idx_map = {idx: i for i, idx in enumerate(valid_indices_dict[baseline_threshold])}
        
        for thresh in valid_thresholds:
            if thresh == baseline_threshold:
                results['ari_vs_baseline'][thresh] = {'ari': 1.0, 'is_baseline': True}
                continue
            
            # Find common indices
            thresh_indices = set(valid_indices_dict[thresh])
            common_indices = list(baseline_indices.intersection(thresh_indices))
            
            if len(common_indices) < n_clusters:
                print(f"\n{thresh}: Insufficient common samples with baseline")
                continue
            
            thresh_idx_map = {idx: i for i, idx in enumerate(valid_indices_dict[thresh])}
            
            labels_baseline = [baseline_labels[baseline_idx_map[idx]] for idx in common_indices]
            labels_thresh = [cluster_labels_dict[thresh][thresh_idx_map[idx]] for idx in common_indices]
            
            ari = adjusted_rand_score(labels_baseline, labels_thresh)
            
            results['ari_vs_baseline'][thresh] = {
                'ari': ari,
                'n_common_samples': len(common_indices),
                'is_baseline': False
            }
            
            print(f"\n{thresh} vs baseline:")
            print(f"  ARI = {ari:.4f}")
    
    # Assess overall stability
    print("\n" + "=" * 70)
    print("STABILITY ASSESSMENT")
    print("=" * 70)
    
    all_ari_values = [v['ari'] for v in results['pairwise_ari'].values()]
    
    if all_ari_values:
        mean_ari = np.mean(all_ari_values)
        min_ari = np.min(all_ari_values)
        max_ari = np.max(all_ari_values)
        std_ari = np.std(all_ari_values)
        
        results['stability_assessment'] = {
            'mean_ari': mean_ari,
            'min_ari': min_ari,
            'max_ari': max_ari,
            'std_ari': std_ari,
            'n_comparisons': len(all_ari_values)
        }
        
        print(f"\nMean pairwise ARI: {mean_ari:.4f}")
        print(f"Range: [{min_ari:.4f}, {max_ari:.4f}]")
        print(f"Std deviation: {std_ari:.4f}")
        
        # Determine stability level
        if min_ari >= 0.90:
            stability_level = "HIGH"
            stability_msg = (
                "Clustering results are HIGHLY STABLE across all tested thresholds. "
                "The choice of 70%/30% is well-justified as results are threshold-independent."
            )
        elif min_ari >= 0.70:
            stability_level = "MODERATE"
            stability_msg = (
                "Clustering results show MODERATE STABILITY. Most threshold choices "
                "produce similar groupings, though some variability exists at extremes."
            )
        else:
            stability_level = "LOW"
            stability_msg = (
                "⚠️ Clustering results show LOW STABILITY - they are THRESHOLD-DEPENDENT. "
                "Consider conducting further investigation to determine optimal thresholds."
            )
        
        results['stability_assessment']['stability_level'] = stability_level
        results['stability_assessment']['interpretation'] = stability_msg
        
        print(f"\n>>> STABILITY LEVEL: {stability_level}")
        print(f">>> {stability_msg}")
    
    # Generate recommendation
    print("\n" + "=" * 70)
    print("RECOMMENDATION")
    print("=" * 70)
    
    # Find threshold with best silhouette score
    valid_results = [r for r in results['per_threshold_results'] if r.get('clustering_valid')]
    
    if valid_results:
        best_silhouette = max(valid_results, key=lambda x: x.get('silhouette_score', -1))
        best_retention = max(valid_results, key=lambda x: x['impact']['isolate_retention_pct'])
        
        # Balanced recommendation
        # Score = 0.5 * normalized_silhouette + 0.3 * normalized_retention + 0.2 * stability_to_baseline
        
        silhouette_scores = [r['silhouette_score'] for r in valid_results if 'silhouette_score' in r]
        retention_pcts = [r['impact']['isolate_retention_pct'] for r in valid_results]
        
        sil_min, sil_max = min(silhouette_scores), max(silhouette_scores)
        ret_min, ret_max = min(retention_pcts), max(retention_pcts)
        
        recommendations = []
        for r in valid_results:
            thresh = r['threshold']
            sil = r.get('silhouette_score', 0)
            ret = r['impact']['isolate_retention_pct']
            
            # Normalize scores
            norm_sil = (sil - sil_min) / (sil_max - sil_min) if sil_max > sil_min else 0.5
            norm_ret = (ret - ret_min) / (ret_max - ret_min) if ret_max > ret_min else 0.5
            
            # Get ARI vs baseline
            ari_baseline = results['ari_vs_baseline'].get(thresh, {}).get('ari', 0)
            
            # Combined score
            combined_score = 0.5 * norm_sil + 0.3 * norm_ret + 0.2 * ari_baseline
            
            recommendations.append({
                'threshold': thresh,
                'silhouette': sil,
                'retention_pct': ret,
                'ari_vs_baseline': ari_baseline,
                'combined_score': combined_score
            })
        
        recommendations.sort(key=lambda x: x['combined_score'], reverse=True)
        
        results['recommendation'] = {
            'ranked_thresholds': recommendations,
            'optimal_threshold': recommendations[0]['threshold'],
            'rationale': (
                f"Based on combined analysis of silhouette score (clustering quality), "
                f"data retention, and stability vs. baseline, the optimal threshold is "
                f"{recommendations[0]['threshold']}."
            )
        }
        
        print(f"\n>>> OPTIMAL THRESHOLD: {recommendations[0]['threshold']}")
        print(f"    Silhouette: {recommendations[0]['silhouette']:.4f}")
        print(f"    Retention: {recommendations[0]['retention_pct']:.1f}%")
        print(f"    Combined Score: {recommendations[0]['combined_score']:.4f}")
        
        print("\n>>> ALL THRESHOLDS RANKED:")
        for i, rec in enumerate(recommendations, 1):
            print(f"    {i}. {rec['threshold']}: score={rec['combined_score']:.4f}")
    
    return results


def generate_sensitivity_report(results: Dict, output_path: str = None) -> str:
    """
    Generate a formatted markdown report from sensitivity analysis results.
    
    Parameters
    ----------
    results : dict
        Results from run_sensitivity_analysis()
    output_path : str, optional
        Path to save the report
    
    Returns
    -------
    str
        Formatted markdown report
    """
    lines = [
        "# Threshold Sensitivity Analysis Report",
        "",
        f"**Generated:** {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "## Executive Summary",
        "",
    ]
    
    if results.get('stability_assessment'):
        sa = results['stability_assessment']
        lines.extend([
            f"**Stability Level:** {sa['stability_level']}",
            "",
            f"**Mean Pairwise ARI:** {sa['mean_ari']:.4f}",
            f"**ARI Range:** [{sa['min_ari']:.4f}, {sa['max_ari']:.4f}]",
            "",
            f"**Interpretation:** {sa['interpretation']}",
            "",
        ])
    
    lines.extend([
        "---",
        "",
        "## Threshold Combinations Tested",
        "",
        "| AB Coverage (%) | Isolate Missing (%) | Isolates Retained | Antibiotics Retained | Silhouette |",
        "|-----------------|---------------------|-------------------|---------------------|------------|",
    ])
    
    for r in results['per_threshold_results']:
        if r.get('clustering_valid'):
            lines.append(
                f"| {r['threshold'][0]:.0f} | {r['threshold'][1]:.0f} | "
                f"{r['impact']['final_isolates']} ({r['impact']['isolate_retention_pct']:.1f}%) | "
                f"{r['impact']['final_antibiotics']} ({r['impact']['antibiotic_retention_pct']:.1f}%) | "
                f"{r['silhouette_score']:.4f} |"
            )
        else:
            lines.append(
                f"| {r['threshold'][0]:.0f} | {r['threshold'][1]:.0f} | "
                f"FAILED: {r.get('error', 'Unknown')} | - | - |"
            )
    
    lines.extend([
        "",
        "---",
        "",
        "## Pairwise ARI Matrix",
        "",
    ])
    
    # Create ARI matrix table
    valid_thresholds = sorted(set(
        t for pair in results['pairwise_ari'].keys() for t in pair
    ))
    
    if valid_thresholds:
        header = "| Threshold | " + " | ".join([f"{t}" for t in valid_thresholds]) + " |"
        separator = "|" + "|".join(["---"] * (len(valid_thresholds) + 1)) + "|"
        lines.extend([header, separator])
        
        for t1 in valid_thresholds:
            row = [f"| {t1}"]
            for t2 in valid_thresholds:
                if t1 == t2:
                    row.append(" 1.000 ")
                else:
                    key = (t1, t2) if (t1, t2) in results['pairwise_ari'] else (t2, t1)
                    if key in results['pairwise_ari']:
                        row.append(f" {results['pairwise_ari'][key]['ari']:.3f} ")
                    else:
                        row.append(" - ")
            row.append("|")
            lines.append("|".join(row))
    
    lines.extend([
        "",
        "---",
        "",
        "## Recommendation",
        "",
    ])
    
    if results.get('recommendation'):
        rec = results['recommendation']
        lines.extend([
            f"**Optimal Threshold:** {rec['optimal_threshold']}",
            "",
            f"**Rationale:** {rec['rationale']}",
            "",
            "### Ranked Thresholds",
            "",
            "| Rank | Threshold | Silhouette | Retention (%) | Combined Score |",
            "|------|-----------|------------|---------------|----------------|",
        ])
        
        for i, r in enumerate(rec['ranked_thresholds'], 1):
            lines.append(
                f"| {i} | {r['threshold']} | {r['silhouette']:.4f} | "
                f"{r['retention_pct']:.1f} | {r['combined_score']:.4f} |"
            )
    
    lines.extend([
        "",
        "---",
        "",
        "## Conclusion",
        "",
        "This sensitivity analysis validates the choice of data cleaning thresholds by demonstrating ",
        "that clustering results are robust (or sensitive) to threshold variations. ",
        "",
        "For thesis defense, this analysis provides evidence-based justification for the selected ",
        "thresholds and documents the stability of the identified resistance phenotypes.",
    ])
    
    report = "\n".join(lines)
    
    if output_path:
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(report)
        print(f"\nReport saved to: {output_path}")
    
    return report


def main(args=None):
    """
    Main function to run sensitivity analysis on the AMR dataset.
    """
    import argparse
    
    parser = argparse.ArgumentParser(description='Threshold Sensitivity Analysis')
    parser.add_argument('--data', type=str, help='Path to encoded dataset CSV')
    parser.add_argument('--output', type=str, help='Path to save report')
    parser.add_argument('--k', type=int, default=5, help='Number of clusters')
    
    if args is None:
        parsed_args = parser.parse_args()
    else:
        parsed_args = parser.parse_args(args)
    
    # Default paths
    project_root = Path(__file__).parent.parent
    data_path = parsed_args.data or project_root / 'data' / 'processed' / 'encoded_dataset.csv'
    output_path = parsed_args.output or project_root / 'docs' / 'results' / 'threshold_sensitivity_report.md'
    
    print(f"Loading data from: {data_path}")
    
    if not Path(data_path).exists():
        print(f"Error: Data file not found at {data_path}")
        print("Please run the main pipeline first to generate the encoded dataset.")
        return
    
    df = pd.read_csv(data_path)
    print(f"Loaded {len(df)} isolates")
    
    # Detect antibiotic columns (those ending with _encoded)
    antibiotic_cols = [c.replace('_encoded', '') for c in df.columns if c.endswith('_encoded')]
    print(f"Detected {len(antibiotic_cols)} antibiotic columns")
    
    # Run sensitivity analysis
    results = run_sensitivity_analysis(
        df, antibiotic_cols, 
        threshold_pairs=THRESHOLD_PAIRS,
        n_clusters=parsed_args.k,
        baseline_threshold=BASELINE_THRESHOLD
    )
    
    # Generate report
    output_path.parent.mkdir(parents=True, exist_ok=True)
    report = generate_sensitivity_report(results, str(output_path))
    
    print("\n" + "=" * 70)
    print("SENSITIVITY ANALYSIS COMPLETE")
    print("=" * 70)


if __name__ == '__main__':
    main()
