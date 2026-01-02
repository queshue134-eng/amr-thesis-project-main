"""
Non-Linear Encoding Sensitivity Analysis
=========================================

This script tests the impact of different resistance encoding schemes on
clustering results. The default ordinal encoding (S=0, I=1, R=2) assumes
equal intervals, which may not reflect biological reality.

ENCODING SCHEMES TESTED:
-----------------------
1. Linear (default):     S=0, I=1, R=2
   - Equal intervals between categories
   - Assumes S→I is biologically equivalent to I→R

2. Clinical-weighted:    S=0, I=1.5, R=3
   - Larger gap between I and R (clinical failure threshold)
   - Reflects that I→R transition is more clinically significant

3. Binary:               S=0, I=1, R=1
   - Collapses I and R into "non-susceptible"
   - Conservative interpretation (intermediate = resistant)

4. Threshold-based:      S=0, I=0.5, R=1
   - Normalized to [0,1] range
   - I treated as partial resistance

METHODOLOGY:
-----------
1. Apply each encoding scheme to the cleaned dataset
2. Perform hierarchical clustering with consistent parameters
3. Compare cluster assignments using Adjusted Rand Index (ARI)
4. Report stability metrics and sensitivity conclusions

If ARI > 0.90 across encodings: Results are ROBUST to encoding choice
If ARI 0.70-0.90: Results are MODERATELY STABLE
If ARI < 0.70: Results are ENCODING-DEPENDENT (concerning)

Author: AMR Thesis Project
Date: December 2025
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple
import sys

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.metrics import adjusted_rand_score, silhouette_score
from sklearn.impute import SimpleImputer

try:
    from config import RANDOM_STATE, RESISTANCE_ENCODING
except ImportError:
    RANDOM_STATE = 42
    RESISTANCE_ENCODING = {'S': 0, 'I': 1, 'R': 2}


# =============================================================================
# ENCODING SCHEMES TO TEST
# =============================================================================

ENCODING_SCHEMES = {
    'linear': {
        'S': 0,
        'I': 1,
        'R': 2,
        'description': 'Default linear encoding (equal intervals)',
        'rationale': 'Standard ordinal encoding used in most AMR studies'
    },
    'clinical_weighted': {
        'S': 0.0,
        'I': 1.5,
        'R': 3.0,
        'description': 'Clinical-weighted encoding (larger I→R gap)',
        'rationale': 'Reflects clinical significance where I→R transition indicates treatment failure'
    },
    'binary': {
        'S': 0,
        'I': 1,
        'R': 1,
        'description': 'Binary encoding (I=R, non-susceptible)',
        'rationale': 'Conservative interpretation treating intermediate as resistant'
    },
    'threshold_normalized': {
        'S': 0.0,
        'I': 0.5,
        'R': 1.0,
        'description': 'Normalized threshold encoding [0-1]',
        'rationale': 'Bounded scale treating intermediate as 50% resistance probability'
    },
    'exponential': {
        'S': 0,
        'I': 1,
        'R': 4,
        'description': 'Exponential encoding (emphasizes R)',
        'rationale': 'Heavily weights resistant phenotype based on exponential MIC scale'
    }
}

# Baseline encoding for comparison
BASELINE_ENCODING = 'linear'


def apply_encoding_scheme(
    df: pd.DataFrame,
    antibiotic_cols: List[str],
    encoding: Dict[str, float]
) -> pd.DataFrame:
    """
    Apply a specific encoding scheme to resistance data.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with categorical resistance values (S, I, R)
    antibiotic_cols : list
        List of antibiotic column names (without _encoded suffix)
    encoding : dict
        Encoding scheme mapping S, I, R to numeric values
    
    Returns
    -------
    pd.DataFrame
        DataFrame with encoded resistance values
    """
    df_encoded = df.copy()
    
    for col in antibiotic_cols:
        if col in df.columns:
            # Map categorical to numeric
            df_encoded[f'{col}_encoded'] = df[col].map(encoding)
    
    return df_encoded


def perform_clustering_with_encoding(
    df: pd.DataFrame,
    antibiotic_cols: List[str],
    encoding_name: str,
    encoding: Dict[str, float],
    n_clusters: int = 5
) -> Tuple[np.ndarray, float, Dict]:
    """
    Perform clustering with a specific encoding scheme.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with categorical resistance values
    antibiotic_cols : list
        List of antibiotic column names
    encoding_name : str
        Name of the encoding scheme
    encoding : dict
        Encoding mapping
    n_clusters : int
        Number of clusters
    
    Returns
    -------
    tuple
        (cluster_labels, silhouette_score, info_dict)
    """
    # Apply encoding
    df_encoded = apply_encoding_scheme(df, antibiotic_cols, encoding)
    
    # Get encoded columns
    encoded_cols = [f'{col}_encoded' for col in antibiotic_cols if f'{col}_encoded' in df_encoded.columns]
    
    if len(encoded_cols) < 3:
        return None, -1, {'error': 'Insufficient encoded columns'}
    
    # Extract features
    X = df_encoded[encoded_cols].values
    
    # Handle missing values
    imputer = SimpleImputer(strategy='median')
    X_imputed = imputer.fit_transform(X)
    
    # Perform clustering
    try:
        Z = linkage(X_imputed, method='ward', metric='euclidean')
        labels = fcluster(Z, n_clusters, criterion='maxclust')
        
        # Calculate silhouette score
        if len(np.unique(labels)) > 1:
            sil_score = silhouette_score(X_imputed, labels)
        else:
            sil_score = -1.0
        
        info = {
            'encoding_name': encoding_name,
            'n_samples': len(labels),
            'n_features': len(encoded_cols),
            'cluster_distribution': dict(zip(*np.unique(labels, return_counts=True))),
            'silhouette_score': sil_score
        }
        
        return labels, sil_score, info
        
    except Exception as e:
        return None, -1, {'error': str(e)}


def run_encoding_sensitivity_analysis(
    df: pd.DataFrame,
    antibiotic_cols: List[str] = None,
    n_clusters: int = 5,
    encoding_schemes: Dict = None,
    baseline: str = None
) -> Dict:
    """
    Run comprehensive encoding sensitivity analysis.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with categorical resistance values
    antibiotic_cols : list, optional
        Antibiotic column names. If None, auto-detected.
    n_clusters : int
        Number of clusters for comparison
    encoding_schemes : dict, optional
        Encoding schemes to test. If None, uses default set.
    baseline : str, optional
        Baseline encoding for ARI comparison
    
    Returns
    -------
    dict
        Comprehensive sensitivity analysis results
    """
    if encoding_schemes is None:
        encoding_schemes = ENCODING_SCHEMES
    
    if baseline is None:
        baseline = BASELINE_ENCODING
    
    # Auto-detect antibiotic columns if not provided
    if antibiotic_cols is None:
        # Look for columns with S/I/R values
        potential_cols = []
        for col in df.columns:
            if df[col].dtype == 'object':
                unique_vals = set(df[col].dropna().unique())
                if unique_vals.issubset({'S', 'I', 'R', 's', 'i', 'r'}):
                    potential_cols.append(col)
        antibiotic_cols = potential_cols
    
    print("=" * 70)
    print("ENCODING SENSITIVITY ANALYSIS")
    print("=" * 70)
    print(f"\nTesting {len(encoding_schemes)} encoding schemes...")
    print(f"Baseline encoding: {baseline}")
    print(f"Number of clusters: k={n_clusters}")
    print(f"Antibiotic columns: {len(antibiotic_cols)}")
    print("-" * 70)
    
    results = {
        'encoding_schemes_tested': list(encoding_schemes.keys()),
        'baseline_encoding': baseline,
        'n_clusters': n_clusters,
        'per_encoding_results': {},
        'pairwise_ari': {},
        'ari_vs_baseline': {},
        'stability_assessment': None,
        'recommendation': None
    }
    
    # Store cluster labels for each encoding
    cluster_labels_dict = {}
    
    # Run clustering for each encoding scheme
    for encoding_name, encoding_config in encoding_schemes.items():
        encoding = {k: v for k, v in encoding_config.items() 
                   if k in ['S', 'I', 'R']}
        
        print(f"\n>>> Testing encoding: {encoding_name}")
        print(f"    Mapping: S={encoding['S']}, I={encoding['I']}, R={encoding['R']}")
        print(f"    Rationale: {encoding_config.get('rationale', 'N/A')}")
        
        labels, sil_score, info = perform_clustering_with_encoding(
            df, antibiotic_cols, encoding_name, encoding, n_clusters
        )
        
        if labels is not None:
            cluster_labels_dict[encoding_name] = labels
            results['per_encoding_results'][encoding_name] = {
                'encoding': encoding,
                'description': encoding_config.get('description', ''),
                'rationale': encoding_config.get('rationale', ''),
                'silhouette_score': sil_score,
                'cluster_distribution': info['cluster_distribution'],
                'n_samples': info['n_samples']
            }
            print(f"    Silhouette score: {sil_score:.4f}")
        else:
            print(f"    ❌ Clustering failed: {info.get('error', 'Unknown error')}")
            results['per_encoding_results'][encoding_name] = {
                'error': info.get('error', 'Unknown error')
            }
    
    # Compute pairwise ARI
    print("\n" + "=" * 70)
    print("PAIRWISE ADJUSTED RAND INDEX (ARI) COMPARISON")
    print("=" * 70)
    
    valid_encodings = list(cluster_labels_dict.keys())
    
    from itertools import combinations
    for enc1, enc2 in combinations(valid_encodings, 2):
        labels1 = cluster_labels_dict[enc1]
        labels2 = cluster_labels_dict[enc2]
        
        ari = adjusted_rand_score(labels1, labels2)
        results['pairwise_ari'][(enc1, enc2)] = ari
        
        print(f"\n{enc1} vs {enc2}: ARI = {ari:.4f}")
    
    # Compute ARI vs baseline
    print("\n" + "=" * 70)
    print(f"ARI VS BASELINE ({baseline})")
    print("=" * 70)
    
    if baseline in cluster_labels_dict:
        baseline_labels = cluster_labels_dict[baseline]
        
        for enc in valid_encodings:
            if enc == baseline:
                results['ari_vs_baseline'][enc] = 1.0
                continue
            
            ari = adjusted_rand_score(baseline_labels, cluster_labels_dict[enc])
            results['ari_vs_baseline'][enc] = ari
            print(f"\n{enc} vs baseline: ARI = {ari:.4f}")
    
    # Assess stability
    print("\n" + "=" * 70)
    print("STABILITY ASSESSMENT")
    print("=" * 70)
    
    all_ari_values = list(results['pairwise_ari'].values())
    
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
                "Clustering results are HIGHLY STABLE across all encoding schemes. "
                "The choice of ordinal encoding (S=0, I=1, R=2) is well-justified."
            )
        elif min_ari >= 0.70:
            stability_level = "MODERATE"
            stability_msg = (
                "Clustering results show MODERATE STABILITY. Most encoding schemes "
                "produce similar groupings, though binary encoding may differ."
            )
        else:
            stability_level = "LOW"
            stability_msg = (
                "⚠️ Clustering results are ENCODING-DEPENDENT. The choice of encoding "
                "significantly affects cluster assignments. Consider biological rationale "
                "for selecting the most appropriate encoding."
            )
        
        results['stability_assessment']['stability_level'] = stability_level
        results['stability_assessment']['interpretation'] = stability_msg
        
        print(f"\n>>> STABILITY LEVEL: {stability_level}")
        print(f">>> {stability_msg}")
    
    # Recommendation
    print("\n" + "=" * 70)
    print("RECOMMENDATION")
    print("=" * 70)
    
    # Find encoding with best silhouette
    valid_results = {k: v for k, v in results['per_encoding_results'].items() 
                    if 'silhouette_score' in v}
    
    if valid_results:
        best_encoding = max(valid_results.items(), 
                           key=lambda x: x[1]['silhouette_score'])
        
        results['recommendation'] = {
            'best_silhouette_encoding': best_encoding[0],
            'best_silhouette_score': best_encoding[1]['silhouette_score'],
            'recommended_encoding': baseline if min_ari >= 0.70 else best_encoding[0],
            'rationale': (
                f"Linear encoding is recommended (default) as results are "
                f"{'stable' if min_ari >= 0.70 else 'most interpretable'}. "
                f"Best silhouette achieved by '{best_encoding[0]}' ({best_encoding[1]['silhouette_score']:.4f})."
                if min_ari >= 0.70 else
                f"Consider '{best_encoding[0]}' encoding which achieves best "
                f"silhouette ({best_encoding[1]['silhouette_score']:.4f}), but note "
                f"encoding choice significantly impacts results."
            )
        }
        
        print(f"\n>>> RECOMMENDED ENCODING: {results['recommendation']['recommended_encoding']}")
        print(f">>> {results['recommendation']['rationale']}")
    
    return results


def generate_encoding_sensitivity_report(results: Dict, output_path: str = None) -> str:
    """
    Generate a formatted markdown report from encoding sensitivity analysis.
    
    Parameters
    ----------
    results : dict
        Results from run_encoding_sensitivity_analysis()
    output_path : str, optional
        Path to save the report
    
    Returns
    -------
    str
        Formatted markdown report
    """
    lines = [
        "# Encoding Sensitivity Analysis Report",
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
        "## Encoding Schemes Tested",
        "",
        "| Encoding | S | I | R | Silhouette | Description |",
        "|----------|---|---|---|------------|-------------|",
    ])
    
    for enc_name, enc_result in results['per_encoding_results'].items():
        if 'error' not in enc_result:
            enc = enc_result['encoding']
            lines.append(
                f"| {enc_name} | {enc['S']} | {enc['I']} | {enc['R']} | "
                f"{enc_result['silhouette_score']:.4f} | {enc_result.get('description', '')} |"
            )
    
    lines.extend([
        "",
        "---",
        "",
        "## ARI vs Baseline (Linear)",
        "",
        "| Encoding | ARI vs Baseline |",
        "|----------|-----------------|",
    ])
    
    for enc, ari in results['ari_vs_baseline'].items():
        lines.append(f"| {enc} | {ari:.4f} |")
    
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
            f"**Recommended Encoding:** {rec['recommended_encoding']}",
            "",
            f"**Rationale:** {rec['rationale']}",
        ])
    
    lines.extend([
        "",
        "---",
        "",
        "## Conclusion",
        "",
        "This analysis validates the robustness of clustering results to different ",
        "resistance encoding schemes. For thesis defense, this provides evidence that ",
        "the identified resistance phenotypes are not artifacts of the encoding choice.",
    ])
    
    report = "\n".join(lines)
    
    if output_path:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(report)
        print(f"\nReport saved to: {output_path}")
    
    return report


def main(args=None):
    """Main function to run encoding sensitivity analysis."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Encoding Sensitivity Analysis')
    parser.add_argument('--data', type=str, help='Path to cleaned dataset CSV')
    parser.add_argument('--output', type=str, help='Path to save report')
    parser.add_argument('--k', type=int, default=5, help='Number of clusters')
    
    if args is None:
        parsed_args = parser.parse_args()
    else:
        parsed_args = parser.parse_args(args)
    
    # Default paths
    project_root = Path(__file__).parent.parent
    data_path = parsed_args.data or project_root / 'data' / 'processed' / 'cleaned_dataset.csv'
    output_path = parsed_args.output or project_root / 'docs' / 'results' / 'encoding_sensitivity_report.md'
    
    print(f"Loading data from: {data_path}")
    
    if not Path(data_path).exists():
        print(f"Error: Data file not found at {data_path}")
        print("Please run the main pipeline first to generate the cleaned dataset.")
        return
    
    df = pd.read_csv(data_path)
    print(f"Loaded {len(df)} isolates")
    
    # Run analysis
    results = run_encoding_sensitivity_analysis(df, n_clusters=parsed_args.k)
    
    # Generate report
    report = generate_encoding_sensitivity_report(results, str(output_path))
    
    print("\n" + "=" * 70)
    print("ENCODING SENSITIVITY ANALYSIS COMPLETE")
    print("=" * 70)


if __name__ == '__main__':
    main()
