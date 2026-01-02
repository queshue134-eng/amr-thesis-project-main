"""
Antibiotic Clustering Analysis for AMR Thesis Project
======================================================

This script clusters ANTIBIOTICS (not isolates) based on their co-resistance
patterns using the phi coefficient matrix from co-resistance analysis.

Purpose:
- Identify groups of antibiotics that tend to be co-resisted together
- Discover potential plasmid-linked resistance gene clusters
- Provide biological interpretation of antibiotic "families" by resistance pattern

Output files:
- data/processed/figures/antibiotic_clusters.csv (cluster assignments)
- data/processed/figures/antibiotic_dendrogram.png (dendrogram visualization)
- data/processed/figures/antibiotic_heatmap.png (clustered heatmap)
"""

import pandas as pd
import numpy as np
import sys
from pathlib import Path
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import squareform
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for pipeline execution
import matplotlib.pyplot as plt
import seaborn as sns

# Add src to path for console import
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))
try:
    from utils.console import console, Colors
except ImportError:
    class FallbackConsole:
        def header(self, t, s=None): print(f"\n{'='*70}\n{t}\n{'='*70}")
        def subheader(self, t): print(f"\n  -- {t} --")
        def step(self, c, t, m, d=None): print(f"[{c}/{t}] {m}")
        def success(self, m): print(f"  ✓ {m}")
        def warning(self, m): print(f"  ⚠ {m}")
        def info(self, m): print(f"  → {m}")
        def kv(self, k, v, i=2): print(f"{'  '*i}{k}: {v}")
        def bullet(self, m, i=2): print(f"{'  '*i}• {m}")
        def complete(self, m=None, s=None): print(f"\n{'='*70}\n{m or 'Complete'}\n{'='*70}")
    console = FallbackConsole()

# Antibiotic class mapping for biological interpretation
ANTIBIOTIC_CLASSES = {
    'AM': 'Penicillins',
    'AMC': 'Penicillins',
    'AN': 'Aminoglycosides',
    'C': 'Phenicols',
    'CF': 'Cephalosporins-1st',
    'CFO': 'Cephalosporins-Vet',
    'CFT': 'Cephalosporins-Vet',
    'CN': 'Cephalosporins-1st',
    'CPD': 'Cephalosporins-3rd',
    'CPT': 'Cephalosporins-5th',
    'CTX': 'Cephalosporins-3rd',
    'CZA': 'Cephalosporins-3rd',
    'DO': 'Tetracyclines',
    'ENR': 'Fluoroquinolones',
    'FT': 'Nitrofurans',
    'GM': 'Aminoglycosides',
    'IPM': 'Carbapenems',
    'MRB': 'Fluoroquinolones',
    'N': 'Aminoglycosides',
    'PRA': 'Fluoroquinolones',
    'SXT': 'Folate-Inhibitors',
    'TE': 'Tetracyclines'
}


def load_coresistance_matrix(data_path: str = None) -> pd.DataFrame:
    """Load the phi coefficient co-resistance matrix."""
    if data_path is None:
        project_root = Path(__file__).parent.parent
        data_path = project_root / "data" / "processed" / "figures" / "coresistance_matrix.csv"
    
    phi_matrix = pd.read_csv(data_path, index_col=0)
    console.info(f"Loaded co-resistance matrix: {phi_matrix.shape[0]} antibiotics")
    return phi_matrix


def convert_similarity_to_distance(phi_matrix: pd.DataFrame) -> np.ndarray:
    """
    Convert phi coefficient (similarity) to distance for clustering.
    
    Distance = 1 - phi (since phi is a correlation-like measure from 0 to 1)
    """
    # Ensure matrix is symmetric and diagonal is 1
    distance_matrix = 1 - phi_matrix.values
    
    # Set diagonal to 0 (self-distance)
    np.fill_diagonal(distance_matrix, 0)
    
    # Ensure non-negative (phi should be 0-1, but handle edge cases)
    distance_matrix = np.maximum(distance_matrix, 0)
    
    # Convert to condensed form for linkage
    condensed_distance = squareform(distance_matrix)
    
    return condensed_distance, distance_matrix


def perform_antibiotic_clustering(phi_matrix: pd.DataFrame, 
                                   n_clusters: int = None,
                                   distance_threshold: float = 0.85,
                                   method: str = 'average') -> tuple:
    """
    Cluster antibiotics based on co-resistance patterns.
    
    Parameters:
    -----------
    phi_matrix : pd.DataFrame
        Phi coefficient matrix (n_antibiotics x n_antibiotics)
    n_clusters : int
        Number of clusters to create (if None, uses distance threshold)
    distance_threshold : float
        Distance threshold for cluster assignment (1 - phi, so 0.85 means phi > 0.15)
    method : str
        Linkage method ('average', 'complete', 'single', 'ward')
    
    Returns:
    --------
    tuple: (cluster_labels, linkage_matrix, antibiotic_names)
    """
    console.subheader("Antibiotic Clustering")
    console.kv("Antibiotics", len(phi_matrix))
    console.kv("Linkage method", method)
    
    antibiotics = list(phi_matrix.columns)
    
    # Convert similarity to distance (distance = 1 - phi)
    condensed_dist, dist_matrix = convert_similarity_to_distance(phi_matrix)
    
    # Perform hierarchical clustering
    Z = linkage(condensed_dist, method=method)
    
    # Assign clusters - use distance threshold approach to separate
    # antibiotics with actual co-resistance from independent ones
    if n_clusters is not None:
        cluster_labels = fcluster(Z, n_clusters, criterion='maxclust')
        console.kv("Using fixed clusters", n_clusters)
    else:
        # Use distance threshold: antibiotics with phi < 0.15 are considered independent
        cluster_labels = fcluster(Z, distance_threshold, criterion='distance')
        console.kv("Distance threshold", f"{distance_threshold} (phi > {1-distance_threshold:.2f} for clustering)")
    
    n_actual_clusters = len(set(cluster_labels))
    console.success(f"Clusters formed: {n_actual_clusters}")
    
    return cluster_labels, Z, antibiotics, dist_matrix


def analyze_antibiotic_clusters(antibiotics: list, 
                                 cluster_labels: np.ndarray,
                                 phi_matrix: pd.DataFrame) -> pd.DataFrame:
    """
    Analyze and characterize each antibiotic cluster.
    """
    console.header("ANTIBIOTIC CLUSTER ANALYSIS")
    
    # Create results dataframe
    results = pd.DataFrame({
        'Antibiotic': antibiotics,
        'Cluster': cluster_labels,
        'Class': [ANTIBIOTIC_CLASSES.get(ab, 'Unknown') for ab in antibiotics]
    })
    
    # Calculate mean phi within each cluster
    mean_phi_within = []
    for ab in antibiotics:
        cluster = results[results['Antibiotic'] == ab]['Cluster'].values[0]
        cluster_abs = results[results['Cluster'] == cluster]['Antibiotic'].tolist()
        cluster_abs.remove(ab)
        if cluster_abs:
            mean_phi = phi_matrix.loc[ab, cluster_abs].mean()
        else:
            mean_phi = 0
        mean_phi_within.append(mean_phi)
    
    results['Mean_Phi_Within_Cluster'] = mean_phi_within
    results = results.sort_values(['Cluster', 'Mean_Phi_Within_Cluster'], 
                                   ascending=[True, False])
    
    # Print cluster summaries
    for cluster_id in sorted(results['Cluster'].unique()):
        cluster_data = results[results['Cluster'] == cluster_id]
        abs_list = cluster_data['Antibiotic'].tolist()
        classes = cluster_data['Class'].value_counts()
        
        console.subheader(f"CLUSTER {cluster_id}: {len(abs_list)} antibiotics")
        console.kv("Antibiotics", ', '.join(abs_list))
        console.kv("Primary class(es)", str(classes.head(2).to_dict()))
        
        # Identify strongest co-resistance pairs within cluster
        if len(abs_list) > 1:
            pairs = []
            for i, ab1 in enumerate(abs_list):
                for ab2 in abs_list[i+1:]:
                    phi = phi_matrix.loc[ab1, ab2]
                    if phi > 0:
                        pairs.append((ab1, ab2, phi))
            
            pairs.sort(key=lambda x: x[2], reverse=True)
            if pairs:
                console.info("Strongest co-resistance pairs:")
                for ab1, ab2, phi in pairs[:3]:
                    console.bullet(f"{ab1} ↔ {ab2}: φ = {phi:.3f}", 4)
        
        console.info("Interpretation:")
        if 'Tetracyclines' in classes.index and 'Folate-Inhibitors' in classes.index:
            console.bullet("Aquaculture-associated resistance cluster", 4)
            console.bullet("TE-SXT linkage suggests mobile genetic elements", 4)
        elif 'Cephalosporins' in ' '.join(classes.index):
            console.bullet("Beta-lactam resistance cluster", 4)
            console.bullet("May indicate ESBL or AmpC production", 4)
        elif 'Aminoglycosides' in classes.index:
            console.bullet("Aminoglycoside resistance cluster", 4)
            console.bullet("Often chromosomally encoded or independent", 4)
        elif all(phi_matrix.loc[ab, [x for x in abs_list if x != ab]].max() == 0 
                 for ab in abs_list if len([x for x in abs_list if x != ab]) > 0):
            console.bullet("Independent/isolated antibiotics (no co-resistance)", 4)
    
    return results


def visualize_antibiotic_dendrogram(Z: np.ndarray, 
                                     antibiotics: list,
                                     cluster_labels: np.ndarray,
                                     output_path: str):
    """Create dendrogram visualization of antibiotic clustering."""
    plt.figure(figsize=(14, 8))
    
    # Color mapping for clusters
    n_clusters = len(set(cluster_labels))
    colors = plt.cm.tab10(np.linspace(0, 1, n_clusters))
    
    # Create dendrogram
    dendrogram_result = dendrogram(
        Z,
        labels=antibiotics,
        leaf_rotation=45,
        leaf_font_size=12,
        color_threshold=0
    )
    
    plt.title('Antibiotic Clustering Dendrogram\n(Based on Co-Resistance Patterns)', 
              fontsize=14, fontweight='bold')
    plt.xlabel('Antibiotic', fontsize=12)
    plt.ylabel('Distance (1 - Phi Coefficient)', fontsize=12)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    console.success(f"Dendrogram saved to: {Path(output_path).name}")


def visualize_clustered_heatmap(phi_matrix: pd.DataFrame,
                                 Z: np.ndarray,
                                 output_path: str):
    """Create clustered heatmap of co-resistance matrix."""
    plt.figure(figsize=(12, 10))
    
    # Create clustermap (automatically reorders by hierarchical clustering)
    g = sns.clustermap(
        phi_matrix,
        cmap='RdYlBu_r',
        linewidths=0.5,
        figsize=(12, 10),
        dendrogram_ratio=(0.15, 0.15),
        cbar_pos=(0.02, 0.8, 0.03, 0.15),
        vmin=0, vmax=1
    )
    
    g.ax_heatmap.set_title('Antibiotic Co-Resistance Matrix\n(Clustered by Phi Coefficient)', 
                           fontsize=12, fontweight='bold', pad=20)
    
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    console.success(f"Clustered heatmap saved to: {Path(output_path).name}")


def save_results(results: pd.DataFrame, output_dir: str):
    """Save antibiotic clustering results."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Save cluster assignments
    csv_path = output_path / 'antibiotic_clusters.csv'
    results.to_csv(csv_path, index=False)
    console.success(f"Cluster assignments saved to: {csv_path.name}")
    
    return csv_path


def main(args=None):
    """Main entry point for antibiotic clustering analysis."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Antibiotic Clustering Analysis")
    parser.add_argument('--threshold', type=float, default=0.70,
                      help='Distance threshold for clustering (1-phi). Default: 0.70')
    parser.add_argument('--clusters', type=int, default=None,
                      help='Force specific number of clusters (overrides threshold)')
    
    # Parse args (handling both CLI and internal calls)
    if args is None:
        parsed_args = parser.parse_args()
    else:
        parsed_args = parser.parse_args(args)

    console.header("ANTIBIOTIC CLUSTERING ANALYSIS", "Clustering antibiotics by co-resistance patterns")
    
    # Load co-resistance matrix
    phi_matrix = load_coresistance_matrix()
    
    # Use distance threshold from args
    distance_threshold = parsed_args.threshold
    n_clusters = parsed_args.clusters
    
    # Perform clustering
    cluster_labels, Z, antibiotics, dist_matrix = perform_antibiotic_clustering(
        phi_matrix, 
        n_clusters=n_clusters,
        distance_threshold=distance_threshold,
        method='average'
    )
    
    # Analyze clusters
    results = analyze_antibiotic_clusters(antibiotics, cluster_labels, phi_matrix)
    
    # Save results
    project_root = Path(__file__).parent.parent
    output_dir = project_root / "data" / "processed" / "figures"
    
    console.header("SAVING RESULTS")
    
    save_results(results, str(output_dir))
    
    # Visualizations
    visualize_antibiotic_dendrogram(
        Z, antibiotics, cluster_labels,
        str(output_dir / 'antibiotic_dendrogram.png')
    )
    
    visualize_clustered_heatmap(
        phi_matrix, Z,
        str(output_dir / 'antibiotic_clustered_heatmap.png')
    )
    
    # Summary
    n_clusters_found = len(results['Cluster'].unique())
    console.complete("ANTIBIOTIC CLUSTERING COMPLETE", {
        'Antibiotics clustered': len(antibiotics),
        'Clusters identified': n_clusters_found,
        'Threshold': distance_threshold
    })
    
    return {
        'results': results,
        'linkage_matrix': Z,
        'phi_matrix': phi_matrix
    }


if __name__ == "__main__":
    main()
