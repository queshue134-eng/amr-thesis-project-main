"""
Verification Script for Clustering Results
Checks all numerical claims in clustering-results.typ against actual data
"""

import pickle
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.metrics import silhouette_score
from pathlib import Path

# Paths
PROJECT_ROOT = Path(r'c:\Users\quesh\Downloads\amr-thesis-project-main-v3.3\amr_thesis_project_code')
DATA_DIR = PROJECT_ROOT / 'data' / 'processed'
ARTIFACTS_DIR = DATA_DIR / 'clustering_artifacts'

def compute_wcss(X, labels):
    """Compute Within-Cluster Sum of Squares."""
    wcss = 0
    for label in np.unique(labels):
        cluster_points = X[labels == label]
        centroid = cluster_points.mean(axis=0)
        wcss += np.sum((cluster_points - centroid) ** 2)
    return wcss

def main():
    print("=" * 70)
    print("CLUSTERING RESULTS VERIFICATION")
    print("=" * 70)
    
    # Load clustered dataset
    clustered_path = DATA_DIR / 'clustered_dataset.csv'
    df = pd.read_csv(clustered_path)
    print(f"\nLoaded: {len(df)} isolates from clustered_dataset.csv")
    
    # Load linkage matrix
    with open(ARTIFACTS_DIR / 'linkage_matrix.pkl', 'rb') as f:
        Z = pickle.load(f)
    print(f"Loaded linkage matrix: shape {Z.shape}")
    
    # Get feature columns
    feature_cols = [c for c in df.columns if c.endswith('_encoded')]
    X = df[feature_cols].fillna(df[feature_cols].median()).values
    print(f"Feature matrix: {X.shape[0]} samples x {X.shape[1]} features")
    
    # =====================================================================
    # 1. VERIFY WCSS AND SILHOUETTE SCORES
    # =====================================================================
    print("\n" + "=" * 70)
    print("1. WCSS AND SILHOUETTE VERIFICATION (k=2 to k=10)")
    print("=" * 70)
    print(f"{'k':<5} {'Silhouette':<15} {'WCSS':<15}")
    print("-" * 35)
    
    for k in range(2, 11):
        labels = fcluster(Z, t=k, criterion='maxclust')
        sil = silhouette_score(X, labels)
        wcss = compute_wcss(X, labels)
        print(f"{k:<5} {sil:<15.3f} {wcss:<15.2f}")
    
    # =====================================================================
    # 2. VERIFY DISTANCE THRESHOLDS
    # =====================================================================
    print("\n" + "=" * 70)
    print("2. EUCLIDEAN DISTANCE THRESHOLDS")
    print("=" * 70)
    
    distances = np.sort(Z[:, 2])
    # dist[-1] = merge 2->1
    # dist[-2] = merge 3->2
    # dist[-3] = merge 4->3
    # dist[-4] = merge 5->4
    
    print(f"Merge 5->4 (lower bound for k=4): {distances[-4]:.2f}")
    print(f"Merge 4->3 (upper bound for k=4): {distances[-3]:.2f}")
    print(f"Merge 3->2 (upper bound for k=3): {distances[-2]:.2f}")
    print(f"Merge 2->1 (upper bound for k=2): {distances[-1]:.2f}")
    
    # =====================================================================
    # 3. VERIFY CLUSTER SIZES
    # =====================================================================
    print("\n" + "=" * 70)
    print("3. CLUSTER SIZE VERIFICATION")
    print("=" * 70)
    
    cluster_counts = df['CLUSTER'].value_counts().sort_index()
    total = len(df)
    for cluster_id, count in cluster_counts.items():
        pct = (count / total) * 100
        print(f"Cluster {cluster_id}: n={count} ({pct:.1f}%)")
    
    # =====================================================================
    # 4. VERIFY MDR PERCENTAGES
    # =====================================================================
    print("\n" + "=" * 70)
    print("4. MDR PERCENTAGE VERIFICATION")
    print("=" * 70)
    
    for cluster_id in sorted(df['CLUSTER'].unique()):
        cluster_df = df[df['CLUSTER'] == cluster_id]
        n = len(cluster_df)
        mdr_count = cluster_df['MDR_FLAG'].sum()
        mdr_pct = (mdr_count / n) * 100 if n > 0 else 0
        print(f"Cluster {cluster_id}: {mdr_count}/{n} MDR = {mdr_pct:.1f}%")
    
    # Total MDR
    total_mdr = df['MDR_FLAG'].sum()
    print(f"\nTotal MDR isolates: {total_mdr}")
    
    # =====================================================================
    # 5. VERIFY SPECIES DISTRIBUTIONS
    # =====================================================================
    print("\n" + "=" * 70)
    print("5. SPECIES DISTRIBUTION VERIFICATION")
    print("=" * 70)
    
    for cluster_id in sorted(df['CLUSTER'].unique()):
        cluster_df = df[df['CLUSTER'] == cluster_id]
        n = len(cluster_df)
        print(f"\nCluster {cluster_id} (n={n}):")
        species_counts = cluster_df['ISOLATE_ID'].value_counts()
        for species, count in species_counts.head(3).items():
            pct = (count / n) * 100
            print(f"  {species}: {count}/{n} = {pct:.1f}%")
    
    # =====================================================================
    # 6. VERIFY REGIONAL DISTRIBUTIONS
    # =====================================================================
    print("\n" + "=" * 70)
    print("6. REGIONAL DISTRIBUTION (Cluster 1 and 3)")
    print("=" * 70)
    
    # Cluster 1
    c1 = df[df['CLUSTER'] == 1]
    c1_region = c1['REGION'].value_counts()
    print(f"\nCluster 1 Regional Distribution:")
    for region, count in c1_region.items():
        pct = (count / len(c1)) * 100
        print(f"  {region}: {count}/{len(c1)} = {pct:.1f}%")
    
    # Cluster 3
    c3 = df[df['CLUSTER'] == 3]
    c3_region = c3['REGION'].value_counts()
    print(f"\nCluster 3 Regional Distribution:")
    for region, count in c3_region.items():
        pct = (count / len(c3)) * 100
        print(f"  {region}: {count}/{len(c3)} = {pct:.1f}%")
    
    # =====================================================================
    # 7. VERIFY ENVIRONMENT DISTRIBUTIONS
    # =====================================================================
    print("\n" + "=" * 70)
    print("7. ENVIRONMENT DISTRIBUTION (Cluster 1 and 3)")
    print("=" * 70)
    
    # Cluster 1
    if 'ENVIRONMENT' in c1.columns:
        c1_env = c1['ENVIRONMENT'].value_counts()
        print(f"\nCluster 1 Environment Distribution:")
        for env, count in c1_env.items():
            pct = (count / len(c1)) * 100
            print(f"  {env}: {count}/{len(c1)} = {pct:.1f}%")
    
    # Cluster 3
    if 'ENVIRONMENT' in c3.columns:
        c3_env = c3['ENVIRONMENT'].value_counts()
        print(f"\nCluster 3 Environment Distribution:")
        for env, count in c3_env.items():
            pct = (count / len(c3)) * 100
            print(f"  {env}: {count}/{len(c3)} = {pct:.1f}%")
    
    print("\n" + "=" * 70)
    print("VERIFICATION COMPLETE")
    print("=" * 70)

if __name__ == "__main__":
    main()
