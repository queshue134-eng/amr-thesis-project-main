"""Compute additional validation metrics without fancy Unicode output."""

import pandas as pd
import numpy as np
from pathlib import Path
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score

def compute_wcss(X, labels):
    """Compute Within-Cluster Sum of Squares."""
    wcss = 0
    for label in np.unique(labels):
        cluster_points = X[labels == label]
        centroid = cluster_points.mean(axis=0)
        wcss += np.sum((cluster_points - centroid) ** 2)
    return wcss

# Load data
project_root = Path(__file__).parent.parent
data_path = project_root / "data" / "processed" / "encoded_dataset.csv"
df = pd.read_csv(data_path)

# Get feature columns
feature_cols = [c for c in df.columns if c.endswith('_encoded')]
X = df[feature_cols].fillna(df[feature_cols].median()).values

# Compute linkage
Z = linkage(X, method='ward', metric='euclidean')

# Compute all metrics
results = []
for k in range(2, 11):
    labels = fcluster(Z, t=k, criterion='maxclust')
    
    # Metrics
    sil = silhouette_score(X, labels) if len(np.unique(labels)) > 1 else 0.0
    wcss = compute_wcss(X, labels)
    ch = calinski_harabasz_score(X, labels) if len(np.unique(labels)) > 1 else 0.0
    db = davies_bouldin_score(X, labels) if len(np.unique(labels)) > 1 else float('inf')
    
    results.append({'k': k, 'silhouette_score': sil, 'wcss': wcss, 'calinski_harabasz': ch, 'davies_bouldin': db})
    print(f"k={k}: Sil={sil:.4f}, WCSS={wcss:.1f}, CH={ch:.1f}, DB={db:.4f}")

results_df = pd.DataFrame(results)

# Save
output_dir = project_root / "data" / "processed" / "figures"
output_dir.mkdir(parents=True, exist_ok=True)
results_df.to_csv(output_dir / "validation_metrics_comparison.csv", index=False)

print("\nSaved to:", output_dir / "validation_metrics_comparison.csv")
