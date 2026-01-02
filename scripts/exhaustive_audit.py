#!/usr/bin/env python3
"""
Exhaustive Audit Script for AMR Thesis Project
Verifies consistency between code outputs and manuscript claims
"""

import pandas as pd
import pickle
import os
from pathlib import Path

def print_section(title):
    print(f"\n{'='*70}")
    print(f" {title}")
    print('='*70)

def main():
    root = Path(__file__).parent
    
    print_section("1. DATA FILE VERIFICATION")
    
    # Check all expected data files exist
    expected_files = [
        "data/processed/clustered_dataset.csv",
        "data/processed/analysis_ready_dataset.csv",
        "data/processed/cleaned_dataset.csv",
        "data/processed/encoded_dataset.csv",
        "data/processed/feature_matrix_X.csv",
        "data/processed/unified_raw_dataset.csv",
        "data/processed/clustering_artifacts/clustering_info.pkl",
        "data/processed/clustering_artifacts/linkage_matrix.pkl",
    ]
    
    for f in expected_files:
        fp = root / f
        status = "✓ EXISTS" if fp.exists() else "✗ MISSING"
        print(f"  {status}: {f}")
    
    print_section("2. ISOLATE COUNT VERIFICATION")
    
    # Load clustered dataset
    clustered_df = pd.read_csv(root / "data/processed/clustered_dataset.csv")
    print(f"  Total isolates: {len(clustered_df)}")
    print(f"  Expected (from manuscript): 491")
    print(f"  Match: {'✓ YES' if len(clustered_df) == 491 else '✗ NO'}")
    
    print_section("3. CLUSTER DISTRIBUTION")
    
    cluster_counts = clustered_df['CLUSTER'].value_counts().sort_index()
    print("  Actual cluster sizes:")
    for cluster, count in cluster_counts.items():
        pct = 100 * count / len(clustered_df)
        print(f"    C{cluster}: {count} ({pct:.1f}%)")
    
    # Expected from manuscript (clustering-results.typ)
    expected = {1: 23, 2: 93, 3: 123, 4: 252}
    print("\n  Expected from manuscript:")
    for cluster, count in expected.items():
        print(f"    C{cluster}: {count}")
    
    print("\n  Comparison:")
    all_match = True
    for cluster, expected_count in expected.items():
        actual = cluster_counts.get(cluster, 0)
        match = "✓" if actual == expected_count else "✗"
        if actual != expected_count:
            all_match = False
        print(f"    C{cluster}: Expected={expected_count}, Actual={actual} {match}")
    
    print_section("4. CLUSTERING INFO FROM PICKLE")
    
    try:
        with open(root / "data/processed/clustering_artifacts/clustering_info.pkl", 'rb') as f:
            info = pickle.load(f)
        
        print(f"  Optimal k: {info.get('optimal_k', 'N/A')}")
        print(f"  Silhouette Score: {info.get('silhouette_score', 'N/A')}")
        print(f"  WCSS: {info.get('wcss', 'N/A')}")
        print(f"  Method: {info.get('method', 'N/A')}")
        print(f"  Metric: {info.get('metric', 'N/A')}")
        
        # Keys available
        print(f"\n  Available keys: {list(info.keys())}")
    except Exception as e:
        print(f"  Error loading pickle: {e}")
    
    print_section("5. SPECIES DISTRIBUTION BY CLUSTER")
    
    if 'ORIGINAL_SPECIES' in clustered_df.columns:
        for cluster in sorted(clustered_df['CLUSTER'].unique()):
            cluster_df = clustered_df[clustered_df['CLUSTER'] == cluster]
            species_dist = cluster_df['ORIGINAL_SPECIES'].value_counts().head(3)
            print(f"\n  Cluster {cluster} (n={len(cluster_df)}):")
            for species, count in species_dist.items():
                pct = 100 * count / len(cluster_df)
                print(f"    {species}: {count} ({pct:.1f}%)")
    
    print_section("6. MDR PREVALENCE BY CLUSTER")
    
    # Check if MDR column exists
    mdr_cols = [c for c in clustered_df.columns if 'MDR' in c.upper()]
    print(f"  MDR-related columns: {mdr_cols}")
    
    if mdr_cols:
        mdr_col = mdr_cols[0]
        for cluster in sorted(clustered_df['CLUSTER'].unique()):
            cluster_df = clustered_df[clustered_df['CLUSTER'] == cluster]
            mdr_count = cluster_df[mdr_col].sum() if cluster_df[mdr_col].dtype in ['int64', 'float64', 'bool'] else 0
            mdr_pct = 100 * mdr_count / len(cluster_df) if len(cluster_df) > 0 else 0
            print(f"  C{cluster}: {mdr_count}/{len(cluster_df)} MDR ({mdr_pct:.1f}%)")
    
    print_section("7. PCA VARIANCE CHECK")
    
    try:
        from sklearn.preprocessing import StandardScaler
        from sklearn.decomposition import PCA
        
        # Load feature matrix
        X = pd.read_csv(root / "data/processed/feature_matrix_X.csv")
        if 'Unnamed: 0' in X.columns:
            X = X.drop('Unnamed: 0', axis=1)
        
        # Run PCA
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
        pca = PCA()
        pca.fit(X_scaled)
        
        print("  First 5 components variance explained:")
        cumulative = 0
        for i, var in enumerate(pca.explained_variance_ratio_[:5]):
            cumulative += var
            print(f"    PC{i+1}: {100*var:.2f}% (Cumulative: {100*cumulative:.2f}%)")
        
        print("\n  Expected from manuscript (validation-results.typ):")
        print("    PC1: 23.53% (Cumulative: 23.53%)")
        print("    PC2: 16.40% (Cumulative: 39.92%)")
    except Exception as e:
        print(f"  Error running PCA: {e}")
    
    print_section("8. SOURCE CODE FILE CHECK")
    
    # Check all Python files in src/
    src_path = root / "src"
    py_files = list(src_path.rglob("*.py"))
    print(f"  Total Python files in src/: {len(py_files)}")
    
    # Check for __pycache__ directories
    pycache_dirs = list(src_path.rglob("__pycache__"))
    print(f"  __pycache__ directories: {len(pycache_dirs)}")
    
    # List all modules
    for subdir in ['clustering', 'preprocessing', 'analysis', 'supervised', 'visualization', 'utils']:
        subpath = src_path / subdir
        if subpath.exists():
            files = [f.name for f in subpath.glob("*.py") if f.name != "__init__.py" and not f.name.startswith("__")]
            print(f"\n  {subdir}/: {files}")
    
    print_section("9. SCRIPTS DIRECTORY CHECK")
    
    scripts_path = root / "scripts"
    if scripts_path.exists():
        scripts = [f.name for f in scripts_path.glob("*.py")]
        print(f"  Scripts: {scripts}")
    
    print_section("10. DOCUMENTATION FILES CHECK")
    
    docs_check = [
        "README.md",
        "USER_MANUAL.md",
        "FAQ.md",
        "DOCUMENTATION_GUIDE.md",
        "REFACTOR_LOG.md",
        "docs/TECHNICAL_REFERENCE.md",
    ]
    
    for doc in docs_check:
        fp = root / doc
        if fp.exists():
            size = fp.stat().st_size
            print(f"  ✓ {doc} ({size} bytes)")
        else:
            print(f"  ✗ {doc} MISSING")
    
    print_section("AUDIT COMPLETE")
    print("  Review the above output for any discrepancies.")

if __name__ == "__main__":
    main()
