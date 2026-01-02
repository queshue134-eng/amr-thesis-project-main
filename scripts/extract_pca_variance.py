"""
PCA Variance Extraction Script for AMR Thesis Project
======================================================

This script extracts and documents PCA explained variance ratios,
addressing a critical gap in the thesis methodology.

Output files:
- data/processed/figures/pca_variance_explained.csv
- data/processed/figures/pca_variance_documentation.md
"""

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')


def load_data() -> tuple:
    """Load the encoded dataset and prepare for PCA."""
    project_root = Path(__file__).parent.parent
    data_path = project_root / "data" / "processed" / "encoded_dataset.csv"
    
    df = pd.read_csv(data_path)
    print(f"Loaded {len(df)} isolates from {data_path}")
    
    # Get encoded resistance columns (features for PCA)
    feature_cols = [c for c in df.columns if c.endswith('_encoded')]
    print(f"Found {len(feature_cols)} resistance features")
    
    return df, feature_cols


def extract_pca_variance(df: pd.DataFrame, feature_cols: list) -> dict:
    """
    Perform PCA and extract explained variance ratios.
    
    Returns:
        Dictionary with PCA variance information
    """
    print("\n" + "="*70)
    print("PCA EXPLAINED VARIANCE EXTRACTION")
    print("="*70)
    
    # Prepare data
    X = df[feature_cols].fillna(df[feature_cols].median())
    
    # Standardize (important for PCA)
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # Full PCA to get all components
    pca_full = PCA()
    pca_full.fit(X_scaled)
    
    # Extract variance ratios
    variance_ratios = pca_full.explained_variance_ratio_
    cumulative_variance = np.cumsum(variance_ratios)
    
    print(f"\n1. PCA Results (n_components = {len(variance_ratios)}):")
    print("-" * 60)
    print(f"{'PC':<8}{'Variance %':<15}{'Cumulative %':<15}{'Interpretation'}")
    print("-" * 60)
    
    results = []
    for i, (var, cum) in enumerate(zip(variance_ratios, cumulative_variance), 1):
        if i <= 10:  # Show first 10
            if i <= 2:
                interp = "Primary axis"
            elif cum < 0.5:
                interp = "Important"
            elif cum < 0.75:
                interp = "Moderate"
            else:
                interp = "Minor"
            
            print(f"PC{i:<6}{var*100:<15.2f}{cum*100:<15.2f}{interp}")
        
        results.append({
            'Component': f'PC{i}',
            'Variance_Explained_%': var * 100,
            'Cumulative_%': cum * 100
        })
    
    results_df = pd.DataFrame(results)
    
    # Key summary statistics
    pc1_var = variance_ratios[0] * 100
    pc2_var = variance_ratios[1] * 100
    pc1_pc2_cum = cumulative_variance[1] * 100
    
    # How many PCs to reach 50%, 75%, 90%?
    pcs_50 = np.argmax(cumulative_variance >= 0.50) + 1
    pcs_75 = np.argmax(cumulative_variance >= 0.75) + 1
    pcs_90 = np.argmax(cumulative_variance >= 0.90) + 1
    
    print("-" * 60)
    print(f"\n2. Summary Statistics:")
    print(f"   PC1 explains: {pc1_var:.1f}%")
    print(f"   PC2 explains: {pc2_var:.1f}%")
    print(f"   PC1+PC2 cumulative: {pc1_pc2_cum:.1f}%")
    print(f"   PCs needed for 50%: {pcs_50}")
    print(f"   PCs needed for 75%: {pcs_75}")
    print(f"   PCs needed for 90%: {pcs_90}")
    
    # Top loadings for PC1 and PC2
    pca_2d = PCA(n_components=2)
    pca_2d.fit(X_scaled)
    
    loadings = pd.DataFrame(
        pca_2d.components_.T,
        columns=['PC1', 'PC2'],
        index=[c.replace('_encoded', '') for c in feature_cols]
    )
    
    print(f"\n3. Top PC1 Loadings (antibiotics driving PC1 axis):")
    pc1_loadings = loadings['PC1'].abs().sort_values(ascending=False).head(5)
    for ab, load in pc1_loadings.items():
        direction = "+" if loadings.loc[ab, 'PC1'] > 0 else "-"
        print(f"   {ab}: {direction}{load:.3f}")
    
    print(f"\n4. Top PC2 Loadings (antibiotics driving PC2 axis):")
    pc2_loadings = loadings['PC2'].abs().sort_values(ascending=False).head(5)
    for ab, load in pc2_loadings.items():
        direction = "+" if loadings.loc[ab, 'PC2'] > 0 else "-"
        print(f"   {ab}: {direction}{load:.3f}")
    
    return {
        'results_df': results_df,
        'pc1_var': pc1_var,
        'pc2_var': pc2_var,
        'cumulative': pc1_pc2_cum,
        'pcs_50': pcs_50,
        'pcs_75': pcs_75,
        'pcs_90': pcs_90,
        'loadings': loadings
    }


def generate_documentation(pca_info: dict) -> str:
    """Generate markdown documentation for PCA variance."""
    
    doc = f"""
## PCA Explained Variance Analysis

### Summary Statistics

| Metric | Value |
|--------|-------|
| PC1 Explained Variance | {pca_info['pc1_var']:.1f}% |
| PC2 Explained Variance | {pca_info['pc2_var']:.1f}% |
| PC1 + PC2 Cumulative | {pca_info['cumulative']:.1f}% |
| PCs for 50% variance | {pca_info['pcs_50']} |
| PCs for 75% variance | {pca_info['pcs_75']} |
| PCs for 90% variance | {pca_info['pcs_90']} |

### Interpretation

"""
    
    if pca_info['cumulative'] >= 60:
        doc += f"""The first two principal components capture **{pca_info['cumulative']:.1f}%** of total 
resistance variation. This substantial proportion indicates the 2D projections 
provide representative views of resistance structure.
"""
    elif pca_info['cumulative'] >= 40:
        doc += f"""The first two principal components capture **{pca_info['cumulative']:.1f}%** of total 
resistance variation. This moderate proportion means the 2D plots represent 
simplified views that emphasize the dominant axes of variation.
"""
    else:
        doc += f"""The first two principal components capture only **{pca_info['cumulative']:.1f}%** of total 
resistance variation. This limited proportion suggests 2D plots should be 
interpreted cautiously as they represent only a fraction of the full 
{len(pca_info['loadings'])}-dimensional resistance space.

**Limitation Acknowledgment:** Full resistance space is multi-dimensional; 
these projections emphasize the two dominant axes of variation but may not 
fully represent cluster separation in higher dimensions.
"""
    
    # Add loadings interpretation
    top_pc1 = pca_info['loadings']['PC1'].abs().idxmax()
    top_pc2 = pca_info['loadings']['PC2'].abs().idxmax()
    
    doc += f"""
### PC Loadings Interpretation

**PC1 Axis:** Primarily driven by **{top_pc1}** resistance, representing 
the dominant axis of resistance variation.

**PC2 Axis:** Primarily driven by **{top_pc2}** resistance, capturing 
secondary patterns orthogonal to PC1.

### Figure Caption Template

> Principal Component Analysis of resistance profiles (n=492 isolates, 
> {len(pca_info['loadings'])} antibiotics). **PC1 explains {pca_info['pc1_var']:.1f}% 
> of variance, PC2 explains {pca_info['pc2_var']:.1f}% (cumulative {pca_info['cumulative']:.1f}%)**. 
> Points colored by [cluster/region/environment].
"""
    
    return doc


def main():
    """Main entry point for PCA variance extraction."""
    print("="*70)
    print("PCA VARIANCE EXTRACTION: Documenting Explained Variance")
    print("="*70)
    
    # Load data
    df, feature_cols = load_data()
    
    # Extract PCA variance
    pca_info = extract_pca_variance(df, feature_cols)
    
    # Save results
    project_root = Path(__file__).parent.parent
    output_dir = project_root / "data" / "processed" / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Save variance table
    variance_path = output_dir / 'pca_variance_explained.csv'
    pca_info['results_df'].to_csv(variance_path, index=False)
    print(f"\n5. Variance table saved to: {variance_path}")
    
    # Save loadings
    loadings_path = output_dir / 'pca_loadings.csv'
    pca_info['loadings'].to_csv(loadings_path)
    print(f"   Loadings saved to: {loadings_path}")
    
    # Generate and save documentation
    doc = generate_documentation(pca_info)
    doc_path = output_dir / 'pca_variance_documentation.md'
    with open(doc_path, 'w') as f:
        f.write(doc)
    print(f"   Documentation saved to: {doc_path}")
    
    print("\n" + "="*70)
    print("PCA VARIANCE EXTRACTION COMPLETE")
    print("="*70)
    
    # Print documentation excerpt
    print("\n" + "="*70)
    print("DOCUMENTATION EXCERPT (for thesis figures):")
    print("="*70)
    print(doc)
    
    return pca_info


if __name__ == "__main__":
    main()
