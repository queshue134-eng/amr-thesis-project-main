"""
Cluster Validation Script for AMR Thesis Project
=================================================

This script validates cluster selection using data-driven analysis:
1. Silhouette analysis (maximize silhouette score)
2. Elbow method (WCSS - Within-Cluster Sum of Squares)
3. Combined multi-criteria selection

The optimal k is determined automatically based on statistical criteria.

Output files:
- data/processed/figures/cluster_validation.png
- data/processed/figures/cluster_validation_metrics.csv

References:
- Silhouette: Rousseeuw, P.J. (1987). Silhouettes: a graphical aid
- Elbow: Thorndike, R.L. (1953). Who belongs in the family?
"""

import pandas as pd
import numpy as np
import sys
from pathlib import Path
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.metrics import silhouette_score, silhouette_samples
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# Add src to path for console import
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))
try:
    from utils.console import console, Colors
except ImportError:
    # Fallback if console not available
    class FallbackConsole:
        def header(self, t, s=None): print(f"\n{'='*70}\n{t}\n{'='*70}")
        def step(self, c, t, m, d=None): print(f"[{c}/{t}] {m}")
        def success(self, m): print(f"  ✓ {m}")
        def warning(self, m): print(f"  ⚠ {m}")
        def error(self, m): print(f"  ✗ {m}")
        def info(self, m): print(f"  → {m}")
        def separator(self, c='-', w=60): print(c * w)
        def table_header(self, cols, widths=None): 
            widths = widths or [15]*len(cols)
            print('  ' + ''.join(f"{c:<{w}}" for c, w in zip(cols, widths)))
            print('  ' + '-' * sum(widths))
            return widths
        def table_row(self, vals, widths, h=False): print('  ' + ''.join(f"{str(v):<{w}}" for v, w in zip(vals, widths)))
        def complete(self, m=None, s=None): print(f"\n{'='*70}\n{m or 'Complete'}\n{'='*70}")
        def kv(self, k, v, i=2): print(f"{'  '*i}{k}: {v}")
        def bullet(self, m, i=2): print(f"{'  '*i}• {m}")
    console = FallbackConsole()
    class Colors:
        @classmethod
        def c(cls, t, *a): return t


def load_data(data_path: str = None) -> tuple:
    """Load the encoded dataset and prepare for clustering."""
    if data_path is None:
        project_root = Path(__file__).parent.parent
        data_path = project_root / "data" / "processed" / "encoded_dataset.csv"
    
    df = pd.read_csv(data_path)
    console.info(f"Loaded {len(df):,} isolates from {Path(data_path).name}")
    
    # Get encoded resistance columns (features for clustering)
    feature_cols = [c for c in df.columns if c.endswith('_encoded')]
    console.info(f"Found {len(feature_cols)} resistance features")
    
    # Prepare feature matrix (impute missing with median)
    X = df[feature_cols].fillna(df[feature_cols].median())
    
    return df, X, feature_cols


def compute_wcss(X: np.ndarray, labels: np.ndarray) -> float:
    """Compute Within-Cluster Sum of Squares (WCSS)."""
    wcss = 0
    for label in np.unique(labels):
        cluster_points = X[labels == label]
        centroid = cluster_points.mean(axis=0)
        wcss += np.sum((cluster_points - centroid) ** 2)
    return wcss


def find_elbow_point(wcss_values: dict) -> int:
    """
    Find the elbow point in WCSS curve using the kneedle algorithm.
    
    The elbow point represents diminishing returns in cluster quality,
    indicating the optimal trade-off between cluster count and cohesion.
    
    Uses vector-based angle calculation to find maximum curvature.
    """
    k_values = sorted(wcss_values.keys())
    wcss_list = [wcss_values[k] for k in k_values]
    
    if len(k_values) < 3:
        return k_values[0]
    
    # Normalize to [0,1] range for fair comparison
    k_norm = np.array([(k - min(k_values)) / (max(k_values) - min(k_values)) for k in k_values])
    wcss_norm = np.array([(w - min(wcss_list)) / (max(wcss_list) - min(wcss_list) + 1e-10) for w in wcss_list])
    
    # Calculate perpendicular distance from each point to line from first to last point
    # This is the standard kneedle algorithm approach
    p1 = np.array([k_norm[0], wcss_norm[0]])
    p2 = np.array([k_norm[-1], wcss_norm[-1]])
    
    distances = []
    for i in range(len(k_values)):
        p = np.array([k_norm[i], wcss_norm[i]])
        # Distance from point to line
        d = np.abs(np.cross(p2 - p1, p1 - p)) / (np.linalg.norm(p2 - p1) + 1e-10)
        distances.append(d)
    
    # Elbow is the point with maximum distance (within valid range)
    elbow_idx = np.argmax(distances)
    return k_values[elbow_idx]


def find_optimal_k(df: pd.DataFrame, feature_cols: list, 
                   k_range: range = range(2, 11),
                   min_k: int = 3, max_k: int = 8,
                   method: str = 'combined') -> dict:
    """
    Find optimal number of clusters using COMBINED elbow + silhouette analysis.
    
    This function provides DATA-DRIVEN cluster selection by:
    1. Computing silhouette scores for each k (cluster cohesion/separation)
    2. Computing WCSS for each k (within-cluster variance)
    3. Finding the elbow point (diminishing returns in WCSS)
    4. Selecting k based on combined evidence
    
    The combined method balances:
    - Silhouette analysis: Measures how well-separated clusters are
    - Elbow method: Identifies where adding clusters yields diminishing returns
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe with encoded resistance data
    feature_cols : list
        List of feature column names to use for clustering
    k_range : range
        Range of k values to test (default 2-10)
    min_k : int
        Minimum acceptable k (to avoid trivial solutions)
    max_k : int
        Maximum acceptable k (to avoid overfitting/small clusters)
    method : str
        Selection method: 'silhouette', 'elbow', or 'combined' (default)
    
    Returns:
    --------
    dict
        Dictionary with:
        - optimal_k: Best k based on combined analysis
        - silhouette_score: Score at optimal k
        - wcss: Within-cluster sum of squares at optimal k
        - elbow_k: Elbow point from WCSS analysis
        - all_scores: Dict of k -> metrics
        - justification: Text explaining why k was selected
    """
    from scipy.cluster.hierarchy import linkage, fcluster
    from sklearn.metrics import silhouette_score
    
    console.header("DATA-DRIVEN CLUSTER SELECTION", "Combined Elbow + Silhouette Analysis")
    
    # Prepare feature matrix
    X = df[feature_cols].fillna(df[feature_cols].median()).values
    
    # Compute linkage matrix
    Z = linkage(X, method='ward', metric='euclidean')
    
    # Test each k value - compute both silhouette and WCSS
    silhouette_scores = {}
    wcss_scores = {}
    
    console.info(f"Testing k={k_range.start} to k={k_range.stop-1}...")
    widths = console.table_header(['k', 'Silhouette', 'WCSS', 'Delta Sil.', 'Status'], [6, 14, 14, 14, 20])
    
    prev_sil = None
    for k in k_range:
        labels = fcluster(Z, t=k, criterion='maxclust')
        
        # Silhouette score
        if len(np.unique(labels)) > 1:
            sil = silhouette_score(X, labels)
        else:
            sil = 0.0
        silhouette_scores[k] = sil
        
        # WCSS (within-cluster sum of squares)
        wcss = compute_wcss(X, labels)
        wcss_scores[k] = wcss
        
        # Calculate improvement from previous k
        delta_sil = (sil - prev_sil) if prev_sil is not None else 0.0
        prev_sil = sil
        
        # Status indicator
        if k < min_k:
            status = "(below min_k)"
        elif k > max_k:
            status = "(above max_k)"
        else:
            status = "✓ candidate"
        
        console.table_row([k, f"{sil:.4f}", f"{wcss:.1f}", f"{delta_sil:+.4f}", status], widths)
    
    console.separator()
    
    # STEP 1: Find elbow point (diminishing returns)
    # Use k=2 to max_k for elbow detection (k=2 anchors the curve,
    # but we don't include k values beyond max_k to avoid curve distortion)
    elbow_wcss = {k: w for k, w in wcss_scores.items() if k <= max_k}
    elbow_k = find_elbow_point(elbow_wcss)
    
    # Ensure elbow is within valid selection range
    if elbow_k < min_k:
        elbow_k = min_k
        console.bullet(f"Elbow Analysis: k={elbow_k} (adjusted to min_k)")
    else:
        console.bullet(f"Elbow Analysis: k={elbow_k} (diminishing returns point)")
    
    # STEP 2: Find max silhouette within bounds
    valid_sil = {k: s for k, s in silhouette_scores.items() if min_k <= k <= max_k}
    max_sil_k = max(valid_sil, key=valid_sil.get)
    console.bullet(f"Max Silhouette: k={max_sil_k} (score={valid_sil[max_sil_k]:.4f})")
    
    # STEP 3: Combined selection logic
    if method == 'elbow':
        optimal_k = elbow_k
        selection_reason = "elbow point (diminishing returns)"
    elif method == 'silhouette':
        optimal_k = max_sil_k
        selection_reason = "maximum silhouette score"
    else:  # combined
        # MULTI-CRITERIA SELECTION with ELBOW PRIORITY
        # ==============================================
        # 
        # Selection criteria (per Rousseeuw 1987, Dolnicar 2002):
        # 1. Silhouette score >= 0.40 (moderate-strong cluster structure)
        # 2. Cluster stability: minimum cluster size >= 20 for reliable estimation
        # 3. ELBOW PRIORITY: If elbow point meets both criteria, use it
        # 4. PARSIMONY FALLBACK: Otherwise, select LOWEST k meeting criteria
        #
        # RATIONALE: The elbow point represents the optimal trade-off between
        # cluster compactness and cluster count. It should be preferred when
        # it satisfies quality criteria.
        
        MIN_CLUSTER_SIZE = 20  # Minimum isolates per cluster for stability
        SILHOUETTE_THRESHOLD = 0.40  # Minimum for "strong" cluster structure
        
        def check_stability(k_val):
            """Check if k produces clusters with minimum size >= MIN_CLUSTER_SIZE."""
            labels = fcluster(Z, t=k_val, criterion='maxclust')
            sizes = pd.Series(labels).value_counts()
            return sizes.min() >= MIN_CLUSTER_SIZE, sizes.min()
        
        console.separator()
        console.info("Multi-criteria k-selection (elbow priority + stability + silhouette):")
        
        # Evaluate all candidates
        candidates_info = []
        qualified_candidates = []
        
        for k in sorted(valid_sil.keys()):
            stable, min_size = check_stability(k)
            strong_sil = valid_sil[k] >= SILHOUETTE_THRESHOLD
            qualified = stable and strong_sil
            
            status = "✓ QUALIFIED" if qualified else "✗"
            console.bullet(f"k={k}: silhouette={valid_sil[k]:.3f}, min_cluster={min_size}, stable={stable} → {status}")
            
            candidates_info.append({
                'k': k, 'silhouette': valid_sil[k], 'stable': stable, 
                'min_size': min_size, 'qualified': qualified
            })
            
            if qualified:
                qualified_candidates.append(k)
        
        console.separator()
        
        # Check if elbow point qualifies
        elbow_stable, elbow_min_size = check_stability(elbow_k)
        elbow_strong_sil = silhouette_scores[elbow_k] >= SILHOUETTE_THRESHOLD
        elbow_qualified = elbow_stable and elbow_strong_sil
        
        console.bullet(f"Elbow point k={elbow_k}: qualified={elbow_qualified}")
        
        if elbow_qualified:
            # ELBOW PRIORITY: Elbow point meets all criteria - use it
            optimal_k = elbow_k
            selection_reason = f"elbow point meeting all criteria (silhouette={silhouette_scores[elbow_k]:.3f}≥{SILHOUETTE_THRESHOLD}, min_cluster={elbow_min_size}≥{MIN_CLUSTER_SIZE})"
            console.bullet(f"Selected k={optimal_k} (elbow priority)")
        elif qualified_candidates:
            # PARSIMONY FALLBACK: Select the LOWEST k that meets all criteria
            optimal_k = min(qualified_candidates)
            selection_reason = f"lowest k meeting all criteria (elbow k={elbow_k} did not qualify)"
            console.bullet(f"Qualified candidates: {qualified_candidates}")
            console.bullet(f"Selected k={optimal_k} (parsimony fallback)")
        else:
            # Fallback: Any stable candidate (relax silhouette requirement)
            stable_only = [c['k'] for c in candidates_info if c['stable']]
            if stable_only:
                # Select the stable k closest to elbow
                optimal_k = min(stable_only, key=lambda k: abs(k - elbow_k))
                selection_reason = f"closest stable k to elbow (no k met silhouette≥{SILHOUETTE_THRESHOLD})"
            else:
                # Ultimate fallback: use elbow point
                optimal_k = elbow_k
                selection_reason = "elbow point (fallback, no stable k found)"
    
    optimal_score = silhouette_scores[optimal_k]
    optimal_wcss = wcss_scores[optimal_k]
    
    print()
    console.success(f"OPTIMAL k = {optimal_k}")
    console.kv("Method", method, 3)
    console.kv("Reason", selection_reason, 3)
    console.kv("Silhouette", f"{optimal_score:.4f}", 3)
    console.kv("WCSS", f"{optimal_wcss:.1f}", 3)
    
    # Check for small clusters warning
    labels = fcluster(Z, t=optimal_k, criterion='maxclust')
    cluster_sizes = pd.Series(labels).value_counts()
    min_cluster_size = cluster_sizes.min()
    if min_cluster_size < 20:
        console.warning(f"Smallest cluster has only {min_cluster_size} isolates. Consider lower k for stability.")
    
    # Generate comprehensive justification
    justification = (
        f"Optimal k={optimal_k} selected using combined elbow + silhouette analysis. "
        f"Elbow point at k={elbow_k} indicates diminishing returns in cluster cohesion beyond this point. "
        f"At k={optimal_k}, silhouette score is {optimal_score:.4f} (indicating "
        f"{'strong' if optimal_score >= 0.4 else 'moderate'} cluster structure) "
        f"with WCSS={optimal_wcss:.1f}. "
        f"This balances statistical cluster quality with interpretability and sample size stability."
    )
    
    console.info(f"Justification: {justification[:80]}...")
    
    return {
        'optimal_k': optimal_k,
        'silhouette_score': optimal_score,
        'wcss': optimal_wcss,
        'elbow_k': elbow_k,
        'max_silhouette_k': max_sil_k,
        'all_silhouette_scores': silhouette_scores,
        'all_wcss_scores': wcss_scores,
        'justification': justification,
        'linkage_matrix': Z,
        'method': method,
        'selection_reason': selection_reason
    }


def validate_clustering(X: pd.DataFrame, k_range: range = range(2, 11)) -> pd.DataFrame:
    """
    Validate cluster quality for different values of k.
    
    Args:
        X: Feature matrix
        k_range: Range of k values to test
        
    Returns:
        DataFrame with metrics for each k
    """
    console.header("CLUSTER VALIDATION ANALYSIS")
    
    # Convert to numpy for scipy
    X_array = X.values
    
    # Compute hierarchical clustering linkage matrix
    console.step(1, 4, "Computing hierarchical clustering", "Ward linkage method")
    Z = linkage(X_array, method='ward', metric='euclidean')
    console.success("Linkage matrix computed")
    
    results = []
    
    console.step(2, 4, "Evaluating cluster quality", "Testing k=2 to k=10")
    widths = console.table_header(['k', 'Silhouette', 'WCSS', 'Interpretation'], [6, 14, 14, 28])
    
    for k in k_range:
        # Get cluster labels
        labels = fcluster(Z, t=k, criterion='maxclust')
        
        # Compute silhouette score
        if len(np.unique(labels)) > 1:
            sil_score = silhouette_score(X_array, labels)
        else:
            sil_score = 0.0
        
        # Compute WCSS
        wcss = compute_wcss(X_array, labels)
        
        # Determine interpretation based on silhouette thresholds
        if sil_score > 0.5:
            interpretation = "Strong clustering"
        elif sil_score > 0.4:
            interpretation = "Moderate-strong clustering"
        elif sil_score > 0.25:
            interpretation = "Moderate clustering"
        else:
            interpretation = "Weak clustering"
        
        console.table_row([k, f"{sil_score:.3f}", f"{wcss:.0f}", interpretation], widths)
        
        results.append({
            'k': k,
            'silhouette_score': sil_score,
            'wcss': wcss,
            'interpretation': interpretation
        })
    
    results_df = pd.DataFrame(results)
    
    # Find optimal k
    best_k_silhouette = results_df.loc[results_df['silhouette_score'].idxmax(), 'k']
    
    console.separator()
    console.step(3, 4, "Analysis Results")
    console.kv("Optimal k (max silhouette)", f"k={best_k_silhouette}")
    best_score = results_df.loc[results_df['silhouette_score'].idxmax(), 'silhouette_score']
    console.kv("Best silhouette score", f"{best_score:.3f}")
    
    return results_df, Z


def create_validation_plot(results_df: pd.DataFrame, output_path: str, 
                           selected_k: int = 5, elbow_k: int = None):
    """Create combined elbow and silhouette plot with dynamic k selection markers.
    
    Parameters:
    -----------
    results_df : pd.DataFrame
        DataFrame with k, silhouette_score, wcss columns
    output_path : str
        Path to save the plot
    selected_k : int
        The k value that was selected (shown in green)
    elbow_k : int, optional
        The elbow point k value (shown with dashed orange line)
    """
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Plot 1: Elbow Method (WCSS)
    ax1 = axes[0]
    ax1.plot(results_df['k'], results_df['wcss'], 'b-o', linewidth=2, markersize=8)
    
    # Mark selected k with red dashed line
    ax1.axvline(x=selected_k, color='red', linestyle='--', alpha=0.7, 
                label=f'k={selected_k} (selected)')
    
    # Mark elbow point if different from selected
    if elbow_k is not None and elbow_k != selected_k:
        ax1.axvline(x=elbow_k, color='orange', linestyle=':', alpha=0.7, 
                    label=f'k={elbow_k} (elbow)')
    
    ax1.set_xlabel('Number of Clusters (k)', fontsize=12)
    ax1.set_ylabel('Within-Cluster Sum of Squares (WCSS)', fontsize=12)
    ax1.set_title('Elbow Method', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    ax1.set_xticks(results_df['k'])
    
    # Plot 2: Silhouette Score
    ax2 = axes[1]
    
    # Color selected k in green, elbow in orange (if different), others in blue
    colors = []
    for k in results_df['k']:
        if k == selected_k:
            colors.append('green')
        elif elbow_k is not None and k == elbow_k:
            colors.append('orange')
        else:
            colors.append('steelblue')
    
    bars = ax2.bar(results_df['k'], results_df['silhouette_score'], color=colors, edgecolor='black')
    ax2.axhline(y=0.25, color='orange', linestyle='--', alpha=0.5, label='Moderate threshold (0.25)')
    ax2.axhline(y=0.40, color='green', linestyle='--', alpha=0.5, label='Strong threshold (0.40)')
    ax2.set_xlabel('Number of Clusters (k)', fontsize=12)
    ax2.set_ylabel('Silhouette Score', fontsize=12)
    ax2.set_title('Silhouette Analysis', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3, axis='y')
    ax2.legend(loc='lower right')
    ax2.set_xticks(results_df['k'])
    
    # Annotate selected k
    selected_sil = results_df[results_df['k']==selected_k]['silhouette_score'].values[0]
    ax2.annotate(f'k={selected_k}: {selected_sil:.3f}\n(SELECTED)', 
                xy=(selected_k, selected_sil), 
                xytext=(selected_k + 1.5, selected_sil + 0.03),
                arrowprops=dict(arrowstyle='->', color='red'),
                fontsize=10, color='red', fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7))
    
    # Annotate elbow if different
    if elbow_k is not None and elbow_k != selected_k:
        elbow_sil = results_df[results_df['k']==elbow_k]['silhouette_score'].values[0]
        ax2.annotate(f'k={elbow_k}: {elbow_sil:.3f}\n(elbow)', 
                    xy=(elbow_k, elbow_sil), 
                    xytext=(elbow_k - 1.5, elbow_sil - 0.05),
                    arrowprops=dict(arrowstyle='->', color='orange'),
                    fontsize=9, color='darkorange', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    console.success(f"Validation plot saved to: {Path(output_path).name}")


def create_silhouette_plot(X: pd.DataFrame, Z: np.ndarray, k: int, output_path: str):
    """Create detailed silhouette plot for the selected k."""
    
    X_array = X.values
    labels = fcluster(Z, t=k, criterion='maxclust')
    
    # Get silhouette values for each sample
    sil_values = silhouette_samples(X_array, labels)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    y_lower = 10
    colors = plt.cm.nipy_spectral(np.linspace(0, 1, k))
    
    for i in range(1, k + 1):
        # Get silhouette values for cluster i
        cluster_sil_values = sil_values[labels == i]
        cluster_sil_values.sort()
        
        size = len(cluster_sil_values)
        y_upper = y_lower + size
        
        ax.fill_betweenx(np.arange(y_lower, y_upper),
                         0, cluster_sil_values,
                         facecolor=colors[i-1], edgecolor=colors[i-1], alpha=0.7)
        
        # Label the clusters
        ax.text(-0.05, y_lower + 0.5 * size, f'C{i}', fontsize=10, fontweight='bold')
        
        y_lower = y_upper + 10
    
    # Overall silhouette score
    avg_sil = silhouette_score(X_array, labels)
    ax.axvline(x=avg_sil, color='red', linestyle='--', linewidth=2, 
               label=f'Average: {avg_sil:.3f}')
    
    ax.set_xlabel('Silhouette Coefficient', fontsize=12)
    ax.set_ylabel('Cluster', fontsize=12)
    ax.set_title(f'Silhouette Plot for k={k} Clusters', fontsize=14, fontweight='bold')
    ax.legend(loc='lower right')
    ax.grid(True, alpha=0.3, axis='x')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    console.success(f"Silhouette detail plot saved to: {Path(output_path).name}")


def generate_documentation(results_df: pd.DataFrame, selected_k: int = None) -> str:
    """Generate markdown documentation for cluster validation."""
    
    best_k_row = results_df.loc[results_df['silhouette_score'].idxmax()]
    optimal_k = selected_k if selected_k else int(best_k_row['k'])
    optimal_row = results_df[results_df['k'] == optimal_k].iloc[0]
    
    doc = f"""
### Cluster Number Selection

Optimal cluster count was determined through combined silhouette analysis and elbow method:

![Cluster Validation](data/processed/figures/cluster_validation.png)

| k | Silhouette Score | WCSS | Interpretation |
|---|------------------|------|----------------|
"""
    
    for _, row in results_df.iterrows():
        k = int(row['k'])
        marker = "**" if k == optimal_k else ""
        doc += f"| {marker}{k}{marker} | {marker}{row['silhouette_score']:.3f}{marker} | {marker}{row['wcss']:.0f}{marker} | {marker}{row['interpretation']}{marker} |\n"
    
    doc += f"""
**Justification for k={optimal_k}:**

Based on convergent evidence from silhouette analysis (score = {optimal_row['silhouette_score']:.3f}) 
and elbow curve inspection, k={optimal_k} was selected as optimal. The silhouette score 
indicates {'moderate' if optimal_row['silhouette_score'] < 0.4 else 'moderate-strong'} clustering structure, 
with distinct resistance phenotypes.

Key observations:
1. Silhouette score at k={optimal_k}: {optimal_row['silhouette_score']:.3f}
2. Maximum silhouette at k={int(best_k_row['k'])}: {best_k_row['silhouette_score']:.3f}
3. WCSS elbow analysis supports this selection
"""
    
    return doc


def main():
    """Main entry point for cluster validation."""
    console.header("CLUSTER VALIDATION", "Dynamic K Selection")
    
    # Load data
    project_root = Path(__file__).parent.parent
    df, X, feature_cols = load_data()
    
    # Run robust combined analysis (Same as pipeline)
    # Using range 2-10 to be comprehensive but respecting logic constraints
    k_result = find_optimal_k(df, feature_cols, k_range=range(2, 11), min_k=3, max_k=8)
    
    selected_k = k_result['optimal_k']
    elbow_k = k_result['elbow_k']
    results_df = pd.DataFrame([
        {'k': k, 'silhouette_score': s, 'wcss': k_result['all_wcss_scores'][k]}
        for k, s in k_result['all_silhouette_scores'].items()
    ]).sort_values('k')
    
    # Create plots
    output_dir = project_root / "data" / "processed" / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    console.step(4, 4, "Creating validation plots")
    
    create_validation_plot(results_df, str(output_dir / 'cluster_validation.png'), 
                          selected_k=selected_k, elbow_k=elbow_k)
    
    # Generate silhouette plot for the SELECTED k
    # We need to re-compute Z for this
    Z = k_result['linkage_matrix']
    create_silhouette_plot(X, Z, k=selected_k, output_path=str(output_dir / f'silhouette_detail_k{selected_k}.png'))
    
    # Also generate k=4 silhouette plot if selected_k != 4 (for comparison)
    if selected_k != 4:
        create_silhouette_plot(X, Z, k=4, output_path=str(output_dir / 'silhouette_detail_k4.png'))
    
    # Save metrics
    metrics_path = output_dir / 'cluster_validation_metrics.csv'
    results_df.to_csv(metrics_path, index=False)
    console.success(f"Metrics saved to: {metrics_path.name}")
    
    # Generate documentation (simple version)
    doc_path = output_dir / 'cluster_validation_documentation.md'
    with open(doc_path, 'w') as f:
        f.write(k_result['justification'])
    console.success(f"Documentation saved to: {doc_path.name}")
    
    console.complete("CLUSTER VALIDATION COMPLETE", {
        'Optimal k': selected_k,
        'Method': k_result['method'],
        'Output': str(output_dir)
    })
    
    return results_df


if __name__ == "__main__":
    main()
