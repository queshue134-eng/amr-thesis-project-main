"""
Regional and Environmental Analysis Module for AMR Thesis Project
Phase 4 - Cluster Distribution Analysis and PCA

PHASE 4 IMPLEMENTATION (Concrete Improvements):
-----------------------------------------------
1. Cluster Distribution Analysis (CRITICAL)
   - Cluster × Region contingency table with counts AND row-wise percentages
   - Cluster × Environment contingency table with counts AND row-wise percentages
   - Stacked bar plots for distribution visualization
   - Chi-square test with Cramér's V effect size

2. One Health Interpretation Discipline (CRITICAL)
   - Controlled language: "over-representation", "association", "co-occurrence"
   - Minimum evidence rule: ≥2 environments with ≥10-15% threshold

3. Multivariate Pattern Analysis - PCA (CRITICAL)
   - Input: resistance fingerprints only (no metadata)
   - Standardization before PCA
   - PCA scatter plots colored by region and environment
   - Explained variance reporting with scree plot

4. Link PCA Back to Clusters (RECOMMENDED)
   - Overlay cluster centroids in PCA space
   - Cluster label annotations
"""

import pandas as pd
import numpy as np
from typing import List, Tuple, Dict, Optional
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from scipy.stats import chi2_contingency


# =============================================================================
# ONE HEALTH INTERPRETATION DISCIPLINE CONSTANTS
# =============================================================================

# Minimum threshold for claiming One Health relevance (proportion within cluster)
ONE_HEALTH_MIN_PROPORTION = 0.10  # 10% minimum proportion

# Minimum number of environments for One Health relevance claim
ONE_HEALTH_MIN_ENVIRONMENTS = 2


def apply_multiple_testing_correction(
    p_values: Dict[str, float],
    method: str = 'bonferroni',
    alpha: float = 0.05
) -> Dict[str, Dict]:
    """
    Apply multiple testing correction to a set of p-values.
    
    STATISTICAL RIGOR ENHANCEMENT:
    When performing multiple chi-square tests (e.g., cluster×region, cluster×environment,
    cluster×species), the probability of at least one false positive increases.
    Multiple testing correction controls the family-wise error rate (FWER) or
    false discovery rate (FDR).
    
    Parameters
    ----------
    p_values : dict
        Dictionary of test_name: p_value pairs
    method : str
        Correction method: 'bonferroni' (conservative) or 'bh' (Benjamini-Hochberg FDR)
    alpha : float
        Significance threshold (default 0.05)
    
    Returns
    -------
    dict
        Dictionary with corrected results for each test:
        {test_name: {'raw_p': float, 'adjusted_p': float, 'significant': bool}}
    
    References
    ----------
    - Bonferroni: Dunn, O.J. (1961). Multiple comparisons among means. JASA.
    - Benjamini-Hochberg: Benjamini & Hochberg (1995). JRSS-B, 57(1):289-300.
    """
    n_tests = len(p_values)
    
    if n_tests == 0:
        return {}
    
    results = {}
    test_names = list(p_values.keys())
    raw_p_values = [p_values[name] for name in test_names]
    
    if method.lower() == 'bonferroni':
        # Bonferroni correction: multiply p-values by number of tests
        # Controls Family-Wise Error Rate (FWER)
        adjusted_alpha = alpha / n_tests
        
        for name, raw_p in zip(test_names, raw_p_values):
            adjusted_p = min(raw_p * n_tests, 1.0)  # Cap at 1.0
            results[name] = {
                'raw_p': raw_p,
                'adjusted_p': adjusted_p,
                'significant_raw': raw_p < alpha,
                'significant_adjusted': raw_p < adjusted_alpha,
                'correction_method': 'Bonferroni',
                'adjusted_alpha': adjusted_alpha,
                'n_tests': n_tests
            }
    
    elif method.lower() in ['bh', 'fdr', 'benjamini-hochberg']:
        # Benjamini-Hochberg FDR correction
        # Less conservative, controls False Discovery Rate
        import numpy as np
        
        # Sort p-values
        sorted_indices = np.argsort(raw_p_values)
        sorted_p = [raw_p_values[i] for i in sorted_indices]
        sorted_names = [test_names[i] for i in sorted_indices]
        
        # Calculate adjusted p-values
        adjusted_p_values = []
        for i, p in enumerate(sorted_p):
            rank = i + 1
            adjusted = min(p * n_tests / rank, 1.0)
            adjusted_p_values.append(adjusted)
        
        # Ensure monotonicity (adjusted p-values should be non-decreasing)
        for i in range(len(adjusted_p_values) - 2, -1, -1):
            adjusted_p_values[i] = min(adjusted_p_values[i], adjusted_p_values[i + 1])
        
        for name, raw_p, adj_p in zip(sorted_names, sorted_p, adjusted_p_values):
            results[name] = {
                'raw_p': raw_p,
                'adjusted_p': adj_p,
                'significant_raw': raw_p < alpha,
                'significant_adjusted': adj_p < alpha,
                'correction_method': 'Benjamini-Hochberg (FDR)',
                'n_tests': n_tests
            }
    
    else:
        raise ValueError(f"Unknown correction method: {method}. Use 'bonferroni' or 'bh'.")
    
    return results


def summarize_multiple_testing_results(corrected_results: Dict[str, Dict]) -> str:
    """
    Generate a summary of multiple testing correction results.
    
    Parameters
    ----------
    corrected_results : dict
        Results from apply_multiple_testing_correction()
    
    Returns
    -------
    str
        Formatted summary string
    """
    if not corrected_results:
        return "No tests to summarize."
    
    first_result = list(corrected_results.values())[0]
    method = first_result.get('correction_method', 'Unknown')
    n_tests = first_result.get('n_tests', len(corrected_results))
    
    lines = [
        f"MULTIPLE TESTING CORRECTION SUMMARY ({method})",
        f"Number of tests: {n_tests}",
        "",
        "| Test | Raw p-value | Adjusted p-value | Significant (raw) | Significant (adjusted) |",
        "|------|-------------|------------------|-------------------|------------------------|"
    ]
    
    n_sig_raw = 0
    n_sig_adjusted = 0
    
    for test_name, result in corrected_results.items():
        raw_sig = "Yes" if result['significant_raw'] else "No"
        adj_sig = "Yes" if result['significant_adjusted'] else "No"
        lines.append(
            f"| {test_name} | {result['raw_p']:.6f} | {result['adjusted_p']:.6f} | {raw_sig} | {adj_sig} |"
        )
        n_sig_raw += result['significant_raw']
        n_sig_adjusted += result['significant_adjusted']
    
    lines.extend([
        "",
        f"Significant tests (α=0.05): {n_sig_raw} (raw), {n_sig_adjusted} (adjusted)",
        "",
        "NOTE: Multiple testing correction reduces false positives but may miss true effects.",
        "      Report BOTH raw and adjusted p-values in thesis.",
    ])
    
    return "\n".join(lines)


def cross_tabulate_clusters(df: pd.DataFrame,
                           cluster_col: str = 'CLUSTER',
                           group_col: str = 'REGION') -> pd.DataFrame:
    """
    Create cross-tabulation of clusters by grouping variable.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe with cluster labels
    cluster_col : str
        Column name for cluster labels
    group_col : str
        Column name for grouping variable
    
    Returns:
    --------
    pd.DataFrame
        Cross-tabulation table
    """
    if group_col not in df.columns or cluster_col not in df.columns:
        print(f"Required columns not found: {cluster_col}, {group_col}")
        return pd.DataFrame()
    
    crosstab = pd.crosstab(df[cluster_col], df[group_col], margins=True)
    
    return crosstab


def create_contingency_table_with_percentages(df: pd.DataFrame,
                                              cluster_col: str = 'CLUSTER',
                                              group_col: str = 'REGION') -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Create contingency table with both absolute counts AND row-wise percentages.
    
    PHASE 4 REQUIREMENT 1.1:
    - Absolute counts
    - Row-wise percentages (within cluster)
    
    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe with cluster labels and grouping variable
    cluster_col : str
        Column name for cluster labels
    group_col : str
        Column name for grouping variable (REGION or SAMPLE_SOURCE)
    
    Returns:
    --------
    tuple
        (counts_df, percentages_df) - Both DataFrames with margins
    """
    if group_col not in df.columns or cluster_col not in df.columns:
        print(f"Required columns not found: {cluster_col}, {group_col}")
        return pd.DataFrame(), pd.DataFrame()
    
    # Absolute counts with margins
    counts = pd.crosstab(df[cluster_col], df[group_col], margins=True, margins_name='Total')
    
    # Row-wise percentages (within cluster)
    # Calculate percentages excluding the 'Total' row and column for proper calculation
    counts_no_margins = counts.iloc[:-1, :-1]  # Exclude Total row and column
    row_totals = counts_no_margins.sum(axis=1)
    percentages_no_margins = counts_no_margins.div(row_totals, axis=0) * 100
    
    # Add Total column (all should sum to 100%)
    percentages_no_margins['Total'] = 100.0
    
    # Add Total row (overall percentages)
    total_counts = counts_no_margins.sum(axis=0)
    overall_total = total_counts.sum()
    total_row = (total_counts / overall_total * 100) if overall_total > 0 else total_counts * 0
    total_row_with_total = pd.concat([total_row, pd.Series({'Total': 100.0})])
    percentages = pd.concat([percentages_no_margins, 
                            pd.DataFrame([total_row_with_total], index=['Total'])])
    
    return counts, percentages.round(1)


def calculate_cramers_v(contingency_table: pd.DataFrame) -> float:
    """
    Calculate Cramér's V effect size for a contingency table.
    
    Cramér's V is a measure of association between two categorical variables.
    It ranges from 0 (no association) to 1 (perfect association).
    
    Interpretation guidelines:
    - V < 0.1: Negligible association
    - 0.1 ≤ V < 0.3: Small association
    - 0.3 ≤ V < 0.5: Medium association
    - V ≥ 0.5: Large association
    
    Parameters:
    -----------
    contingency_table : pd.DataFrame
        Contingency table (counts, without margins)
    
    Returns:
    --------
    float
        Cramér's V value
    """
    # Remove margins if present
    if 'Total' in contingency_table.columns:
        contingency_table = contingency_table.drop('Total', axis=1)
    if 'Total' in contingency_table.index:
        contingency_table = contingency_table.drop('Total', axis=0)
    
    # Convert to numpy array
    observed = contingency_table.values
    
    # Perform chi-square test
    try:
        chi2, p_value, dof, expected = chi2_contingency(observed)
        
        # Calculate Cramér's V
        n = observed.sum()
        min_dim = min(observed.shape[0] - 1, observed.shape[1] - 1)
        
        if min_dim == 0 or n == 0:
            return 0.0
        
        cramers_v = np.sqrt(chi2 / (n * min_dim))
        return cramers_v
    except Exception:
        return 0.0


def perform_chi_square_test(df: pd.DataFrame,
                            cluster_col: str = 'CLUSTER',
                            group_col: str = 'REGION') -> Dict:
    """
    Perform Chi-square test of independence with Cramér's V effect size.
    
    PHASE 4 REQUIREMENT 1.3:
    - χ² statistic
    - p-value
    - Effect size (Cramér's V)
    
    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe with cluster labels and grouping variable
    cluster_col : str
        Column name for cluster labels
    group_col : str
        Column name for grouping variable
    
    Returns:
    --------
    dict
        Chi-square test results with Cramér's V
    """
    if group_col not in df.columns or cluster_col not in df.columns:
        return {'error': f"Required columns not found: {cluster_col}, {group_col}"}
    
    # Create contingency table (without margins)
    contingency = pd.crosstab(df[cluster_col], df[group_col])
    
    # Perform chi-square test
    try:
        chi2, p_value, dof, expected = chi2_contingency(contingency)
        
        # Calculate Cramér's V
        cramers_v = calculate_cramers_v(contingency)
        
        # Interpret Cramér's V
        if cramers_v < 0.1:
            effect_interpretation = "Negligible association"
        elif cramers_v < 0.3:
            effect_interpretation = "Small association"
        elif cramers_v < 0.5:
            effect_interpretation = "Medium association"
        else:
            effect_interpretation = "Large association"
        
        # LANGUAGE DISCIPLINE: Use controlled language
        result = {
            'chi_square': float(chi2),
            'p_value': float(p_value),
            'degrees_of_freedom': int(dof),
            'cramers_v': float(cramers_v),
            'effect_interpretation': effect_interpretation,
            'significant': p_value < 0.05,
            # CONTROLLED LANGUAGE - say "statistical association", NOT "causal relationship"
            'interpretation': (
                f"Statistical association detected (Chi-sq={chi2:.2f}, p={p_value:.4f}, "
                f"Cramér's V={cramers_v:.3f} - {effect_interpretation}). "
                "NOTE: This indicates association, NOT causal relationship."
            ) if p_value < 0.05 else (
                f"No significant statistical association (Chi-sq={chi2:.2f}, p={p_value:.4f}). "
                "Clusters appear independent of this grouping variable."
            )
        }
        
        return result
        
    except Exception as e:
        return {'error': str(e)}


def identify_one_health_relevant_clusters(df: pd.DataFrame,
                                          cluster_col: str = 'CLUSTER',
                                          environment_col: str = 'SAMPLE_SOURCE',
                                          min_environments: int = ONE_HEALTH_MIN_ENVIRONMENTS,
                                          min_proportion: float = ONE_HEALTH_MIN_PROPORTION) -> Dict:
    """
    Identify clusters with potential One Health relevance.
    
    PHASE 4 REQUIREMENT 2.2 - Minimum Evidence Rule:
    Only claim One Health relevance if:
    - A cluster appears in ≥2 environments
    - With meaningful proportion (≥10-15%)
    
    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe with cluster labels and environment data
    cluster_col : str
        Column name for cluster labels
    environment_col : str
        Column name for environment/sample source
    min_environments : int
        Minimum number of environments for relevance claim
    min_proportion : float
        Minimum proportion threshold (0.0 to 1.0)
    
    Returns:
    --------
    dict
        One Health relevance analysis for each cluster
    """
    if environment_col not in df.columns or cluster_col not in df.columns:
        return {'error': f"Required columns not found: {cluster_col}, {environment_col}"}
    
    results = {
        'cluster_analysis': {},
        'one_health_relevant_clusters': [],
        'methodology_note': (
            f"Minimum evidence rule applied: Clusters must appear in >={min_environments} "
            f"environments with >={min_proportion*100:.0f}% proportion to claim One Health relevance."
        )
    }
    
    # Analyze each cluster
    for cluster_id in sorted(df[cluster_col].unique()):
        cluster_df = df[df[cluster_col] == cluster_id]
        total_in_cluster = len(cluster_df)
        
        # Count environments with meaningful presence
        env_counts = cluster_df[environment_col].value_counts()
        env_proportions = env_counts / total_in_cluster
        
        # Filter environments meeting minimum proportion threshold
        significant_envs = env_proportions[env_proportions >= min_proportion]
        n_significant_envs = len(significant_envs)
        
        cluster_analysis = {
            'total_isolates': total_in_cluster,
            'all_environments': env_counts.to_dict(),
            'environment_proportions': (env_proportions * 100).round(1).to_dict(),
            'significant_environments': significant_envs.index.tolist(),
            'n_significant_environments': n_significant_envs,
            'meets_one_health_criteria': n_significant_envs >= min_environments,
            # CONTROLLED LANGUAGE: Use "over-representation" and "co-occurrence"
            'interpretation': ""
        }
        
        # Generate interpretation with controlled language
        if n_significant_envs >= min_environments:
            env_list = ", ".join(significant_envs.index.tolist())
            cluster_analysis['interpretation'] = (
                f"C{cluster_id} shows CO-OCCURRENCE across {n_significant_envs} environments "
                f"({env_list}) with meaningful proportions (>={min_proportion*100:.0f}%). "
                "This suggests potential One Health relevance. "
                "NOTE: This is an ASSOCIATION, not evidence of transmission."
            )
            results['one_health_relevant_clusters'].append(int(cluster_id))
        else:
            cluster_analysis['interpretation'] = (
                f"C{cluster_id} does not meet minimum evidence criteria for One Health relevance "
                f"(present in only {n_significant_envs} environment(s) with >={min_proportion*100:.0f}% proportion)."
            )
        
        results['cluster_analysis'][int(cluster_id)] = cluster_analysis
    
    return results


def _get_environment_column(df: pd.DataFrame) -> Optional[str]:
    """
    Get the environment/sample source column name.
    
    Different datasets may use different column names for environment data.
    This function checks for common variations.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe to check
    
    Returns:
    --------
    str or None
        The column name if found, otherwise None
    """
    # Check for common environment column names in order of preference
    possible_names = ['SAMPLE_SOURCE', 'ENVIRONMENT', 'SAMPLING_SOURCE', 'ENV', 'SOURCE']
    
    for name in possible_names:
        if name in df.columns:
            return name
    
    return None


def analyze_cluster_distribution(df: pd.DataFrame,
                                 cluster_col: str = 'CLUSTER') -> Dict:
    """
    Comprehensive cluster distribution analysis.
    
    PHASE 4 REQUIREMENTS IMPLEMENTED:
    - Cluster × Region contingency table with counts AND percentages
    - Cluster × Environment contingency table with counts AND percentages
    - Chi-square tests with Cramér's V effect size
    - One Health relevance assessment using minimum evidence rule
    
    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe with cluster labels and metadata
    cluster_col : str
        Column name for cluster labels
    
    Returns:
    --------
    dict
        Distribution analysis results including counts, percentages, and statistical tests
    """
    print("=" * 70)
    print("PHASE 4.1: Cluster Distribution Analysis (Enhanced)")
    print("=" * 70)
    
    # Detect environment column
    env_col = _get_environment_column(df)
    if env_col:
        print(f"\nDetected environment column: {env_col}")
    
    analysis = {
        'by_region': {
            'counts': None,
            'percentages': None
        },
        'environment_column': env_col,
        'by_environment': {
            'counts': None,
            'percentages': None
        },
        'by_species': {
            'counts': None,
            'percentages': None
        },
        'chi_square_tests': {},
        'one_health_assessment': None
    }
    
    # 1. Cross-tabulation by REGION with counts AND percentages
    if 'REGION' in df.columns:
        print("\n1. Cluster × Region Distribution:")
        print("-" * 50)
        counts, percentages = create_contingency_table_with_percentages(df, cluster_col, 'REGION')
        analysis['by_region']['counts'] = counts
        analysis['by_region']['percentages'] = percentages
        
        print("\n   ABSOLUTE COUNTS:")
        print(counts.to_string())
        print("\n   ROW-WISE PERCENTAGES (within cluster):")
        print(percentages.to_string())
        
        # Chi-square test with Cramér's V
        chi_result = perform_chi_square_test(df, cluster_col, 'REGION')
        if 'error' not in chi_result:
            analysis['chi_square_tests']['region'] = chi_result
            print(f"\n   Chi-square test (Clusters vs Region):")
            print(f"   Chi-sq = {chi_result['chi_square']:.4f}")
            print(f"   p-value = {chi_result['p_value']:.6f}")
            print(f"   Cramér's V = {chi_result['cramers_v']:.4f} ({chi_result['effect_interpretation']})")
            print(f"   {chi_result['interpretation']}")
    
    # 2. Cross-tabulation by ENVIRONMENT with counts AND percentages
    # Use detected environment column (SAMPLE_SOURCE, ENVIRONMENT, SAMPLING_SOURCE, etc.)
    if env_col:
        print(f"\n2. Cluster × Environment Distribution (using '{env_col}'):")
        print("-" * 50)
        counts, percentages = create_contingency_table_with_percentages(df, cluster_col, env_col)
        analysis['by_environment']['counts'] = counts
        analysis['by_environment']['percentages'] = percentages
        
        print("\n   ABSOLUTE COUNTS:")
        print(counts.to_string())
        print("\n   ROW-WISE PERCENTAGES (within cluster):")
        print(percentages.to_string())
        
        # Chi-square test with Cramér's V
        chi_result = perform_chi_square_test(df, cluster_col, env_col)
        if 'error' not in chi_result:
            analysis['chi_square_tests']['environment'] = chi_result
            print(f"\n   Chi-square test (Clusters vs Environment):")
            print(f"   Chi-sq = {chi_result['chi_square']:.4f}")
            print(f"   p-value = {chi_result['p_value']:.6f}")
            print(f"   Cramér's V = {chi_result['cramers_v']:.4f} ({chi_result['effect_interpretation']})")
            print(f"   {chi_result['interpretation']}")
        
        # One Health relevance assessment
        print("\n3. One Health Relevance Assessment:")
        print("-" * 50)
        one_health = identify_one_health_relevant_clusters(df, cluster_col, env_col)
        analysis['one_health_assessment'] = one_health
        
        print(f"   {one_health['methodology_note']}")
        if one_health.get('one_health_relevant_clusters'):
            print(f"\n   Clusters with potential One Health relevance: "
                  f"{one_health['one_health_relevant_clusters']}")
            for cluster_id in one_health['one_health_relevant_clusters']:
                cluster_info = one_health['cluster_analysis'].get(cluster_id, {})
                print(f"\n   {cluster_info.get('interpretation', '')}")
        else:
            print("\n   No clusters meet the minimum evidence criteria for One Health relevance.")
    
    # 3. Cross-tabulation by SPECIES (for reference)
    if 'ISOLATE_ID' in df.columns:
        print("\n4. Cluster × Species Distribution (Reference):")
        print("-" * 50)
        counts, percentages = create_contingency_table_with_percentages(df, cluster_col, 'ISOLATE_ID')
        analysis['by_species']['counts'] = counts
        analysis['by_species']['percentages'] = percentages
        
        print("\n   ABSOLUTE COUNTS:")
        print(counts.to_string())
        
        # Chi-square test with Cramér's V
        chi_result = perform_chi_square_test(df, cluster_col, 'ISOLATE_ID')
        if 'error' not in chi_result:
            analysis['chi_square_tests']['species'] = chi_result
            print(f"\n   Chi-square test (Clusters vs Species):")
            print(f"   Chi-sq = {chi_result['chi_square']:.4f}")
            print(f"   p-value = {chi_result['p_value']:.6f}")
            print(f"   Cramér's V = {chi_result['cramers_v']:.4f} ({chi_result['effect_interpretation']})")
    
    return analysis


def perform_pca(df: pd.DataFrame,
                feature_cols: List[str],
                n_components: int = None) -> Tuple[np.ndarray, PCA, Dict]:
    """
    Perform Principal Component Analysis on resistance fingerprints only.
    
    PHASE 4 REQUIREMENT 3.1 - PCA Input Discipline:
    - Input: resistance fingerprints ONLY (no metadata)
    - Scale features (standardization) BEFORE PCA
    - PCA is UNSUPERVISED - labels used only for coloring
    
    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe with encoded resistance values
    feature_cols : list
        List of feature column names (resistance data ONLY)
    n_components : int, optional
        Number of components to extract. If None, extracts all components.
    
    Returns:
    --------
    tuple
        (Transformed data, PCA object, PCA info dict)
    """
    print("\n" + "=" * 70)
    print("PHASE 4.3: Principal Component Analysis (PCA)")
    print("=" * 70)
    
    # PHASE 4 REQUIREMENT: Explicit statement about PCA methodology
    print("\n>>> PCA INPUT DISCIPLINE (Phase 4 Requirement 3.1):")
    print("    - Input: RESISTANCE FINGERPRINTS ONLY (no metadata)")
    print("    - Standardization applied BEFORE PCA")
    print("    - PCA is UNSUPERVISED - labels used only for coloring")
    
    # Prepare data - ONLY resistance features
    existing_cols = [c for c in feature_cols if c in df.columns]
    X = df[existing_cols].copy()
    
    print(f"\n1. Input data: {X.shape[0]} isolates × {X.shape[1]} resistance features")
    
    # Handle missing values
    imputer = SimpleImputer(strategy='median')
    X_imputed = imputer.fit_transform(X)
    
    # PHASE 4 REQUIREMENT: Standardize features BEFORE PCA
    print("2. Standardizing features (z-score normalization)...")
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_imputed)
    
    # Perform PCA - extract all components for scree plot if n_components not specified
    if n_components is None:
        n_components_fit = min(X_scaled.shape[0], X_scaled.shape[1])
    else:
        n_components_fit = n_components
    
    pca_full = PCA(n_components=n_components_fit)
    X_pca_full = pca_full.fit_transform(X_scaled)
    
    # For return, use only first 2 components (or requested number)
    n_return = min(2, n_components_fit) if n_components is None else n_components
    X_pca = X_pca_full[:, :n_return]
    
    # Create PCA object with requested components for compatibility
    pca = PCA(n_components=n_return)
    pca.fit(X_scaled)
    
    # Analysis info with full variance information
    pca_info = {
        'n_components_returned': n_return,
        'n_components_total': n_components_fit,
        'explained_variance_ratio': pca_full.explained_variance_ratio_.tolist(),
        'cumulative_variance': np.cumsum(pca_full.explained_variance_ratio_).tolist(),
        'feature_names': [c.replace('_encoded', '') for c in existing_cols],
        'loadings': {},
        # PHASE 4 REQUIREMENT: Explicit variance reporting
        'pc1_variance': float(pca_full.explained_variance_ratio_[0]) if len(pca_full.explained_variance_ratio_) > 0 else 0,
        'pc2_variance': float(pca_full.explained_variance_ratio_[1]) if len(pca_full.explained_variance_ratio_) > 1 else 0,
        'methodology_statement': (
            "PCA performed on standardized resistance fingerprints only. "
            "No metadata was included in the analysis. "
            "PCA is an unsupervised method - any observed patterns are data-driven, "
            "not influenced by category labels."
        )
    }
    
    # Component loadings (for first 2 components)
    for i in range(min(2, n_components_fit)):
        loadings = dict(zip(pca_info['feature_names'], pca_full.components_[i].tolist()))
        loadings = dict(sorted(loadings.items(), key=lambda x: abs(x[1]), reverse=True))
        pca_info['loadings'][f'PC{i+1}'] = loadings
    
    # PHASE 4 REQUIREMENT 3.3: Report explained variance
    print(f"\n3. EXPLAINED VARIANCE REPORTING (Phase 4 Requirement 3.3):")
    print(f"   PC1: {pca_info['pc1_variance']*100:.2f}% variance explained")
    print(f"   PC2: {pca_info['pc2_variance']*100:.2f}% variance explained")
    print(f"   PC1+PC2 cumulative: {(pca_info['pc1_variance'] + pca_info['pc2_variance'])*100:.2f}%")
    
    # Show first few components
    print(f"\n4. Variance explained by components:")
    for i in range(min(5, n_components_fit)):
        var_ratio = pca_full.explained_variance_ratio_[i]
        cum_var = pca_info['cumulative_variance'][i]
        print(f"   PC{i+1}: {var_ratio*100:.2f}% (cumulative: {cum_var*100:.2f}%)")
    
    print("\n5. Top loadings per component (PC1, PC2):")
    for pc, loadings in pca_info['loadings'].items():
        top_loadings = list(loadings.items())[:5]
        print(f"   {pc}: {', '.join([f'{k}({v:.3f})' for k, v in top_loadings])}")
    
    return X_pca, pca_full, pca_info


def create_scree_plot(pca_info: Dict,
                      figsize: Tuple[int, int] = (10, 6),
                      save_path: str = None) -> Figure:
    """
    Create scree plot showing variance explained by each principal component.
    
    PHASE 4 REQUIREMENT 3.3 (Recommended):
    Include scree plot for explained variance visualization.
    
    Parameters:
    -----------
    pca_info : dict
        PCA information dictionary with explained_variance_ratio
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save the figure
    
    Returns:
    --------
    matplotlib.Figure
        Scree plot figure
    """
    fig, ax1 = plt.subplots(figsize=figsize)
    
    n_components = len(pca_info['explained_variance_ratio'])
    components = range(1, n_components + 1)
    variance_ratio = np.array(pca_info['explained_variance_ratio']) * 100
    cumulative_variance = np.array(pca_info['cumulative_variance']) * 100
    
    # Limit to first 10 components for clarity
    max_display = min(10, n_components)
    components = range(1, max_display + 1)
    variance_ratio = variance_ratio[:max_display]
    cumulative_variance = cumulative_variance[:max_display]
    
    # Bar plot for individual variance
    color_bar = '#1f77b4'
    ax1.bar(components, variance_ratio, alpha=0.7, color=color_bar, 
            edgecolor='black', label='Individual')
    ax1.set_xlabel('Principal Component', fontsize=12)
    ax1.set_ylabel('Variance Explained (%)', color=color_bar, fontsize=12)
    ax1.tick_params(axis='y', labelcolor=color_bar)
    ax1.set_xticks(components)
    
    # Line plot for cumulative variance (secondary y-axis)
    ax2 = ax1.twinx()
    color_line = '#d62728'
    ax2.plot(components, cumulative_variance, 'o-', color=color_line, 
             linewidth=2, markersize=8, label='Cumulative')
    ax2.set_ylabel('Cumulative Variance (%)', color=color_line, fontsize=12)
    ax2.tick_params(axis='y', labelcolor=color_line)
    ax2.set_ylim(0, 105)
    
    # Add value labels on bars
    for i, (comp, var) in enumerate(zip(components, variance_ratio)):
        ax1.text(comp, var + 1, f'{var:.1f}%', ha='center', va='bottom', fontsize=9)
    
    # Add 80% threshold line
    ax2.axhline(y=80, color='gray', linestyle='--', alpha=0.7)
    ax2.text(max_display * 0.85, 82, '80% threshold', fontsize=9, color='gray')
    
    # Title and legend
    plt.title('Scree Plot: Variance Explained by Principal Components\n'
              '(PCA performed on resistance fingerprints only)', fontsize=12, fontweight='bold')
    
    # Combine legends
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='center right')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Scree plot saved to: {save_path}")
    
    return fig


def create_pca_plot(X_pca: np.ndarray,
                    df: pd.DataFrame,
                    color_col: str = 'CLUSTER',
                    pca_info: Dict = None,
                    figsize: Tuple[int, int] = (10, 8),
                    save_path: str = None,
                    show_cluster_centroids: bool = False,
                    cluster_col: str = 'CLUSTER') -> Figure:
    """
    Create PCA scatter plot with optional cluster centroids overlay.
    
    PHASE 4 REQUIREMENTS:
    - 3.2: PCA scatter plot colored by region/environment
    - 5: Link PCA back to clusters (overlay cluster centroids)
    
    INTERPRETATION DISCIPLINE (Phase 4 Requirement 4):
    - ✅ "Partial separation"
    - ✅ "Substantial overlap"
    - ✅ "No clear stratification"
    - ❌ "Clear regional clustering" (unless extremely obvious)
    - ❌ "Environmental origin"
    
    Parameters:
    -----------
    X_pca : np.ndarray
        PCA-transformed data
    df : pd.DataFrame
        Original dataframe for coloring
    color_col : str
        Column name for color coding (REGION, SAMPLE_SOURCE, CLUSTER, etc.)
    pca_info : dict, optional
        PCA information for axis labels
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save the figure
    show_cluster_centroids : bool
        If True, overlay cluster centroids with annotations
    cluster_col : str
        Column name for cluster labels (used for centroid calculation)
    
    Returns:
    --------
    matplotlib.Figure
        PCA plot figure
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # Get colors
    if color_col in df.columns:
        categories = sorted(df[color_col].dropna().unique())
        colors = plt.cm.Set1(np.linspace(0, 1, len(categories)))
        color_map = dict(zip(categories, colors))
        
        for category in categories:
            mask = df[color_col] == category
            ax.scatter(X_pca[mask, 0], X_pca[mask, 1], 
                      c=[color_map[category]], label=str(category),
                      alpha=0.6, s=50, edgecolors='none')
        
        ax.legend(title=color_col, bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)
    else:
        ax.scatter(X_pca[:, 0], X_pca[:, 1], alpha=0.6, s=50, c='steelblue')
    
    # PHASE 4 REQUIREMENT 5: Overlay cluster centroids
    if show_cluster_centroids and cluster_col in df.columns:
        clusters = sorted(df[cluster_col].unique())
        centroid_colors = plt.cm.tab10(np.linspace(0, 1, len(clusters)))
        
        for i, cluster_id in enumerate(clusters):
            mask = df[cluster_col] == cluster_id
            centroid_x = np.mean(X_pca[mask, 0])
            centroid_y = np.mean(X_pca[mask, 1])
            
            # Plot centroid as a larger marker with black edge
            ax.scatter(centroid_x, centroid_y, 
                      c=[centroid_colors[i]], s=300, marker='*',
                      edgecolors='black', linewidths=1.5, zorder=10)
            
            # Add cluster label annotation
            ax.annotate(f'C{cluster_id}', (centroid_x, centroid_y),
                       textcoords="offset points", xytext=(8, 8),
                       fontsize=10, fontweight='bold',
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
    
    # Axis labels with variance explained
    if pca_info:
        var_ratio = pca_info['explained_variance_ratio']
        ax.set_xlabel(f'PC1 ({var_ratio[0]*100:.1f}% variance)', fontsize=11)
        ax.set_ylabel(f'PC2 ({var_ratio[1]*100:.1f}% variance)', fontsize=11)
    else:
        ax.set_xlabel('PC1', fontsize=11)
        ax.set_ylabel('PC2', fontsize=11)
    
    # Title with methodology note
    centroid_note = " (★ = cluster centroids)" if show_cluster_centroids else ""
    ax.set_title(f'PCA of Resistance Profiles (colored by {color_col}){centroid_note}\n'
                 'NOTE: PCA is UNSUPERVISED - labels used only for visualization',
                 fontsize=11)
    
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.4)
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.4)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"PCA plot saved to: {save_path}")
    
    return fig


def create_stacked_bar_plot(df: pd.DataFrame,
                            cluster_col: str = 'CLUSTER',
                            group_col: str = 'REGION',
                            figsize: Tuple[int, int] = (12, 7),
                            save_path: str = None) -> Figure:
    """
    Create stacked bar plot showing cluster composition by region or environment.
    
    PHASE 4 REQUIREMENT 1.2:
    - Clusters on x-axis
    - Region/environment as color-coded categories
    - Use PROPORTIONS, not raw counts
    
    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe with cluster labels and grouping variable
    cluster_col : str
        Column name for cluster labels
    group_col : str
        Column name for grouping variable (REGION or SAMPLE_SOURCE)
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save the figure
    
    Returns:
    --------
    matplotlib.Figure
        Stacked bar plot figure
    """
    if cluster_col not in df.columns or group_col not in df.columns:
        print(f"Required columns not found: {cluster_col}, {group_col}")
        return None
    
    # Create cross-tabulation (proportions within cluster)
    crosstab = pd.crosstab(df[cluster_col], df[group_col], normalize='index') * 100
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create stacked bar plot
    bottom = np.zeros(len(crosstab))
    
    # Use a colormap with distinct colors
    n_categories = len(crosstab.columns)
    colors = plt.cm.Set2(np.linspace(0, 1, n_categories))
    
    cluster_labels = [f'C{c}' for c in crosstab.index]
    x_pos = np.arange(len(cluster_labels))
    
    for i, col in enumerate(crosstab.columns):
        bars = ax.bar(x_pos, crosstab[col], bottom=bottom, label=col, 
                     color=colors[i], edgecolor='white', linewidth=0.5)
        
        # Add percentage labels for significant portions (>10%)
        for j, (val, b) in enumerate(zip(crosstab[col], bottom)):
            if val > 10:  # Only show labels for proportions > 10%
                ax.text(j, b + val/2, f'{val:.0f}%', ha='center', va='center', 
                       fontsize=8, fontweight='bold', color='black')
        
        bottom += crosstab[col].values
    
    ax.set_xlabel(f'Cluster (Resistance Phenotype)', fontsize=11)
    ax.set_ylabel('Proportion (%)', fontsize=11)
    ax.set_title(f'Cluster Composition by {group_col}\n'
                 '(Proportions within each cluster - Row-wise normalization)',
                 fontsize=12, fontweight='bold')
    
    ax.set_xticks(x_pos)
    ax.set_xticklabels(cluster_labels, fontsize=10)
    ax.set_ylim(0, 100)
    
    # Add legend outside plot
    ax.legend(title=group_col, bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=9)
    
    # Add grid for readability
    ax.yaxis.grid(True, linestyle='--', alpha=0.3)
    ax.set_axisbelow(True)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Stacked bar plot saved to: {save_path}")
    
    return fig


def create_pca_biplot(X_pca: np.ndarray,
                      pca: PCA,
                      feature_names: List[str],
                      df: pd.DataFrame = None,
                      color_col: str = 'CLUSTER',
                      figsize: Tuple[int, int] = (12, 10),
                      save_path: str = None) -> Figure:
    """
    Create PCA biplot with loadings.
    
    Parameters:
    -----------
    X_pca : np.ndarray
        PCA-transformed data
    pca : PCA
        Fitted PCA object
    feature_names : list
        List of feature names
    df : pd.DataFrame, optional
        Original dataframe for coloring
    color_col : str
        Column name for color coding
    figsize : tuple
        Figure size
    save_path : str, optional
        Path to save the figure
    
    Returns:
    --------
    matplotlib.Figure
        Biplot figure
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # Scale for visualization
    scale_factor = np.max(np.abs(X_pca)) / np.max(np.abs(pca.components_))
    
    # Plot data points
    if df is not None and color_col in df.columns:
        categories = df[color_col].unique()
        colors = plt.cm.Set1(np.linspace(0, 1, len(categories)))
        color_map = dict(zip(categories, colors))
        
        for category in categories:
            mask = df[color_col] == category
            ax.scatter(X_pca[mask, 0], X_pca[mask, 1],
                      c=[color_map[category]], label=str(category),
                      alpha=0.5, s=30)
        
        ax.legend(title=color_col, bbox_to_anchor=(1.05, 1), loc='upper left')
    else:
        ax.scatter(X_pca[:, 0], X_pca[:, 1], alpha=0.5, s=30, c='gray')
    
    # Plot loadings
    for i, name in enumerate(feature_names):
        ax.arrow(0, 0,
                pca.components_[0, i] * scale_factor * 0.8,
                pca.components_[1, i] * scale_factor * 0.8,
                head_width=0.05 * scale_factor, head_length=0.02 * scale_factor,
                fc='red', ec='red', alpha=0.8)
        ax.text(pca.components_[0, i] * scale_factor,
               pca.components_[1, i] * scale_factor,
               name, color='red', fontsize=8)
    
    ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
    ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
    ax.set_title('PCA Biplot: Resistance Profiles and Antibiotic Loadings')
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.3)
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Biplot saved to: {save_path}")
    
    return fig


def run_regional_environmental_analysis(df: pd.DataFrame,
                                        feature_cols: List[str],
                                        output_dir: str = None) -> Dict:
    """
    Run complete regional and environmental analysis pipeline.
    
    PHASE 4 DELIVERABLES CHECKLIST:
    ✅ Cluster × Region table (counts + %)
    ✅ Cluster × Environment table (counts + %)
    ✅ Distribution plots (region & environment)
    ✅ Chi-square test results with Cramér's V
    ✅ PCA scatter plots (region & environment colored)
    ✅ Explained variance reported
    ✅ Scree plot
    ✅ Language discipline enforced
    ✅ Cluster centroids overlay (optional)
    
    Parameters:
    -----------
    df : pd.DataFrame
        Clustered dataframe
    feature_cols : list
        List of feature column names
    output_dir : str, optional
        Directory to save outputs
    
    Returns:
    --------
    dict
        Complete analysis results
    """
    import os
    
    print("\n" + "=" * 70)
    print("PHASE 4: REGIONAL & ENVIRONMENTAL PATTERN ANALYSIS")
    print("=" * 70)
    
    results = {
        'cluster_distribution': None,
        'one_health_assessment': None,
        'pca': None,
        'figures': {},
        'tables': {}
    }
    
    # Detect environment column
    env_col = _get_environment_column(df)
    results['environment_column'] = env_col
    
    # =========================================================================
    # 1. CLUSTER DISTRIBUTION ANALYSIS (CRITICAL)
    # =========================================================================
    if 'CLUSTER' in df.columns:
        results['cluster_distribution'] = analyze_cluster_distribution(df)
    
    # =========================================================================
    # 2. PCA ANALYSIS (CRITICAL)
    # =========================================================================
    print("\n" + "-" * 70)
    X_pca, pca, pca_info = perform_pca(df, feature_cols)
    results['pca'] = {
        'transformed_data': X_pca,
        'pca_object': pca,
        'pca_info': pca_info
    }
    
    # =========================================================================
    # 3. GENERATE VISUALIZATIONS
    # =========================================================================
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        
        print("\n" + "=" * 70)
        print("PHASE 4: Generating Visualizations")
        print("=" * 70)
        
        # -----------------------------------------------------------------
        # 3.1 Stacked bar plots (CRITICAL - Phase 4 Requirement 1.2)
        # -----------------------------------------------------------------
        print("\n1. Creating stacked bar plots...")
        
        # Cluster composition by REGION
        if 'REGION' in df.columns and 'CLUSTER' in df.columns:
            region_bar_path = os.path.join(output_dir, 'cluster_composition_by_region.png')
            fig = create_stacked_bar_plot(df, 'CLUSTER', 'REGION', save_path=region_bar_path)
            if fig:
                results['figures']['stacked_bar_region'] = region_bar_path
                plt.close(fig)
        
        # Cluster composition by ENVIRONMENT (use detected column)
        if env_col and 'CLUSTER' in df.columns:
            env_bar_path = os.path.join(output_dir, 'cluster_composition_by_environment.png')
            fig = create_stacked_bar_plot(df, 'CLUSTER', env_col, save_path=env_bar_path)
            if fig:
                results['figures']['stacked_bar_environment'] = env_bar_path
                plt.close(fig)
        
        # -----------------------------------------------------------------
        # 3.2 PCA scatter plots (CRITICAL - Phase 4 Requirement 3.2)
        # -----------------------------------------------------------------
        print("\n2. Creating PCA scatter plots...")
        
        # PCA plot by REGION (separate plot)
        if 'REGION' in df.columns:
            pca_region_path = os.path.join(output_dir, 'pca_by_region.png')
            fig = create_pca_plot(X_pca, df, 'REGION', pca_info, save_path=pca_region_path,
                                 show_cluster_centroids=True, cluster_col='CLUSTER')
            results['figures']['pca_region'] = pca_region_path
            plt.close(fig)
        
        # PCA plot by ENVIRONMENT (use detected column)
        if env_col:
            pca_env_path = os.path.join(output_dir, 'pca_by_environment.png')
            fig = create_pca_plot(X_pca, df, env_col, pca_info, save_path=pca_env_path,
                                 show_cluster_centroids=True, cluster_col='CLUSTER')
            results['figures']['pca_environment'] = pca_env_path
            plt.close(fig)
        
        # PCA plot by CLUSTER (with centroids - Phase 4 Requirement 5)
        if 'CLUSTER' in df.columns:
            pca_cluster_path = os.path.join(output_dir, 'pca_by_cluster.png')
            fig = create_pca_plot(X_pca, df, 'CLUSTER', pca_info, save_path=pca_cluster_path,
                                 show_cluster_centroids=True, cluster_col='CLUSTER')
            results['figures']['pca_cluster'] = pca_cluster_path
            plt.close(fig)
        
        # PCA plot by MDR status
        if 'MDR_CATEGORY' in df.columns:
            pca_mdr_path = os.path.join(output_dir, 'pca_by_mdr.png')
            fig = create_pca_plot(X_pca, df, 'MDR_CATEGORY', pca_info, save_path=pca_mdr_path,
                                 show_cluster_centroids=True, cluster_col='CLUSTER')
            results['figures']['pca_mdr'] = pca_mdr_path
            plt.close(fig)
        
        # -----------------------------------------------------------------
        # 3.3 Scree plot (RECOMMENDED - Phase 4 Requirement 3.3)
        # -----------------------------------------------------------------
        print("\n3. Creating scree plot...")
        scree_path = os.path.join(output_dir, 'pca_scree_plot.png')
        fig = create_scree_plot(pca_info, save_path=scree_path)
        results['figures']['scree_plot'] = scree_path
        plt.close(fig)
        
        # -----------------------------------------------------------------
        # 3.4 PCA biplot
        # -----------------------------------------------------------------
        print("\n4. Creating PCA biplot...")
        feature_names = [c.replace('_encoded', '') for c in feature_cols if c in df.columns]
        biplot_path = os.path.join(output_dir, 'pca_biplot.png')
        fig = create_pca_biplot(X_pca, pca, feature_names, df, 'CLUSTER', save_path=biplot_path)
        results['figures']['biplot'] = biplot_path
        plt.close(fig)
        
        # -----------------------------------------------------------------
        # 3.5 Save contingency tables as CSV
        # -----------------------------------------------------------------
        print("\n5. Saving contingency tables...")
        
        cluster_dist = results.get('cluster_distribution', {})
        
        # Region table with counts and percentages
        if cluster_dist.get('by_region', {}).get('counts') is not None:
            region_counts = cluster_dist['by_region']['counts']
            region_pct = cluster_dist['by_region']['percentages']
            
            region_counts_path = os.path.join(output_dir, 'cluster_region_counts.csv')
            region_counts.to_csv(region_counts_path)
            results['tables']['region_counts'] = region_counts_path
            print(f"   Saved: {region_counts_path}")
            
            region_pct_path = os.path.join(output_dir, 'cluster_region_percentages.csv')
            region_pct.to_csv(region_pct_path)
            results['tables']['region_percentages'] = region_pct_path
            print(f"   Saved: {region_pct_path}")
        
        # Environment table with counts and percentages
        if cluster_dist.get('by_environment', {}).get('counts') is not None:
            env_counts = cluster_dist['by_environment']['counts']
            env_pct = cluster_dist['by_environment']['percentages']
            
            env_counts_path = os.path.join(output_dir, 'cluster_environment_counts.csv')
            env_counts.to_csv(env_counts_path)
            results['tables']['environment_counts'] = env_counts_path
            print(f"   Saved: {env_counts_path}")
            
            env_pct_path = os.path.join(output_dir, 'cluster_environment_percentages.csv')
            env_pct.to_csv(env_pct_path)
            results['tables']['environment_percentages'] = env_pct_path
            print(f"   Saved: {env_pct_path}")
        
        # -----------------------------------------------------------------
        # 3.6 Save statistical test results
        # -----------------------------------------------------------------
        print("\n6. Saving statistical test results...")
        
        chi_square_results = cluster_dist.get('chi_square_tests', {})
        if chi_square_results:
            chi_summary = []
            for test_name, result in chi_square_results.items():
                if 'error' not in result:
                    chi_summary.append({
                        'Test': f'Cluster vs {test_name.capitalize()}',
                        'Chi-Square': result.get('chi_square', 'N/A'),
                        'p-value': result.get('p_value', 'N/A'),
                        'Degrees of Freedom': result.get('degrees_of_freedom', 'N/A'),
                        'Cramers V': result.get('cramers_v', 'N/A'),
                        'Effect Size': result.get('effect_interpretation', 'N/A'),
                        'Significant': result.get('significant', 'N/A')
                    })
            
            if chi_summary:
                chi_df = pd.DataFrame(chi_summary)
                chi_path = os.path.join(output_dir, 'chi_square_test_results.csv')
                chi_df.to_csv(chi_path, index=False)
                results['tables']['chi_square_results'] = chi_path
                print(f"   Saved: {chi_path}")
        
        plt.close('all')
        
        print("\n" + "=" * 70)
        print("PHASE 4 ANALYSIS COMPLETE")
        print("=" * 70)
        print(f"\nOutputs saved to: {output_dir}")
        print("\nGenerated figures:")
        for name, path in results['figures'].items():
            print(f"  - {name}: {os.path.basename(path)}")
        print("\nGenerated tables:")
        for name, path in results['tables'].items():
            print(f"  - {name}: {os.path.basename(path)}")
    
    return results


if __name__ == "__main__":
    from pathlib import Path
    
    project_root = Path(__file__).parent.parent.parent
    clustered_path = project_root / "data" / "processed" / "clustered_dataset.csv"
    
    if clustered_path.exists():
        df = pd.read_csv(clustered_path)
        feature_cols = [c for c in df.columns if c.endswith('_encoded')]
        
        output_dir = project_root / "data" / "processed" / "figures"
        results = run_regional_environmental_analysis(df, feature_cols, str(output_dir))
        
        print("\nAnalysis complete!")
    else:
        print(f"Clustered dataset not found at {clustered_path}")
        print("Run hierarchical_clustering.py first.")
