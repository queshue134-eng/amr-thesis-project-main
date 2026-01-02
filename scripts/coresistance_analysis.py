"""
Co-Resistance Network Analysis for AMR Thesis Project
======================================================

This script replaces the circular MDR discrimination task with scientifically 
rigorous co-resistance analysis following the academic review recommendations.

Two main analyses:
1. Co-Resistance Network Construction - statistically significant antibiotic 
   co-resistance pairs with Bonferroni correction
2. Predictive Modeling - predict resistance to key antibiotics (TE, NAL, IPM) 
   from other antibiotics using Random Forest

Output files:
- data/processed/figures/coresistance_network.graphml
- data/processed/figures/coresistance_matrix.csv
- data/processed/figures/coresistance_predictions.csv
- data/processed/figures/coresistance_network.png
"""

import pandas as pd
import numpy as np
import sys
from pathlib import Path
from scipy import stats
from itertools import combinations
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_predict, StratifiedKFold
from sklearn.metrics import roc_auc_score, classification_report, confusion_matrix
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# Add src to path for console import
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

# Import RANDOM_STATE from config for reproducibility
try:
    from config import RANDOM_STATE
except ImportError:
    RANDOM_STATE = 42  # Fallback if config not available

try:
    from utils.console import console, Colors
except ImportError:
    class FallbackConsole:
        def header(self, t, s=None): print(f"\n{'='*70}\n{t}\n{'='*70}")
        def subheader(self, t): print(f"\n  -- {t} --")
        def step(self, c, t, m, d=None): print(f"[{c}/{t}] {m}")
        def success(self, m): print(f"  ✓ {m}")
        def warning(self, m): print(f"  ⚠ {m}")
        def error(self, m): print(f"  ✗ {m}")
        def info(self, m): print(f"  → {m}")
        def kv(self, k, v, i=2): print(f"{'  '*i}{k}: {v}")
        def bullet(self, m, i=2): print(f"{'  '*i}• {m}")
        def complete(self, m=None, s=None): print(f"\n{'='*70}\n{m or 'Complete'}\n{'='*70}")
        def separator(self, c='-', w=60): print(c * w)
    console = FallbackConsole()

# Try to import networkx, if available
try:
    import networkx as nx
    NETWORKX_AVAILABLE = True
except ImportError:
    NETWORKX_AVAILABLE = False
    console.warning("networkx not installed. Network export will be skipped.")


def load_data(data_path: str = None) -> pd.DataFrame:
    """Load the clustered dataset with resistance data."""
    if data_path is None:
        project_root = Path(__file__).parent.parent
        data_path = project_root / "data" / "processed" / "clustered_dataset.csv"
    
    df = pd.read_csv(data_path)
    console.info(f"Loaded {len(df):,} isolates from {Path(data_path).name}")
    return df


def get_resistance_columns(df: pd.DataFrame) -> list:
    """Get columns with binary resistance status (R=1, S/I=0)."""
    return [col for col in df.columns if col.endswith('_RESISTANT')]


def compute_coresistance_matrix(df: pd.DataFrame, 
                                 resistance_cols: list,
                                 alpha: float = 0.05) -> tuple:
    """
    Compute pairwise co-resistance matrix with Bonferroni-corrected chi-square tests.
    
    Returns:
        phi_matrix: DataFrame of phi coefficients (effect sizes)
        pvalue_matrix: DataFrame of p-values
        significant_pairs: List of tuples (ab1, ab2, phi, p_adj) for significant pairs
    """
    n_antibiotics = len(resistance_cols)
    n_tests = n_antibiotics * (n_antibiotics - 1) // 2
    bonferroni_alpha = alpha / n_tests
    
    console.subheader("Co-resistance Matrix Computation")
    console.kv("Antibiotics", n_antibiotics)
    console.kv("Pairwise tests", n_tests)
    console.kv("Bonferroni-corrected alpha", f"{bonferroni_alpha:.6f}")
    
    # Initialize matrices
    antibiotics = [col.replace('_RESISTANT', '') for col in resistance_cols]
    phi_matrix = pd.DataFrame(0.0, index=antibiotics, columns=antibiotics)
    pvalue_matrix = pd.DataFrame(1.0, index=antibiotics, columns=antibiotics)
    significant_pairs = []
    
    for col1, col2 in combinations(resistance_cols, 2):
        ab1 = col1.replace('_RESISTANT', '')
        ab2 = col2.replace('_RESISTANT', '')
        
        # Get binary resistance data
        r1 = df[col1].fillna(0).astype(int)
        r2 = df[col2].fillna(0).astype(int)
        
        # Create contingency table
        contingency = pd.crosstab(r1, r2)
        
        # Chi-square test
        if contingency.shape == (2, 2) and contingency.values.min() >= 5:
            chi2, p_value, dof, expected = stats.chi2_contingency(contingency)
            
            # Phi coefficient (effect size for 2x2 tables)
            n = len(df)
            phi = np.sqrt(chi2 / n)
            
            phi_matrix.loc[ab1, ab2] = phi
            phi_matrix.loc[ab2, ab1] = phi
            pvalue_matrix.loc[ab1, ab2] = p_value
            pvalue_matrix.loc[ab2, ab1] = p_value
            
            # Check significance after Bonferroni correction
            if p_value < bonferroni_alpha:
                significant_pairs.append((ab1, ab2, phi, p_value))
    
    # Set diagonal to 1 (perfect self-correlation)
    np.fill_diagonal(phi_matrix.values, 1.0)
    np.fill_diagonal(pvalue_matrix.values, 0.0)
    
    console.success(f"Significant pairs (after Bonferroni): {len(significant_pairs)}")
    
    return phi_matrix, pvalue_matrix, significant_pairs


def build_coresistance_network(significant_pairs: list,
                               phi_threshold: float = 0.2) -> 'nx.Graph':
    """
    Build a NetworkX graph from significant co-resistance pairs.
    
    Args:
        significant_pairs: List of (ab1, ab2, phi, p_value) tuples
        phi_threshold: Minimum phi coefficient to include edge
    
    Returns:
        NetworkX Graph object
    """
    if not NETWORKX_AVAILABLE:
        print("Warning: NetworkX not available, skipping network construction")
        return None
    
    G = nx.Graph()
    
    # Add edges for significant pairs above threshold
    for ab1, ab2, phi, pvalue in significant_pairs:
        if phi >= phi_threshold:
            G.add_edge(ab1, ab2, weight=phi, pvalue=pvalue)
    
    # Add node attributes (degree)
    for node in G.nodes():
        G.nodes[node]['degree'] = G.degree(node)
        G.nodes[node]['weighted_degree'] = sum(
            G[node][neighbor]['weight'] for neighbor in G.neighbors(node)
        )
    
    console.subheader("Network Statistics")
    console.kv("Nodes (antibiotics)", G.number_of_nodes())
    console.kv("Edges (co-resistance pairs)", G.number_of_edges())
    if G.number_of_nodes() > 0:
        console.kv("Density", f"{nx.density(G):.3f}")
        
        # Find connected components
        components = list(nx.connected_components(G))
        console.kv("Connected components", len(components))
        
        # Identify hub antibiotics (high degree)
        degrees = dict(G.degree())
        if degrees:
            top_hubs = sorted(degrees.items(), key=lambda x: x[1], reverse=True)[:5]
            console.info("Top hub antibiotics:")
            for ab, deg in top_hubs:
                console.bullet(f"{ab}: {deg} connections", 4)
    
    return G


def train_resistance_predictors(df: pd.DataFrame,
                                 encoded_cols: list,
                                 target_antibiotics: list = ['TE', 'NAL', 'IPM']) -> dict:
    """
    Train Random Forest models to predict resistance to key antibiotics.
    
    This provides scientifically valid supervised learning (not circular like MDR).
    
    Args:
        df: Dataframe with encoded resistance data
        encoded_cols: List of encoded resistance columns
        target_antibiotics: Antibiotics to predict
    
    Returns:
        Dictionary with model results per target antibiotic
    """
    console.header("PREDICTIVE MODELING", "Resistance Prediction from Other Antibiotics")
    
    results = {}
    
    for target_ab in target_antibiotics:
        target_col = f'{target_ab}_encoded'
        
        if target_col not in df.columns:
            console.warning(f"{target_col} not found in data, skipping...")
            continue
        
        # Features: all OTHER antibiotics (not the target)
        feature_cols = [c for c in encoded_cols if c != target_col]
        
        # Create binary target (R=1 vs S/I=0)
        y = (df[target_col] == 2).astype(int)
        X = df[feature_cols].fillna(df[feature_cols].median())
        
        # Check class balance
        n_resistant = y.sum()
        n_susceptible = len(y) - n_resistant
        prevalence = n_resistant / len(y) * 100
        
        console.subheader(f"Target: {target_ab} Resistance")
        console.kv("Resistant", f"{n_resistant} ({prevalence:.1f}%)")
        console.kv("Susceptible", f"{n_susceptible} ({100-prevalence:.1f}%)")
        console.kv("Features", f"{len(feature_cols)} other antibiotics")
        
        if n_resistant < 10 or n_susceptible < 10:
            console.warning("Skipping: insufficient class balance")
            continue
        
        # Train Random Forest with cross-validation
        rf = RandomForestClassifier(
            n_estimators=100,
            max_depth=5,
            min_samples_leaf=10,
            class_weight='balanced',
            random_state=RANDOM_STATE,
            n_jobs=-1
        )
        
        # 5-fold stratified cross-validation
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_STATE)
        
        # Get predicted probabilities for AUC
        y_pred_proba = cross_val_predict(rf, X, y, cv=cv, method='predict_proba')[:, 1]
        y_pred = (y_pred_proba >= 0.5).astype(int)
        
        # Calculate metrics
        auc = roc_auc_score(y, y_pred_proba)
        
        console.success(f"AUC-ROC: {auc:.3f}")
        
        # Classification report
        report = classification_report(y, y_pred, output_dict=True, zero_division=0)
        console.kv("Accuracy", f"{report['accuracy']:.3f}")
        console.kv("Precision (R)", f"{report['1']['precision']:.3f}")
        console.kv("Recall (R)", f"{report['1']['recall']:.3f}")
        console.kv("F1 (R)", f"{report['1']['f1-score']:.3f}")
        
        # Feature importance
        rf.fit(X, y)
        importance = pd.Series(rf.feature_importances_, index=feature_cols)
        importance = importance.sort_values(ascending=False)
        top_predictors = importance.head(5)
        
        console.info("Top 5 Predictive Antibiotics:")
        for ab, imp in top_predictors.items():
            ab_name = ab.replace('_encoded', '')
            console.bullet(f"{ab_name}: {imp:.3f}", 4)
        
        # Biological interpretation
        print(f"\n  Biological Interpretation:")
        top_ab = top_predictors.index[0].replace('_encoded', '')
        if target_ab == 'TE':
            print(f"    {target_ab} resistance is best predicted by {top_ab},")
            print(f"    suggesting potential co-carriage on mobile genetic elements.")
        elif target_ab in ['NAL', 'ENR']:
            print(f"    Fluoroquinolone ({target_ab}) resistance prediction suggests")
            print(f"    cross-resistance or co-selection patterns.")
        elif target_ab == 'IPM':
            print(f"    Carbapenem ({target_ab}) resistance prediction is clinically")
            print(f"    significant for identifying high-risk isolates.")
        
        results[target_ab] = {
            'auc': auc,
            'accuracy': report['accuracy'],
            'precision': report['1']['precision'],
            'recall': report['1']['recall'],
            'f1': report['1']['f1-score'],
            'prevalence': prevalence,
            'n_resistant': n_resistant,
            'n_susceptible': n_susceptible,
            'top_predictors': top_predictors.to_dict(),
            'confusion_matrix': confusion_matrix(y, y_pred).tolist()
        }
    
    return results


def visualize_coresistance_network(G, output_path: str):
    """Create visualization of the co-resistance network."""
    if G is None or G.number_of_nodes() == 0:
        print("No network to visualize")
        return
    
    plt.figure(figsize=(12, 10))
    
    # Position nodes using spring layout
    pos = nx.spring_layout(G, seed=42, k=2)
    
    # Node sizes based on degree
    node_sizes = [300 + G.degree(node) * 100 for node in G.nodes()]
    
    # Edge widths based on phi coefficient
    edge_weights = [G[u][v]['weight'] * 3 for u, v in G.edges()]
    
    # Draw network
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, 
                          node_color='lightblue', alpha=0.8)
    nx.draw_networkx_edges(G, pos, width=edge_weights, alpha=0.5, 
                          edge_color='gray')
    nx.draw_networkx_labels(G, pos, font_size=10, font_weight='bold')
    
    plt.title('Co-Resistance Network\n(Bonferroni-corrected significant associations)', 
             fontsize=14)
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    console.success(f"Network visualization saved to: {Path(output_path).name}")


def save_results(phi_matrix: pd.DataFrame,
                 significant_pairs: list,
                 prediction_results: dict,
                 G,
                 output_dir: str):
    """Save all analysis results to files."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    console.header("SAVING RESULTS")
    
    # Save phi coefficient matrix
    matrix_path = output_path / 'coresistance_matrix.csv'
    phi_matrix.to_csv(matrix_path)
    console.success(f"Co-resistance matrix: {matrix_path.name}")
    
    # Save significant pairs
    pairs_df = pd.DataFrame(significant_pairs, 
                           columns=['Antibiotic_1', 'Antibiotic_2', 'Phi', 'P_value'])
    pairs_df = pairs_df.sort_values('Phi', ascending=False)
    pairs_path = output_path / 'coresistance_significant_pairs.csv'
    pairs_df.to_csv(pairs_path, index=False)
    console.success(f"Significant pairs: {pairs_path.name}")
    
    # Save prediction results
    if prediction_results:
        pred_rows = []
        for ab, res in prediction_results.items():
            pred_rows.append({
                'Target_Antibiotic': ab,
                'AUC': res['auc'],
                'Accuracy': res['accuracy'],
                'Precision': res['precision'],
                'Recall': res['recall'],
                'F1': res['f1'],
                'Prevalence_%': res['prevalence'],
                'N_Resistant': res['n_resistant'],
                'N_Susceptible': res['n_susceptible']
            })
        pred_df = pd.DataFrame(pred_rows)
        pred_path = output_path / 'coresistance_predictions.csv'
        pred_df.to_csv(pred_path, index=False)
        console.success(f"Prediction results: {pred_path.name}")
    
    # Save network as GraphML
    if G is not None and NETWORKX_AVAILABLE:
        network_path = output_path / 'coresistance_network.graphml'
        nx.write_graphml(G, str(network_path))
        console.success(f"Network (GraphML): {network_path.name}")
        
        # Save visualization
        viz_path = output_path / 'coresistance_network.png'
        visualize_coresistance_network(G, str(viz_path))


def main():
    """Main entry point for co-resistance analysis."""
    console.header("CO-RESISTANCE NETWORK ANALYSIS", "Alternative to circular MDR discrimination")
    
    # Load data
    project_root = Path(__file__).parent.parent
    df = load_data()
    
    # Get resistance columns
    resistance_cols = get_resistance_columns(df)
    encoded_cols = [c for c in df.columns if c.endswith('_encoded')]
    
    console.kv("Binary resistance columns", len(resistance_cols))
    console.kv("Encoded resistance columns", len(encoded_cols))
    
    # Compute co-resistance matrix with Bonferroni correction
    phi_matrix, pvalue_matrix, significant_pairs = compute_coresistance_matrix(
        df, resistance_cols, alpha=0.05
    )
    
    # Build network
    G = build_coresistance_network(significant_pairs, phi_threshold=0.15)
    
    # Train predictive models
    # Target clinically important antibiotics
    target_antibiotics = ['TE', 'AM', 'SXT', 'CN', 'IPM']
    prediction_results = train_resistance_predictors(
        df, encoded_cols, target_antibiotics
    )
    
    # Save all results
    output_dir = project_root / "data" / "processed" / "figures"
    save_results(phi_matrix, significant_pairs, prediction_results, G, str(output_dir))
    
    # Summary
    console.complete("CO-RESISTANCE ANALYSIS COMPLETE", {
        'Significant pairs': len(significant_pairs),
        'Network nodes': G.number_of_nodes() if G else 0,
        'Models trained': len(prediction_results)
    })
    
    return {
        'phi_matrix': phi_matrix,
        'significant_pairs': significant_pairs,
        'prediction_results': prediction_results,
        'network': G
    }


if __name__ == "__main__":
    main()
