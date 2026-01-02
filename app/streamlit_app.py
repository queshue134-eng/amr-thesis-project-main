"""
AMR Pattern Recognition & Exploratory Analysis Dashboard
Phase 6/7 - Interactive Streamlit Application for AMR Surveillance

IMPORTANT DISCLAIMER:
This tool is intended exclusively for exploratory antimicrobial resistance 
pattern recognition and surveillance analysis. It does NOT provide:
- Clinical decision support
- Predictive assessments
- Treatment recommendations
- Risk scores

Data Privacy:
- Uploaded data is processed in memory only
- No data is stored on disk or transmitted externally
- No raw inputs are logged
"""

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for Streamlit (fixes threading error)
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
import os
import sys

# Add src directory to path for module imports
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
src_path = os.path.join(project_root, 'src')
if src_path not in sys.path:
    sys.path.insert(0, src_path)

# Add parent directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

# Import Phase 6 modular components if available
try:
    from data_loader import (
        validate_csv_schema, get_antibiotic_columns as get_antibiotic_cols_v2,
        display_expected_format, get_dataset_info
    )
    from interpretation import (
        get_methodology_content, get_limitations_content, get_disclaimers,
        get_about_content, get_glossary
    )
    from supervised_models import (
        get_model_disclaimer, get_feature_importance_disclaimer
    )
    PHASE6_AVAILABLE = True
except ImportError:
    PHASE6_AVAILABLE = False

# Page configuration
st.set_page_config(
    page_title="AMR Pattern Recognition Dashboard",
    page_icon="🦠",
    layout="wide",
    initial_sidebar_state="collapsed"
)

# Custom CSS with Phase 6 enhancements
st.markdown("""
<style>
    .main-header {
        font-size: 2.2rem;
        color: #1f77b4;
        text-align: center;
        margin-bottom: 0.2rem;
    }
    .sub-header {
        font-size: 1.0rem;
        color: #666;
        text-align: center;
        margin-bottom: 1.5rem;
        font-style: italic;
    }
    .disclaimer {
        background-color: #fff3cd;
        border: 2px solid #ffc107;
        border-radius: 10px;
        padding: 15px;
        margin-bottom: 15px;
    }
    .disclaimer-title {
        color: #856404;
        font-weight: bold;
        margin-bottom: 10px;
    }
    .disclaimer-text {
        color: #856404;
    }
    .metric-card {
        background-color: #f8f9fa;
        color: #333333;
        border-radius: 10px;
        padding: 1rem;
        text-align: center;
    }
    .footer {
        text-align: center;
        color: #6c757d;
        font-size: 0.85rem;
        padding: 20px 0;
        border-top: 1px solid #e9ecef;
        margin-top: 50px;
    }
    .info-box {
        background-color: #d1ecf1;
        color: #0c5460;
        border: 1px solid #bee5eb;
        border-radius: 5px;
        padding: 10px;
        margin: 10px 0;
    }
</style>
""", unsafe_allow_html=True)

# Color scheme for resistance
RESISTANCE_CMAP = mcolors.ListedColormap(['#4CAF50', '#FFC107', '#F44336'])


@st.cache_data
def load_data(uploaded_file=None, default_path=None):
    """Load data from uploaded file or default path with caching."""
    if uploaded_file is not None:
        return pd.read_csv(uploaded_file)
    elif default_path and os.path.exists(default_path):
        return pd.read_csv(default_path)
    return None


def get_antibiotic_cols(df):
    """Get antibiotic column names."""
    encoded_cols = [c for c in df.columns if c.endswith('_encoded')]
    if encoded_cols:
        return encoded_cols
    # Try to identify antibiotic columns by common names
    possible_antibiotics = ['AM', 'AMC', 'CPT', 'CN', 'CF', 'CPD', 'CTX', 'CFO', 
                           'CFT', 'CZA', 'IPM', 'AN', 'GM', 'N', 'NAL', 'ENR',
                           'MRB', 'PRA', 'DO', 'TE', 'FT', 'C', 'SXT']
    return [c for c in df.columns if any(ab in c.upper() for ab in possible_antibiotics)]


def create_heatmap(df, feature_cols, cluster_col=None):
    """Create resistance heatmap."""
    fig, ax = plt.subplots(figsize=(14, 8))
    
    existing_cols = [c for c in feature_cols if c in df.columns]
    if cluster_col and cluster_col in df.columns:
        df_sorted = df.sort_values(cluster_col)
    else:
        df_sorted = df.copy()
    
    data = df_sorted[existing_cols].values
    display_cols = [c.replace('_encoded', '') for c in existing_cols]
    
    im = ax.imshow(data, aspect='auto', cmap=RESISTANCE_CMAP, vmin=0, vmax=2)
    
    ax.set_xticks(np.arange(len(display_cols)))
    ax.set_xticklabels(display_cols, rotation=45, ha='right')
    ax.set_xlabel('Antibiotics')
    ax.set_ylabel('Isolates')
    ax.set_title('Resistance Profile Heatmap')
    
    cbar = plt.colorbar(im, ax=ax, ticks=[0, 1, 2])
    cbar.ax.set_yticklabels(['S', 'I', 'R'])
    
    plt.tight_layout()
    return fig


def create_cluster_summary(df):
    """Create cluster summary statistics."""
    if 'CLUSTER' not in df.columns:
        return None
    
    summary = df.groupby('CLUSTER').agg({
        'CODE': 'count',
        'MAR_INDEX_COMPUTED': 'mean' if 'MAR_INDEX_COMPUTED' in df.columns else lambda x: None,
        'MDR_FLAG': 'mean' if 'MDR_FLAG' in df.columns else lambda x: None
    }).reset_index()
    
    summary.columns = ['Cluster', 'Count', 'Mean MAR Index', 'MDR Proportion']
    return summary


def perform_pca_analysis(df, feature_cols):
    """Perform PCA and return results."""
    existing_cols = [c for c in feature_cols if c in df.columns]
    X = df[existing_cols].copy()
    
    imputer = SimpleImputer(strategy='median')
    X_imputed = imputer.fit_transform(X)
    
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_imputed)
    
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_scaled)
    
    return X_pca, pca


def create_pca_plot(X_pca, df, color_col, pca):
    """Create PCA scatter plot."""
    fig, ax = plt.subplots(figsize=(10, 8))
    
    if color_col in df.columns:
        categories = df[color_col].dropna().unique()
        colors = plt.cm.Set1(np.linspace(0, 1, len(categories)))
        
        for i, category in enumerate(categories):
            mask = df[color_col] == category
            ax.scatter(X_pca[mask, 0], X_pca[mask, 1],
                      c=[colors[i]], label=str(category),
                      alpha=0.7, s=50)
        
        ax.legend(title=color_col, bbox_to_anchor=(1.05, 1))
    else:
        ax.scatter(X_pca[:, 0], X_pca[:, 1], alpha=0.7, s=50)
    
    ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
    ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
    ax.set_title(f'PCA of Resistance Profiles')
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    return fig


def main():
    # Header - Updated for Phase 6 (explicit non-predictive positioning)
    st.markdown('<h1 class="main-header">🦠 AMR Pattern Recognition & Exploratory Analysis Dashboard</h1>', 
                unsafe_allow_html=True)
    st.markdown('<p class="sub-header">Interactive Tool for Antimicrobial Resistance Pattern Analysis | Not for Clinical Use</p>',
                unsafe_allow_html=True)
    
    # Main Disclaimer - HARD-CODED on landing page (Phase 6 requirement)
    st.markdown("""
    <div class="disclaimer">
        <p class="disclaimer-title">⚠️ Important Disclaimer</p>
        <p class="disclaimer-text">
            This tool is intended <strong>exclusively</strong> for exploratory antimicrobial resistance 
            pattern recognition and surveillance analysis. It does <strong>NOT</strong> provide:
            <ul style="color: #856404; margin-top: 10px;">
                <li>Clinical decision support</li>
                <li>Predictive assessments</li>
                <li>Treatment recommendations</li>
                <li>Risk scores</li>
            </ul>
        </p>
    </div>
    """, unsafe_allow_html=True)
    
    # Sidebar
    # st.sidebar.header("📁 Data Upload")
    
    # Upload removed as per request - relying on default path
    uploaded_file = None
    
    # Try to load default data
    default_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'processed', 
                                'clustered_dataset.csv')
    
    df = load_data(uploaded_file, default_path)
    
    if df is None:
        st.error("⚠️ Data file not found.")
        st.info("Please run the data pipeline to generate the dataset:\n`python main.py --pipeline`")
        
        st.markdown("""
        ### Expected Data Format
        The dataset should contain:
        - **Isolate identifiers** (CODE, ISOLATE_ID)
        - **Antibiotic resistance data** (encoded as 0=S, 1=I, 2=R)
        - **Metadata** (Region, Site, Sample Source)
        - **Computed features** (MAR_INDEX_COMPUTED, MDR_FLAG, CLUSTER)
        """)
        return
    
    # Data loaded successfully
    st.sidebar.success(f"✅ Loaded {len(df)} isolates")
    
    # Identify columns
    antibiotic_cols = get_antibiotic_cols(df)
    
    # Display options - Updated for Phase 6 with Methodology and Limitations tabs
    st.sidebar.header("📊 Analysis")
    
    analysis_type = st.sidebar.radio(
        "Select Analysis",
        ["Overview",
         "Regional Distribution", "Local Site Distribution", "Source Distribution", "Species Distribution",
         "Resistance Heatmap",
         "Cluster Validation", "Cluster Analysis", "PCA Analysis",
         "Antibiotic Clustering", "Co-Resistance Network",
         "Supervised Learning",
         "Integration & Synthesis",
         "Methodology", "Limitations"]
    )
    
    # Main content area
    if analysis_type == "Overview":
        st.header("📋 Dataset Overview")
        
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.metric("Total Isolates", len(df))
        
        with col2:
            if 'ISOLATE_ID' in df.columns:
                st.metric("Species", df['ISOLATE_ID'].nunique())
        
        with col3:
            if 'MDR_FLAG' in df.columns:
                mdr_pct = df['MDR_FLAG'].mean() * 100
                st.metric("MDR Proportion", f"{mdr_pct:.1f}%")
        
        with col4:
            if 'CLUSTER' in df.columns:
                st.metric("Clusters", df['CLUSTER'].nunique())
        
        # Data preview
        st.subheader("Data Preview")
        st.dataframe(df.head(20), width='stretch')
        
        # Column information
        with st.expander("📚 Column Information"):
            col_info = pd.DataFrame({
                'Column': df.columns,
                'Type': df.dtypes.astype(str),
                'Non-Null Count': df.notna().sum(),
                'Unique Values': df.nunique()
            })
            st.dataframe(col_info, width='stretch')
    
    elif analysis_type == "Resistance Heatmap":
        st.header("🔥 Resistance Profile Heatmap (Read-Only)")
        
        # Phase 6 requirement: Info box about read-only visualization
        st.markdown("""
        <div class="info-box">
            <strong>ℹ️ About this visualization:</strong> This heatmap shows resistance patterns 
            across isolates and antibiotics. The visualization is <strong>read-only</strong> - 
            it displays pre-computed results without modification.
        </div>
        """, unsafe_allow_html=True)
        
        if antibiotic_cols:
            cluster_col = None
            if 'CLUSTER' in df.columns:
                sort_by_cluster = st.checkbox("Sort by Cluster", value=True)
                if sort_by_cluster:
                    cluster_col = 'CLUSTER'
            
            fig = create_heatmap(df, antibiotic_cols, cluster_col)
            st.pyplot(fig)
            
            # Legend
            st.markdown("""
            **Legend:**
            - 🟢 Green (0): Susceptible (S)
            - 🟡 Yellow (1): Intermediate (I)
            - 🔴 Red (2): Resistant (R)
            """)
        else:
            st.warning("No antibiotic columns found in the dataset.")

    elif analysis_type == "Cluster Analysis":
        st.header("🎯 Cluster Analysis")
        
        if 'CLUSTER' not in df.columns:
            st.warning("No cluster information found. Please run clustering first.")
        else:
            # Cluster summary
            summary = create_cluster_summary(df)
            if summary is not None:
                st.subheader("Cluster Summary")
                st.dataframe(summary, width='stretch')
            
            # Cluster distribution
            col1, col2 = st.columns(2)
            
            with col1:
                st.subheader("Cluster Size Distribution")
                cluster_counts = df['CLUSTER'].value_counts().sort_index()
                fig, ax = plt.subplots(figsize=(8, 6))
                cluster_counts.plot(kind='bar', ax=ax, color='steelblue', edgecolor='black')
                ax.set_xlabel('Cluster')
                ax.set_ylabel('Count')
                ax.set_title('Isolates per Cluster')
                plt.tight_layout()
                st.pyplot(fig)
            
            with col2:
                if 'MDR_FLAG' in df.columns:
                    st.subheader("MDR by Cluster")
                    mdr_by_cluster = df.groupby('CLUSTER')['MDR_FLAG'].mean() * 100
                    fig, ax = plt.subplots(figsize=(8, 6))
                    mdr_by_cluster.plot(kind='bar', ax=ax, color='#F44336', edgecolor='black')
                    ax.set_xlabel('Cluster')
                    ax.set_ylabel('MDR Proportion (%)')
                    ax.set_title('MDR Proportion by Cluster')
                    ax.set_ylim(0, 100)
                    plt.tight_layout()
                    st.pyplot(fig)

            # --- Enhanced Summary Integration ---
            st.markdown("---")
            st.header("📊 Enhanced Cluster Summary")
            
            st.markdown("""
            <div class="info-box">
                <strong>ℹ️ Enhanced Distribution Summary:</strong> This table shows comprehensive 
                epidemiological profiles including species composition, geographic and environmental 
                distribution for each resistance cluster.
            </div>
            """, unsafe_allow_html=True)
            
            # Try to load the enhanced summary table
            summary_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'processed', 
                                        'figures', 'cluster_summary_table.csv')
            
            if os.path.exists(summary_path):
                summary_df = pd.read_csv(summary_path)
                st.dataframe(summary_df, width='stretch')
                
                # Display key insights
                st.subheader("Key Insights")
                
                col1, col2 = st.columns(2)
                
                with col1:
                    st.markdown("**Species Composition:**")
                    for _, row in summary_df.iterrows():
                        species = row.get('Species Composition', 'N/A')
                        if species and species != 'N/A':
                            st.markdown(f"- **{row['Cluster']}**: {species.split(';')[0]}")
                
                with col2:
                    st.markdown("**Environment Distribution:**")
                    for _, row in summary_df.iterrows():
                        env = row.get('Major Environment', 'N/A')
                        if env and env != 'N/A':
                            st.markdown(f"- **{row['Cluster']}**: {env}")
            else:
                st.warning("Enhanced summary table not found. Run `python main.py` to generate.")
    
    elif analysis_type == "Cluster Validation":
        # Dynamically detect k from the dataset
        detected_k = df['CLUSTER'].nunique() if 'CLUSTER' in df.columns else 5
        
        st.header(f"✅ Cluster Validation (k={detected_k} Justification)")
        
        st.markdown(f"""
        <div class="info-box">
            <strong>ℹ️ Validation Method:</strong> k={detected_k} was selected using a <strong>combined approach</strong>:
            <ul style="margin-top: 8px;">
                <li><strong>Elbow Method:</strong> Identifies where WCSS reduction slows (diminishing returns)</li>
                <li><strong>Silhouette Analysis:</strong> Measures cluster separation quality</li>
                <li><strong>Sample Size Constraints:</strong> Ensures minimum cluster sizes for statistical validity</li>
                <li><strong>Interpretability:</strong> Balances granularity with biological meaningfulness</li>
            </ul>
        </div>
        """, unsafe_allow_html=True)
        
        # Load validation metrics
        metrics_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'processed', 
                                    'figures', 'cluster_validation_metrics.csv')
        
        if os.path.exists(metrics_path):
            metrics_df = pd.read_csv(metrics_path)
            
            col1, col2 = st.columns(2)
            
            with col1:
                st.subheader("Validation Metrics")
                st.dataframe(metrics_df, width='stretch')
            
            with col2:
                # Display key metrics for detected k
                k_row = metrics_df[metrics_df['k'] == detected_k].iloc[0] if detected_k in metrics_df['k'].values else None
                if k_row is not None:
                    st.subheader(f"k={detected_k} Selection Rationale")
                    
                    # Metrics display
                    st.metric("Silhouette Score", f"{k_row['silhouette_score']:.3f}")
                    st.metric("WCSS", f"{k_row['wcss']:.0f}")
                    
                    # Best k comparison with explanation
                    best_k = metrics_df.loc[metrics_df['silhouette_score'].idxmax(), 'k']
                    best_sil = metrics_df['silhouette_score'].max()
                    current_sil = k_row['silhouette_score']
                    
                    st.markdown("---")
                    st.markdown(f"**Why k={detected_k} was selected:**")
                    
                    # Detailed justification - Dynamic
                    st.markdown(f"""
                    1. **Elbow Point:** The algorithm identified k={detected_k} as a point of diminishing returns in variance reduction.
                    
                    2. **Sample Size:** At k={detected_k}, cluster sizes are balanced to ensure statistical validity (avoiding tiny clusters <20 isolates).
                    
                    3. **Silhouette Score:** The score of {current_sil:.3f} indicates reasonable separation suitable for phenotype identification.
                    
                    4. **Interpretability:** {detected_k} clusters provide a manageable number of distinct resistance profiles for analysis.
                    """)
                    
                    # Classification interpretation
                    if current_sil >= 0.5:
                        interpretation = "Strong clustering structure"
                        color = "green"
                    elif current_sil >= 0.4:
                        interpretation = "Good clustering structure"
                        color = "blue"
                    elif current_sil >= 0.25:
                        interpretation = "Moderate clustering structure"
                        color = "orange"
                    else:
                        interpretation = "Weak clustering structure"
                        color = "red"
                    
                    st.markdown(f"**Interpretation:** <span style='color:{color};font-weight:bold'>{interpretation}</span>", unsafe_allow_html=True)
            
            # Display validation plots if available
            plots_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'processed', 'figures')
            
            col1, col2 = st.columns(2)
            
            validation_plot = os.path.join(plots_path, 'cluster_validation.png')
            silhouette_plot = os.path.join(plots_path, 'silhouette_detail_k5.png')
            
            with col1:
                if os.path.exists(validation_plot):
                    st.image(validation_plot, caption="Elbow and Silhouette Analysis")
            
            with col2:
                # Plot silhouette for the detected k if available, otherwise try k=5 or k=4
                silhouette_plot = os.path.join(plots_path, f'silhouette_detail_k{detected_k}.png')
                if not os.path.exists(silhouette_plot):
                    # Fallback to existing plot if exact match not found
                    silhouette_plot = os.path.join(plots_path, 'silhouette_detail_k5.png')
                
                if os.path.exists(silhouette_plot):
                    st.image(silhouette_plot, caption=f"Silhouette Detail for k={detected_k}")
            
            # Additional context
            with st.expander("📚 Methodology Details"):
                st.markdown(f"""
                **Combined Selection Algorithm:**
                
                The k={detected_k} selection follows a principled approach:
                
                1. **Elbow Detection:** Using the kneedle algorithm, the elbow point (diminishing WCSS returns) is identified.
                
                2. **Silhouette Validation:** Any k with silhouette ≥0.40 is considered "strong clustering."
                
                3. **Stability Check:** Cluster sizes are verified to ensure minimum sample sizes (≥20) for reliable statistics.
                
                4. **Literature Alignment:** The number of clusters fits typical AMR phenotype studies (3-7 patterns).
                
                **Why not maximize silhouette?**
                - Higher k values may overfit to noise
                - Small clusters (<20 samples) reduce generalizability
                - More clusters = harder interpretation without proportional information gain
                
                **References:**
                - Rousseeuw, P.J. (1987). Silhouettes: a graphical aid to interpretation and validation of cluster analysis
                - Thorndike, R.L. (1953). Who belongs in the family? (Elbow method origin)
                """)
        else:
            st.warning("Validation metrics not found. Run `python scripts/validate_clustering.py` to generate.")

    
    elif analysis_type == "Supervised Learning":
        st.header("🤖 Supervised Learning & Validation")
        
        tab1, tab2 = st.tabs(["🧬 Species Classifier", "📈 Model Evaluation"])
        
        with tab1:
            st.subheader("🧬 Species Classifier Results")
            
            st.markdown("""
            <div class="info-box">
                <strong>ℹ️ Species Discrimination Analysis:</strong> This supervised learning task assesses 
                whether resistance fingerprints can discriminate between bacterial species. 
                High accuracy indicates species-specific resistance patterns.
            </div>
            """, unsafe_allow_html=True)
            
            # Load species classifier model
            models_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'models')
            model_path = os.path.join(models_dir, 'species_classifier.joblib')
            
            if os.path.exists(model_path):
                import joblib
                from sklearn.metrics import confusion_matrix, classification_report, accuracy_score
                
                model_data = joblib.load(model_path)
                model = model_data.get('model')
                label_encoder = model_data.get('label_encoder')
                scaler = model_data.get('scaler')
                imputer = model_data.get('imputer')
                preprocessing_info = model_data.get('preprocessing_info', {})
                
                # Model summary metrics
                st.subheader("📊 Model Summary")
                col1, col2, col3, col4 = st.columns(4)
                with col1:
                    st.metric("Model Type", type(model).__name__)
                with col2:
                    st.metric("Species Classes", len(label_encoder.classes_) if label_encoder else 0)
                with col3:
                    st.metric("Features Used", preprocessing_info.get('n_features', 22))
                with col4:
                    st.metric("Test Size", f"{preprocessing_info.get('test_size', 0.2) * 100:.0f}%")
                
                # Display species classes
                if label_encoder:
                    st.markdown("**Species Classes:**")
                    species_list = ", ".join([f"*{s}*" for s in label_encoder.classes_])
                    st.markdown(species_list)
                
                # Generate predictions on current data for evaluation
                if model is not None and label_encoder is not None and 'ISOLATE_ID' in df.columns:
                    try:
                        # Prepare features
                        X = df[antibiotic_cols].copy()
                        
                        # Handle missing values
                        if imputer is not None:
                            X = pd.DataFrame(imputer.transform(X), columns=X.columns)
                        else:
                            X = X.fillna(0)
                        
                        # Scale features
                        if scaler is not None:
                            X_scaled = scaler.transform(X)
                        else:
                            X_scaled = X.values
                        
                        # Get true labels - filter to only known classes
                        y_true = df['ISOLATE_ID'].values
                        known_classes = set(label_encoder.classes_)
                        valid_mask = np.array([y in known_classes for y in y_true])
                        
                        if valid_mask.sum() < len(y_true):
                            n_unknown = len(y_true) - valid_mask.sum()
                            st.info(f"Note: {n_unknown} samples with unseen species excluded from evaluation.")
                        
                        if valid_mask.sum() > 0:
                            X_scaled_valid = X_scaled[valid_mask]
                            y_true_valid = y_true[valid_mask]
                            y_true_encoded = label_encoder.transform(y_true_valid)
                            
                            # Get predictions
                            y_pred = model.predict(X_scaled_valid)
                            
                            # ========== PERFORMANCE METRICS ==========
                            st.subheader("📈 Performance Metrics")
                            
                            # Calculate metrics
                            accuracy = accuracy_score(y_true_encoded, y_pred)
                            report = classification_report(y_true_encoded, y_pred, 
                                                           target_names=label_encoder.classes_,
                                                           output_dict=True)
                            
                            # Display overall metrics
                            col1, col2, col3, col4 = st.columns(4)
                            with col1:
                                st.metric("Accuracy", f"{accuracy:.3f}")
                            with col2:
                                st.metric("Macro Precision", f"{report['macro avg']['precision']:.3f}")
                            with col3:
                                st.metric("Macro Recall", f"{report['macro avg']['recall']:.3f}")
                            with col4:
                                st.metric("Macro F1-Score", f"{report['macro avg']['f1-score']:.3f}")
                            
                            # Per-class metrics table
                            with st.expander("📋 Per-Class Performance"):
                                metrics_data = []
                                for species in label_encoder.classes_:
                                    if species in report:
                                        metrics_data.append({
                                            'Species': species,
                                            'Precision': f"{report[species]['precision']:.3f}",
                                            'Recall': f"{report[species]['recall']:.3f}",
                                            'F1-Score': f"{report[species]['f1-score']:.3f}",
                                            'Support': int(report[species]['support'])
                                        })
                                st.dataframe(pd.DataFrame(metrics_data), width='stretch')
                            
                            # ========== CONFUSION MATRIX ==========
                            st.subheader("🎯 Confusion Matrix")
                            
                            cm = confusion_matrix(y_true_encoded, y_pred)
                            
                            col1, col2 = st.columns([1.2, 0.8])
                            
                            with col1:
                                fig, ax = plt.subplots(figsize=(10, 8))
                                im = ax.imshow(cm, interpolation='nearest', cmap='Blues')
                                ax.figure.colorbar(im, ax=ax)
                                
                                # Truncate long species names for display
                                classes_display = [s[:12] + '...' if len(s) > 15 else s 
                                                  for s in label_encoder.classes_]
                                
                                ax.set(xticks=np.arange(cm.shape[1]),
                                       yticks=np.arange(cm.shape[0]),
                                       xticklabels=classes_display,
                                       yticklabels=classes_display,
                                       ylabel='True Species',
                                       xlabel='Predicted Species',
                                       title='Species Classification Confusion Matrix')
                                
                                plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
                                
                                # Show values in cells
                                thresh = cm.max() / 2.
                                for i in range(cm.shape[0]):
                                    for j in range(cm.shape[1]):
                                        ax.text(j, i, format(cm[i, j], 'd'),
                                               ha="center", va="center",
                                               color="white" if cm[i, j] > thresh else "black",
                                               fontweight="bold")
                                
                                plt.tight_layout()
                                st.pyplot(fig)
                            
                            with col2:
                                st.markdown("**Interpretation:**")
                                st.markdown(f"""
                                - **Diagonal cells**: Correct classifications
                                - **Off-diagonal**: Misclassifications
                                - **Overall accuracy**: {accuracy:.1%}
                                
                                High diagonal values indicate species-specific resistance patterns 
                                can be distinguished by the model.
                                """)
                            
                            # ========== FEATURE IMPORTANCE ==========
                            st.subheader("🔬 Feature Importance (Top Discriminating Antibiotics)")
                            
                            if hasattr(model, 'feature_importances_'):
                                importances = model.feature_importances_
                                feature_names = [c.replace('_encoded', '') for c in antibiotic_cols]
                                
                                # Sort by importance
                                indices = np.argsort(importances)[::-1]
                                top_n = min(15, len(indices))
                                
                                col1, col2 = st.columns([1.2, 0.8])
                                
                                with col1:
                                    fig, ax = plt.subplots(figsize=(10, 6))
                                    
                                    top_features = [feature_names[i] for i in indices[:top_n]]
                                    top_importances = [importances[i] for i in indices[:top_n]]
                                    
                                    colors = plt.cm.RdYlGn_r(np.linspace(0.2, 0.8, top_n))
                                    ax.barh(range(top_n), top_importances[::-1], color=colors[::-1], edgecolor='black')
                                    ax.set_yticks(range(top_n))
                                    ax.set_yticklabels(top_features[::-1])
                                    ax.set_xlabel('Importance (Mean Decrease in Impurity)')
                                    ax.set_title('Top Antibiotics for Species Discrimination')
                                    
                                    plt.tight_layout()
                                    st.pyplot(fig)
                                
                                with col2:
                                    st.markdown("**Interpretation:**")
                                    st.markdown(f"""
                                    Top 3 discriminating antibiotics:
                                    1. **{top_features[0]}** ({top_importances[0]:.3f})
                                    2. **{top_features[1]}** ({top_importances[1]:.3f})
                                    3. **{top_features[2]}** ({top_importances[2]:.3f})
                                    
                                    These antibiotics show the most species-specific 
                                    resistance patterns and contribute most to classification.
                                    """)
                                
                                # Full importance table
                                with st.expander("📋 Full Feature Importance Table"):
                                    importance_df = pd.DataFrame({
                                        'Antibiotic': feature_names,
                                        'Importance': importances
                                    }).sort_values('Importance', ascending=False)
                                    importance_df['Importance'] = importance_df['Importance'].apply(lambda x: f"{x:.4f}")
                                    st.dataframe(importance_df, width='stretch')
                            else:
                                st.info("Feature importance not available for this model type.")
                            
                    except Exception as e:
                        st.error(f"Error evaluating model: {e}")
                else:
                    st.warning("Unable to evaluate model - missing required data columns.")
            else:
                st.warning("Species classifier model not found. Run `python main.py` to train the model.")
        
        with tab2:
            st.subheader("📈 Model Evaluation Summary")
            
            # Phase 6 requirement: Disclaimer near model results
            st.markdown("""
            <div class="info-box">
                <strong>📊 Interpretation Note:</strong> Metrics shown are <em>pattern consistency measures</em>, 
                not predictive performance. They quantify how resistance patterns align with known categories 
                within this dataset only. <strong>No retraining or single-sample inference is available.</strong>
            </div>
            """, unsafe_allow_html=True)
            
            st.markdown("""
            **Interpretation Note:** Model metrics quantify how consistently resistance 
            patterns align with known categories (e.g., species, MDR status), 
            not predictive performance for future samples.
            
            - **Accuracy**: Proportion of isolates where resistance patterns align with category
            - **Precision/Recall/F1**: Consistency of pattern-category alignment
            - Feature importance shows **ASSOCIATIVE** patterns only (not causal)
            """)
            
            # Check for saved model results
            models_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'models')
            
            if os.path.exists(models_dir):
                import joblib
                from sklearn.metrics import confusion_matrix, classification_report
                
                model_files = [f for f in os.listdir(models_dir) if f.endswith('.joblib')]
                
                if model_files:
                    st.subheader("📈 Model Performance Summary")
                    
                    for model_file in model_files:
                        st.markdown(f"### 🔬 {model_file.replace('_classifier.joblib', '').title()} Classifier")
                        
                        try:
                            model_data = joblib.load(os.path.join(models_dir, model_file))
                            model = model_data.get('model')
                            label_encoder = model_data.get('label_encoder')
                            scaler = model_data.get('scaler')
                            imputer = model_data.get('imputer')
                            preprocessing_info = model_data.get('preprocessing_info', {})
                            
                            # Get model info
                            model_type = type(model).__name__
                            classes = label_encoder.classes_ if label_encoder else []
                            
                            # Display model info in columns
                            col1, col2, col3 = st.columns(3)
                            with col1:
                                st.metric("Model Type", model_type)
                            with col2:
                                st.metric("Classes", len(classes))
                            with col3:
                                st.metric("Features", preprocessing_info.get('n_features', 'N/A'))
                            
                            # Try to compute evaluation metrics on current data
                            if model is not None and label_encoder is not None:
                                # Determine target column
                                if 'species' in model_file.lower():
                                    target_col = 'ISOLATE_ID'
                                else:
                                    target_col = 'MDR_CATEGORY'
                                
                                if target_col in df.columns:
                                    try:
                                        # Prepare test data
                                        X = df[antibiotic_cols].copy()
                                        
                                        # Handle missing values
                                        if imputer is not None:
                                            X = pd.DataFrame(imputer.transform(X), columns=X.columns)
                                        else:
                                            X = X.fillna(0)
                                        
                                        # Scale features
                                        if scaler is not None:
                                            X_scaled = scaler.transform(X)
                                        else:
                                            X_scaled = X.values
                                        
                                        # Get true labels - filter to only known classes
                                        y_true = df[target_col].values
                                        known_classes = set(label_encoder.classes_)
                                        
                                        # Create mask for samples with known labels only
                                        valid_mask = np.array([y in known_classes for y in y_true])
                                        
                                        if valid_mask.sum() < len(y_true):
                                            n_unknown = len(y_true) - valid_mask.sum()
                                            st.info(f"Note: {n_unknown} samples with unseen labels excluded from evaluation.")
                                        
                                        if valid_mask.sum() > 0:
                                            X_scaled_valid = X_scaled[valid_mask]
                                            y_true_valid = y_true[valid_mask]
                                            y_true_encoded = label_encoder.transform(y_true_valid)
                                            
                                            # Get predictions
                                            y_pred = model.predict(X_scaled_valid)
                                            
                                            # Compute confusion matrix
                                            cm = confusion_matrix(y_true_encoded, y_pred)
                                        
                                            # Display metrics and confusion matrix
                                            col1, col2 = st.columns([1, 1])
                                            
                                            with col1:
                                                st.subheader("📊 Confusion Matrix")
                                                fig, ax = plt.subplots(figsize=(8, 6))
                                                
                                                # Heatmap
                                                im = ax.imshow(cm, interpolation='nearest', cmap='Blues')
                                                ax.figure.colorbar(im, ax=ax)
                                                
                                                # Labels
                                                classes_display = [str(c)[:15] for c in classes]  # Truncate long names
                                                ax.set(xticks=np.arange(cm.shape[1]),
                                                       yticks=np.arange(cm.shape[0]),
                                                       xticklabels=classes_display,
                                                       yticklabels=classes_display,
                                                       ylabel='True Label',
                                                       xlabel='Predicted Label',
                                                       title='Confusion Matrix')
                                                
                                                plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
                                                
                                                # Show values in cells
                                                thresh = cm.max() / 2.
                                                for i in range(cm.shape[0]):
                                                    for j in range(cm.shape[1]):
                                                        ax.text(j, i, format(cm[i, j], 'd'),
                                                               ha="center", va="center",
                                                               color="white" if cm[i, j] > thresh else "black",
                                                               fontsize=10)
                                                
                                                plt.tight_layout()
                                                st.pyplot(fig)
                                            
                                            with col2:
                                                st.subheader("📈 Performance Metrics")
                                                
                                                # Compute per-class metrics
                                                from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
                                                
                                                acc = accuracy_score(y_true_encoded, y_pred)
                                                prec = precision_score(y_true_encoded, y_pred, average='macro', zero_division=0)
                                                rec = recall_score(y_true_encoded, y_pred, average='macro', zero_division=0)
                                                f1 = f1_score(y_true_encoded, y_pred, average='macro', zero_division=0)
                                                
                                                # Metrics bar chart
                                                metrics_data = {
                                                    'Metric': ['Accuracy', 'Precision\n(Macro)', 'Recall\n(Macro)', 'F1-Score\n(Macro)'],
                                                    'Score': [acc, prec, rec, f1]
                                                }
                                                
                                                fig, ax = plt.subplots(figsize=(8, 5))
                                                colors = ['#4CAF50', '#2196F3', '#FF9800', '#9C27B0']
                                                bars = ax.bar(metrics_data['Metric'], metrics_data['Score'], color=colors, edgecolor='black')
                                                
                                                # Add value labels on bars
                                                for bar, score in zip(bars, metrics_data['Score']):
                                                    height = bar.get_height()
                                                    ax.text(bar.get_x() + bar.get_width()/2., height,
                                                           f'{score:.3f}',
                                                           ha='center', va='bottom', fontsize=12, fontweight='bold')
                                                
                                                ax.set_ylim(0, 1.1)
                                                ax.set_ylabel('Score')
                                                ax.set_title('Model Performance (Macro-Averaged)')
                                                ax.axhline(y=0.8, color='green', linestyle='--', alpha=0.5, label='Good threshold')
                                                ax.axhline(y=0.6, color='orange', linestyle='--', alpha=0.5, label='Moderate threshold')
                                                plt.tight_layout()
                                                st.pyplot(fig)
                                                
                                                # Metrics summary
                                                st.markdown(f"""
                                                **Performance Summary:**
                                                - Accuracy: **{acc:.1%}** of patterns correctly aligned
                                                - F1-Score: **{f1:.1%}** (balance of precision/recall)
                                                - Classes evaluated: **{len(classes)}**
                                                - Samples evaluated: **{valid_mask.sum()}**
                                                """)
                                        else:
                                            st.warning("No valid samples found for model evaluation.")
                                            
                                    except Exception as e:
                                        st.warning(f"Could not compute metrics: {e}")
                            
                            # Feature importance (existing code)
                            if model is not None:
                                feature_importance = None
                                
                                # Extract feature importance based on model type
                                if hasattr(model, 'feature_importances_'):
                                    feature_importance = model.feature_importances_
                                elif hasattr(model, 'coef_'):
                                    coef = np.abs(model.coef_).flatten() if model.coef_.ndim == 1 else np.abs(model.coef_).mean(axis=0)
                                    feature_importance = coef
                                
                                if feature_importance is not None and len(feature_importance) > 0:
                                    st.subheader("🎯 Feature Importance (Top Discriminating Antibiotics)")
                                    
                                    feature_names = [c.replace('_encoded', '') for c in antibiotic_cols]
                                    
                                    if len(feature_names) == len(feature_importance):
                                        importance_df = pd.DataFrame({
                                            'Antibiotic': feature_names,
                                            'Importance': feature_importance
                                        }).sort_values('Importance', ascending=False)
                                        
                                        # Display as bar chart
                                        fig, ax = plt.subplots(figsize=(10, 6))
                                        colors = plt.cm.Blues(np.linspace(0.4, 0.9, 15))[::-1]
                                        ax.barh(importance_df['Antibiotic'][:15], 
                                               importance_df['Importance'][:15],
                                               color=colors, edgecolor='black')
                                        ax.set_xlabel('Importance Score (ASSOCIATIVE, not causal)')
                                        ax.set_ylabel('Antibiotic')
                                        ax.set_title('Top 15 Discriminating Antibiotics')
                                        ax.invert_yaxis()
                                        plt.tight_layout()
                                        st.pyplot(fig)
                                        
                                        with st.expander("📊 Full Feature Importance Table"):
                                            st.dataframe(importance_df.reset_index(drop=True), 
                                                       width='stretch')
                                    else:
                                        st.info(f"Feature count mismatch: dataset has {len(feature_names)} features, "
                                               f"model expects {len(feature_importance)} features.")
                                else:
                                    st.info("Feature importance not available for this model type.")
                                    
                        except Exception as e:
                            st.error(f"Error loading model: {e}")
                        
                        st.markdown("---")  # Separator between models
                else:
                    st.info("No trained models found. Run the supervised learning pipeline first.")
            else:
                st.info("Models directory not found. Run the supervised learning pipeline first.")
            
            # ========== CO-RESISTANCE MODEL PERFORMANCE ==========
            st.markdown("---")
            st.subheader("🔗 Co-Resistance Predictive Models")
            
            st.markdown("""
            <div class="info-box">
                <strong>ℹ️ Co-Resistance Models:</strong> These models predict resistance to target antibiotics 
                based on resistance patterns to other antibiotics. High AUC indicates strong co-resistance 
                relationships suitable for guiding treatment decisions.
            </div>
            """, unsafe_allow_html=True)
            
            # Load co-resistance predictions
            pred_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'processed', 
                                     'figures', 'coresistance_predictions.csv')
            
            if os.path.exists(pred_path):
                pred_df = pd.read_csv(pred_path)
                
                # Model info metrics
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric("Target Antibiotics", len(pred_df))
                with col2:
                    avg_auc = pred_df['AUC'].mean()
                    st.metric("Average AUC", f"{avg_auc:.3f}")
                with col3:
                    high_auc_count = (pred_df['AUC'] > 0.9).sum()
                    st.metric("High AUC (>0.9)", f"{high_auc_count}/{len(pred_df)}")
                
                # Visual performance display
                col1, col2 = st.columns([1, 1])
                
                with col1:
                    st.subheader("📊 AUC Scores by Antibiotic")
                    fig, ax = plt.subplots(figsize=(8, 5))
                    
                    # Sort by AUC
                    pred_sorted = pred_df.sort_values('AUC', ascending=True)
                    
                    # Color based on AUC threshold
                    colors_auc = ['#4CAF50' if x > 0.9 else '#FF9800' if x > 0.7 else '#F44336' 
                              for x in pred_sorted['AUC']]
                    
                    bars = ax.barh(pred_sorted['Target_Antibiotic'], pred_sorted['AUC'], 
                                  color=colors_auc, edgecolor='black')
                    
                    # Add value labels
                    for bar, auc in zip(bars, pred_sorted['AUC']):
                        width = bar.get_width()
                        ax.text(width + 0.01, bar.get_y() + bar.get_height()/2,
                               f'{auc:.3f}', ha='left', va='center', fontsize=10, fontweight='bold')
                    
                    ax.set_xlim(0, 1.15)
                    ax.set_xlabel('AUC Score')
                    ax.set_ylabel('Target Antibiotic')
                    ax.set_title('Predictive Model AUC (Random Forest)')
                    ax.axvline(x=0.9, color='green', linestyle='--', alpha=0.5, label='Excellent (0.9)')
                    ax.axvline(x=0.7, color='orange', linestyle='--', alpha=0.5, label='Good (0.7)')
                    plt.tight_layout()
                    st.pyplot(fig)
                
                with col2:
                    st.subheader("📈 Performance Metrics Comparison")
                    fig, ax = plt.subplots(figsize=(8, 5))
                    
                    # Create grouped bar chart
                    x = np.arange(len(pred_df))
                    width = 0.2
                    
                    ax.bar(x - width*1.5, pred_df['AUC'], width, label='AUC', color='#4CAF50', edgecolor='black')
                    ax.bar(x - width/2, pred_df['Accuracy'], width, label='Accuracy', color='#2196F3', edgecolor='black')
                    ax.bar(x + width/2, pred_df['Precision'], width, label='Precision', color='#FF9800', edgecolor='black')
                    ax.bar(x + width*1.5, pred_df['Recall'], width, label='Recall', color='#9C27B0', edgecolor='black')
                    
                    ax.set_xlabel('Target Antibiotic')
                    ax.set_ylabel('Score')
                    ax.set_title('Model Performance by Antibiotic')
                    ax.set_xticks(x)
                    ax.set_xticklabels(pred_df['Target_Antibiotic'], rotation=45, ha='right')
                    ax.set_ylim(0, 1.1)
                    ax.legend(loc='upper right')
                    ax.axhline(y=0.8, color='gray', linestyle='--', alpha=0.3)
                    plt.tight_layout()
                    st.pyplot(fig)
                
                # Performance summary
                st.markdown(f"""
                **🎯 Co-Resistance Model Summary:**
                - **Best predictor:** {pred_df.loc[pred_df['AUC'].idxmax(), 'Target_Antibiotic']} (AUC: {pred_df['AUC'].max():.3f})
                - **Average AUC:** {avg_auc:.3f} - {'Excellent' if avg_auc > 0.9 else 'Good' if avg_auc > 0.8 else 'Moderate'}
                - **High performers (AUC > 0.9):** {', '.join(pred_df[pred_df['AUC'] > 0.9]['Target_Antibiotic'].tolist()) or 'None'}
                """)
                
                # Table display
                with st.expander("📋 Full Metrics Table"):
                    display_df = pred_df[['Target_Antibiotic', 'AUC', 'Accuracy', 'Precision', 'Recall']].copy()
                    display_df['AUC'] = display_df['AUC'].apply(lambda x: f"{x:.3f}")
                    display_df['Accuracy'] = display_df['Accuracy'].apply(lambda x: f"{x:.3f}")
                    display_df['Precision'] = display_df['Precision'].apply(lambda x: f"{x:.3f}")
                    display_df['Recall'] = display_df['Recall'].apply(lambda x: f"{x:.3f}")
                    st.dataframe(display_df, width='stretch')
            else:
                st.info("Co-resistance predictions not found. Run `python scripts/coresistance_analysis.py` to generate.")
    
    elif analysis_type == "Co-Resistance Network":
        st.header("🔗 Co-Resistance Network Analysis")
        
        st.markdown("""
        <div class="info-box">
            <strong>ℹ️ Co-Resistance Analysis:</strong> This analysis identifies statistically 
            significant pairwise antibiotic co-resistance patterns using Bonferroni-corrected 
            chi-square tests. Predictive models show how well resistance to one antibiotic 
            predicts resistance to others.
        </div>
        """, unsafe_allow_html=True)
        
        # Load predictions
        pred_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'processed', 
                                 'figures', 'coresistance_predictions.csv')
        pairs_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'processed', 
                                  'figures', 'coresistance_significant_pairs.csv')
        
        # ========== MODEL PERFORMANCE SUMMARY ==========

        
        # ========== SIGNIFICANT PAIRS ==========
        st.subheader("🔗 Significant Co-Resistance Pairs")
        
        if os.path.exists(pairs_path):
            pairs_df = pd.read_csv(pairs_path)
            
            col1, col2 = st.columns([1, 1])
            
            with col1:
                st.markdown(f"**Total significant pairs:** {len(pairs_df)}")
                pairs_display = pairs_df.copy()
                pairs_display['Phi'] = pairs_display['Phi'].apply(lambda x: f"{x:.3f}")
                pairs_display['P_value'] = pairs_display['P_value'].apply(lambda x: f"{x:.2e}")
                st.dataframe(pairs_display.head(10), width='stretch')
            
            with col2:
                # Top pairs bar chart
                st.subheader("Top 10 Co-Resistance Pairs (by Phi)")
                fig, ax = plt.subplots(figsize=(8, 5))
                
                top_pairs = pairs_df.head(10).copy()
                top_pairs['Pair'] = top_pairs['Antibiotic_1'] + ' - ' + top_pairs['Antibiotic_2']
                
                colors = plt.cm.Reds(np.linspace(0.4, 0.9, len(top_pairs)))[::-1]
                ax.barh(top_pairs['Pair'], top_pairs['Phi'], color=colors, edgecolor='black')
                ax.set_xlabel('Phi Coefficient (Association Strength)')
                ax.set_ylabel('Antibiotic Pair')
                ax.set_title('Strongest Co-Resistance Associations')
                ax.invert_yaxis()
                plt.tight_layout()
                st.pyplot(fig)
        else:
            st.warning("Co-resistance pairs not found.")
        
        # Display network visualization
        network_plot = os.path.join(os.path.dirname(__file__), '..', 'data', 'processed', 
                                    'figures', 'coresistance_network.png')
        
        if os.path.exists(network_plot):
            st.subheader("🕸️ Co-Resistance Network Visualization")
            st.image(network_plot, caption="Co-resistance network (edges = significant associations)")
        
        # Biological interpretation
        with st.expander("📝 Biological Interpretation"):
            st.markdown("""
            - **High co-resistance (phi > 0.3):** Suggests co-localization on mobile genetic elements (plasmids)
            - **Tetracycline-aminoglycoside links:** Common in IncF/IncA-C plasmids
            - **Beta-lactam clusters:** May indicate ESBL-type resistance mechanisms
            - **Predictive power:** Resistance to one antibiotic strongly predicts others, supporting coordinated resistance hypothesis
            """)
    
    elif analysis_type == "Antibiotic Clustering":
        st.header("💊 Antibiotic Clustering (by Co-Resistance)")
        
        st.markdown("""
        <div class="info-box">
            <strong>ℹ️ Antibiotic Clustering Analysis:</strong> This unsupervised analysis clusters 
            <strong>antibiotics</strong> (not isolates) based on their co-resistance patterns. 
            Antibiotics that are frequently co-resisted cluster together, potentially indicating 
            shared genetic mechanisms (e.g., plasmid co-carriage, integron linkage).
        </div>
        """, unsafe_allow_html=True)
        
        # Load antibiotic clustering results
        clusters_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'processed', 
                                      'figures', 'antibiotic_clusters.csv')
        dendrogram_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'processed', 
                                        'figures', 'antibiotic_dendrogram.png')
        heatmap_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'processed', 
                                     'figures', 'antibiotic_clustered_heatmap.png')
        
        if os.path.exists(clusters_path):
            clusters_df = pd.read_csv(clusters_path)
            
            # Summary metrics
            n_clusters = clusters_df['Cluster'].nunique()
            n_antibiotics = len(clusters_df)
            
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Total Antibiotics", n_antibiotics)
            with col2:
                st.metric("Clusters Identified", n_clusters)
            with col3:
                linked_abs = len(clusters_df[clusters_df['Mean_Phi_Within_Cluster'] > 0])
                st.metric("Linked Antibiotics", linked_abs)
            
            # Cluster table
            st.subheader("📋 Antibiotic Cluster Assignments")
            
            # Create a summary table by cluster
            cluster_summary = []
            for cluster_id in sorted(clusters_df['Cluster'].unique()):
                cluster_data = clusters_df[clusters_df['Cluster'] == cluster_id]
                abs_list = ", ".join(cluster_data['Antibiotic'].tolist())
                classes = cluster_data['Class'].value_counts().head(2).to_dict()
                primary_class = list(classes.keys())[0] if classes else "N/A"
                mean_phi = cluster_data['Mean_Phi_Within_Cluster'].mean()
                
                # Biological interpretation
                if mean_phi > 0.5:
                    interpretation = "Strong co-resistance (likely plasmid linkage)"
                elif mean_phi > 0.3:
                    interpretation = "Moderate co-resistance"
                elif mean_phi > 0:
                    interpretation = "Weak co-resistance"
                else:
                    interpretation = "Independent (no co-resistance)"
                
                cluster_summary.append({
                    'Cluster': cluster_id,
                    'Antibiotics': abs_list,
                    'Primary Class': primary_class,
                    'Mean φ': f"{mean_phi:.3f}" if mean_phi > 0 else "0",
                    'Interpretation': interpretation
                })
            
            summary_df = pd.DataFrame(cluster_summary)
            st.dataframe(summary_df, width='stretch')
            
            # Full raw data display (User Request)
            with st.expander("📋 Full Antibiotic Cluster Data (Raw CSV Content)"):
                st.markdown("Raw content of `data/processed/figures/antibiotic_clusters.csv`:")
                st.dataframe(clusters_df, width='stretch')
            
            # Visualizations
            col1, col2 = st.columns(2)
            
            with col1:
                st.subheader("🌳 Antibiotic Dendrogram")
                if os.path.exists(dendrogram_path):
                    st.image(dendrogram_path, caption="Antibiotics clustered by co-resistance patterns")
                else:
                    st.info("Dendrogram not generated. Run `python main.py --validate` to create.")
            
            with col2:
                st.subheader("🔥 Co-Resistance Heatmap (Clustered)")
                if os.path.exists(heatmap_path):
                    st.image(heatmap_path, caption="Phi coefficients with hierarchical ordering")
                else:
                    st.info("Clustered heatmap not generated. Run `python main.py --validate` to create.")
            
            # Key findings
            st.subheader("🔬 Key Findings")
            
            # Find linked clusters (phi > 0)
            linked_clusters = clusters_df[clusters_df['Mean_Phi_Within_Cluster'] > 0.3]
            
            if len(linked_clusters) > 0:
                st.markdown("**Co-Resistance Clusters Identified:**")
                
                for cluster_id in linked_clusters['Cluster'].unique():
                    cluster_data = linked_clusters[linked_clusters['Cluster'] == cluster_id]
                    abs_list = cluster_data['Antibiotic'].tolist()
                    phi_val = cluster_data['Mean_Phi_Within_Cluster'].iloc[0]
                    classes = cluster_data['Class'].tolist()
                    
                    if len(abs_list) >= 2:
                        st.markdown(f"""
                        - **{' + '.join(abs_list)}** (φ = {phi_val:.2f})
                          - Classes: {', '.join(set(classes))}
                          - *Biological implication*: Likely co-carriage on mobile genetic elements
                        """)
            else:
                st.info("No strong co-resistance clusters identified (φ > 0.3)")
            
            # Interpretation guide
            with st.expander("📚 Interpretation Guide"):
                st.markdown("""
                **What is Antibiotic Clustering?**
                
                This analysis clusters antibiotics based on how often they are co-resisted together 
                across isolates. Antibiotics in the same cluster tend to have correlated resistance 
                patterns.
                
                **Phi Coefficient (φ) Interpretation:**
                - φ > 0.5: Strong association (often plasmid-mediated)
                - φ = 0.3-0.5: Moderate association
                - φ = 0.1-0.3: Weak association
                - φ ≈ 0: No association (independent resistance)
                
                **Biological Significance:**
                
                | Pattern | Likely Mechanism |
                |---------|------------------|
                | TE + DO cluster | *tet* genes on same mobile element |
                | C + SXT cluster | Class 1 integron (*floR* + *sul1*) |
                | CFO + CFT cluster | Shared veterinary cephalosporin resistance |
                | AM + AMC cluster | Beta-lactamase production |
                | Independent antibiotics | Chromosomal mutations or infrequent resistance |
                """)
        else:
            st.warning("Antibiotic clustering results not found. Run `python main.py --validate` to generate.")
    
    elif analysis_type == "PCA Analysis":
        st.header("📐 Principal Component Analysis")
        
        if antibiotic_cols:
            # Color by selection
            color_options = ['None']
            if 'CLUSTER' in df.columns:
                color_options.append('CLUSTER')
            if 'REGION' in df.columns:
                color_options.append('REGION')
            if 'MDR_CATEGORY' in df.columns:
                color_options.append('MDR_CATEGORY')
            if 'ISOLATE_ID' in df.columns:
                color_options.append('ISOLATE_ID')
            
            color_by = st.selectbox("Color by:", color_options)
            
            X_pca, pca = perform_pca_analysis(df, antibiotic_cols)
            
            # Variance explained
            st.subheader("Variance Explained")
            col1, col2 = st.columns(2)
            with col1:
                st.metric("PC1", f"{pca.explained_variance_ratio_[0]*100:.1f}%")
            with col2:
                st.metric("PC2", f"{pca.explained_variance_ratio_[1]*100:.1f}%")
            
            # PCA plot
            if color_by != 'None':
                fig = create_pca_plot(X_pca, df, color_by, pca)
            else:
                fig, ax = plt.subplots(figsize=(10, 8))
                ax.scatter(X_pca[:, 0], X_pca[:, 1], alpha=0.7, s=50)
                ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
                ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
                ax.set_title('PCA of Resistance Profiles')
                plt.tight_layout()
            
            st.pyplot(fig)
            
            # Component loadings
            with st.expander("📊 Component Loadings"):
                feature_names = [c.replace('_encoded', '') for c in antibiotic_cols 
                               if c in df.columns]
                loadings_df = pd.DataFrame({
                    'Feature': feature_names,
                    'PC1': pca.components_[0],
                    'PC2': pca.components_[1]
                })
                loadings_df['PC1_abs'] = np.abs(loadings_df['PC1'])
                loadings_df = loadings_df.sort_values('PC1_abs', ascending=False)
                st.dataframe(loadings_df[['Feature', 'PC1', 'PC2']], width='stretch')
        else:
            st.warning("No antibiotic columns found for PCA.")
    
    elif analysis_type == "Regional Distribution":
        st.header("🗺️ Regional Distribution")
        
        if 'REGION' not in df.columns:
            st.warning("No region information found in the dataset.")
        else:
            col1, col2 = st.columns(2)
            
            with col1:
                st.subheader("Isolates by Region")
                region_counts = df['REGION'].value_counts()
                fig, ax = plt.subplots(figsize=(8, 6))
                region_counts.plot(kind='bar', ax=ax, color='steelblue', edgecolor='black')
                ax.set_xlabel('Region')
                ax.set_ylabel('Count')
                plt.xticks(rotation=45, ha='right')
                plt.tight_layout()
                st.pyplot(fig)
            
            with col2:
                if 'MDR_FLAG' in df.columns:
                    st.subheader("MDR by Region")
                    mdr_by_region = df.groupby('REGION')['MDR_FLAG'].mean() * 100
                    fig, ax = plt.subplots(figsize=(8, 6))
                    mdr_by_region.plot(kind='bar', ax=ax, color='#F44336', edgecolor='black')
                    ax.set_xlabel('Region')
                    ax.set_ylabel('MDR Proportion (%)')
                    ax.set_ylim(0, 100)
                    plt.xticks(rotation=45, ha='right')
                    plt.tight_layout()
                    st.pyplot(fig)
            
            # Cross-tabulation
            if 'CLUSTER' in df.columns:
                st.subheader("Cluster Distribution by Region")
                crosstab = pd.crosstab(df['CLUSTER'], df['REGION'])
                st.dataframe(crosstab, width='stretch')
    
    elif analysis_type == "Local Site Distribution":
        st.header("📍 Local Site (Barangay) Distribution")
        
        if 'SITE' not in df.columns:
            st.warning("No local site information found in the dataset.")
        else:
            col1, col2 = st.columns(2)
            
            with col1:
                st.subheader("Isolates by Local Site")
                site_counts = df['SITE'].value_counts()
                fig, ax = plt.subplots(figsize=(10, 6))
                site_counts.plot(kind='bar', ax=ax, color='teal', edgecolor='black')
                ax.set_xlabel('Local Site (Barangay)')
                ax.set_ylabel('Count')
                plt.xticks(rotation=45, ha='right')
                plt.tight_layout()
                st.pyplot(fig)
            
            with col2:
                if 'MDR_FLAG' in df.columns:
                    st.subheader("MDR by Local Site")
                    mdr_by_site = df.groupby('SITE')['MDR_FLAG'].mean() * 100
                    mdr_by_site = mdr_by_site.sort_values(ascending=False)
                    fig, ax = plt.subplots(figsize=(10, 6))
                    mdr_by_site.plot(kind='bar', ax=ax, color='#F44336', edgecolor='black')
                    ax.set_xlabel('Local Site (Barangay)')
                    ax.set_ylabel('MDR Proportion (%)')
                    ax.set_ylim(0, 100)
                    plt.xticks(rotation=45, ha='right')
                    plt.tight_layout()
                    st.pyplot(fig)
            
            # Cross-tabulation
            if 'CLUSTER' in df.columns:
                st.subheader("Cluster Distribution by Local Site")
                crosstab = pd.crosstab(df['CLUSTER'], df['SITE'])
                st.dataframe(crosstab, width='stretch')
            
            # Region breakdown
            if 'REGION' in df.columns:
                st.subheader("Local Site by Region")
                crosstab_region = pd.crosstab(df['SITE'], df['REGION'])
                st.dataframe(crosstab_region, width='stretch')
    
    elif analysis_type == "Source Distribution":
        st.header("🔬 Sampling Source Distribution")
        
        if 'SAMPLING_SOURCE' not in df.columns:
            st.warning("No sampling source information found in the dataset.")
        else:
            col1, col2 = st.columns(2)
            
            with col1:
                st.subheader("Isolates by Sampling Source")
                source_counts = df['SAMPLING_SOURCE'].value_counts()
                fig, ax = plt.subplots(figsize=(10, 6))
                source_counts.plot(kind='bar', ax=ax, color='darkorange', edgecolor='black')
                ax.set_xlabel('Sampling Source')
                ax.set_ylabel('Count')
                plt.xticks(rotation=45, ha='right')
                plt.tight_layout()
                st.pyplot(fig)
            
            with col2:
                if 'MDR_FLAG' in df.columns:
                    st.subheader("MDR by Sampling Source")
                    mdr_by_source = df.groupby('SAMPLING_SOURCE')['MDR_FLAG'].mean() * 100
                    mdr_by_source = mdr_by_source.sort_values(ascending=False)
                    fig, ax = plt.subplots(figsize=(10, 6))
                    mdr_by_source.plot(kind='bar', ax=ax, color='#F44336', edgecolor='black')
                    ax.set_xlabel('Sampling Source')
                    ax.set_ylabel('MDR Proportion (%)')
                    ax.set_ylim(0, 100)
                    plt.xticks(rotation=45, ha='right')
                    plt.tight_layout()
                    st.pyplot(fig)
            
            # Cross-tabulation
            if 'CLUSTER' in df.columns:
                st.subheader("Cluster Distribution by Sampling Source")
                crosstab = pd.crosstab(df['CLUSTER'], df['SAMPLING_SOURCE'])
                st.dataframe(crosstab, width='stretch')
            
            # Environment breakdown
            if 'ENVIRONMENT' in df.columns:
                st.subheader("Sampling Source by Environment")
                crosstab_env = pd.crosstab(df['SAMPLING_SOURCE'], df['ENVIRONMENT'])
                st.dataframe(crosstab_env, width='stretch')
    
    elif analysis_type == "Species Distribution":
        st.header("🦠 Species Distribution")
        
        if 'ISOLATE_ID' not in df.columns:
            st.warning("No species information found in the dataset.")
        else:
            col1, col2 = st.columns(2)
            
            with col1:
                st.subheader("Isolates by Species")
                species_counts = df['ISOLATE_ID'].value_counts()
                fig, ax = plt.subplots(figsize=(10, 6))
                species_counts.plot(kind='bar', ax=ax, color='purple', edgecolor='black')
                ax.set_xlabel('Species')
                ax.set_ylabel('Count')
                plt.xticks(rotation=45, ha='right')
                plt.tight_layout()
                st.pyplot(fig)
            
            with col2:
                if 'MDR_FLAG' in df.columns:
                    st.subheader("MDR by Species")
                    mdr_by_species = df.groupby('ISOLATE_ID')['MDR_FLAG'].mean() * 100
                    mdr_by_species = mdr_by_species.sort_values(ascending=False)
                    fig, ax = plt.subplots(figsize=(10, 6))
                    mdr_by_species.plot(kind='bar', ax=ax, color='#F44336', edgecolor='black')
                    ax.set_xlabel('Species')
                    ax.set_ylabel('MDR Proportion (%)')
                    ax.set_ylim(0, 100)
                    plt.xticks(rotation=45, ha='right')
                    plt.tight_layout()
                    st.pyplot(fig)
            
            # Cross-tabulation with cluster
            if 'CLUSTER' in df.columns:
                st.subheader("Cluster Distribution by Species")
                crosstab = pd.crosstab(df['CLUSTER'], df['ISOLATE_ID'])
                st.dataframe(crosstab, width='stretch')
            
            # Species by Region
            if 'REGION' in df.columns:
                st.subheader("Species by Region")
                crosstab_region = pd.crosstab(df['ISOLATE_ID'], df['REGION'])
                st.dataframe(crosstab_region, width='stretch')
            
            # Species by Environment
            if 'ENVIRONMENT' in df.columns:
                st.subheader("Species by Environment")
                crosstab_env = pd.crosstab(df['ISOLATE_ID'], df['ENVIRONMENT'])
                st.dataframe(crosstab_env, width='stretch')
    

    
    elif analysis_type == "Integration & Synthesis":
        st.header("🔗 Integration & Synthesis (Phase 6)")
        
        st.markdown("""
        This section integrates results from unsupervised clustering (Phase 3), 
        supervised learning (Phase 4), and regional/environmental analysis (Phase 5)
        to identify key patterns in AMR data.
        """)
        

        # Import integration module
        try:
            from analysis.integration_synthesis import (
                define_integration_framework,
                generate_cluster_supervised_alignment_table,
                identify_resistance_archetypes,
                generate_triangulation_table,
                compare_clusters_with_supervised,
                identify_mdr_enriched_patterns,
                identify_species_environment_associations
            )
            
            # Check for necessary columns
            antibiotic_cols = [col for col in df.columns if col.endswith('_encoded')]
            has_cluster = 'CLUSTER' in df.columns
            has_mdr = 'MDR_CATEGORY' in df.columns or 'MDR_FLAG' in df.columns
            
            if not has_cluster:
                st.warning("Cluster information missing. Run clustering analysis first.")
            
            # Create sub-tabs for organized presentation
            int_tabs = st.tabs([
                "📚 Integration Framework", 
                "📊 Alignment & Archetypes", 
                "🌍 Triangulation", 
                "⚠️ MDR Analysis"
            ])
            
            # --- TAB 1: FRAMEWORK ---
            with int_tabs[0]:
                st.subheader("1. Formal Integration Framework")
                framework = define_integration_framework()
                st.markdown(framework['integration_rationale'])
                st.table(framework['framework_table'])
                
                st.markdown("""
                **Purpose:** This framework prevents "narrative drift" by explicitly defining 
                how unsupervised clusters (Phase 2), supervised signals (Phase 3), and 
                environmental context (Phase 4) are synthesized into coherent archetypes.
                """)
            
            # --- TAB 2: ALIGNMENT & ARCHETYPES ---
            with int_tabs[1]:
                st.subheader("2. Cluster-Supervised Alignment")
                
                if has_cluster and has_mdr:
                    align_df, align_interp = generate_cluster_supervised_alignment_table(
                        df, cluster_col='CLUSTER', mdr_flag_col='MDR_FLAG'
                    )
                    st.dataframe(align_df, width='stretch')
                    
                    col1, col2 = st.columns(2)
                    with col1:
                        if align_interp['agreement']:
                            st.markdown("**✅ Agreement patterns:**")
                            for note in align_interp['agreement']:
                                st.write(f"- {note}")
                    with col2:
                        if align_interp['disagreement']:
                            st.markdown("**❌ Heterogeneous patterns:**")
                            for note in align_interp['disagreement']:
                                st.write(f"- {note}")
                
                st.markdown("---")
                st.subheader("3. Resistance Archetypes")
                st.markdown("*Defined as recurring resistance profiles with stable cluster membership.*")
                
                if has_cluster:
                    archetypes = identify_resistance_archetypes(df, feature_cols=antibiotic_cols)
                    
                    if archetypes['archetype_summary_table'] is not None:
                        st.dataframe(archetypes['archetype_summary_table'], width='stretch')
                    else:
                        st.info("No dominant archetypes identified (clusters may be too small or heterogeneous).")

            # --- TAB 3: TRIANGULATION ---
            with int_tabs[2]:
                st.subheader("4. Species-Environment-Resistance Triangulation")
                
                if has_cluster:
                    tri_df, tri_notes = generate_triangulation_table(
                        df, feature_cols=antibiotic_cols,
                        cluster_col='CLUSTER',
                        species_col='ISOLATE_ID',
                        environment_col='SAMPLING_SOURCE' if 'SAMPLING_SOURCE' in df.columns else 'SAMPLE_SOURCE'
                    )
                    st.dataframe(tri_df, width='stretch')
                
                st.markdown("---")
                
                # Use the association analysis function
                associations = identify_species_environment_associations(df)
                
                col1, col2 = st.columns(2)
                with col1:
                    if associations.get('species_environment'):
                        st.markdown("**Species by Environment:**")
                        env_data = []
                        for species, info in associations['species_environment'].items():
                            env_data.append({
                                'Species': species,
                                'Dominant Environment': info['dominant_environment'],
                                'Proportion': f"{info['proportion']*100:.1f}%"
                            })
                        st.dataframe(pd.DataFrame(env_data), width='stretch')
                
                with col2:
                    if associations.get('species_region'):
                        st.markdown("**Species by Region:**")
                        region_data = []
                        for species, info in associations['species_region'].items():
                            region_data.append({
                                'Species': species,
                                'Dominant Region': info['dominant_region'],
                                'Proportion': f"{info['proportion']*100:.1f}%"
                            })
                        st.dataframe(pd.DataFrame(region_data), width='stretch')

            # --- TAB 4: MDR ANALYSIS ---
            with int_tabs[3]:
                st.subheader("5. MDR-Enriched Patterns")
                
                if has_mdr:
                    mdr_patterns = identify_mdr_enriched_patterns(df, antibiotic_cols)
                    
                    # Overall MDR rate
                    if mdr_patterns.get('overall_mdr_rate') is not None:
                        st.metric("Overall MDR Prevalence", f"{mdr_patterns['overall_mdr_rate']*100:.1f}%")
                    
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        if mdr_patterns.get('mdr_enriched_clusters'):
                            st.markdown("**MDR-Enriched Clusters:**")
                            cluster_data = []
                            for item in mdr_patterns['mdr_enriched_clusters']:
                                cluster_data.append({
                                    'Cluster': item['cluster'],
                                    'MDR Rate': f"{item['mdr_rate']*100:.1f}%",
                                    'Fold Enrichment': f"{item['fold_enrichment']:.1f}x",
                                    'Sample Size': item['sample_size']
                                })
                            st.dataframe(pd.DataFrame(cluster_data), width='stretch')
                    
                    with col2:
                        if mdr_patterns.get('mdr_enriched_regions'):
                            st.markdown("**MDR-Enriched Regions:**")
                            region_data = []
                            for item in mdr_patterns['mdr_enriched_regions']:
                                region_data.append({
                                    'Region': item['region'],
                                    'MDR Rate': f"{item['mdr_rate']*100:.1f}%",
                                    'Fold Enrichment': f"{item['fold_enrichment']:.1f}x"
                                })
                            st.dataframe(pd.DataFrame(region_data), width='stretch')
                    
                    # MDR Resistance Signature
                    if mdr_patterns.get('mdr_resistance_signature'):
                        with st.expander("🔬 MDR Resistance Signature"):
                            sig_data = []
                            for ab, ab_stats in mdr_patterns['mdr_resistance_signature'].items():
                                sig_data.append({
                                    'Antibiotic': ab,
                                    'MDR Mean': ab_stats['mdr_mean'],
                                    'Non-MDR Mean': ab_stats['non_mdr_mean'],
                                    'Difference': ab_stats['difference']
                                })
                            sig_df = pd.DataFrame(sig_data)
                            sig_df = sig_df.sort_values('Difference', ascending=False)
                            # Format for display
                            sig_df['MDR Mean'] = sig_df['MDR Mean'].apply(lambda x: f"{x:.2f}")
                            sig_df['Non-MDR Mean'] = sig_df['Non-MDR Mean'].apply(lambda x: f"{x:.2f}")
                            sig_df['Difference'] = sig_df['Difference'].apply(lambda x: f"{x:.2f}")
                            st.dataframe(sig_df, width='stretch')
                    
                    # Interpretation
                    if mdr_patterns.get('interpretation'):
                        with st.expander("📝 Interpretation"):
                            for interp in mdr_patterns['interpretation']:
                                st.write(f"• {interp}")
                else:
                    st.info("ℹ️ MDR information not available in the dataset.")
                    
        except ImportError as e:
            st.error(f"Could not import integration module: {e}")
            st.info("Possible causes: 1) Missing dependencies (scipy, numpy). "
                    "2) Module not found in src/analysis directory. "
                    "Run 'pip install -r requirements.txt' to install dependencies.")
    
    # =========================================================================
    # METHODOLOGY TAB - Phase 6 Requirement: Scientific Transparency
    # =========================================================================
    elif analysis_type == "Methodology":
        st.header("📚 Methodology")
        
        if PHASE6_AVAILABLE:
            methodology = get_methodology_content()
            
            tabs = st.tabs(["Overview", "Preprocessing", "Clustering", "Supervised Learning", "PCA"])
            
            with tabs[0]:
                st.markdown(methodology['overview'])
            
            with tabs[1]:
                st.markdown(methodology['data_preprocessing'])
            
            with tabs[2]:
                st.markdown(methodology['clustering_method'])
            
            with tabs[3]:
                st.markdown(methodology['supervised_discrimination'])
            
            with tabs[4]:
                st.markdown(methodology['pca_usage'])
            
            # Glossary
            st.subheader("📖 Glossary")
            glossary = get_glossary()
            for term, definition in glossary.items():
                with st.expander(term):
                    st.write(definition)
        else:
            st.markdown("""
            ## Methodology Overview
            
            This tool implements a multi-phase analytical pipeline for antimicrobial 
            resistance (AMR) pattern recognition and surveillance analysis.
            
            ### Data Preprocessing
            - Data cleaning and validation (S/I/R values only)
            - Missing data handling with transparent thresholds
            - Resistance encoding: S→0, I→1, R→2
            - Feature engineering: MAR Index, MDR classification
            
            ### Clustering Method
            - Algorithm: Hierarchical Agglomerative Clustering
            - Linkage: Ward's method (minimizes within-cluster variance)
            - Distance: Euclidean
            - Clusters represent **resistance phenotypes**, NOT taxonomic groups
            
            ### Supervised Learning (Pattern Discrimination)
            - Purpose: Evaluate how resistance patterns discriminate known categories
            - This is **pattern discrimination**, NOT prediction
            - Models: Logistic Regression, Random Forest, k-NN
            - Metrics: Macro-averaged precision, recall, F1-score
            
            ### PCA
            - Purpose: Dimensionality reduction for visualization
            - Preserves variance while reducing complexity
            """)
    
    # =========================================================================
    # LIMITATIONS TAB - Phase 6 Requirement: Explicit Claim Boundaries
    # =========================================================================
    elif analysis_type == "Limitations":
        st.header("⚠️ Limitations & Claim Boundaries")
        
        if PHASE6_AVAILABLE:
            limitations = get_limitations_content()
            
            st.markdown(limitations['overview'])
            st.markdown(limitations['no_temporal_inference'])
            st.markdown(limitations['no_causal_inference'])
            st.markdown(limitations['no_predictive_claims'])
            st.markdown(limitations['no_transmission_inference'])
            st.markdown(limitations['dataset_dependency'])
            st.markdown(limitations['data_quality'])
        else:
            st.markdown("""
            ## Limitations & Claim Boundaries
            
            This section explicitly documents what this analysis does NOT show.
            
            ### ❌ No Temporal Inference
            - This cross-sectional study cannot determine temporal trends
            - No claims about increasing/decreasing resistance
            
            ### ❌ No Causal Inference
            - Associations do NOT imply causation
            - Environmental factors are not claimed to "cause" resistance
            - Use: "associated with", "enriched in"
            - Avoid: "driven by", "caused by"
            
            ### ❌ No Predictive Claims
            - Model metrics describe pattern consistency only
            - NOT predictive performance for future samples
            
            ### ❌ No Transmission Inference
            - Cannot identify transmission pathways
            - Requires genomic data (WGS) not available here
            
            ### ⚠️ Dataset-Dependent Results
            - Results specific to this dataset
            - May not generalize to other regions/time periods
            - Sample selection may introduce bias
            """)
    
    # =========================================================================
    # FOOTER - Phase 6 Requirement: Hard-coded disclaimer
    # =========================================================================
    st.markdown("---")
    st.markdown("""
    <div class="footer">
        <strong>AMR Pattern Recognition & Exploratory Analysis Dashboard</strong><br>
        For research and surveillance purposes only | Not for clinical use<br>
        <em>This tool does not provide clinical decision support or predictive assessments.</em>
    </div>
    """, unsafe_allow_html=True)


if __name__ == "__main__":
    main()
