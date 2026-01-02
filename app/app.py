"""
AMR Pattern Recognition & Exploratory Analysis Dashboard
Phase 6 - Concrete Implementation Improvements

This tool is intended exclusively for exploratory antimicrobial resistance 
pattern recognition and surveillance analysis. It does not provide clinical 
decision support or predictive assessments.

DEPLOYMENT SCOPE:
- Pattern Recognition & Exploratory Analysis Dashboard
- NOT a predictive tool
- NOT for clinical decision support
"""

import streamlit as st
import pandas as pd
import numpy as np
import os
import sys

# Add app directory to path
app_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, app_dir)

# Import modular components
from data_loader import (
    load_uploaded_file, load_default_dataset, get_dataset_info,
    get_antibiotic_columns, display_upload_status, display_expected_format,
    validate_csv_schema
)
from visualization import (
    create_resistance_heatmap, create_heatmap_with_dendrogram,
    perform_pca_analysis, create_pca_plot, create_cluster_summary_table,
    create_feature_importance_plot, create_mdr_distribution_plot,
    create_regional_distribution_plot
)
from supervised_models import (
    load_pretrained_model, get_available_models, get_model_info,
    extract_feature_importance, create_confusion_matrix_plot,
    format_metrics_table, get_model_disclaimer, get_feature_importance_disclaimer
)
from interpretation import (
    get_methodology_content, get_limitations_content, get_disclaimers,
    get_about_content, get_glossary
)


# Page configuration
st.set_page_config(
    page_title="AMR Pattern Recognition Dashboard",
    page_icon="🦠",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS
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
    .disclaimer-box {
        background-color: #fff3cd;
        border: 2px solid #ffc107;
        border-radius: 10px;
        padding: 15px;
        margin: 15px 0;
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
        border: 1px solid #bee5eb;
        border-radius: 5px;
        padding: 10px;
        margin: 10px 0;
    }
</style>
""", unsafe_allow_html=True)


def show_main_disclaimer():
    """Display the main disclaimer prominently."""
    disclaimers = get_disclaimers()
    st.markdown(disclaimers['landing_page'], unsafe_allow_html=True)


def show_footer():
    """Display the footer disclaimer."""
    disclaimers = get_disclaimers()
    st.markdown(disclaimers['footer'], unsafe_allow_html=True)


def main():
    # =========================================================================
    # HEADER
    # =========================================================================
    st.markdown('<h1 class="main-header">🦠 AMR Pattern Recognition & Exploratory Analysis Dashboard</h1>', 
                unsafe_allow_html=True)
    st.markdown('<p class="sub-header">Interactive Tool for Antimicrobial Resistance Pattern Analysis</p>',
                unsafe_allow_html=True)
    
    # Main disclaimer on landing page
    show_main_disclaimer()
    
    # =========================================================================
    # SIDEBAR - Data Upload
    # =========================================================================
    st.sidebar.header("📁 Data Upload")
    
    uploaded_file = st.sidebar.file_uploader(
        "Upload CSV Dataset",
        type=['csv'],
        help="Upload your AMR dataset in CSV format. See 'Expected Format' in the main panel."
    )
    
    # Data loading
    df = None
    default_path = os.path.join(app_dir, '..', 'data', 'processed', 'analysis_ready_dataset.csv')
    
    if uploaded_file is not None:
        df, errors, warnings = load_uploaded_file(uploaded_file)
        if errors:
            display_upload_status(errors, warnings)
        elif warnings:
            display_upload_status([], warnings)
            st.sidebar.success(f"✅ Loaded {len(df)} isolates from uploaded file")
        else:
            st.sidebar.success(f"✅ Loaded {len(df)} isolates from uploaded file")
    else:
        df, errors, warnings = load_default_dataset(default_path)
        if df is not None:
            st.sidebar.info(f"📊 Using default dataset ({len(df)} isolates)")
        else:
            st.sidebar.info("👆 Upload a CSV file to begin analysis")
    
    # =========================================================================
    # SIDEBAR - Navigation
    # =========================================================================
    st.sidebar.header("📊 Analysis")
    
    analysis_options = [
        "Overview",
        "Resistance Heatmap",
        "Cluster Analysis",
        "PCA Analysis",
        "Regional Distribution",
        "Model Evaluation",
        "Methodology",
        "Limitations"
    ]
    
    analysis_type = st.sidebar.selectbox(
        "Select Analysis",
        analysis_options,
        help="Choose an analysis to view"
    )
    
    # =========================================================================
    # MAIN CONTENT
    # =========================================================================
    
    if df is None:
        st.header("📋 Getting Started")
        display_expected_format()
        st.markdown(get_about_content())
        show_footer()
        return
    
    # Get dataset info
    dataset_info = get_dataset_info(df)
    antibiotic_cols = dataset_info['antibiotic_columns']
    
    # =========================================================================
    # OVERVIEW
    # =========================================================================
    if analysis_type == "Overview":
        st.header("📋 Dataset Overview")
        
        # Metrics row
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.metric("Total Isolates", dataset_info['n_isolates'])
        
        with col2:
            st.metric("Antibiotics", dataset_info['n_antibiotics'])
        
        with col3:
            if dataset_info.get('n_species'):
                st.metric("Species", dataset_info['n_species'])
            elif dataset_info.get('n_clusters'):
                st.metric("Clusters", dataset_info['n_clusters'])
        
        with col4:
            if dataset_info.get('mdr_proportion') is not None:
                st.metric("MDR Proportion", f"{dataset_info['mdr_proportion']:.1f}%")
        
        # Data preview
        st.subheader("Data Preview")
        st.dataframe(df.head(20), use_container_width=True)
        
        # Column information
        with st.expander("📚 Column Information"):
            col_info = pd.DataFrame({
                'Column': df.columns,
                'Type': df.dtypes.astype(str),
                'Non-Null': df.notna().sum(),
                'Unique': df.nunique()
            })
            st.dataframe(col_info, use_container_width=True)
    
    # =========================================================================
    # RESISTANCE HEATMAP
    # =========================================================================
    elif analysis_type == "Resistance Heatmap":
        st.header("🔥 Resistance Profile Heatmap")
        
        st.markdown("""
        <div class="info-box">
            <strong>ℹ️ About this visualization:</strong> This heatmap shows resistance patterns 
            across isolates and antibiotics. The visualization is <strong>read-only</strong> - 
            it displays pre-computed results without modification.
        </div>
        """, unsafe_allow_html=True)
        
        if antibiotic_cols:
            # Options
            col1, col2 = st.columns(2)
            with col1:
                show_dendrogram = st.checkbox("Show with Dendrogram", value=True,
                                             help="Display heatmap with hierarchical clustering dendrogram")
            with col2:
                if 'CLUSTER' in df.columns:
                    sort_by_cluster = st.checkbox("Sort by Cluster", value=True)
                else:
                    sort_by_cluster = False
            
            if show_dendrogram:
                fig = create_heatmap_with_dendrogram(df, antibiotic_cols)
            else:
                cluster_col = 'CLUSTER' if sort_by_cluster and 'CLUSTER' in df.columns else None
                fig = create_resistance_heatmap(df, antibiotic_cols, cluster_col)
            
            st.pyplot(fig)
            
            # Legend
            st.markdown("""
            **Legend:**
            - 🟢 **Green (0)**: Susceptible (S)
            - 🟡 **Yellow (1)**: Intermediate (I)
            - 🔴 **Red (2)**: Resistant (R)
            """)
        else:
            st.warning("No antibiotic columns found in the dataset.")
    
    # =========================================================================
    # CLUSTER ANALYSIS
    # =========================================================================
    elif analysis_type == "Cluster Analysis":
        st.header("🎯 Cluster Analysis")
        
        if 'CLUSTER' not in df.columns:
            st.warning("⚠️ No cluster information found. Run the analysis pipeline first.")
        else:
            # Cluster selection dropdown
            clusters = sorted(df['CLUSTER'].unique())
            selected_clusters = st.multiselect(
                "Select Clusters to View",
                options=['All'] + [f'C{int(c)}' for c in clusters],
                default=['All'],
                help="Filter to specific clusters"
            )
            
            # Filter data
            if 'All' not in selected_clusters:
                selected_nums = [int(c.replace('C', '')) for c in selected_clusters]
                df_filtered = df[df['CLUSTER'].isin(selected_nums)]
            else:
                df_filtered = df
            
            # Cluster summary table
            st.subheader("Cluster Summary")
            summary_table = create_cluster_summary_table(df_filtered, antibiotic_cols)
            if not summary_table.empty:
                st.dataframe(summary_table, use_container_width=True)
            
            # Visualizations
            col1, col2 = st.columns(2)
            
            with col1:
                st.subheader("Cluster Size Distribution")
                cluster_counts = df_filtered['CLUSTER'].value_counts().sort_index()
                import matplotlib.pyplot as plt
                fig, ax = plt.subplots(figsize=(8, 6))
                ax.bar([f'C{int(c)}' for c in cluster_counts.index], cluster_counts.values, 
                      color='steelblue', edgecolor='black')
                ax.set_xlabel('Cluster')
                ax.set_ylabel('Count')
                ax.set_title('Isolates per Cluster')
                plt.tight_layout()
                st.pyplot(fig)
            
            with col2:
                fig = create_mdr_distribution_plot(df_filtered)
                if fig:
                    st.subheader("MDR by Cluster")
                    st.pyplot(fig)
    
    # =========================================================================
    # PCA ANALYSIS
    # =========================================================================
    elif analysis_type == "PCA Analysis":
        st.header("📐 Principal Component Analysis")
        
        st.markdown("""
        <div class="info-box">
            <strong>ℹ️ About PCA:</strong> Principal Component Analysis reduces the 
            dimensionality of resistance data for visualization. Points closer together 
            have more similar resistance profiles.
        </div>
        """, unsafe_allow_html=True)
        
        if antibiotic_cols:
            # Color by selection
            color_options = ['None']
            if 'CLUSTER' in df.columns:
                color_options.append('CLUSTER')
            if 'REGION' in df.columns:
                color_options.append('REGION')
            if 'MDR_CATEGORY' in df.columns:
                color_options.append('MDR_CATEGORY')
            if 'SAMPLE_SOURCE' in df.columns:
                color_options.append('SAMPLE_SOURCE (Environment)')
            if 'ISOLATE_ID' in df.columns:
                color_options.append('ISOLATE_ID (Species)')
            
            col1, col2 = st.columns(2)
            with col1:
                color_by = st.selectbox("Color by:", color_options, index=1 if len(color_options) > 1 else 0,
                                       help="Select variable for coloring points")
            
            # Map display name to actual column
            color_col_map = {
                'SAMPLE_SOURCE (Environment)': 'SAMPLE_SOURCE',
                'ISOLATE_ID (Species)': 'ISOLATE_ID'
            }
            actual_color_col = color_col_map.get(color_by, color_by)
            
            # Perform PCA
            X_pca, pca, used_cols = perform_pca_analysis(df, antibiotic_cols)
            
            # Variance explained
            st.subheader("Variance Explained")
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("PC1", f"{pca.explained_variance_ratio_[0]*100:.1f}%")
            with col2:
                st.metric("PC2", f"{pca.explained_variance_ratio_[1]*100:.1f}%")
            with col3:
                st.metric("Total", f"{sum(pca.explained_variance_ratio_)*100:.1f}%")
            
            # PCA plot
            if actual_color_col != 'None' and actual_color_col in df.columns:
                fig = create_pca_plot(X_pca, df, actual_color_col, pca)
            else:
                import matplotlib.pyplot as plt
                fig, ax = plt.subplots(figsize=(10, 8))
                ax.scatter(X_pca[:, 0], X_pca[:, 1], alpha=0.7, s=50, c='steelblue')
                ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
                ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
                ax.set_title('PCA of Resistance Profiles')
                ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
                ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
            
            st.pyplot(fig)
            
            # Component loadings
            with st.expander("📊 Component Loadings (Top Contributors)"):
                feature_names = [c.replace('_encoded', '') for c in used_cols]
                loadings_df = pd.DataFrame({
                    'Antibiotic': feature_names,
                    'PC1 Loading': pca.components_[0],
                    'PC2 Loading': pca.components_[1]
                })
                loadings_df['PC1 |Loading|'] = np.abs(loadings_df['PC1 Loading'])
                loadings_df = loadings_df.sort_values('PC1 |Loading|', ascending=False)
                st.dataframe(loadings_df[['Antibiotic', 'PC1 Loading', 'PC2 Loading']].head(10), 
                           use_container_width=True)
        else:
            st.warning("No antibiotic columns found for PCA.")
    
    # =========================================================================
    # REGIONAL DISTRIBUTION
    # =========================================================================
    elif analysis_type == "Regional Distribution":
        st.header("🗺️ Regional Distribution")
        
        if 'REGION' not in df.columns:
            st.warning("⚠️ No region information found in the dataset.")
        else:
            # Filter controls
            if 'REGION' in df.columns:
                regions = df['REGION'].unique().tolist()
                selected_regions = st.multiselect(
                    "Filter by Region",
                    options=['All'] + regions,
                    default=['All']
                )
                
                if 'All' not in selected_regions:
                    df_filtered = df[df['REGION'].isin(selected_regions)]
                else:
                    df_filtered = df
            else:
                df_filtered = df
            
            col1, col2 = st.columns(2)
            
            with col1:
                fig = create_regional_distribution_plot(df_filtered)
                if fig:
                    st.subheader("Isolates by Region")
                    st.pyplot(fig)
            
            with col2:
                if 'MDR_FLAG' in df_filtered.columns:
                    st.subheader("MDR by Region")
                    fig = create_mdr_distribution_plot(df_filtered, group_col='REGION')
                    if fig:
                        st.pyplot(fig)
            
            # Cross-tabulation
            if 'CLUSTER' in df.columns:
                st.subheader("Cluster Distribution by Region")
                crosstab = pd.crosstab(df_filtered['CLUSTER'], df_filtered['REGION'])
                crosstab.index = [f'C{int(i)}' for i in crosstab.index]
                st.dataframe(crosstab, use_container_width=True)
    
    # =========================================================================
    # MODEL EVALUATION (Read-Only)
    # =========================================================================
    elif analysis_type == "Model Evaluation":
        st.header("🤖 Model Evaluation Summary")
        
        # Disclaimer near model results
        st.markdown("""
        <div class="info-box">
            <strong>📊 Interpretation Note:</strong> Metrics shown are <em>pattern consistency measures</em>, 
            not predictive performance. They quantify how resistance patterns align with known categories 
            within this dataset only. <strong>No retraining or single-sample inference is available.</strong>
        </div>
        """, unsafe_allow_html=True)
        
        st.markdown(get_model_disclaimer())
        
        # Check for saved models
        models_dir = os.path.join(app_dir, '..', 'data', 'models')
        available_models = get_available_models(models_dir)
        
        if not available_models:
            st.info("ℹ️ No pretrained models found. Run the supervised learning pipeline first.")
        else:
            st.subheader("Available Pretrained Models")
            
            for model_info in available_models:
                with st.expander(f"📁 {model_info['name']} ({model_info['task_type']})"):
                    model_data = load_pretrained_model(model_info['path'])
                    
                    if model_data:
                        info = get_model_info(model_data)
                        
                        col1, col2 = st.columns(2)
                        with col1:
                            st.write(f"**Model Type:** {info['model_type']}")
                            st.write(f"**Task:** {model_info['model_type']}")
                        with col2:
                            if info['classes']:
                                st.write(f"**Classes:** {', '.join(str(c) for c in info['classes'])}")
                            if info['n_features']:
                                st.write(f"**Features:** {info['n_features']}")
                        
                        # Feature importance (read-only)
                        if info['has_feature_importance']:
                            st.subheader("Feature Importance (Read-Only)")
                            st.markdown(get_feature_importance_disclaimer())
                            
                            importance_scores = extract_feature_importance(model_data, antibiotic_cols)
                            
                            if importance_scores:
                                fig = create_feature_importance_plot(
                                    importance_scores,
                                    title=f'{model_info["name"]} - Feature Importance'
                                )
                                st.pyplot(fig)
                                
                                # Table view
                                with st.expander("📊 Full Feature Importance Table"):
                                    imp_df = pd.DataFrame([
                                        {'Antibiotic': k.replace('_encoded', ''), 
                                         'Importance': round(v, 4)}
                                        for k, v in importance_scores.items()
                                    ])
                                    st.dataframe(imp_df, use_container_width=True)
                    else:
                        st.error(f"Error loading model from {model_info['path']}")
    
    # =========================================================================
    # METHODOLOGY
    # =========================================================================
    elif analysis_type == "Methodology":
        st.header("📚 Methodology")
        
        methodology = get_methodology_content()
        
        # Tabs for different methodology sections
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
    
    # =========================================================================
    # LIMITATIONS
    # =========================================================================
    elif analysis_type == "Limitations":
        st.header("⚠️ Limitations & Claim Boundaries")
        
        limitations = get_limitations_content()
        
        st.markdown(limitations['overview'])
        
        # Display all limitations
        st.markdown(limitations['no_temporal_inference'])
        st.markdown(limitations['no_causal_inference'])
        st.markdown(limitations['no_predictive_claims'])
        st.markdown(limitations['no_transmission_inference'])
        st.markdown(limitations['dataset_dependency'])
        st.markdown(limitations['data_quality'])
    
    # =========================================================================
    # FOOTER
    # =========================================================================
    show_footer()


if __name__ == "__main__":
    main()
