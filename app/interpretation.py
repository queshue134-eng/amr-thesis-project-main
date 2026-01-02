"""
Interpretation Module for AMR Pattern Recognition Dashboard
Phase 6 - Scientific Transparency and Interpretation

This module provides:
- Methodology explanation
- Limitations documentation
- Interpretation guidance
- Claim boundaries
"""


def get_methodology_content() -> dict:
    """
    Get methodology explanation content.
    
    Returns:
    --------
    dict
        Dictionary with methodology sections
    """
    return {
        'overview': """
        ## Methodology Overview
        
        This tool implements a multi-phase analytical pipeline for antimicrobial 
        resistance (AMR) pattern recognition and surveillance analysis. The methodology 
        follows established scientific practices for exploratory data analysis.
        """,
        
        'data_preprocessing': """
        ### Data Preprocessing
        
        **Step 1: Data Ingestion**
        - CSV files containing antimicrobial susceptibility testing (AST) results
        - Metadata extraction (region, site, sample source)
        - Standardization of species and antibiotic names
        
        **Step 2: Data Cleaning**
        - Validation of resistance values (only S, I, R allowed)
        - Duplicate detection and removal
        - Missing data handling with transparent thresholds
        
        **Step 3: Resistance Encoding**
        - Susceptible (S) → 0
        - Intermediate (I) → 1
        - Resistant (R) → 2
        
        This ordinal encoding preserves the biological meaning of resistance levels 
        while enabling numerical analysis.
        
        **Step 4: Feature Engineering**
        - MAR Index calculation (Krumperman, 1983)
        - MDR classification based on Magiorakos et al. (2012) criteria
        - Per-antibiotic binary resistance indicators
        """,
        
        'clustering_method': """
        ### Clustering Method
        
        **Algorithm: Hierarchical Agglomerative Clustering**
        
        | Parameter | Value | Rationale |
        |-----------|-------|-----------|
        | Linkage | Ward's method | Minimizes within-cluster variance |
        | Distance | Euclidean | Standard for numerical resistance data |
        | Imputation | Median | Robust to outliers |
        
        **Important Notes:**
        - Clustering uses ONLY resistance data (no metadata)
        - Metadata (species, region) is excluded to prevent data leakage
        - Clusters represent **resistance phenotypes**, NOT taxonomic groups
        - Post-hoc associations with metadata are correlational only
        
        **Cluster Cut Rule:**
        - Number of clusters determined by combined elbow + silhouette analysis
        - Optimal k is selected dynamically based on statistical criteria
        """,
        
        'supervised_discrimination': """
        ### Supervised Learning (Pattern Discrimination)
        
        **Purpose:** Evaluate how well resistance fingerprints discriminate known categories
        
        This is **pattern discrimination**, NOT prediction. The models assess:
        - How consistently resistance patterns align with species identity
        - How consistently resistance patterns align with MDR status
        
        **Important Distinction:**
        - ✅ Pattern discrimination: Assesses pattern consistency within existing data
        - ❌ Prediction: Forecasting outcomes for new, unseen samples (NOT our goal)
        
        **Model Set (Rationalized):**
        1. **Logistic Regression** - Linear baseline with coefficient interpretation
        2. **Random Forest** - Nonlinear model with Gini feature importance
        3. **k-Nearest Neighbors** - Distance-based consistency check
        
        **Data Splitting:**
        - 80% training / 20% testing (stratified)
        - Split performed BEFORE any preprocessing to prevent leakage
        - Scaling and imputation fit on training data only
        
        **Metrics:**
        - Macro-averaged precision, recall, F1-score
        - Macro averaging treats all classes equally
        
        **MDR Task Transparency:**
        MDR labels are derived from the SAME AST features used as input. This is 
        explicitly a **self-consistency check**, not prediction.
        """,
        
        'pca_usage': """
        ### Principal Component Analysis (PCA)
        
        **Purpose:** Dimensionality reduction for visualization
        
        PCA transforms the high-dimensional resistance data into 2D for visualization, 
        preserving as much variance as possible.
        
        **Procedure:**
        1. Missing value imputation (median strategy)
        2. Standardization (StandardScaler)
        3. PCA transformation to 2 components
        
        **Interpretation:**
        - PC1 and PC2 capture the most variance in resistance patterns
        - Component loadings indicate which antibiotics contribute most
        - Clustering in PCA space reflects resistance profile similarity
        
        **Limitations:**
        - Only 2 components shown (may lose information)
        - Linear transformation only
        - Should not be over-interpreted
        """
    }


def get_limitations_content() -> dict:
    """
    Get limitations documentation.
    
    Returns:
    --------
    dict
        Dictionary with limitations sections
    """
    return {
        'overview': """
        ## Limitations & Claim Boundaries
        
        This section explicitly documents what this analysis does NOT show 
        and the boundaries of claims that can be made from these results.
        """,
        
        'no_temporal_inference': """
        ### ❌ No Temporal Inference
        
        **This analysis does NOT:**
        - Determine if resistance patterns are increasing or decreasing over time
        - Identify trends or temporal evolution
        - Predict future resistance rates
        
        **Why:** This is cross-sectional data from a single time point. Temporal 
        inference requires longitudinal sampling with consistent methodology.
        """,
        
        'no_causal_inference': """
        ### ❌ No Causal Inference
        
        **This analysis does NOT:**
        - Determine what causes resistance patterns
        - Identify mechanisms of resistance acquisition
        - Establish that environmental factors "drive" resistance
        
        **Why:** Observational data can identify associations, not causation. 
        Experimental studies are needed to establish causal relationships.
        
        **Language Discipline:**
        - ✅ Use: "associated with", "enriched in", "observed across"
        - ❌ Avoid: "driven by", "caused by", "transmission pathway"
        """,
        
        'no_predictive_claims': """
        ### ❌ No Predictive Claims
        
        **This analysis does NOT:**
        - Predict resistance for new/future isolates
        - Provide clinical decision support
        - Forecast outbreak risk
        
        **Why:** Model metrics describe pattern consistency within the current 
        dataset only. Predictive validity requires prospective validation studies.
        """,
        
        'no_transmission_inference': """
        ### ❌ No Transmission Inference
        
        **This analysis does NOT:**
        - Identify transmission pathways
        - Determine how resistance genes spread
        - Trace the origin of resistant isolates
        
        **Why:** Transmission inference requires genomic data (e.g., WGS) and 
        epidemiological investigation. AST data alone cannot establish transmission.
        """,
        
        'dataset_dependency': """
        ### ⚠️ Dataset-Dependent Results
        
        **Important Considerations:**
        - Results are specific to this dataset and sampling design
        - Findings may not generalize to other regions or time periods
        - Sample selection may introduce bias (e.g., clinical vs. environmental)
        
        **Generalizability Constraints:**
        - Geographic: Limited to sampled regions
        - Temporal: Cross-sectional snapshot only
        - Species: Limited to identified bacterial species
        - Antibiotics: Limited to tested panel
        """,
        
        'data_quality': """
        ### ⚠️ Data Quality Considerations
        
        **Potential Issues:**
        - Missing AST data may affect pattern detection
        - Sample size variation across sites affects statistical power
        - Laboratory methodology variation may introduce noise
        
        **Mitigations Applied:**
        - Threshold-based exclusion of low-coverage antibiotics
        - Median imputation for remaining missing values
        - Chi-square tests with appropriate significance levels
        """
    }


def get_disclaimers() -> dict:
    """
    Get all disclaimer texts.
    
    Returns:
    --------
    dict
        Dictionary with different disclaimer variants
    """
    return {
        'main_disclaimer': """
        > **⚠️ Important Disclaimer**
        >
        > This tool is intended **exclusively** for exploratory antimicrobial resistance 
        > pattern recognition and surveillance analysis. It does **NOT** provide:
        > - Clinical decision support
        > - Predictive assessments
        > - Treatment recommendations
        > - Risk scores
        >
        > Results should be interpreted as exploratory findings only.
        """,
        
        'landing_page': """
        <div style="background-color: #fff3cd; border: 2px solid #ffc107; border-radius: 10px; padding: 20px; margin: 20px 0;">
            <h4 style="color: #856404; margin-top: 0;">⚠️ Important Disclaimer</h4>
            <p style="color: #856404; margin-bottom: 0;">
                This tool is intended <strong>exclusively</strong> for exploratory antimicrobial 
                resistance pattern recognition and surveillance analysis. It does <strong>NOT</strong> 
                provide clinical decision support or predictive assessments.
            </p>
        </div>
        """,
        
        'model_results': """
        <div style="background-color: #d1ecf1; border: 1px solid #bee5eb; border-radius: 5px; padding: 10px; margin: 10px 0;">
            <strong>📊 Interpretation Note:</strong> Metrics shown are <em>pattern consistency measures</em>, 
            not predictive performance. They quantify how resistance patterns align with known categories 
            within this dataset only.
        </div>
        """,
        
        'footer': """
        ---
        <div style="text-align: center; color: #6c757d; font-size: 0.85rem; padding: 20px 0;">
            <strong>AMR Pattern Recognition & Exploratory Analysis Dashboard</strong><br>
            For research and surveillance purposes only | Not for clinical use<br>
            <em>This tool does not provide clinical decision support or predictive assessments.</em>
        </div>
        """
    }


def get_about_content() -> str:
    """
    Get the about/introduction content.
    
    Returns:
    --------
    str
        About section content
    """
    return """
    ## About This Tool
    
    **AMR Pattern Recognition & Exploratory Analysis Dashboard**
    
    This interactive dashboard provides tools for exploring antimicrobial resistance 
    (AMR) patterns in bacterial isolate data. It is designed for:
    
    - 🔬 **Surveillance Analysis**: Visualize resistance patterns across regions and environments
    - 📊 **Pattern Recognition**: Identify natural groupings in resistance profiles
    - 📈 **Exploratory Analysis**: Examine associations between resistance and metadata
    
    ### What This Tool Does
    
    ✅ Visualizes resistance patterns through heatmaps and PCA plots  
    ✅ Displays cluster analysis results from hierarchical clustering  
    ✅ Shows model consistency metrics (how patterns align with known categories)  
    ✅ Provides regional and environmental distribution analysis  
    
    ### What This Tool Does NOT Do
    
    ❌ Predict resistance for new samples  
    ❌ Provide clinical decision support  
    ❌ Calculate risk scores  
    ❌ Make treatment recommendations  
    
    ### Data Privacy
    
    - Uploaded data is processed **in memory only**
    - No data is stored on disk or transmitted externally
    - No raw inputs are logged
    
    ---
    
    *This tool was developed as part of a thesis project on AMR pattern recognition.*
    """


def get_glossary() -> dict:
    """
    Get glossary of terms.
    
    Returns:
    --------
    dict
        Dictionary of term definitions
    """
    return {
        'MDR (Multi-Drug Resistant)': 'Isolate showing resistance to ≥3 antibiotic classes (Magiorakos et al., 2012)',
        'MAR Index': 'Multiple Antibiotic Resistance Index: ratio of resistant antibiotics to total tested',
        'AST': 'Antimicrobial Susceptibility Testing - laboratory method to determine resistance',
        'S (Susceptible)': 'Isolate is killed or inhibited by typical antibiotic concentrations',
        'I (Intermediate)': 'Isolate shows borderline susceptibility; clinical efficacy uncertain',
        'R (Resistant)': 'Isolate is not killed or inhibited by typical antibiotic concentrations',
        'Resistance Fingerprint': 'Pattern of resistance/susceptibility results across all tested antibiotics',
        'Pattern Discrimination': 'Evaluating how well resistance patterns distinguish between known categories',
        'Cluster': 'Group of isolates with similar resistance patterns identified by hierarchical clustering',
        'PCA': 'Principal Component Analysis - dimensionality reduction technique',
        'Feature Importance': 'Measure of how much an antibiotic contributes to group discrimination (associative, not causal)'
    }
