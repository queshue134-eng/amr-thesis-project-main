"""
Supervised Learning Module for AMR Thesis Project
Phase 4 - Supervised Learning for Pattern Discrimination

CRITICAL IMPLEMENTATION NOTES (Phase 3 Improvements):
=====================================================
1. LEAKAGE-SAFE PIPELINE:
   - Input matrix contains ONLY resistance fingerprints (no metadata like region, site, environment)
   - Train-test split performed BEFORE scaling, imputation, and feature selection
   - Fixed random_state for reproducibility

2. MDR TARGET TRANSPARENCY:
   - MDR label is derived from the SAME AST features used as input
   - MDR task is treated as SELF-CONSISTENCY DISCRIMINATION
   - This neutralizes the "predicted MDR from MDR" objection by making it explicit

3. TASK SEPARATION:
   - Task A: Species discrimination (multi-class) - resistance fingerprints → species
   - Task B: MDR discrimination (binary) - resistance fingerprints → MDR flag
   - Each task has separate pipelines, confusion matrices, and metric tables

4. RATIONALIZED MODEL SET (3-4 models max):
   - Logistic Regression → Linear baseline (coefficient magnitude for importance)
   - Random Forest → Nonlinear + feature importance (Gini importance)
   - k-Nearest Neighbors → Distance-based consistency check
   - (Optional) SVM → Margin-based separation

5. EVALUATION METRICS DISCIPLINE:
   - Always report: Accuracy, Macro-averaged Precision/Recall/F1
   - Macro averaging prevents class imbalance bias
   - Confusion matrix per task

6. INTERPRETATION LANGUAGE:
   - Use: "Shows consistent alignment", "Demonstrates discriminative capacity"
   - Avoid: "Performs well", "Predicts accurately"
   - Feature importance is ASSOCIATIVE, not CAUSAL

OBJECTIVE (4.1):
    Evaluate how well resistance fingerprints discriminate known categories:
    - Task A: Bacterial species (multi-class classification)
    - Task B: MDR vs non-MDR groups (binary classification)
    Note: This is pattern discrimination, NOT forecasting/prediction of future outcomes.

DATA SPLITTING (4.2):
    - 80%-20% train-test split BEFORE any preprocessing
    - Scaling applied separately to train/test after split (prevents leakage)
    - Purpose: Assess model generalization, avoid overfitting, support robustness
    - Framed as model validation, not prediction

MODEL SELECTION (4.3) - RATIONALIZED SET:
    - Logistic Regression → Linear baseline
    - Random Forest → Tree-based, feature importance
    - k-Nearest Neighbors → Distance-based
    Model groups: Linear, Tree-based, Distance-based

MODEL TRAINING (4.4):
    - Inputs: Resistance fingerprints ONLY (encoded antibiotic susceptibility)
    - Targets: Known labels (species or MDR category)
    - NO metadata enters the model (region, site, environment excluded)

MODEL EVALUATION (4.5):
    - Accuracy, Macro Precision, Macro Recall, Macro F1-score, Confusion matrix
    - Interpretation: Metrics quantify how consistently resistance patterns
      align with known categories, NOT predictive performance for future samples.

MODEL INTERPRETATION (4.6):
    - Feature importance analysis identifies antibiotics contributing most
      to group separation (ASSOCIATIVE importance, not causal)
    - Random Forest: Gini importance
    - Logistic Regression: Absolute coefficient magnitude
"""

import pandas as pd
import numpy as np
from typing import List, Tuple, Dict, Optional
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, f1_score,
    confusion_matrix, classification_report
)
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
import warnings

# Import from centralized configuration
try:
    from config import ANTIBIOTIC_CLASSES, RANDOM_STATE, SUPERVISED_CONFIG
except ImportError:
    # Fallback for standalone execution
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).parent.parent))
    from config import ANTIBIOTIC_CLASSES, RANDOM_STATE, SUPERVISED_CONFIG

warnings.filterwarnings('ignore')


# ============================================================================
# RATIONALIZED MODEL SET (Phase 3 Requirement 3.1)
# ============================================================================
# Reduced to 3 models with clear methodological justification:
# - Linear: Logistic Regression (coefficient interpretation)
# - Tree-based: Random Forest (Gini importance, nonlinear)
# - Distance-based: k-NN (local pattern consistency)
# ============================================================================

# Model configurations with conceptual grouping
# Uses RANDOM_STATE from centralized config for reproducibility
MODELS = {
    'Logistic Regression': {
        'instance': LogisticRegression(random_state=RANDOM_STATE, max_iter=1000, solver='lbfgs'),
        'category': 'Linear',
        'description': 'Linear baseline with coefficient interpretation'
    },
    'Random Forest': {
        'instance': RandomForestClassifier(n_estimators=100, random_state=RANDOM_STATE, n_jobs=-1),
        'category': 'Tree-based',
        'description': 'Nonlinear model with Gini feature importance'
    },
    'k-Nearest Neighbors': {
        'instance': KNeighborsClassifier(n_neighbors=5),
        'category': 'Distance-based',
        'description': 'Distance-based consistency check'
    }
}

# Model category descriptions for methodology section
MODEL_CATEGORIES = {
    'Linear': 'Models that learn linear decision boundaries. Interpretable through coefficient magnitudes.',
    'Tree-based': 'Models that partition feature space using decision rules. Provide Gini-based importance.',
    'Distance-based': 'Models that classify based on similarity to training instances. No parameters to interpret.'
}

# Antibiotic class mapping imported from centralized config
# See src/config.py for the complete ANTIBIOTIC_CLASSES dictionary


def prepare_data_for_classification(df: pd.DataFrame,
                                    feature_cols: List[str],
                                    target_col: str,
                                    test_size: float = 0.2,
                                    random_state: int = 42) -> Tuple:
    """
    Prepare data for supervised learning with LEAKAGE-SAFE preprocessing (Phase 4.2 & 4.4).
    
    CRITICAL: Train-test split is performed BEFORE scaling and imputation to prevent
    data leakage. This ensures no information from test set influences preprocessing.
    
    LEAKAGE-SAFE PIPELINE ORDER:
    1. Validate target column and filter invalid samples
    2. Extract ONLY resistance fingerprints (no metadata)
    3. Perform train-test split (80/20)
    4. Fit imputer on TRAIN only, transform both
    5. Fit scaler on TRAIN only, transform both
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe with features and target
    feature_cols : list
        List of feature column names (resistance fingerprints ONLY - no metadata)
    target_col : str
        Name of target column (known labels: species or MDR category)
    test_size : float
        Proportion of data for testing (default 0.2 = 20%)
    random_state : int
        Random seed for reproducibility (MUST be fixed and reported)
    
    Returns:
    --------
    tuple
        (X_train, X_test, y_train, y_test, label_encoder, scaler, imputer, preprocessing_info)
    """
    # Remove rows with missing target
    df_valid = df[df[target_col].notna()].copy()
    
    # Filter out classes with fewer than 2 samples (required for stratified split)
    class_counts = df_valid[target_col].value_counts()
    valid_classes = class_counts[class_counts >= 2].index
    
    if (class_counts < 2).any():
        removed_classes = class_counts[class_counts < 2].index.tolist()
        warnings.warn(
            f"Removed {len(removed_classes)} class(es) with fewer than 2 samples "
            f"(required for stratified split): {removed_classes}",
            category=UserWarning
        )
    
    df_valid = df_valid[df_valid[target_col].isin(valid_classes)].copy()
    
    if len(df_valid) == 0:
        raise ValueError("No samples remaining after filtering classes with fewer than 2 samples.")
    
    # ===========================================================================
    # STRICT FEATURE-LABEL SEPARATION (Phase 3 Requirement 1.1)
    # Only resistance fingerprints enter the model - NO metadata
    # ===========================================================================
    existing_cols = [c for c in feature_cols if c in df_valid.columns]
    
    # Verify all columns are resistance fingerprints (must end with _encoded)
    non_fingerprint_cols = [c for c in existing_cols if not c.endswith('_encoded')]
    if non_fingerprint_cols:
        warnings.warn(
            f"Non-resistance columns detected and will be excluded: {non_fingerprint_cols}. "
            "Only '_encoded' columns are allowed as features.",
            category=UserWarning
        )
        existing_cols = [c for c in existing_cols if c.endswith('_encoded')]
    
    X = df_valid[existing_cols].values  # Convert to numpy immediately
    y = df_valid[target_col].copy()
    
    # Encode target if categorical
    label_encoder = None
    if y.dtype == 'object' or y.dtype.name == 'category':
        label_encoder = LabelEncoder()
        y = label_encoder.fit_transform(y)
    else:
        y = y.values
    
    # ===========================================================================
    # TRAIN-TEST SPLIT DISCIPLINE (Phase 3 Requirement 1.2)
    # Split BEFORE scaling, feature selection, and imputation
    # ===========================================================================
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=test_size, random_state=random_state, stratify=y
    )
    
    # ===========================================================================
    # FIT PREPROCESSING ON TRAIN ONLY (prevents leakage)
    # ===========================================================================
    # Imputation: fit on train, transform both
    imputer = SimpleImputer(strategy='median')
    X_train_imputed = imputer.fit_transform(X_train)
    X_test_imputed = imputer.transform(X_test)
    
    # Scaling: fit on train, transform both
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train_imputed)
    X_test_scaled = scaler.transform(X_test_imputed)
    
    # Preprocessing info for transparency
    preprocessing_info = {
        'random_state': random_state,
        'test_size': test_size,
        'n_features': len(existing_cols),
        'feature_columns': existing_cols,
        'n_train_samples': len(X_train),
        'n_test_samples': len(X_test),
        'preprocessing_order': [
            '1. Train-test split (stratified)',
            '2. Median imputation (fit on train)',
            '3. Standard scaling (fit on train)'
        ],
        'leakage_safe': True
    }
    
    return X_train_scaled, X_test_scaled, y_train, y_test, label_encoder, scaler, imputer, preprocessing_info


def train_model(model, X_train: np.ndarray, y_train: np.ndarray):
    """
    Train a single model.
    
    Parameters:
    -----------
    model : sklearn estimator
        Model to train
    X_train : np.ndarray
        Training features
    y_train : np.ndarray
        Training labels
    
    Returns:
    --------
    sklearn estimator
        Trained model
    """
    model.fit(X_train, y_train)
    return model


def evaluate_model(model, X_test: np.ndarray, y_test: np.ndarray,
                   label_encoder=None) -> Dict:
    """
    Evaluate model performance with MACRO-AVERAGED metrics (Phase 4.5).
    
    METRICS DISCIPLINE (Phase 3 Requirement 4.1):
    - Always report: Accuracy, Macro-averaged Precision, Recall, F1
    - Macro averaging treats all classes equally, preventing class imbalance bias
    - Include confusion matrix for detailed analysis
    
    INTERPRETATION LANGUAGE (Phase 3 Requirement 4.2):
    - These metrics quantify how consistently resistance patterns align 
      with known categories (NOT predictive performance)
    - Use: "Shows consistent alignment", "Demonstrates discriminative capacity"
    - Avoid: "Performs well", "Predicts accurately"
    
    Parameters:
    -----------
    model : sklearn estimator
        Trained model
    X_test : np.ndarray
        Test features
    y_test : np.ndarray
        Test labels
    label_encoder : LabelEncoder, optional
        Label encoder for decoding classes
    
    Returns:
    --------
    dict
        Evaluation metrics including macro-averaged scores and confusion matrix
    """
    # Predictions
    y_pred = model.predict(X_test)
    
    # Get unique labels present in either y_test or y_pred
    present_labels = np.unique(np.concatenate([y_test, y_pred]))
    
    # ===========================================================================
    # MACRO-AVERAGED METRICS (Phase 3 Requirement 4.1)
    # Macro averaging prevents class imbalance bias
    # ===========================================================================
    results = {
        'accuracy': accuracy_score(y_test, y_pred),
        # MACRO-averaged metrics (treats all classes equally)
        'precision_macro': precision_score(y_test, y_pred, average='macro', zero_division=0),
        'recall_macro': recall_score(y_test, y_pred, average='macro', zero_division=0),
        'f1_score_macro': f1_score(y_test, y_pred, average='macro', zero_division=0),
        # Keep weighted for reference but primary is macro
        'precision_weighted': precision_score(y_test, y_pred, average='weighted', zero_division=0),
        'recall_weighted': recall_score(y_test, y_pred, average='weighted', zero_division=0),
        'f1_score_weighted': f1_score(y_test, y_pred, average='weighted', zero_division=0),
        # Confusion matrix
        'confusion_matrix': confusion_matrix(y_test, y_pred, labels=present_labels).tolist(),
        'class_labels': present_labels.tolist()
    }
    
    # Classification report (provides per-class breakdown)
    if label_encoder is not None:
        target_names = label_encoder.inverse_transform(present_labels).tolist()
        results['class_names'] = target_names
    else:
        target_names = [str(c) for c in present_labels]
        results['class_names'] = target_names
    
    results['classification_report'] = classification_report(
        y_test, y_pred, labels=present_labels, target_names=target_names, output_dict=True
    )
    
    return results


def bootstrap_confidence_interval(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    metric_func,
    n_iterations: int = 1000,
    confidence_level: float = 0.95,
    random_state: int = None,
    **metric_kwargs
) -> Dict:
    """
    Compute bootstrap confidence interval for a classification metric.
    
    STATISTICAL RIGOR ENHANCEMENT:
    This function addresses the gap of reporting point estimates without
    confidence intervals. Bootstrap resampling provides non-parametric
    CI estimation that works for any metric.
    
    Parameters
    ----------
    y_true : np.ndarray
        True labels
    y_pred : np.ndarray
        Predicted labels
    metric_func : callable
        Metric function (e.g., accuracy_score, f1_score)
    n_iterations : int
        Number of bootstrap iterations (default 1000)
    confidence_level : float
        Confidence level for interval (default 0.95 = 95% CI)
    random_state : int, optional
        Random seed for reproducibility
    **metric_kwargs : dict
        Additional arguments passed to metric_func
    
    Returns
    -------
    dict
        Dictionary with point_estimate, ci_lower, ci_upper, std_error
    
    Reference
    ---------
    Efron, B., & Tibshirani, R. J. (1993). An Introduction to the Bootstrap.
    Chapman & Hall/CRC.
    """
    rng = np.random.RandomState(random_state)
    n_samples = len(y_true)
    
    # Collect bootstrap samples
    bootstrap_scores = []
    
    for _ in range(n_iterations):
        # Resample with replacement
        indices = rng.choice(n_samples, size=n_samples, replace=True)
        y_true_boot = y_true[indices]
        y_pred_boot = y_pred[indices]
        
        # Skip if only one class present (can't compute some metrics)
        if len(np.unique(y_true_boot)) < 2:
            continue
        
        try:
            score = metric_func(y_true_boot, y_pred_boot, **metric_kwargs)
            bootstrap_scores.append(score)
        except Exception:
            # Skip failed iterations
            continue
    
    if len(bootstrap_scores) < 100:
        # Not enough valid bootstrap samples
        return {
            'point_estimate': metric_func(y_true, y_pred, **metric_kwargs),
            'ci_lower': np.nan,
            'ci_upper': np.nan,
            'std_error': np.nan,
            'n_valid_iterations': len(bootstrap_scores),
            'warning': 'Insufficient valid bootstrap samples for reliable CI'
        }
    
    # Calculate percentile-based confidence interval
    alpha = 1 - confidence_level
    ci_lower = np.percentile(bootstrap_scores, (alpha / 2) * 100)
    ci_upper = np.percentile(bootstrap_scores, (1 - alpha / 2) * 100)
    
    return {
        'point_estimate': metric_func(y_true, y_pred, **metric_kwargs),
        'ci_lower': ci_lower,
        'ci_upper': ci_upper,
        'std_error': np.std(bootstrap_scores),
        'n_valid_iterations': len(bootstrap_scores),
        'confidence_level': confidence_level
    }


def evaluate_model_with_ci(
    model,
    X_test: np.ndarray,
    y_test: np.ndarray,
    label_encoder=None,
    n_bootstrap: int = 1000,
    confidence_level: float = 0.95,
    random_state: int = None
) -> Dict:
    """
    Evaluate model performance with bootstrap confidence intervals.
    
    This function extends evaluate_model() by adding 95% CI for all metrics,
    addressing the statistical rigor requirement for thesis defense.
    
    Parameters
    ----------
    model : sklearn estimator
        Trained model
    X_test : np.ndarray
        Test features
    y_test : np.ndarray
        Test labels
    label_encoder : LabelEncoder, optional
        Label encoder for decoding classes
    n_bootstrap : int
        Number of bootstrap iterations (default 1000)
    confidence_level : float
        Confidence level (default 0.95)
    random_state : int, optional
        Random seed for reproducibility
    
    Returns
    -------
    dict
        Evaluation metrics with confidence intervals
    """
    # Get base evaluation
    results = evaluate_model(model, X_test, y_test, label_encoder)
    
    # Get predictions for CI calculation
    y_pred = model.predict(X_test)
    
    # Metrics to compute CI for
    metrics_with_ci = {
        'accuracy': (accuracy_score, {}),
        'precision_macro': (precision_score, {'average': 'macro', 'zero_division': 0}),
        'recall_macro': (recall_score, {'average': 'macro', 'zero_division': 0}),
        'f1_score_macro': (f1_score, {'average': 'macro', 'zero_division': 0}),
    }
    
    results['confidence_intervals'] = {}
    
    print(f"    Computing {confidence_level*100:.0f}% bootstrap confidence intervals...")
    
    for metric_name, (metric_func, kwargs) in metrics_with_ci.items():
        ci_result = bootstrap_confidence_interval(
            y_test, y_pred, metric_func,
            n_iterations=n_bootstrap,
            confidence_level=confidence_level,
            random_state=random_state,
            **kwargs
        )
        results['confidence_intervals'][metric_name] = ci_result
        
        # Add formatted CI string for easy reporting
        if not np.isnan(ci_result['ci_lower']):
            results[f'{metric_name}_ci'] = (
                f"{ci_result['point_estimate']:.4f} "
                f"({ci_result['ci_lower']:.4f}-{ci_result['ci_upper']:.4f})"
            )
        else:
            results[f'{metric_name}_ci'] = f"{ci_result['point_estimate']:.4f} (CI unavailable)"
    
    results['bootstrap_config'] = {
        'n_iterations': n_bootstrap,
        'confidence_level': confidence_level,
        'random_state': random_state
    }
    
    return results


def get_feature_importance(model, feature_names: List[str], model_name: str = None) -> Dict:
    """
    Extract ASSOCIATIVE feature importance from model (Phase 4.6).
    
    FEATURE IMPORTANCE EXTRACTION (Phase 3 Requirement 5.1):
    - Random Forest: Gini importance (feature_importances_)
    - Logistic Regression: Absolute coefficient magnitude
    - k-NN: No native importance (returns empty dict)
    
    BIOLOGICAL RESTRAINT (Phase 3 Requirement 5.3):
    - Importance scores indicate ASSOCIATIVE patterns, NOT causal mechanisms
    - High importance means the antibiotic helps discriminate groups
    - Does NOT imply the antibiotic causes the group membership
    
    Parameters:
    -----------
    model : sklearn estimator
        Trained model
    feature_names : list
        List of feature names (antibiotic names)
    model_name : str, optional
        Name of the model for importance method annotation
    
    Returns:
    --------
    dict
        Feature importance scores with method annotation
    """
    importance_data = {
        'scores': {},
        'method': None,
        'model_type': model_name,
        'interpretation_note': 'ASSOCIATIVE importance - does not imply causation'
    }
    
    if hasattr(model, 'feature_importances_'):
        # Tree-based models: Gini importance
        importance_data['method'] = 'Gini importance (mean decrease in impurity)'
        for name, imp in zip(feature_names, model.feature_importances_):
            importance_data['scores'][name] = float(imp)
    elif hasattr(model, 'coef_'):
        # Linear models: Absolute coefficient magnitude
        importance_data['method'] = 'Absolute coefficient magnitude'
        coef = np.abs(model.coef_)
        if len(coef.shape) > 1:
            coef = coef.mean(axis=0)
        for name, imp in zip(feature_names, coef):
            importance_data['scores'][name] = float(imp)
    else:
        # Models without native importance (e.g., k-NN)
        importance_data['method'] = 'No native feature importance available'
        importance_data['scores'] = {}
    
    # Sort by importance
    importance_data['scores'] = dict(
        sorted(importance_data['scores'].items(), key=lambda x: x[1], reverse=True)
    )
    
    return importance_data


def _get_resistance_interpretation(resistance_rate: float) -> str:
    """
    Get human-readable interpretation of resistance rate.
    
    Parameters:
    -----------
    resistance_rate : float
        Proportion of isolates showing resistance (0.0-1.0)
    
    Returns:
    --------
    str
        Interpretation category (High/Moderate/Low resistance)
    """
    if resistance_rate > 0.5:
        return 'High resistance'
    elif resistance_rate > 0.2:
        return 'Moderate resistance'
    else:
        return 'Low resistance'


def interpret_feature_importance(feature_importance: Dict,
                                  df: pd.DataFrame = None,
                                  target_col: str = None,
                                  task_type: str = None) -> Dict:
    """
    Interpret feature importance findings with BIOLOGICAL RESTRAINT (Phase 4.6).
    
    INTERPRETATION LANGUAGE (Phase 3 Requirement 5.3):
    - Use: "Associative importance", "Shows consistent alignment"
    - Avoid: Causal or mechanistic claims
    
    Parameters:
    -----------
    feature_importance : dict
        Feature importance data from get_feature_importance()
        Expected format: {'scores': {...}, 'method': str, ...}
    df : pd.DataFrame, optional
        Original dataframe for computing additional statistics
    target_col : str, optional
        Target column name for context
    task_type : str, optional
        'species' or 'mdr' for task-specific interpretation
    
    Returns:
    --------
    dict
        Interpretation results with biological restraint language
    """
    interpretation = {
        'top_discriminators': [],
        'antibiotic_classes_involved': set(),
        'interpretation_notes': [],
        'biological_restraint_note': (
            "IMPORTANT: Feature importance indicates ASSOCIATIVE patterns only. "
            "High importance means the antibiotic helps discriminate groups, "
            "but does NOT imply causal or mechanistic relationships."
        )
    }
    
    # Handle both old format (dict of scores) and new format (nested dict)
    if isinstance(feature_importance, dict) and 'scores' in feature_importance:
        scores = feature_importance['scores']
        interpretation['importance_method'] = feature_importance.get('method', 'Unknown')
    else:
        scores = feature_importance
        interpretation['importance_method'] = 'Unknown'
    
    # Get top 5 discriminators
    top_features = list(scores.items())[:5]
    
    for antibiotic, score in top_features:
        ab_class = ANTIBIOTIC_CLASSES.get(antibiotic, 'Unknown class')
        interpretation['top_discriminators'].append({
            'antibiotic': antibiotic,
            'importance_score': score,
            'antibiotic_class': ab_class,
            'interpretation': 'Shows associative importance for group discrimination'
        })
        interpretation['antibiotic_classes_involved'].add(ab_class)
    
    # Convert set to list for JSON serialization
    interpretation['antibiotic_classes_involved'] = list(
        interpretation['antibiotic_classes_involved']
    )
    
    # ===========================================================================
    # DISCIPLINED INTERPRETATION LANGUAGE (Phase 3 Requirement 4.2)
    # ===========================================================================
    interpretation['interpretation_notes'].append(
        "Feature importance scores indicate antibiotics showing consistent alignment "
        "with group membership (ASSOCIATIVE, not causal)."
    )
    interpretation['interpretation_notes'].append(
        "Higher scores suggest these antibiotics demonstrate stronger discriminative "
        "capacity within the analyzed categories."
    )
    
    # Task-specific interpretation
    if task_type == 'mdr' or (target_col and 'MDR' in target_col.upper()):
        interpretation['interpretation_notes'].append(
            "MDR TASK NOTE: For MDR discrimination, top features indicate antibiotics "
            "whose resistance patterns show consistent alignment with MDR status. "
            "This is a self-consistency check since MDR is derived from the same AST features."
        )
    elif task_type == 'species' or (target_col and target_col.upper() in ['ISOLATE_ID', 'SPECIES']):
        interpretation['interpretation_notes'].append(
            "SPECIES TASK NOTE: For species discrimination, top features indicate antibiotics "
            "whose resistance patterns show consistent alignment with species identity."
        )
    
    # Compute resistance rates if dataframe is provided
    if df is not None:
        resistance_stats = {}
        for ab, _ in top_features:
            col_name = f"{ab}_encoded" if f"{ab}_encoded" in df.columns else ab
            if col_name in df.columns:
                resistance_rate = (df[col_name] == 2).mean()
                resistance_stats[ab] = {
                    'resistance_rate': float(resistance_rate),
                    'interpretation': _get_resistance_interpretation(resistance_rate)
                }
        interpretation['ast_results'] = resistance_stats
    
    return interpretation


def run_all_models(X_train: np.ndarray, X_test: np.ndarray,
                   y_train: np.ndarray, y_test: np.ndarray,
                   feature_names: List[str],
                   label_encoder=None) -> Dict:
    """
    Train and evaluate RATIONALIZED MODEL SET (Phase 3 Requirement 3.1).
    
    Models are grouped conceptually:
    - Linear: Logistic Regression
    - Tree-based: Random Forest  
    - Distance-based: k-NN
    
    Parameters:
    -----------
    X_train, X_test : np.ndarray
        Training and test features
    y_train, y_test : np.ndarray
        Training and test labels
    feature_names : list
        List of feature names
    label_encoder : LabelEncoder, optional
        Label encoder
    
    Returns:
    --------
    dict
        Results for all models with conceptual grouping
    """
    results = {}
    
    for name, model_config in MODELS.items():
        model_template = model_config['instance']
        category = model_config['category']
        description = model_config['description']
        
        print(f"  Training {name} ({category})...")
        
        # Create fresh model instance
        model_instance = model_template.__class__(**model_template.get_params())
        
        # Train
        trained_model = train_model(model_instance, X_train, y_train)
        
        # Evaluate
        metrics = evaluate_model(trained_model, X_test, y_test, label_encoder)
        metrics['model_name'] = name
        metrics['model_category'] = category
        metrics['model_description'] = description
        
        # Feature importance with method annotation
        metrics['feature_importance'] = get_feature_importance(trained_model, feature_names, name)
        
        results[name] = {
            'model': trained_model,
            'metrics': metrics,
            'category': category,
            'description': description
        }
    
    return results


def compare_models(results: Dict) -> pd.DataFrame:
    """
    Compare performance of all models using MACRO-AVERAGED metrics.
    
    Phase 3 Requirement 4.1: Use macro-averaged metrics to prevent class imbalance bias.
    
    Parameters:
    -----------
    results : dict
        Results from run_all_models
    
    Returns:
    --------
    pd.DataFrame
        Comparison table with macro-averaged metrics and model categories
    """
    comparison = []
    
    for name, result in results.items():
        metrics = result['metrics']
        comparison.append({
            'Model': name,
            'Category': result.get('category', 'Unknown'),
            'Accuracy': metrics['accuracy'],
            'Precision (Macro)': metrics['precision_macro'],
            'Recall (Macro)': metrics['recall_macro'],
            'F1-Score (Macro)': metrics['f1_score_macro']
        })
    
    df_comparison = pd.DataFrame(comparison)
    df_comparison = df_comparison.sort_values('F1-Score (Macro)', ascending=False)
    
    return df_comparison


def compare_feature_importance_across_models(results: Dict, top_n: int = 10) -> Dict:
    """
    Compare top-ranked antibiotics across models (Phase 3 Requirement 6.1).
    
    MODEL AGREEMENT CHECK: Identifies overlapping important features across
    different model types, which strengthens confidence in the patterns.
    
    Parameters:
    -----------
    results : dict
        Results from run_all_models
    top_n : int
        Number of top features to compare
    
    Returns:
    --------
    dict
        Model agreement analysis with overlapping features
    """
    agreement = {
        'model_rankings': {},
        'overlapping_features': [],
        'agreement_summary': ''
    }
    
    all_top_features = []
    
    for name, result in results.items():
        importance = result['metrics'].get('feature_importance', {})
        if isinstance(importance, dict) and 'scores' in importance:
            scores = importance['scores']
        else:
            scores = importance if importance else {}
        
        top_features = list(scores.keys())[:top_n]
        agreement['model_rankings'][name] = top_features
        all_top_features.extend(top_features)
    
    # Find features that appear in multiple models' top rankings
    from collections import Counter
    feature_counts = Counter(all_top_features)
    n_models = len(results)
    
    # Features appearing in at least 2 models (if 3 models) or all models
    min_agreement = max(2, n_models - 1)
    overlapping = [f for f, count in feature_counts.items() if count >= min_agreement]
    
    agreement['overlapping_features'] = overlapping
    agreement['agreement_summary'] = (
        f"{len(overlapping)} antibiotic(s) appear in top-{top_n} rankings across "
        f"at least {min_agreement} of {n_models} models: {', '.join(overlapping) if overlapping else 'None'}"
    )
    
    return agreement


def check_stability_across_seeds(df: pd.DataFrame,
                                  feature_cols: List[str],
                                  target_col: str,
                                  seeds: List[int] = None,
                                  n_seeds: int = 3) -> Dict:
    """
    Check stability of results across different random seeds (Phase 3 Requirement 6.2).
    
    Re-runs the split with different seeds and reports consistency qualitatively.
    Shows robustness without heavy computation.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe
    feature_cols : list
        List of feature column names
    target_col : str
        Target column name
    seeds : list, optional
        List of random seeds to test
    n_seeds : int
        Number of seeds to test if seeds not provided
    
    Returns:
    --------
    dict
        Stability analysis results
    """
    if seeds is None:
        seeds = [42, 123, 456][:n_seeds]
    
    stability = {
        'seeds_tested': seeds,
        'metrics_per_seed': [],
        'feature_rankings_per_seed': {},
        'stability_summary': ''
    }
    
    accuracies = []
    f1_scores = []
    
    for seed in seeds:
        # Prepare data with this seed
        result = prepare_data_for_classification(
            df, feature_cols, target_col, random_state=seed
        )
        X_train, X_test, y_train, y_test, label_encoder, scaler, imputer, _ = result
        
        # Train Random Forest (representative model)
        model_config = MODELS['Random Forest']
        model = model_config['instance'].__class__(**model_config['instance'].get_params())
        model.fit(X_train, y_train)
        
        # Evaluate
        y_pred = model.predict(X_test)
        acc = accuracy_score(y_test, y_pred)
        f1 = f1_score(y_test, y_pred, average='macro', zero_division=0)
        
        accuracies.append(acc)
        f1_scores.append(f1)
        
        stability['metrics_per_seed'].append({
            'seed': seed,
            'accuracy': acc,
            'f1_macro': f1
        })
        
        # Get feature importance ranking - use feature_cols directly (no need to check against df)
        feature_names = [c.replace('_encoded', '') for c in feature_cols]
        importance = get_feature_importance(model, feature_names, 'Random Forest')
        top_features = list(importance['scores'].keys())[:5]
        stability['feature_rankings_per_seed'][seed] = top_features
    
    # Calculate stability metrics
    acc_std = np.std(accuracies)
    f1_std = np.std(f1_scores)
    
    # Check feature ranking consistency
    all_top5 = [set(feats) for feats in stability['feature_rankings_per_seed'].values()]
    common_features = set.intersection(*all_top5) if all_top5 else set()
    
    if acc_std < 0.05 and f1_std < 0.05:
        stability_level = 'High stability'
    elif acc_std < 0.10 and f1_std < 0.10:
        stability_level = 'Moderate stability'
    else:
        stability_level = 'Variable results'
    
    stability['stability_summary'] = (
        f"{stability_level}: Accuracy std={acc_std:.4f}, F1 std={f1_std:.4f}. "
        f"{len(common_features)} feature(s) consistently in top-5 across all seeds."
    )
    
    return stability


def create_antibiotic_importance_table(results: Dict, df: pd.DataFrame, task_type: str) -> pd.DataFrame:
    """
    Create antibiotic-level summary table (Phase 3 Requirement 5.2).
    
    Combines importance scores from models with biological context.
    
    Parameters:
    -----------
    results : dict
        Model results from run_all_models
    df : pd.DataFrame
        Original dataframe for resistance rates
    task_type : str
        'species' or 'mdr' for interpretation context
    
    Returns:
    --------
    pd.DataFrame
        Antibiotic importance summary table
    """
    importance_data = []
    
    # Get importance from models that have native feature importance
    # Dynamically check which models have importance data
    models_with_importance = []
    for model_name in results:
        importance = results[model_name]['metrics'].get('feature_importance', {})
        if isinstance(importance, dict) and 'scores' in importance:
            if importance['scores']:  # Non-empty scores
                models_with_importance.append(model_name)
        elif importance:
            models_with_importance.append(model_name)
    
    for model_name in models_with_importance:
        importance = results[model_name]['metrics'].get('feature_importance', {})
        if isinstance(importance, dict) and 'scores' in importance:
            scores = importance['scores']
            method = importance.get('method', 'Unknown')
        else:
            scores = importance if importance else {}
            method = 'Unknown'
        
        if not scores:
            continue
        
        for antibiotic, score in scores.items():
            ab_class = ANTIBIOTIC_CLASSES.get(antibiotic, 'Unknown')
            
            # Get resistance rate
            col_name = f"{antibiotic}_encoded"
            if col_name in df.columns:
                resistance_rate = (df[col_name] == 2).mean()
            else:
                resistance_rate = None
            
            # Determine interpretation based on score rank
            all_scores = list(scores.values())
            rank = all_scores.index(score) + 1 if score in all_scores else None
            
            if rank and rank <= 3:
                importance_level = 'High'
            elif rank and rank <= 7:
                importance_level = 'Moderate'
            else:
                importance_level = 'Low'
            
            # Task-specific interpretation
            if task_type == 'mdr':
                if importance_level == 'High':
                    interpretation = 'Common MDR marker (high discriminative capacity)'
                elif importance_level == 'Moderate':
                    interpretation = 'Moderate MDR association'
                else:
                    interpretation = 'Low MDR discrimination'
            else:  # species
                if importance_level == 'High':
                    interpretation = 'Species-specific resistance pattern'
                elif importance_level == 'Moderate':
                    interpretation = 'Moderate species association'
                else:
                    interpretation = 'Low species discrimination'
            
            importance_data.append({
                'Antibiotic': antibiotic,
                'Class': ab_class,
                'Model': model_name,
                'Importance_Score': round(score, 4),
                'Importance_Level': importance_level,
                'Resistance_Rate': f"{resistance_rate*100:.1f}%" if resistance_rate else 'N/A',
                'Task': task_type.upper(),
                'Interpretation': interpretation
            })
    
    table = pd.DataFrame(importance_data)
    if not table.empty:
        table = table.sort_values(['Model', 'Importance_Score'], ascending=[True, False])
    
    return table


def run_supervised_pipeline(df: pd.DataFrame,
                            feature_cols: List[str],
                            target_col: str,
                            test_size: float = 0.2,
                            random_state: int = 42,
                            task_type: str = None) -> Dict:
    """
    Main LEAKAGE-SAFE supervised learning pipeline (Phase 4).
    
    CRITICAL IMPROVEMENTS (Phase 3):
    ================================
    1. LEAKAGE-SAFE: Split BEFORE scaling/imputation
    2. STRICT FEATURE SEPARATION: Only resistance fingerprints
    3. MACRO METRICS: Prevents class imbalance bias
    4. RATIONALIZED MODELS: 3 models with conceptual grouping
    5. MODEL AGREEMENT: Cross-model feature importance comparison
    6. DISCIPLINED LANGUAGE: Associative, not causal interpretation
    
    This pipeline implements supervised learning for pattern discrimination,
    evaluating how well resistance fingerprints discriminate known categories.
    
    Note: This is NOT for forecasting or predicting future outcomes. The metrics
    quantify how consistently resistance patterns align with known categories.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe
    feature_cols : list
        List of feature column names (resistance fingerprints ONLY)
    target_col : str
        Name of target column (species or MDR category)
    test_size : float
        Test set proportion (default 0.2 = 80%-20% split)
    random_state : int
        Random seed for reproducibility (MUST be fixed and reported)
    task_type : str, optional
        'species' or 'mdr' for task-specific interpretation
    
    Returns:
    --------
    dict
        Complete pipeline results with leakage-safe preprocessing info
    """
    print("=" * 70)
    print("PHASE 4: Supervised Learning for Pattern Discrimination")
    print("        (Leakage-Safe Pipeline with Macro-Averaged Metrics)")
    print("=" * 70)
    
    # Determine task type
    if task_type is None:
        if 'MDR' in target_col.upper():
            task_type = 'mdr'
        else:
            task_type = 'species'
    
    print(f"\n1. TASK: {task_type.upper()} discrimination")
    print(f"   Target variable: {target_col}")
    
    # ===========================================================================
    # MDR TARGET TRANSPARENCY (Phase 3 Requirement 1.3)
    # ===========================================================================
    if task_type == 'mdr':
        print("\n   *** MDR TARGET TRANSPARENCY NOTE ***")
        print("   MDR label is derived from the SAME AST features used as input.")
        print("   This task evaluates SELF-CONSISTENCY DISCRIMINATION:")
        print("   - How consistently do resistance fingerprints align with MDR status?")
        print("   - This is NOT prediction; it's a pattern consistency check.")
        print("   **************************************")
    
    # Prepare data with LEAKAGE-SAFE preprocessing
    print(f"\n2. Preparing data (LEAKAGE-SAFE pipeline)...")
    print(f"   - random_state: {random_state} (fixed for reproducibility)")
    print(f"   - test_size: {test_size*100:.0f}%")
    
    result = prepare_data_for_classification(
        df, feature_cols, target_col, test_size, random_state
    )
    X_train, X_test, y_train, y_test, label_encoder, scaler, imputer, preprocessing_info = result
    
    print(f"   - Training samples: {preprocessing_info['n_train_samples']}")
    print(f"   - Test samples: {preprocessing_info['n_test_samples']}")
    print(f"   - Features used: {preprocessing_info['n_features']} (resistance fingerprints only)")
    print(f"   - Preprocessing order: {' -> '.join(preprocessing_info['preprocessing_order'])}")
    print(f"   - Leakage-safe: {preprocessing_info['leakage_safe']}")
    
    if label_encoder:
        print(f"   - Classes: {list(label_encoder.classes_)}")
    
    # Get feature names for interpretation
    feature_names = [c.replace('_encoded', '') for c in preprocessing_info['feature_columns']]
    
    # ===========================================================================
    # TRAIN AND EVALUATE RATIONALIZED MODEL SET
    # ===========================================================================
    print("\n3. Training RATIONALIZED MODEL SET (Linear, Tree-based, Distance-based)...")
    model_results = run_all_models(
        X_train, X_test, y_train, y_test, feature_names, label_encoder
    )
    
    # Compare models with MACRO-averaged metrics
    print("\n4. Model Comparison (MACRO-averaged metrics):")
    comparison_df = compare_models(model_results)
    print(comparison_df.to_string(index=False))
    
    # Get best model based on macro F1
    best_model_name = comparison_df.iloc[0]['Model']
    best_model = model_results[best_model_name]['model']
    best_f1 = comparison_df.iloc[0]['F1-Score (Macro)']
    
    # Dynamic interpretation based on F1 score
    if best_f1 >= 0.8:
        f1_interpretation = "strong discriminative capacity"
    elif best_f1 >= 0.6:
        f1_interpretation = "moderate discriminative capacity"
    elif best_f1 >= 0.4:
        f1_interpretation = "limited discriminative capacity"
    else:
        f1_interpretation = "weak pattern alignment"
    
    print(f"\n5. Best performing model: {best_model_name}")
    print(f"   Category: {model_results[best_model_name]['category']}")
    print(f"   F1-Score (Macro): {best_f1:.4f}")
    print(f"   Interpretation: Model demonstrates {f1_interpretation}")
    
    # ===========================================================================
    # FEATURE IMPORTANCE WITH MODEL AGREEMENT CHECK
    # ===========================================================================
    print("\n6. Feature Importance Analysis (ASSOCIATIVE, not causal):")
    feature_imp = model_results[best_model_name]['metrics']['feature_importance']
    
    if isinstance(feature_imp, dict) and 'scores' in feature_imp:
        print(f"   Method: {feature_imp.get('method', 'Unknown')}")
        scores = feature_imp['scores']
    else:
        scores = feature_imp
    
    print("   Top 10 discriminating antibiotics:")
    for i, (feat, imp) in enumerate(list(scores.items())[:10]):
        ab_class = ANTIBIOTIC_CLASSES.get(feat, 'Unknown')
        print(f"   {i+1}. {feat} ({ab_class}): {imp:.4f}")
    
    # Model agreement check (Phase 3 Requirement 6.1)
    print("\n7. Model Agreement Check (cross-model comparison):")
    agreement = compare_feature_importance_across_models(model_results, top_n=10)
    print(f"   {agreement['agreement_summary']}")
    
    # Interpretation with biological restraint
    print("\n8. Interpretation (with BIOLOGICAL RESTRAINT):")
    interpretation = interpret_feature_importance(feature_imp, df, target_col, task_type)
    
    print("   Top discriminating antibiotics by class:")
    for disc in interpretation['top_discriminators']:
        print(f"   - {disc['antibiotic']} ({disc['antibiotic_class']}): "
              f"{disc['importance_score']:.4f} - {disc.get('interpretation', 'N/A')}")
    
    if 'ast_results' in interpretation:
        print("\n   AST Results Context:")
        for ab, stats in interpretation['ast_results'].items():
            print(f"   - {ab}: {stats['resistance_rate']*100:.1f}% resistance rate ({stats['interpretation']})")
    
    print(f"\n   BIOLOGICAL RESTRAINT NOTE:")
    print(f"   {interpretation.get('biological_restraint_note', 'N/A')}")
    
    # Generate antibiotic importance table
    print("\n9. Generating Antibiotic Importance Summary Table...")
    importance_table = create_antibiotic_importance_table(model_results, df, task_type)
    
    # ===========================================================================
    # CONFUSION MATRIX OUTPUT (Phase 3 Requirement 4.1)
    # ===========================================================================
    print(f"\n10. Confusion Matrix ({task_type.upper()} task):")
    best_metrics = model_results[best_model_name]['metrics']
    cm = best_metrics['confusion_matrix']
    class_names = best_metrics.get('class_names', [])
    
    print(f"   Classes: {class_names}")
    print("   Matrix:")
    for i, row in enumerate(cm):
        class_label = class_names[i] if i < len(class_names) else f"Class {i}"
        print(f"   {class_label}: {row}")
    
    # Complete results
    pipeline_results = {
        'task_type': task_type,
        'target_variable': target_col,
        'preprocessing_info': preprocessing_info,
        'training_samples': preprocessing_info['n_train_samples'],
        'test_samples': preprocessing_info['n_test_samples'],
        'classes': label_encoder.classes_.tolist() if label_encoder else None,
        'model_results': {name: result['metrics'] for name, result in model_results.items()},
        'comparison': comparison_df.to_dict('records'),
        'best_model': {
            'name': best_model_name,
            'category': model_results[best_model_name]['category'],
            'metrics': model_results[best_model_name]['metrics'],
            'model_object': best_model
        },
        'model_agreement': agreement,
        'scaler': scaler,
        'imputer': imputer,
        'label_encoder': label_encoder,
        'interpretation': interpretation,
        'antibiotic_importance_table': importance_table.to_dict('records') if not importance_table.empty else [],
        'mdr_transparency_note': (
            "MDR label is derived from the same AST features used as input. "
            "MDR discrimination evaluates self-consistency, not prediction."
        ) if task_type == 'mdr' else None
    }
    
    return pipeline_results


def run_species_discrimination(df: pd.DataFrame,
                               feature_cols: List[str],
                               random_state: int = 42) -> Dict:
    """
    Run supervised learning for SPECIES discrimination (Task A).
    
    TASK SEPARATION (Phase 3 Requirement 2.1):
    This is an INDEPENDENT experiment separate from MDR discrimination.
    
    | Task | Input | Target | Type |
    | Species discrimination | Resistance fingerprints | Species | Multi-class |
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe
    feature_cols : list
        List of feature column names (resistance fingerprints only)
    random_state : int
        Random seed for reproducibility
    
    Returns:
    --------
    dict
        Pipeline results for species discrimination
    """
    print("\n" + "=" * 70)
    print("TASK A: SPECIES DISCRIMINATION ANALYSIS")
    print("       (Independent Pipeline - Multi-class Classification)")
    print("=" * 70)
    
    return run_supervised_pipeline(
        df, feature_cols, 'ISOLATE_ID', 
        random_state=random_state, task_type='species'
    )


def run_mdr_discrimination(df: pd.DataFrame,
                           feature_cols: List[str],
                           random_state: int = 42) -> Dict:
    """
    Run supervised learning for MDR discrimination (Task B).
    
    TASK SEPARATION (Phase 3 Requirement 2.1):
    This is an INDEPENDENT experiment separate from species discrimination.
    
    | Task | Input | Target | Type |
    | MDR discrimination | Resistance fingerprints | MDR flag | Binary |
    
    MDR TARGET TRANSPARENCY (Phase 3 Requirement 1.3):
    The MDR label is derived from the SAME AST features used as input.
    This task is treated as SELF-CONSISTENCY DISCRIMINATION:
    - Evaluates how consistently resistance fingerprints align with MDR status
    - Neutralizes "predicted MDR from MDR" objection by explicit acknowledgment
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe
    feature_cols : list
        List of feature column names (resistance fingerprints only)
    random_state : int
        Random seed for reproducibility
    
    Returns:
    --------
    dict
        Pipeline results for MDR discrimination
    """
    print("\n" + "=" * 70)
    print("TASK B: MDR DISCRIMINATION ANALYSIS")
    print("       (Independent Pipeline - Binary Classification)")
    print("       NOTE: MDR is derived from same AST features (self-consistency)")
    print("=" * 70)
    
    return run_supervised_pipeline(
        df, feature_cols, 'MDR_CATEGORY',
        random_state=random_state, task_type='mdr'
    )


def save_model(model, scaler, label_encoder, output_path: str, imputer=None, preprocessing_info=None):
    """
    Save trained model and preprocessors.
    
    Parameters:
    -----------
    model : sklearn estimator
        Trained model
    scaler : StandardScaler
        Fitted scaler
    label_encoder : LabelEncoder
        Fitted label encoder
    output_path : str
        Path to save the model
    imputer : SimpleImputer, optional
        Fitted imputer for missing value handling
    preprocessing_info : dict, optional
        Preprocessing metadata for reproducibility
    """
    import joblib
    
    model_data = {
        'model': model,
        'scaler': scaler,
        'label_encoder': label_encoder,
        'imputer': imputer,
        'preprocessing_info': preprocessing_info
    }
    
    joblib.dump(model_data, output_path)
    print(f"Model saved to: {output_path}")


def load_model(model_path: str) -> Dict:
    """
    Load saved model and preprocessors.
    
    Parameters:
    -----------
    model_path : str
        Path to saved model
    
    Returns:
    --------
    dict
        Dictionary with model, scaler, and label_encoder
    """
    import joblib
    return joblib.load(model_path)


if __name__ == "__main__":
    from pathlib import Path
    
    print("=" * 70)
    print("PHASE 3 IMPLEMENTATION: Concrete Supervised Learning Improvements")
    print("=" * 70)
    print("\nImplemented improvements:")
    print("  1. Leakage-safe preprocessing (split BEFORE scaling)")
    print("  2. Strict feature-label separation (resistance fingerprints only)")
    print("  3. MDR target transparency (self-consistency acknowledgment)")
    print("  4. Task separation (independent species and MDR pipelines)")
    print("  5. Rationalized model set (3 models: Linear, Tree-based, Distance-based)")
    print("  6. Macro-averaged metrics (prevents class imbalance bias)")
    print("  7. Model agreement check (cross-model feature comparison)")
    print("  8. Disciplined interpretation language (associative, not causal)")
    print("=" * 70)
    
    project_root = Path(__file__).parent.parent.parent
    analysis_path = project_root / "data" / "processed" / "analysis_ready_dataset.csv"
    
    if analysis_path.exists():
        df = pd.read_csv(analysis_path)
        feature_cols = [c for c in df.columns if c.endswith('_encoded')]
        
        # Create models directory
        models_dir = project_root / "data" / "models"
        models_dir.mkdir(exist_ok=True, parents=True)
        
        # ===========================================================================
        # TASK B: MDR DISCRIMINATION (Binary Classification)
        # ===========================================================================
        if 'MDR_CATEGORY' in df.columns:
            mdr_results = run_mdr_discrimination(df, feature_cols)
            
            # Save best model with all preprocessors
            save_model(
                mdr_results['best_model']['model_object'],
                mdr_results['scaler'],
                mdr_results['label_encoder'],
                str(models_dir / "mdr_classifier.joblib"),
                imputer=mdr_results.get('imputer'),
                preprocessing_info=mdr_results.get('preprocessing_info')
            )
            
            # Save antibiotic importance table
            if mdr_results.get('antibiotic_importance_table'):
                import json
                table_path = models_dir / "mdr_antibiotic_importance.json"
                with open(table_path, 'w') as f:
                    json.dump(mdr_results['antibiotic_importance_table'], f, indent=2)
                print(f"\nMDR antibiotic importance table saved to: {table_path}")
        
        # ===========================================================================
        # TASK A: SPECIES DISCRIMINATION (Multi-class Classification)
        # ===========================================================================
        if 'ISOLATE_ID' in df.columns and df['ISOLATE_ID'].nunique() > 1:
            species_results = run_species_discrimination(df, feature_cols)
            
            save_model(
                species_results['best_model']['model_object'],
                species_results['scaler'],
                species_results['label_encoder'],
                str(models_dir / "species_classifier.joblib"),
                imputer=species_results.get('imputer'),
                preprocessing_info=species_results.get('preprocessing_info')
            )
            
            # Save antibiotic importance table
            if species_results.get('antibiotic_importance_table'):
                import json
                table_path = models_dir / "species_antibiotic_importance.json"
                with open(table_path, 'w') as f:
                    json.dump(species_results['antibiotic_importance_table'], f, indent=2)
                print(f"\nSpecies antibiotic importance table saved to: {table_path}")
        
        print("\n" + "=" * 70)
        print("SUPERVISED LEARNING PIPELINE COMPLETE")
        print("=" * 70)
    else:
        print(f"Analysis-ready dataset not found at {analysis_path}")
