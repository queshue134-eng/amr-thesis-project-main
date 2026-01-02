"""
Model Handling Module for AMR Pattern Recognition Dashboard
Phase 6 - Model Loading and Read-Only Outputs

This module provides model handling with:
- Load pretrained models using joblib
- No retraining in deployment (all training disabled)
- Read-only model outputs (confusion matrices, consistency metrics)
- No single-sample inference to prevent clinical misuse
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Dict, List, Optional, Tuple
import joblib


def load_pretrained_model(model_path: str) -> Optional[Dict]:
    """
    Load a pretrained model from disk.
    
    Parameters:
    -----------
    model_path : str
        Path to the saved model file (.joblib)
    
    Returns:
    --------
    dict or None
        Model data dictionary containing model, scaler, label_encoder, etc.
    """
    if not os.path.exists(model_path):
        return None
    
    try:
        model_data = joblib.load(model_path)
        return model_data
    except Exception as e:
        print(f"Error loading model: {e}")
        return None


def get_available_models(models_dir: str) -> List[Dict]:
    """
    Get list of available pretrained models.
    
    Parameters:
    -----------
    models_dir : str
        Directory containing model files
    
    Returns:
    --------
    list
        List of dictionaries with model info
    """
    available_models = []
    
    if not os.path.exists(models_dir):
        return available_models
    
    for filename in os.listdir(models_dir):
        if filename.endswith('.joblib'):
            model_path = os.path.join(models_dir, filename)
            model_name = filename.replace('.joblib', '').replace('_', ' ').title()
            
            # Determine model type from filename
            if 'mdr' in filename.lower():
                model_type = 'MDR Discrimination'
                task_type = 'Binary Classification'
            elif 'species' in filename.lower():
                model_type = 'Species Discrimination'
                task_type = 'Multi-class Classification'
            else:
                model_type = 'Unknown'
                task_type = 'Unknown'
            
            available_models.append({
                'name': model_name,
                'filename': filename,
                'path': model_path,
                'model_type': model_type,
                'task_type': task_type
            })
    
    return available_models


def get_model_info(model_data: Dict) -> Dict:
    """
    Extract information from loaded model data.
    
    Parameters:
    -----------
    model_data : dict
        Loaded model dictionary
    
    Returns:
    --------
    dict
        Model information
    """
    info = {
        'model_type': None,
        'model_name': None,
        'classes': None,
        'n_features': None,
        'preprocessing_info': None,
        'has_feature_importance': False
    }
    
    model = model_data.get('model')
    if model is not None:
        info['model_type'] = type(model).__name__
        info['model_name'] = type(model).__name__
        
        # Check for feature importance
        if hasattr(model, 'feature_importances_') or hasattr(model, 'coef_'):
            info['has_feature_importance'] = True
    
    # Get classes from label encoder
    label_encoder = model_data.get('label_encoder')
    if label_encoder is not None and hasattr(label_encoder, 'classes_'):
        info['classes'] = list(label_encoder.classes_)
    
    # Get preprocessing info
    preprocessing_info = model_data.get('preprocessing_info')
    if preprocessing_info:
        info['preprocessing_info'] = preprocessing_info
        info['n_features'] = preprocessing_info.get('n_features')
    
    return info


def extract_feature_importance(model_data: Dict, 
                              feature_names: List[str] = None) -> Dict[str, float]:
    """
    Extract feature importance from a trained model.
    
    Note: This shows ASSOCIATIVE importance, not causal relationships.
    
    Parameters:
    -----------
    model_data : dict
        Loaded model dictionary
    feature_names : list, optional
        List of feature names
    
    Returns:
    --------
    dict
        Feature importance scores
    """
    model = model_data.get('model')
    if model is None:
        return {}
    
    importance_scores = {}
    
    if hasattr(model, 'feature_importances_'):
        # Tree-based models (Random Forest, etc.)
        importances = model.feature_importances_
    elif hasattr(model, 'coef_'):
        # Linear models (Logistic Regression, etc.)
        coef = np.abs(model.coef_)
        if len(coef.shape) > 1:
            importances = coef.mean(axis=0)
        else:
            importances = coef
    else:
        return {}
    
    # Match with feature names
    if feature_names is None:
        preprocessing_info = model_data.get('preprocessing_info', {})
        feature_names = preprocessing_info.get('feature_columns', [])
    
    if len(feature_names) == len(importances):
        for name, imp in zip(feature_names, importances):
            importance_scores[name] = float(imp)
    else:
        # Use generic names if mismatch
        for i, imp in enumerate(importances):
            importance_scores[f'Feature_{i}'] = float(imp)
    
    # Sort by importance
    importance_scores = dict(sorted(importance_scores.items(), 
                                   key=lambda x: x[1], reverse=True))
    
    return importance_scores


def create_confusion_matrix_plot(confusion_matrix: np.ndarray,
                                class_names: List[str],
                                title: str = 'Confusion Matrix',
                                figsize: Tuple[int, int] = (8, 6)) -> plt.Figure:
    """
    Create a confusion matrix visualization.
    
    Note: This shows pattern consistency metrics, NOT predictive performance.
    
    Parameters:
    -----------
    confusion_matrix : np.ndarray
        Confusion matrix values
    class_names : list
        Names of classes
    title : str
        Plot title
    figsize : tuple
        Figure size
    
    Returns:
    --------
    matplotlib.Figure
        Confusion matrix figure
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    im = ax.imshow(confusion_matrix, interpolation='nearest', cmap='Blues')
    ax.figure.colorbar(im, ax=ax)
    
    ax.set(xticks=np.arange(len(class_names)),
           yticks=np.arange(len(class_names)),
           xticklabels=class_names,
           yticklabels=class_names,
           ylabel='Actual',
           xlabel='Assigned')
    
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    
    # Add text annotations
    thresh = confusion_matrix.max() / 2.
    for i in range(len(class_names)):
        for j in range(len(class_names)):
            ax.text(j, i, format(confusion_matrix[i, j], 'd'),
                   ha="center", va="center",
                   color="white" if confusion_matrix[i, j] > thresh else "black")
    
    ax.set_title(f'{title}\n(Pattern Consistency Metric - Not Predictive)', fontweight='bold')
    
    plt.tight_layout()
    return fig


def format_metrics_table(metrics: Dict) -> pd.DataFrame:
    """
    Format model metrics as a display table.
    
    Note: Metrics are labeled as 'consistency metrics' not 'performance metrics'.
    
    Parameters:
    -----------
    metrics : dict
        Dictionary of metric values
    
    Returns:
    --------
    pd.DataFrame
        Formatted metrics table
    """
    # Rename metrics to emphasize consistency, not prediction
    metric_names = {
        'accuracy': 'Consistency Accuracy',
        'precision_macro': 'Macro Precision (Consistency)',
        'recall_macro': 'Macro Recall (Consistency)',
        'f1_score_macro': 'Macro F1 (Consistency)',
        'precision_weighted': 'Weighted Precision (Consistency)',
        'recall_weighted': 'Weighted Recall (Consistency)',
        'f1_score_weighted': 'Weighted F1 (Consistency)'
    }
    
    rows = []
    for key, display_name in metric_names.items():
        if key in metrics:
            value = metrics[key]
            if isinstance(value, (int, float)):
                rows.append({
                    'Metric': display_name,
                    'Value': f'{value:.4f}'
                })
    
    return pd.DataFrame(rows)


def get_model_disclaimer() -> str:
    """
    Get the standard model interpretation disclaimer.
    
    Returns:
    --------
    str
        Disclaimer text
    """
    return """
    **⚠️ Important Interpretation Note**
    
    The metrics shown here quantify how consistently resistance patterns align 
    with known categories (species, MDR status). They represent **pattern consistency 
    metrics** within the analyzed dataset, NOT predictive performance for future samples.
    
    - **Accuracy**: Proportion of isolates where resistance patterns align with category
    - **Precision/Recall/F1**: Consistency of pattern-category alignment across groups
    
    These metrics should be interpreted as demonstrating **discriminative capacity** 
    of resistance fingerprints, not clinical prediction accuracy.
    """


def get_feature_importance_disclaimer() -> str:
    """
    Get the standard feature importance disclaimer.
    
    Returns:
    --------
    str
        Disclaimer text
    """
    return """
    **⚠️ Feature Importance Interpretation**
    
    Feature importance scores indicate **ASSOCIATIVE** patterns only:
    - High importance means the antibiotic helps discriminate between groups
    - This does NOT imply the antibiotic **causes** group membership
    - Biological mechanisms should not be inferred from importance scores alone
    
    These scores reflect statistical associations in the training data.
    """
