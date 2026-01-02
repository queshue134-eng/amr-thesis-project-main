"""
Data Loader Module for AMR Pattern Recognition Dashboard
Phase 6 - Controlled Data Upload with Schema Validation

This module handles CSV file loading with strict validation:
- Accept CSV only
- Validate required columns
- Validate S/I/R encoding
- Reject invalid schemas with user feedback
"""

import pandas as pd
import streamlit as st
from typing import Tuple, List, Dict, Optional
import os


# Required columns for valid AMR dataset
REQUIRED_COLUMNS = ['CODE']  # At minimum, CODE is required

# Optional but recommended columns
RECOMMENDED_COLUMNS = ['ISOLATE_ID', 'REGION', 'SITE', 'CLUSTER', 'MDR_FLAG', 'MAR_INDEX_COMPUTED']

# Valid resistance values (CLSI standard)
VALID_RESISTANCE_VALUES = {'S', 'I', 'R', 0, 1, 2, 0.0, 1.0, 2.0}

# Valid encoded resistance values
VALID_ENCODED_VALUES = {0, 1, 2, 0.0, 1.0, 2.0}


def validate_csv_schema(df: pd.DataFrame) -> Tuple[bool, List[str], List[str]]:
    """
    Validate the CSV schema for AMR analysis.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Loaded dataframe to validate
        
    Returns:
    --------
    tuple
        (is_valid, errors_list, warnings_list)
    """
    errors = []
    warnings = []
    
    # Check for required columns
    missing_required = [col for col in REQUIRED_COLUMNS if col not in df.columns]
    if missing_required:
        errors.append(f"Missing required columns: {', '.join(missing_required)}")
    
    # Check for recommended columns
    missing_recommended = [col for col in RECOMMENDED_COLUMNS if col not in df.columns]
    if missing_recommended:
        warnings.append(f"Missing recommended columns: {', '.join(missing_recommended)}")
    
    # Check for encoded antibiotic columns
    encoded_cols = [c for c in df.columns if c.endswith('_encoded')]
    if not encoded_cols:
        # Check for raw antibiotic columns (S/I/R format)
        possible_ab_cols = [c for c in df.columns if c not in REQUIRED_COLUMNS + RECOMMENDED_COLUMNS 
                          and not c.endswith('_encoded') and not c.startswith('_')]
        if possible_ab_cols:
            warnings.append("No encoded columns found. Data may need preprocessing.")
        else:
            errors.append("No antibiotic resistance data found (no '_encoded' columns)")
    
    # Validate encoded values are in valid range
    for col in encoded_cols:
        unique_vals = df[col].dropna().unique()
        invalid_vals = [v for v in unique_vals if v not in VALID_ENCODED_VALUES]
        if invalid_vals:
            errors.append(f"Invalid values in {col}: {invalid_vals}. Expected 0 (S), 1 (I), or 2 (R).")
    
    # Check minimum row count
    if len(df) < 1:
        errors.append("Dataset is empty")
    elif len(df) < 5:
        warnings.append(f"Dataset has only {len(df)} rows. Results may be unreliable.")
    
    is_valid = len(errors) == 0
    return is_valid, errors, warnings


def get_antibiotic_columns(df: pd.DataFrame) -> List[str]:
    """
    Identify antibiotic resistance columns in the dataframe.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe
        
    Returns:
    --------
    list
        List of antibiotic column names (preferring encoded columns)
    """
    # First try encoded columns
    encoded_cols = [c for c in df.columns if c.endswith('_encoded')]
    if encoded_cols:
        return encoded_cols
    
    # Fall back to identifying by common antibiotic abbreviations
    possible_antibiotics = ['AM', 'AMC', 'CPT', 'CN', 'CF', 'CPD', 'CTX', 'CFO', 
                           'CFT', 'CZA', 'IPM', 'AN', 'GM', 'N', 'NAL', 'ENR',
                           'MRB', 'PRA', 'DO', 'TE', 'FT', 'C', 'SXT']
    
    return [c for c in df.columns if any(ab in c.upper() for ab in possible_antibiotics)]


def load_uploaded_file(uploaded_file) -> Tuple[Optional[pd.DataFrame], List[str], List[str]]:
    """
    Load and validate an uploaded CSV file.
    
    Data is processed in memory only - not stored on disk.
    
    Parameters:
    -----------
    uploaded_file : streamlit.UploadedFile
        File uploaded through Streamlit file uploader
        
    Returns:
    --------
    tuple
        (dataframe or None, errors_list, warnings_list)
    """
    errors = []
    warnings = []
    
    if uploaded_file is None:
        return None, ["No file uploaded"], []
    
    # Check file type
    if not uploaded_file.name.endswith('.csv'):
        errors.append("Invalid file type. Please upload a CSV file.")
        return None, errors, warnings
    
    # Try to load the CSV
    try:
        df = pd.read_csv(uploaded_file)
    except Exception as e:
        errors.append(f"Error reading CSV file: {str(e)}")
        return None, errors, warnings
    
    # Validate schema
    is_valid, schema_errors, schema_warnings = validate_csv_schema(df)
    errors.extend(schema_errors)
    warnings.extend(schema_warnings)
    
    if not is_valid:
        return None, errors, warnings
    
    return df, [], warnings


def load_default_dataset(default_path: str) -> Tuple[Optional[pd.DataFrame], List[str], List[str]]:
    """
    Load the default dataset from disk.
    
    Parameters:
    -----------
    default_path : str
        Path to the default dataset
        
    Returns:
    --------
    tuple
        (dataframe or None, errors_list, warnings_list)
    """
    errors = []
    warnings = []
    
    if not os.path.exists(default_path):
        return None, [f"Default dataset not found at {default_path}"], []
    
    try:
        df = pd.read_csv(default_path)
        is_valid, schema_errors, schema_warnings = validate_csv_schema(df)
        errors.extend(schema_errors)
        warnings.extend(schema_warnings)
        
        if not is_valid:
            return None, errors, warnings
        
        return df, [], warnings
    except Exception as e:
        errors.append(f"Error loading default dataset: {str(e)}")
        return None, errors, warnings


def get_dataset_info(df: pd.DataFrame) -> Dict:
    """
    Get summary information about the loaded dataset.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Loaded dataframe
        
    Returns:
    --------
    dict
        Dictionary with dataset summary information
    """
    info = {
        'n_isolates': len(df),
        'n_columns': len(df.columns),
        'antibiotic_columns': get_antibiotic_columns(df),
        'n_antibiotics': len(get_antibiotic_columns(df)),
        'has_clusters': 'CLUSTER' in df.columns,
        'has_mdr': 'MDR_FLAG' in df.columns or 'MDR_CATEGORY' in df.columns,
        'has_region': 'REGION' in df.columns,
        'has_species': 'ISOLATE_ID' in df.columns,
    }
    
    if info['has_clusters']:
        info['n_clusters'] = df['CLUSTER'].nunique()
    
    if info['has_mdr'] and 'MDR_FLAG' in df.columns:
        info['mdr_proportion'] = df['MDR_FLAG'].mean() * 100
    
    if info['has_region']:
        info['regions'] = df['REGION'].unique().tolist()
    
    if info['has_species']:
        info['n_species'] = df['ISOLATE_ID'].nunique()
    
    return info


def display_upload_status(errors: List[str], warnings: List[str]):
    """
    Display upload status messages in Streamlit.
    
    Parameters:
    -----------
    errors : list
        List of error messages
    warnings : list
        List of warning messages
    """
    if errors:
        for error in errors:
            st.error(f"❌ {error}")
    
    if warnings:
        for warning in warnings:
            st.warning(f"⚠️ {warning}")


def display_expected_format():
    """Display the expected data format for user guidance."""
    st.markdown("""
    ### Expected Data Format
    
    The dataset should be a **CSV file** containing:
    
    | Column Type | Required | Description |
    |-------------|----------|-------------|
    | `CODE` | Yes | Unique isolate identifier |
    | `ISOLATE_ID` | Recommended | Bacterial species name |
    | `REGION` | Recommended | Geographic region |
    | `CLUSTER` | Recommended | Cluster assignment (from analysis) |
    | `MDR_FLAG` | Recommended | MDR status (0 or 1) |
    | `MAR_INDEX_COMPUTED` | Recommended | MAR index value |
    | `{AB}_encoded` | Yes | Encoded resistance data (0=S, 1=I, 2=R) |
    
    #### Resistance Encoding
    - **0** = Susceptible (S)
    - **1** = Intermediate (I)
    - **2** = Resistant (R)
    
    #### Example
    | CODE | ISOLATE_ID | AM_encoded | AMC_encoded | ... |
    |------|------------|------------|-------------|-----|
    | EC_001 | Escherichia coli | 2 | 1 | ... |
    | KP_002 | Klebsiella pneumoniae | 0 | 2 | ... |
    """)
