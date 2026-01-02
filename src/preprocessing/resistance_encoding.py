"""
Resistance Encoding Module for AMR Thesis Project
Phase 2.4 - Encode AST results and create resistance fingerprints
"""

import pandas as pd
import numpy as np
from typing import List, Tuple, Dict, Optional


# Resistance encoding mapping
RESISTANCE_ENCODING = {
    'S': 0,  # Susceptible
    'I': 1,  # Intermediate
    'R': 2   # Resistant
}

# Reverse mapping
RESISTANCE_DECODING = {v: k for k, v in RESISTANCE_ENCODING.items()}


def encode_resistance_value(value) -> Optional[int]:
    """
    Encode a single resistance value.
    
    Parameters:
    -----------
    value : str
        Resistance interpretation (S, I, or R)
    
    Returns:
    --------
    int or None
        Encoded value (0, 1, 2) or None for missing
    """
    if pd.isna(value) or value is None:
        return np.nan
    
    value_str = str(value).strip().upper()
    return RESISTANCE_ENCODING.get(value_str, np.nan)


def decode_resistance_value(value) -> Optional[str]:
    """
    Decode an encoded resistance value back to S/I/R.
    
    Parameters:
    -----------
    value : int
        Encoded resistance value (0, 1, 2)
    
    Returns:
    --------
    str or None
        Resistance interpretation (S, I, R) or None
    """
    if pd.isna(value):
        return None
    
    return RESISTANCE_DECODING.get(int(value), None)


def encode_resistance_profile(df: pd.DataFrame, 
                              antibiotic_cols: List[str]) -> pd.DataFrame:
    """
    Encode all resistance values in the dataframe.
    
    Creates a numerical resistance fingerprint for each isolate.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe with S/I/R values
    antibiotic_cols : list
        List of antibiotic column names to encode
    
    Returns:
    --------
    pd.DataFrame
        Dataframe with encoded resistance values
    """
    df_encoded = df.copy()
    
    for col in antibiotic_cols:
        if col in df_encoded.columns:
            df_encoded[col] = df_encoded[col].apply(encode_resistance_value)
    
    return df_encoded


def create_binary_resistance_matrix(df: pd.DataFrame,
                                   antibiotic_cols: List[str]) -> pd.DataFrame:
    """
    Create binary resistance indicators (R vs non-R).
    
    Useful for certain analyses where only resistant/not resistant matters.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe with S/I/R values
    antibiotic_cols : list
        List of antibiotic column names
    
    Returns:
    --------
    pd.DataFrame
        Dataframe with binary resistance values (1=R, 0=S or I)
    """
    df_binary = df.copy()
    
    for col in antibiotic_cols:
        if col in df_binary.columns:
            # R -> 1, S or I -> 0
            df_binary[col] = df_binary[col].apply(
                lambda x: 1 if str(x).strip().upper() == 'R' else (0 if pd.notna(x) else np.nan)
            )
    
    return df_binary


def extract_resistance_matrix(df: pd.DataFrame,
                             antibiotic_cols: List[str]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Separate feature matrix (antibiotics) from metadata.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Full dataframe with encoded resistance values
    antibiotic_cols : list
        List of antibiotic column names
    
    Returns:
    --------
    tuple
        (Feature matrix, Metadata dataframe)
    """
    # Identify existing antibiotic columns
    existing_ab_cols = [c for c in antibiotic_cols if c in df.columns]
    
    # Feature matrix (antibiotics only)
    feature_matrix = df[existing_ab_cols].copy()
    
    # Metadata (everything else)
    metadata_cols = [c for c in df.columns if c not in existing_ab_cols]
    metadata = df[metadata_cols].copy()
    
    return feature_matrix, metadata


def get_resistance_fingerprint(row: pd.Series, antibiotic_cols: List[str]) -> str:
    """
    Generate a string representation of the resistance fingerprint.
    
    Parameters:
    -----------
    row : pd.Series
        Row from dataframe with resistance values
    antibiotic_cols : list
        List of antibiotic column names
    
    Returns:
    --------
    str
        Resistance fingerprint string (e.g., "SSIRSRRSIS")
    """
    fingerprint_parts = []
    
    for col in antibiotic_cols:
        if col in row.index:
            value = row[col]
            if pd.notna(value):
                if isinstance(value, (int, float)):
                    fingerprint_parts.append(RESISTANCE_DECODING.get(int(value), '?'))
                else:
                    fingerprint_parts.append(str(value).upper()[0])
            else:
                fingerprint_parts.append('-')
    
    return ''.join(fingerprint_parts)


def create_encoded_dataset(df: pd.DataFrame,
                          min_antibiotic_coverage: float = 50.0) -> Tuple[pd.DataFrame, Dict]:
    """
    Main function to create fully encoded dataset.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Cleaned dataframe with S/I/R values
    min_antibiotic_coverage : float
        Minimum coverage percentage for antibiotics
    
    Returns:
    --------
    tuple
        (Encoded dataframe, encoding info dictionary)
    """
    print("=" * 50)
    print("PHASE 2.4: Encoding Resistance Outcomes")
    print("=" * 50)
    
    # Identify antibiotic columns (exclude metadata and summary columns)
    metadata_cols = ['CODE', 'ISOLATE_ID', 'REGION', 'SITE', 'ENVIRONMENT',
                     'SAMPLING_SOURCE', 'NATIONAL_SITE', 'LOCAL_SITE', 
                     'REPLICATE', 'COLONY', 'ESBL', 'SOURCE_FILE', 
                     'SCORED_RESISTANCE', 'NUM_ANTIBIOTICS_TESTED', 'MAR_INDEX']
    
    antibiotic_cols = [c for c in df.columns if c not in metadata_cols]
    
    # Filter to antibiotics with sufficient coverage
    ab_coverage = {}
    for col in antibiotic_cols:
        if col in df.columns:
            coverage = (df[col].notna().sum() / len(df)) * 100
            ab_coverage[col] = coverage
    
    valid_antibiotics = [c for c, cov in ab_coverage.items() if cov >= min_antibiotic_coverage]
    
    print(f"\n1. Identified {len(valid_antibiotics)} antibiotics with >={min_antibiotic_coverage}% coverage")
    
    # Encode resistance values
    df_encoded = df.copy()
    
    for col in valid_antibiotics:
        df_encoded[f'{col}_encoded'] = df_encoded[col].apply(encode_resistance_value)
    
    print("2. Encoded resistance values (S=0, I=1, R=2)")
    
    # Create fingerprint column
    encoded_cols = [f'{c}_encoded' for c in valid_antibiotics]
    
    df_encoded['resistance_fingerprint'] = df_encoded.apply(
        lambda row: get_resistance_fingerprint(row, valid_antibiotics),
        axis=1
    )
    
    print("3. Generated resistance fingerprints")
    
    # Encoding info
    encoding_info = {
        'encoding_scheme': RESISTANCE_ENCODING,
        'antibiotics_encoded': valid_antibiotics,
        'encoded_columns': encoded_cols,
        'total_isolates': len(df_encoded),
        'coverage_stats': ab_coverage
    }
    
    print(f"\n4. Final encoded dataset: {len(df_encoded)} isolates")
    print(f"   Encoded columns: {len(encoded_cols)}")
    
    return df_encoded, encoding_info


def prepare_feature_matrix_for_clustering(df: pd.DataFrame,
                                          encoded_cols: List[str],
                                          fillna_value: Optional[float] = None) -> pd.DataFrame:
    """
    Prepare the feature matrix for clustering algorithms.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Encoded dataframe
    encoded_cols : list
        List of encoded column names
    fillna_value : float, optional
        Value to fill NaN (if None, rows with NaN are dropped)
    
    Returns:
    --------
    pd.DataFrame
        Feature matrix ready for clustering
    """
    feature_matrix = df[encoded_cols].copy()
    
    if fillna_value is not None:
        feature_matrix = feature_matrix.fillna(fillna_value)
    else:
        # Drop rows with any NaN values
        feature_matrix = feature_matrix.dropna()
    
    return feature_matrix


if __name__ == "__main__":
    from pathlib import Path
    
    project_root = Path(__file__).parent.parent.parent
    clean_path = project_root / "data" / "processed" / "cleaned_dataset.csv"
    
    if clean_path.exists():
        df = pd.read_csv(clean_path)
        df_encoded, info = create_encoded_dataset(df)
        
        # Save encoded dataset
        encoded_path = project_root / "data" / "processed" / "encoded_dataset.csv"
        df_encoded.to_csv(encoded_path, index=False)
        print(f"\nEncoded dataset saved to: {encoded_path}")
        
        # Print encoding info
        print("\nEncoding Info:")
        print(f"  Encoding scheme: {info['encoding_scheme']}")
        print(f"  Antibiotics encoded: {info['antibiotics_encoded']}")
    else:
        print(f"Cleaned dataset not found at {clean_path}")
        print("Run data_cleaning.py first.")
