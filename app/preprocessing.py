"""
Preprocessing Module for AMR Pattern Recognition Dashboard
Phase 6 - Data preprocessing utilities

This module provides preprocessing functions that align with the main pipeline
but are designed for read-only dashboard use.
"""

import pandas as pd
import numpy as np
from typing import List, Tuple, Optional


def encode_resistance_values(df: pd.DataFrame, 
                            antibiotic_cols: List[str]) -> pd.DataFrame:
    """
    Encode S/I/R resistance values to numeric format.
    
    This is for reference only - the dashboard expects pre-encoded data.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe with S/I/R values
    antibiotic_cols : list
        List of antibiotic column names
    
    Returns:
    --------
    pd.DataFrame
        Dataframe with encoded columns added
    """
    encoding_map = {'S': 0, 'I': 1, 'R': 2}
    
    df_encoded = df.copy()
    
    for col in antibiotic_cols:
        if col in df.columns and not col.endswith('_encoded'):
            encoded_col = f'{col}_encoded'
            df_encoded[encoded_col] = df[col].map(encoding_map)
    
    return df_encoded


def calculate_mar_index(df: pd.DataFrame, 
                       encoded_cols: List[str]) -> pd.Series:
    """
    Calculate Multiple Antibiotic Resistance (MAR) index.
    
    MAR = (Number of antibiotics resistant) / (Number of antibiotics tested)
    Reference: Krumperman PH (1983)
    
    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe with encoded resistance values
    encoded_cols : list
        List of encoded column names
    
    Returns:
    --------
    pd.Series
        MAR index values
    """
    existing_cols = [c for c in encoded_cols if c in df.columns]
    
    # Count resistant (value == 2)
    resistance_counts = (df[existing_cols] == 2).sum(axis=1)
    
    # Count tested (non-null values)
    tested_counts = df[existing_cols].notna().sum(axis=1)
    
    # Calculate MAR
    mar_index = resistance_counts / tested_counts
    mar_index = mar_index.fillna(0)
    
    return mar_index


def calculate_mdr_status(df: pd.DataFrame,
                        encoded_cols: List[str]) -> Tuple[pd.Series, pd.Series]:
    """
    Calculate MDR (Multi-Drug Resistant) status.
    
    MDR is defined as resistance to ≥3 antibiotic classes.
    Reference: Magiorakos et al. (2012)
    
    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe with encoded resistance values
    encoded_cols : list
        List of encoded column names
    
    Returns:
    --------
    tuple
        (MDR flag series, MDR category series)
    """
    # Antibiotic class mapping
    antibiotic_classes = {
        'AM': 'Penicillins', 'AMP': 'Penicillins',
        'AMC': 'BetaLactamBLI', 'PRA': 'BetaLactamBLI',
        'CN': 'Ceph1st', 'CF': 'Ceph1st',
        'CPD': 'Ceph3rd4th', 'CTX': 'Ceph3rd4th', 
        'CFT': 'Ceph3rd4th', 'CPT': 'Ceph3rd4th',
        'CFO': 'Cephamycins',
        'CZA': 'CephBLI',
        'IPM': 'Carbapenems', 'MRB': 'Carbapenems',
        'AN': 'Aminoglycosides', 'GM': 'Aminoglycosides', 'N': 'Aminoglycosides',
        'NAL': 'Quinolones', 'ENR': 'Quinolones',
        'DO': 'Tetracyclines', 'TE': 'Tetracyclines',
        'FT': 'Nitrofurans',
        'C': 'Phenicols',
        'SXT': 'FolateInhibitors'
    }
    
    def count_resistant_classes(row):
        resistant_classes = set()
        for col in encoded_cols:
            if col in row.index and pd.notna(row[col]) and row[col] == 2:
                ab_name = col.replace('_encoded', '')
                if ab_name in antibiotic_classes:
                    resistant_classes.add(antibiotic_classes[ab_name])
        return len(resistant_classes)
    
    existing_cols = [c for c in encoded_cols if c in df.columns]
    resistant_class_counts = df[existing_cols].apply(count_resistant_classes, axis=1)
    
    mdr_flag = (resistant_class_counts >= 3).astype(int)
    mdr_category = mdr_flag.map({0: 'Non-MDR', 1: 'MDR'})
    
    return mdr_flag, mdr_category


def get_resistance_statistics(df: pd.DataFrame,
                             encoded_cols: List[str]) -> pd.DataFrame:
    """
    Calculate resistance statistics for each antibiotic.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe with encoded resistance values
    encoded_cols : list
        List of encoded column names
    
    Returns:
    --------
    pd.DataFrame
        Resistance statistics table
    """
    existing_cols = [c for c in encoded_cols if c in df.columns]
    
    stats = []
    for col in existing_cols:
        ab_name = col.replace('_encoded', '')
        values = df[col].dropna()
        n_tested = len(values)
        
        if n_tested > 0:
            n_susceptible = (values == 0).sum()
            n_intermediate = (values == 1).sum()
            n_resistant = (values == 2).sum()
            
            stats.append({
                'Antibiotic': ab_name,
                'Tested': n_tested,
                'Susceptible': n_susceptible,
                'Intermediate': n_intermediate,
                'Resistant': n_resistant,
                'S%': round(n_susceptible / n_tested * 100, 1),
                'I%': round(n_intermediate / n_tested * 100, 1),
                'R%': round(n_resistant / n_tested * 100, 1)
            })
    
    return pd.DataFrame(stats)
