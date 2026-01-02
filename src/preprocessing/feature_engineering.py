"""
Feature Engineering Module for AMR Thesis Project
Phase 2.5 - Compute MAR index, MDR flag, and other derived features

This module implements formalized feature engineering with explicit definitions:

MAR Index (Multiple Antibiotic Resistance Index):
    Formula: MAR = a / b
    Where:
        a = Number of antibiotics to which the isolate is resistant (R)
        b = Total number of antibiotics tested on the isolate
    Reference: Krumperman PH. (1983). Multiple antibiotic resistance indexing of 
               Escherichia coli to identify high-risk sources of fecal contamination 
               of foods. Applied and Environmental Microbiology, 46(1), 165-170.

MDR (Multi-Drug Resistant) Classification:
    Definition: An isolate is classified as MDR if it exhibits resistance to 
                at least one agent in ≥3 antimicrobial categories.
    Reference: Magiorakos AP, et al. (2012). Multidrug-resistant, extensively 
               drug-resistant and pandrug-resistant bacteria: an international 
               expert proposal for interim standard definitions for acquired 
               resistance. Clinical Microbiology and Infection, 18(3), 268-281.
               DOI: 10.1111/j.1469-0691.2011.03570.x
"""

import pandas as pd
import numpy as np
from typing import List, Tuple, Dict, Optional

# Import from centralized configuration
try:
    from config import (
        ANTIBIOTIC_CLASSES, 
        RESISTANCE_ENCODING, 
        RESISTANCE_THRESHOLD,
        MDR_MIN_CLASSES,
        METADATA_COLUMNS,
        get_mdr_classes_for_species  # Species-specific MDR class mappings
    )
except ImportError:
    # Fallback for standalone execution
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).parent.parent))
    from config import (
        ANTIBIOTIC_CLASSES, 
        RESISTANCE_ENCODING, 
        RESISTANCE_THRESHOLD,
        MDR_MIN_CLASSES,
        METADATA_COLUMNS,
        get_mdr_classes_for_species  # Species-specific MDR class mappings
    )

# Note: ANTIBIOTIC_CLASSES and RESISTANCE_ENCODING are now imported from config.py
# This ensures a single source of truth across all modules


def compute_mar_index(row: pd.Series, 
                      antibiotic_cols: List[str],
                      resistance_threshold: int = 2) -> Optional[float]:
    """
    Compute Multiple Antibiotic Resistance (MAR) Index.
    
    Formula: MAR = a / b
    Where:
        a = Number of antibiotics to which the isolate is resistant (encoded value >= threshold)
        b = Total number of antibiotics tested on the isolate (non-null values)
    
    Reference: Krumperman PH. (1983). Multiple antibiotic resistance indexing of 
               Escherichia coli. Applied and Environmental Microbiology, 46(1), 165-170.
    
    Parameters:
    -----------
    row : pd.Series
        Row containing resistance values (encoded: S=0, I=1, R=2)
    antibiotic_cols : list
        List of antibiotic column names (typically ending in '_encoded')
    resistance_threshold : int
        Encoded value considered resistant (default: 2 for R)
        Note: Using 2 means only 'R' is counted as resistant
        Using 1 would include 'I' (intermediate) as resistant
    
    Returns:
    --------
    float or None
        MAR index value (0.0 to 1.0) or None if no antibiotics were tested
    """
    tested = 0
    resistant = 0
    
    for col in antibiotic_cols:
        if col in row.index:
            value = row[col]
            if pd.notna(value):
                tested += 1
                try:
                    if int(value) >= resistance_threshold:
                        resistant += 1
                except (ValueError, TypeError):
                    # Handle non-numeric values
                    if str(value).strip().upper() == 'R':
                        resistant += 1
    
    if tested == 0:
        return None
    
    return round(resistant / tested, 4)


def compute_resistance_count(row: pd.Series,
                            antibiotic_cols: List[str],
                            resistance_threshold: int = 2) -> int:
    """
    Count number of resistant antibiotics.
    
    Parameters:
    -----------
    row : pd.Series
        Row containing resistance values
    antibiotic_cols : list
        List of antibiotic column names
    resistance_threshold : int
        Encoded value considered resistant (default: 2 for R)
    
    Returns:
    --------
    int
        Number of resistant antibiotics
    """
    resistant = 0
    
    for col in antibiotic_cols:
        if col in row.index:
            value = row[col]
            if pd.notna(value):
                try:
                    if int(value) >= resistance_threshold:
                        resistant += 1
                except (ValueError, TypeError):
                    if str(value).strip().upper() == 'R':
                        resistant += 1
    
    return resistant


def count_resistant_classes(row: pd.Series,
                           antibiotic_cols: List[str],
                           resistance_threshold: int = 2,
                           antibiotic_class_mapping: Dict[str, str] = None) -> int:
    """
    Count number of antibiotic classes with at least one resistant antibiotic.
    
    This is used for MDR classification per Magiorakos et al. (2012) definition.
    
    SPECIES-SPECIFIC SUPPORT:
    When antibiotic_class_mapping is provided (e.g., from get_mdr_classes_for_species()),
    only antibiotics defined in that mapping are counted. This excludes intrinsic
    resistances from MDR calculation (e.g., Pseudomonas ampicillin resistance).
    
    Parameters:
    -----------
    row : pd.Series
        Row containing resistance values (encoded)
    antibiotic_cols : list
        List of antibiotic column names
    resistance_threshold : int
        Encoded value considered resistant (default: 2 for R)
    antibiotic_class_mapping : dict, optional
        Species-specific antibiotic-to-class mapping. If None, uses universal
        ANTIBIOTIC_CLASSES (may overestimate MDR for species with intrinsic resistances).
    
    Returns:
    --------
    int
        Number of distinct antibiotic classes with at least one resistant agent
    """
    # Use provided mapping or fall back to universal
    class_mapping = antibiotic_class_mapping if antibiotic_class_mapping else ANTIBIOTIC_CLASSES
    
    resistant_classes = set()
    
    for col in antibiotic_cols:
        if col in row.index:
            value = row[col]
            if pd.notna(value):
                is_resistant = False
                try:
                    if int(value) >= resistance_threshold:
                        is_resistant = True
                except (ValueError, TypeError):
                    if str(value).strip().upper() == 'R':
                        is_resistant = True
                
                if is_resistant:
                    # Get the base antibiotic name (remove _encoded suffix)
                    ab_name = col.replace('_encoded', '').upper()
                    # Only count if antibiotic is in the species-specific mapping
                    ab_class = class_mapping.get(ab_name, None)
                    if ab_class is not None:
                        resistant_classes.add(ab_class)
    
    return len(resistant_classes)


def determine_mdr_status(row: pd.Series,
                        antibiotic_cols: List[str],
                        min_classes: int = 3,
                        resistance_threshold: int = 2,
                        species: str = None) -> bool:
    """
    Determine if an isolate is Multi-Drug Resistant (MDR).
    
    Definition: An isolate is classified as MDR if it exhibits resistance to 
                at least one agent in ≥3 antimicrobial categories.
    
    SPECIES-SPECIFIC SUPPORT:
    When species is provided, uses get_mdr_classes_for_species() to get
    species-appropriate antibiotic class mappings. This ensures intrinsic
    resistances are excluded from MDR calculation (per EUCAST Expert Rules).
    
    Reference: Magiorakos AP, et al. (2012). Multidrug-resistant, extensively 
               drug-resistant and pandrug-resistant bacteria: an international 
               expert proposal for interim standard definitions for acquired 
               resistance. Clinical Microbiology and Infection, 18(3), 268-281.
               DOI: 10.1111/j.1469-0691.2011.03570.x
    
    Parameters:
    -----------
    row : pd.Series
        Row containing encoded resistance values (S=0, I=1, R=2)
    antibiotic_cols : list
        List of antibiotic column names (typically ending in '_encoded')
    min_classes : int
        Minimum number of resistant classes for MDR classification (default: 3)
        Per Magiorakos et al., MDR requires ≥3 antimicrobial categories
    resistance_threshold : int
        Encoded value considered resistant (default: 2 for R only)
    species : str, optional
        Species name for species-specific MDR calculation. If None, uses
        universal ANTIBIOTIC_CLASSES (may overestimate MDR for species with
        intrinsic resistances like Pseudomonas aeruginosa).
    
    Returns:
    --------
    bool
        True if isolate meets MDR criteria, False otherwise
    """
    # Get species-specific class mapping if species provided
    antibiotic_class_mapping = get_mdr_classes_for_species(species) if species else None
    
    resistant_classes_count = count_resistant_classes(
        row, antibiotic_cols, resistance_threshold, antibiotic_class_mapping
    )
    return resistant_classes_count >= min_classes


def _safe_encode_binary_resistance(value, resistance_threshold: int = 2) -> Optional[int]:
    """
    Safely convert encoded resistance value to binary (R=1, non-R=0).
    
    Parameters:
    -----------
    value : any
        Encoded resistance value
    resistance_threshold : int
        Threshold for resistance (default: 2 for R)
    
    Returns:
    --------
    int or None
        1 if resistant, 0 if susceptible/intermediate, None if missing
    """
    if pd.isna(value):
        return None
    try:
        # Try float first to handle '2.0' strings, then convert to int
        return 1 if int(float(value)) >= resistance_threshold else 0
    except (ValueError, TypeError):
        # Handle non-numeric string values
        value_str = str(value).strip().upper()
        if value_str == 'R':
            return 1
        elif value_str in ('S', 'I'):
            return 0
        return None


def create_binary_resistance_features(df: pd.DataFrame,
                                     antibiotic_cols: List[str]) -> pd.DataFrame:
    """
    Create binary resistance features (R vs non-R) for each antibiotic.
    
    Binary features are useful for:
    - Simplified pattern analysis
    - Feature importance interpretation
    - Direct resistance prevalence calculations
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe with encoded resistance values
    antibiotic_cols : list
        List of encoded antibiotic column names
    
    Returns:
    --------
    pd.DataFrame
        Dataframe with additional binary resistance columns
    """
    df_binary = df.copy()
    
    for col in antibiotic_cols:
        if col in df_binary.columns:
            ab_name = col.replace('_encoded', '')
            # R (encoded as 2) -> 1, S or I -> 0, missing -> NaN
            df_binary[f'{ab_name}_RESISTANT'] = df_binary[col].apply(_safe_encode_binary_resistance)
    
    return df_binary


def add_derived_features(df: pd.DataFrame,
                        antibiotic_cols: List[str] = None) -> pd.DataFrame:
    """
    Add all derived features to the dataframe.
    
    Features added:
    - MAR_INDEX_COMPUTED: Multiple Antibiotic Resistance Index (0-1)
    - RESISTANCE_COUNT: Total number of resistant antibiotics
    - RESISTANT_CLASSES_COUNT: Number of antimicrobial categories with resistance
    - MDR_FLAG: Boolean Multi-Drug Resistant status
    - MDR_CATEGORY: Categorical "MDR" or "Non-MDR"
    - {AB}_RESISTANT: Binary resistance indicators for each antibiotic
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe with encoded resistance values
    antibiotic_cols : list, optional
        List of antibiotic column names (encoded). If None, auto-detected.
    
    Returns:
    --------
    pd.DataFrame
        Dataframe with all derived features added
    """
    print("=" * 70)
    print("PHASE 2.5: Feature Engineering (Formalized Definitions)")
    print("=" * 70)
    
    # Print formal definitions
    print("\nFORMAL FEATURE DEFINITIONS:")
    print("-" * 50)
    print("MAR Index (Multiple Antibiotic Resistance Index):")
    print("  Formula: MAR = a / b")
    print("  Where: a = resistant antibiotics, b = total tested")
    print("  Reference: Krumperman PH. (1983). Appl Environ Microbiol.")
    print("")
    print("MDR (Multi-Drug Resistant) Classification:")
    print("  Definition: Resistance to >=1 agent in >=3 antimicrobial categories")
    print("  Reference: Magiorakos AP, et al. (2012). Clin Microbiol Infect.")
    print("-" * 50)
    
    df_features = df.copy()
    
    # Auto-detect encoded columns if not provided
    if antibiotic_cols is None:
        antibiotic_cols = [c for c in df.columns if c.endswith('_encoded')]
    
    print(f"\n1. Using {len(antibiotic_cols)} antibiotic columns for feature engineering")
    
    # Compute MAR index
    print("2. Computing MAR index (Krumperman, 1983)...")
    df_features['MAR_INDEX_COMPUTED'] = df_features.apply(
        lambda row: compute_mar_index(row, antibiotic_cols),
        axis=1
    )
    
    # Compute resistance count
    print("3. Computing resistance count...")
    df_features['RESISTANCE_COUNT'] = df_features.apply(
        lambda row: compute_resistance_count(row, antibiotic_cols),
        axis=1
    )
    
    # Count resistant classes
    print("4. Counting resistant antibiotic classes...")
    df_features['RESISTANT_CLASSES_COUNT'] = df_features.apply(
        lambda row: count_resistant_classes(row, antibiotic_cols),
        axis=1
    )
    
    # Determine MDR status (Magiorakos et al., 2012) with SPECIES-SPECIFIC support
    print("5. Determining MDR status (Magiorakos et al., 2012)...")
    
    # Detect species column (ISOLATE_ID contains species information)
    species_col = 'ISOLATE_ID' if 'ISOLATE_ID' in df_features.columns else None
    if species_col:
        print(f"   → Using species-specific MDR classification (column: {species_col})")
        df_features['MDR_FLAG'] = df_features.apply(
            lambda row: determine_mdr_status(
                row, antibiotic_cols, 
                species=row.get(species_col, None)
            ),
            axis=1
        )
    else:
        print("   → Using universal MDR classification (no species column detected)")
        df_features['MDR_FLAG'] = df_features.apply(
            lambda row: determine_mdr_status(row, antibiotic_cols),
            axis=1
        )
    df_features['MDR_CATEGORY'] = df_features['MDR_FLAG'].map({True: 'MDR', False: 'Non-MDR'})
    
    # Add binary resistance indicators using the helper function
    print("6. Creating binary resistance indicators (R=1, S/I=0)...")
    df_features = create_binary_resistance_features(df_features, antibiotic_cols)
    
    # Summary statistics
    mdr_count = df_features['MDR_FLAG'].sum()
    mdr_pct = (mdr_count / len(df_features)) * 100
    
    print(f"\n7. FEATURE ENGINEERING SUMMARY:")
    print(f"   Total isolates: {len(df_features)}")
    print(f"   MDR isolates: {mdr_count} ({mdr_pct:.1f}%)")
    print(f"   Non-MDR isolates: {len(df_features) - mdr_count} ({100-mdr_pct:.1f}%)")
    print(f"   Mean MAR index: {df_features['MAR_INDEX_COMPUTED'].mean():.4f}")
    print(f"   Mean resistance count: {df_features['RESISTANCE_COUNT'].mean():.2f}")
    print(f"   Mean resistant classes: {df_features['RESISTANT_CLASSES_COUNT'].mean():.2f}")
    
    return df_features


def prepare_analysis_ready_dataset(df: pd.DataFrame,
                                   antibiotic_cols: List[str] = None,
                                   output_path: str = None) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, Dict]:
    """
    Create the final analysis-ready dataset with structural data separation.
    
    This function physically separates:
    - Feature matrix (X): Encoded resistance values for clustering/modeling
    - Metadata (meta): Sample identification and derived features
    
    This separation improves pipeline clarity and downstream modeling safety
    by preventing accidental use of metadata as model features.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Encoded dataframe
    antibiotic_cols : list, optional
        List of antibiotic columns (encoded). If None, auto-detected.
    output_path : str, optional
        Path to save the full analysis-ready dataset
    
    Returns:
    --------
    tuple
        (Full dataset, Feature matrix X, Metadata dataframe, Feature info dict)
        
    Output Files (when output_path provided):
        - analysis_ready_dataset.csv: Full combined dataset
        - feature_matrix_X.csv: Feature matrix only (for clustering/ML)
        - metadata.csv: Metadata only (for interpretation)
    """
    # Add derived features
    df_full = add_derived_features(df, antibiotic_cols)
    
    # Prepare feature matrix (X) - encoded resistance values only
    if antibiotic_cols is None:
        antibiotic_cols = [c for c in df_full.columns if c.endswith('_encoded')]
    
    feature_matrix = df_full[antibiotic_cols].copy()
    
    # Prepare metadata (meta) - everything except raw antibiotic columns and encoded columns
    metadata_cols = ['CODE', 'ISOLATE_ID', 'REGION', 'SITE', 'ENVIRONMENT',
                     'SAMPLING_SOURCE', 'NATIONAL_SITE', 'LOCAL_SITE', 
                     'REPLICATE', 'COLONY', 'ESBL', 'SOURCE_FILE', 
                     'resistance_fingerprint', 'MAR_INDEX_COMPUTED', 
                     'RESISTANCE_COUNT', 'RESISTANT_CLASSES_COUNT', 
                     'MDR_FLAG', 'MDR_CATEGORY']
    
    existing_metadata = [c for c in metadata_cols if c in df_full.columns]
    
    # Also include binary resistance features in metadata for easier interpretation
    binary_cols = [c for c in df_full.columns if c.endswith('_RESISTANT')]
    existing_metadata.extend(binary_cols)
    
    metadata = df_full[existing_metadata].copy()
    
    # Feature info
    feature_info = {
        'antibiotic_columns': antibiotic_cols,
        'total_antibiotics': len(antibiotic_cols),
        'total_isolates': len(df_full),
        'feature_matrix_shape': feature_matrix.shape,
        'metadata_shape': metadata.shape,
        'mdr_count': int(df_full['MDR_FLAG'].sum()),
        'mdr_percentage': float((df_full['MDR_FLAG'].sum() / len(df_full)) * 100),
        'mean_mar_index': float(df_full['MAR_INDEX_COMPUTED'].mean()),
        'mean_resistance_count': float(df_full['RESISTANCE_COUNT'].mean()),
        'references': {
            'mar_index': 'Krumperman PH. (1983). Appl Environ Microbiol, 46(1):165-170.',
            'mdr_definition': 'Magiorakos AP, et al. (2012). Clin Microbiol Infect, 18(3):268-281.'
        }
    }
    
    # Save outputs if output_path provided
    if output_path:
        import os
        output_dir = os.path.dirname(output_path)
        
        # Save full dataset
        df_full.to_csv(output_path, index=False)
        print(f"\n8. STRUCTURAL DATA SEPARATION:")
        print(f"   Full dataset saved to: {output_path}")
        
        # Save separated feature matrix and metadata
        feature_matrix_path = os.path.join(output_dir, 'feature_matrix_X.csv')
        metadata_path = os.path.join(output_dir, 'metadata.csv')
        
        feature_matrix.to_csv(feature_matrix_path, index=False)
        print(f"   Feature matrix (X) saved to: {feature_matrix_path}")
        print(f"     Shape: {feature_matrix.shape}")
        
        metadata.to_csv(metadata_path, index=False)
        print(f"   Metadata saved to: {metadata_path}")
        print(f"     Shape: {metadata.shape}")
    
    return df_full, feature_matrix, metadata, feature_info


if __name__ == "__main__":
    from pathlib import Path
    import os
    
    project_root = Path(__file__).parent.parent.parent
    encoded_path = project_root / "data" / "processed" / "encoded_dataset.csv"
    
    if encoded_path.exists():
        df = pd.read_csv(encoded_path)
        
        output_path = project_root / "data" / "processed" / "analysis_ready_dataset.csv"
        df_full, feature_matrix, metadata, info = prepare_analysis_ready_dataset(
            df, output_path=str(output_path)
        )
        
        print("\nFeature Info:")
        for key, value in info.items():
            if key != 'references':
                print(f"  {key}: {value}")
        
        print("\nReferences:")
        for ref_key, ref_value in info['references'].items():
            print(f"  {ref_key}: {ref_value}")
    else:
        print(f"Encoded dataset not found at {encoded_path}")
        print("Run resistance_encoding.py first.")
