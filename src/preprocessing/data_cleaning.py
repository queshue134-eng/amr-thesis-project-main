"""
Data Cleaning Module for AMR Thesis Project
Phase 2.2 and 2.3 - Data cleaning and handling missing data

This module implements:
1. Validation rules - Only {S, I, R} allowed, no multi-label entries
2. Duplicate isolate detection and resolution
3. Formal missing data strategy with transparent methodology
4. Controlled vocabularies for species and antibiotic names
5. Logging of all cleaning actions for reproducibility
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional
import re
from datetime import datetime


# Valid resistance values (CLSI standard)
VALID_RESISTANCE_VALUES = {'S', 'I', 'R'}

# Standardized species names mapping (Controlled Vocabulary)
SPECIES_STANDARDIZATION = {
    # Escherichia coli variants
    'escherichia coli': 'Escherichia coli',
    'e. coli': 'Escherichia coli',
    'e.coli': 'Escherichia coli',
    
    # Klebsiella pneumoniae variants
    'klebsiella pneumoniae ssp pneumoniae': 'Klebsiella pneumoniae',
    'klebsiella pneumoniae': 'Klebsiella pneumoniae',
    'k. pneumoniae': 'Klebsiella pneumoniae',
    
    # Enterobacter variants
    'enterobacter cloacae complex': 'Enterobacter cloacae',
    'enterobacter cloacae': 'Enterobacter cloacae',
    'enterobacter aerogenes': 'Enterobacter aerogenes',
    
    # Pseudomonas variants
    'pseudomonas aeruginosa': 'Pseudomonas aeruginosa',
    'p. aeruginosa': 'Pseudomonas aeruginosa',
    
    # Vibrio variants
    'vibrio fluvialis': 'Vibrio fluvialis',
    'v. fluvialis': 'Vibrio fluvialis',
}

# Standardized antibiotic names mapping (Controlled Vocabulary)
ANTIBIOTIC_STANDARDIZATION = {
    # Full names to abbreviations
    'ampicillin': 'AM',
    'amoxicillin-clavulanate': 'AMC',
    'cefepime': 'CPT',
    'cephalothin': 'CN',
    'cefazolin': 'CF',
    'cefpodoxime': 'CPD',
    'cefotaxime': 'CTX',
    'cefoxitin': 'CFO',
    'ceftriaxone': 'CFT',
    'ceftazidime-avibactam': 'CZA',
    'imipenem': 'IPM',
    'amikacin': 'AN',
    'gentamicin': 'GM',
    'neomycin': 'N',
    'nalidixic acid': 'NAL',
    'enrofloxacin': 'ENR',
    'meropenem': 'MRB',
    'piperacillin-tazobactam': 'PRA',
    'doxycycline': 'DO',
    'tetracycline': 'TE',
    'nitrofurantoin': 'FT',
    'chloramphenicol': 'C',
    'trimethoprim-sulfamethoxazole': 'SXT',
    
    # Alternate abbreviations
    'AMP': 'AM',
    'AMO': 'AMC',
    'CFX': 'CN',
    'CFP': 'CF',
    'CFA': 'CTX',
    'CFV': 'CFO',
    'CTF': 'CFT',
    'CFZ': 'CZA',
    'IME': 'IPM',
    'AMI': 'AN',
    'GEN': 'GM',
    'NEO': 'N',
    'NLA': 'NAL',
    'MAR': 'MRB',
    'DOX': 'DO',
    'TET': 'TE',
    'NIT': 'FT',
    'CHL': 'C',
}


def standardize_species_name(name: str) -> str:
    """
    Standardize bacterial species names.
    
    Parameters:
    -----------
    name : str
        Raw species name
    
    Returns:
    --------
    str
        Standardized species name
    """
    if pd.isna(name) or not isinstance(name, str):
        return np.nan
    
    name_lower = name.strip().lower()
    return SPECIES_STANDARDIZATION.get(name_lower, name.strip())


def standardize_antibiotic_name(name: str) -> str:
    """
    Standardize antibiotic names/abbreviations.
    
    Parameters:
    -----------
    name : str
        Raw antibiotic name
    
    Returns:
    --------
    str
        Standardized antibiotic abbreviation
    """
    if pd.isna(name) or not isinstance(name, str):
        return name
    
    name_clean = name.strip().upper()
    return ANTIBIOTIC_STANDARDIZATION.get(name_clean, name_clean)


def validate_resistance_value(value) -> Tuple[bool, Optional[str], str]:
    """
    Validate a single resistance value according to CLSI standards.
    
    Validation rules:
    - Only {S, I, R} allowed
    - No multi-label entries (e.g., "S/R" is invalid)
    
    Parameters:
    -----------
    value : any
        Resistance value to validate
    
    Returns:
    --------
    tuple
        (is_valid: bool, standardized_value: str or None, issue: str)
    """
    if pd.isna(value) or value is None:
        return (True, None, '')
    
    value_str = str(value).strip().upper()
    
    # Check for multi-label entries (invalid)
    if '/' in value_str or ',' in value_str or ';' in value_str:
        return (False, None, f'Multi-label entry: {value}')
    
    # Standardize and validate
    if value_str in ['S', 'SUSCEPTIBLE']:
        return (True, 'S', '')
    elif value_str in ['I', 'INTERMEDIATE']:
        return (True, 'I', '')
    elif value_str in ['R', 'RESISTANT', '*R']:  # *R indicates borderline resistant
        return (True, 'R', '')
    elif value_str in ['', 'NAN', 'NONE', '-']:
        return (True, None, '')
    else:
        return (False, None, f'Invalid value: {value}')


def standardize_resistance_value(value) -> Optional[str]:
    """
    Standardize resistance interpretation values to S, I, or R.
    
    Parameters:
    -----------
    value : any
        Raw resistance value
    
    Returns:
    --------
    str or None
        Standardized value (S, I, R) or None for missing/invalid
    """
    is_valid, std_value, _ = validate_resistance_value(value)
    return std_value


def validate_resistance_data(df: pd.DataFrame, antibiotic_cols: List[str]) -> Dict[str, List]:
    """
    Validate all resistance data in the dataframe.
    
    Detects:
    - Multi-label entries (e.g., "S/R")
    - Invalid values (not S, I, R, or empty)
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe
    antibiotic_cols : list
        List of antibiotic column names
    
    Returns:
    --------
    dict
        Dictionary containing validation issues:
        - invalid_values: list of (index, column, value, issue)
        - total_invalid_count: int
    """
    validation_results = {
        'invalid_values': [],
        'total_invalid_count': 0,
        'validation_summary': {}
    }
    
    for col in antibiotic_cols:
        if col not in df.columns:
            continue
            
        col_invalid = 0
        for idx, value in df[col].items():
            is_valid, _, issue = validate_resistance_value(value)
            if not is_valid:
                validation_results['invalid_values'].append({
                    'index': idx,
                    'column': col,
                    'value': value,
                    'issue': issue
                })
                col_invalid += 1
        
        validation_results['validation_summary'][col] = {
            'invalid_count': col_invalid
        }
        validation_results['total_invalid_count'] += col_invalid
    
    return validation_results


def clean_resistance_data(df: pd.DataFrame, antibiotic_cols: List[str]) -> Tuple[pd.DataFrame, Dict]:
    """
    Clean and standardize resistance data in antibiotic columns.
    
    Validates values before standardization and logs all cleaning actions.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe
    antibiotic_cols : list
        List of antibiotic column names
    
    Returns:
    --------
    tuple
        (Cleaned dataframe, cleaning log dictionary)
    """
    df_clean = df.copy()
    cleaning_log = {
        'values_standardized': 0,
        'invalid_values_nullified': 0,
        'actions': []
    }
    
    for col in antibiotic_cols:
        if col in df_clean.columns:
            for idx in df_clean.index:
                original_value = df_clean.at[idx, col]
                is_valid, std_value, issue = validate_resistance_value(original_value)
                
                if pd.notna(original_value):
                    if not is_valid:
                        cleaning_log['invalid_values_nullified'] += 1
                        cleaning_log['actions'].append({
                            'type': 'invalid_nullified',
                            'index': idx,
                            'column': col,
                            'original': original_value,
                            'new': None,
                            'reason': issue
                        })
                    elif std_value != str(original_value).strip().upper() if original_value else True:
                        cleaning_log['values_standardized'] += 1
                
                df_clean.at[idx, col] = std_value
    
    return df_clean, cleaning_log


def remove_duplicate_isolates(df: pd.DataFrame, key_cols: List[str] = None) -> Tuple[pd.DataFrame, int, List[Dict]]:
    """
    Remove duplicate isolates based on key columns.
    
    Detects and removes duplicate isolates, logging all removals for transparency.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe
    key_cols : list, optional
        Columns to use for identifying duplicates (default: ['CODE'])
    
    Returns:
    --------
    tuple
        (Cleaned dataframe, number of duplicates removed, list of removed duplicates info)
    """
    if key_cols is None:
        key_cols = ['CODE']
    
    initial_count = len(df)
    
    # Identify duplicates before removing
    duplicates_mask = df.duplicated(subset=key_cols, keep='first')
    duplicates_info = []
    
    for idx in df[duplicates_mask].index:
        duplicates_info.append({
            'index': idx,
            'code': df.at[idx, 'CODE'] if 'CODE' in df.columns else None,
            'reason': f'Duplicate based on {key_cols}'
        })
    
    df_clean = df.drop_duplicates(subset=key_cols, keep='first')
    removed_count = initial_count - len(df_clean)
    
    return df_clean, removed_count, duplicates_info


def check_inconsistent_values(df: pd.DataFrame, antibiotic_cols: List[str]) -> Dict[str, List]:
    """
    Check for impossible or inconsistent values in resistance data.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe
    antibiotic_cols : list
        List of antibiotic column names
    
    Returns:
    --------
    dict
        Dictionary of issues found
    """
    issues = {
        'mixed_values': [],  # Cells with multiple S/I/R values
        'invalid_values': [],  # Values that aren't S, I, R, or NaN
    }
    
    for idx, row in df.iterrows():
        for col in antibiotic_cols:
            if col in df.columns:
                value = row[col]
                if pd.notna(value):
                    value_str = str(value).strip().upper()
                    
                    # Check for mixed values (e.g., "S/R")
                    if '/' in value_str or ',' in value_str:
                        issues['mixed_values'].append({
                            'index': idx,
                            'column': col,
                            'value': value
                        })
                    
                    # Check for invalid values after standardization
                    std_value = standardize_resistance_value(value)
                    if std_value not in {'S', 'I', 'R', None}:
                        issues['invalid_values'].append({
                            'index': idx,
                            'column': col,
                            'value': value
                        })
    
    return issues


def compute_antibiotic_test_coverage(df: pd.DataFrame, antibiotic_cols: List[str]) -> Dict[str, Dict]:
    """
    Compute antibiotic test coverage (% of isolates tested for each antibiotic).
    
    This is a key component of the formal missing data strategy.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe
    antibiotic_cols : list
        List of antibiotic column names
    
    Returns:
    --------
    dict
        Dictionary with coverage statistics for each antibiotic:
        - tested_count: Number of isolates tested
        - missing_count: Number of isolates not tested
        - coverage_pct: Percentage of isolates tested
    """
    total_isolates = len(df)
    coverage_stats = {}
    
    for col in antibiotic_cols:
        if col in df.columns:
            # Count non-missing values (tested isolates)
            tested_count = df[col].notna().sum()
            # Also exclude empty strings
            if df[col].dtype == 'object':
                tested_count = tested_count - (df[col] == '').sum()
            
            missing_count = total_isolates - tested_count
            coverage_pct = (tested_count / total_isolates) * 100 if total_isolates > 0 else 0
            
            coverage_stats[col] = {
                'tested_count': int(tested_count),
                'missing_count': int(missing_count),
                'coverage_pct': coverage_pct
            }
    
    return coverage_stats


def analyze_missing_data(df: pd.DataFrame, antibiotic_cols: List[str]) -> Dict[str, float]:
    """
    Analyze missing data patterns in antibiotic columns.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe
    antibiotic_cols : list
        List of antibiotic column names
    
    Returns:
    --------
    dict
        Dictionary with missing percentage for each antibiotic
    """
    missing_stats = {}
    total_isolates = len(df)
    
    for col in antibiotic_cols:
        if col in df.columns:
            # Count missing (None, NaN, empty string)
            missing_count = df[col].isna().sum()
            missing_count += (df[col] == '').sum() if df[col].dtype == 'object' else 0
            missing_stats[col] = (missing_count / total_isolates) * 100
    
    return missing_stats


def filter_antibiotics_by_coverage(df: pd.DataFrame, 
                                   antibiotic_cols: List[str],
                                   min_coverage: float = 70.0) -> Tuple[List[str], List[str], Dict]:
    """
    Filter antibiotics to retain only those tested in sufficient isolates.
    
    Part of the formal missing data strategy - retains antibiotics tested
    in ≥X% of isolates (default 70%) to ensure robust pattern analysis.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe
    antibiotic_cols : list
        List of antibiotic column names
    min_coverage : float
        Minimum percentage of isolates that must have data (default: 70%)
    
    Returns:
    --------
    tuple
        (retained_antibiotics, excluded_antibiotics, coverage_details)
    """
    coverage_stats = compute_antibiotic_test_coverage(df, antibiotic_cols)
    
    retained_antibiotics = []
    excluded_antibiotics = []
    
    for col in antibiotic_cols:
        if col in coverage_stats:
            if coverage_stats[col]['coverage_pct'] >= min_coverage:
                retained_antibiotics.append(col)
            else:
                excluded_antibiotics.append(col)
    
    return retained_antibiotics, excluded_antibiotics, coverage_stats


def remove_isolates_with_excessive_missing(df: pd.DataFrame,
                                            antibiotic_cols: List[str],
                                            max_missing_pct: float = 30.0) -> Tuple[pd.DataFrame, int, List[Dict]]:
    """
    Remove isolates with excessive missing AST values.
    
    Part of the formal missing data strategy - removes isolates where
    more than X% of antibiotic test results are missing.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe
    antibiotic_cols : list
        List of antibiotic column names
    max_missing_pct : float
        Maximum percentage of missing values allowed (default: 30%)
    
    Returns:
    --------
    tuple
        (Cleaned dataframe, number of isolates removed, list of removed isolate info)
    """
    df_clean = df.copy()
    initial_count = len(df_clean)
    
    # Calculate missing percentage for each isolate
    existing_cols = [c for c in antibiotic_cols if c in df_clean.columns]
    
    if not existing_cols:
        return df_clean, 0, []
    
    missing_counts = df_clean[existing_cols].isna().sum(axis=1)
    missing_pcts = (missing_counts / len(existing_cols)) * 100
    
    # Track removed isolates
    removed_isolates_info = []
    for idx in df_clean[missing_pcts > max_missing_pct].index:
        removed_isolates_info.append({
            'index': idx,
            'code': df_clean.at[idx, 'CODE'] if 'CODE' in df_clean.columns else None,
            'missing_pct': missing_pcts[idx],
            'reason': f'Missing {missing_pcts[idx]:.1f}% of antibiotic tests (threshold: {max_missing_pct}%)'
        })
    
    # Filter out isolates with too much missing data
    df_clean = df_clean[missing_pcts <= max_missing_pct]
    removed_count = initial_count - len(df_clean)
    
    return df_clean, removed_count, removed_isolates_info


def clean_dataset(df: pd.DataFrame, 
                  min_antibiotic_coverage: float = 70.0,
                  max_isolate_missing: float = 30.0) -> Tuple[pd.DataFrame, Dict]:
    """
    Main data cleaning function implementing formal missing data strategy.
    
    This function provides transparent, defensible methodology for:
    1. Validation of resistance values (only S, I, R allowed)
    2. Detection and resolution of duplicate isolates
    3. Species and antibiotic name standardization using controlled vocabularies
    4. Formal missing data strategy with explicit thresholds
    5. Comprehensive logging of all cleaning actions
    
    Parameters:
    -----------
    df : pd.DataFrame
        Raw unified dataset
    min_antibiotic_coverage : float
        Minimum test coverage percentage for antibiotics (default: 70%)
        Antibiotics tested in fewer isolates are excluded
    max_isolate_missing : float
        Maximum missing percentage allowed for isolates (default: 30%)
        Isolates exceeding this threshold are excluded
    
    Returns:
    --------
    tuple
        (Cleaned dataframe, comprehensive cleaning report dictionary)
    """
    print("=" * 50)
    print("PHASE 2.2 & 2.3: Data Cleaning and Missing Data Handling")
    print("=" * 50)
    
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    report = {
        'timestamp': timestamp,
        'initial_isolates': len(df),
        'initial_columns': len(df.columns),
        'thresholds': {
            'min_antibiotic_coverage_pct': min_antibiotic_coverage,
            'max_isolate_missing_pct': max_isolate_missing
        },
        'duplicates_removed': 0,
        'duplicates_info': [],
        'isolates_removed_missing': 0,
        'isolates_removed_info': [],
        'antibiotics_retained': [],
        'antibiotics_excluded': [],
        'antibiotic_coverage': {},
        'validation_issues': {},
        'cleaning_log': {},
        'final_isolates': 0,
        'final_columns': 0,
        'missing_data_stats': {}
    }
    
    df_clean = df.copy()
    
    # Identify antibiotic columns (exclude metadata and summary columns)
    metadata_cols = ['CODE', 'ISOLATE_ID', 'REGION', 'SITE', 'ENVIRONMENT',
                     'SAMPLING_SOURCE', 'NATIONAL_SITE', 'LOCAL_SITE', 
                     'REPLICATE', 'COLONY', 'ESBL', 'SOURCE_FILE', 
                     'SCORED_RESISTANCE', 'NUM_ANTIBIOTICS_TESTED', 'MAR_INDEX']
    
    antibiotic_cols = [c for c in df_clean.columns if c not in metadata_cols]
    
    print(f"\n1. Initial dataset: {len(df_clean)} isolates, {len(antibiotic_cols)} antibiotic columns")
    
    # Step 1: Validate resistance data
    print("\n2. Validating resistance data...")
    validation_results = validate_resistance_data(df_clean, antibiotic_cols)
    report['validation_issues'] = validation_results
    
    if validation_results['total_invalid_count'] > 0:
        print(f"   Found {validation_results['total_invalid_count']} invalid values")
    else:
        print("   All values valid")
    
    # Step 2: Standardize species names using controlled vocabulary
    if 'ISOLATE_ID' in df_clean.columns:
        df_clean['ISOLATE_ID'] = df_clean['ISOLATE_ID'].apply(standardize_species_name)
        print("3. Standardized species names (controlled vocabulary)")
    
    # Step 3: Clean and standardize resistance values
    df_clean, cleaning_log = clean_resistance_data(df_clean, antibiotic_cols)
    report['cleaning_log'] = cleaning_log
    print(f"4. Standardized resistance values (S, I, R)")
    if cleaning_log['invalid_values_nullified'] > 0:
        print(f"   Nullified {cleaning_log['invalid_values_nullified']} invalid values")
    
    # Step 4: Remove duplicates
    df_clean, dup_removed, dup_info = remove_duplicate_isolates(df_clean)
    report['duplicates_removed'] = dup_removed
    report['duplicates_info'] = dup_info
    print(f"5. Removed {dup_removed} duplicate isolates")
    
    # Step 5: Compute antibiotic test coverage
    print(f"\n6. Computing antibiotic test coverage...")
    retained_antibiotics, excluded_antibiotics, coverage_stats = filter_antibiotics_by_coverage(
        df_clean, antibiotic_cols, min_antibiotic_coverage
    )
    
    report['antibiotic_coverage'] = coverage_stats
    report['antibiotics_retained'] = retained_antibiotics
    report['antibiotics_excluded'] = excluded_antibiotics
    
    print(f"   Threshold: >={min_antibiotic_coverage}% coverage")
    print(f"   Retained: {len(retained_antibiotics)} antibiotics")
    print(f"   Excluded: {len(excluded_antibiotics)} antibiotics")
    if excluded_antibiotics:
        for ab in excluded_antibiotics:
            cov = coverage_stats.get(ab, {}).get('coverage_pct', 0)
            print(f"     - {ab}: {cov:.1f}% coverage")
    
    # Step 6: Analyze missing data patterns
    missing_stats = analyze_missing_data(df_clean, retained_antibiotics)
    report['missing_data_stats'] = missing_stats
    
    # Step 7: Remove isolates with excessive missing data
    print(f"\n7. Removing isolates with excessive missing data...")
    print(f"   Threshold: >{max_isolate_missing}% missing")
    df_clean, iso_removed, iso_removed_info = remove_isolates_with_excessive_missing(
        df_clean, retained_antibiotics, max_isolate_missing
    )
    report['isolates_removed_missing'] = iso_removed
    report['isolates_removed_info'] = iso_removed_info
    print(f"   Removed: {iso_removed} isolates")
    
    # Final stats
    report['final_isolates'] = len(df_clean)
    report['final_columns'] = len(df_clean.columns)
    
    print(f"\n8. Final dataset: {len(df_clean)} isolates")
    print(f"   Data retention: {(len(df_clean) / report['initial_isolates'] * 100):.1f}%")
    
    return df_clean, report


def generate_cleaning_report(report: Dict, output_path: str = None) -> str:
    """
    Generate a comprehensive text report of the cleaning process.
    
    Produces a transparent, defensible documentation of all cleaning
    decisions and actions for reproducibility.
    
    Parameters:
    -----------
    report : dict
        Cleaning report dictionary
    output_path : str, optional
        Path to save the report
    
    Returns:
    --------
    str
        Report text
    """
    lines = [
        "=" * 70,
        "DATA CLEANING REPORT - FORMAL MISSING DATA STRATEGY",
        "=" * 70,
        "",
        f"Generated: {report.get('timestamp', 'N/A')}",
        "",
        "METHODOLOGY SUMMARY",
        "-" * 50,
        "This report documents all cleaning actions performed as part of a",
        "transparent, defensible missing data methodology.",
        "",
        "THRESHOLDS APPLIED",
        "-" * 50,
    ]
    
    thresholds = report.get('thresholds', {})
    lines.append(f"  Minimum antibiotic test coverage: {thresholds.get('min_antibiotic_coverage_pct', 'N/A')}%")
    lines.append(f"  Maximum isolate missing data: {thresholds.get('max_isolate_missing_pct', 'N/A')}%")
    
    lines.extend([
        "",
        "DATA RETENTION SUMMARY",
        "-" * 50,
        f"  Initial isolates: {report['initial_isolates']}",
        f"  Duplicates removed: {report['duplicates_removed']}",
        f"  Isolates removed (missing data): {report['isolates_removed_missing']}",
        f"  Final isolates: {report['final_isolates']}",
        f"  Retention rate: {(report['final_isolates'] / report['initial_isolates'] * 100):.1f}%",
        "",
        "VALIDATION SUMMARY",
        "-" * 50,
    ])
    
    validation = report.get('validation_issues', {})
    lines.append(f"  Total invalid values found: {validation.get('total_invalid_count', 0)}")
    
    lines.extend([
        "",
        "ANTIBIOTIC TEST COVERAGE",
        "-" * 50,
        f"  Retained ({len(report['antibiotics_retained'])}): {', '.join(report['antibiotics_retained'])}",
        f"  Excluded ({len(report['antibiotics_excluded'])}): {', '.join(report['antibiotics_excluded']) if report['antibiotics_excluded'] else 'None'}",
        "",
        "  Coverage by Antibiotic:",
    ])
    
    coverage = report.get('antibiotic_coverage', {})
    for ab in sorted(coverage.keys(), key=lambda x: coverage[x]['coverage_pct'], reverse=True):
        stats = coverage[ab]
        status = "RETAINED" if ab in report['antibiotics_retained'] else "EXCLUDED"
        lines.append(f"    {ab}: {stats['coverage_pct']:.1f}% ({stats['tested_count']}/{stats['tested_count'] + stats['missing_count']}) [{status}]")
    
    lines.extend([
        "",
        "EXCLUSION SUMMARY TABLE",
        "-" * 50,
    ])
    
    # Excluded antibiotics table
    if report['antibiotics_excluded']:
        lines.append("  Excluded Antibiotics (insufficient coverage):")
        for ab in report['antibiotics_excluded']:
            if ab in coverage:
                lines.append(f"    - {ab}: {coverage[ab]['coverage_pct']:.1f}% coverage")
    else:
        lines.append("  No antibiotics excluded.")
    
    lines.append("")
    
    # Excluded isolates summary
    if report.get('isolates_removed_info'):
        lines.append("  Excluded Isolates (excessive missing data):")
        lines.append(f"    Total: {len(report['isolates_removed_info'])} isolates")
        # Show first 10 as examples
        for i, info in enumerate(report['isolates_removed_info'][:10]):
            lines.append(f"    - {info.get('code', 'Unknown')}: {info.get('missing_pct', 0):.1f}% missing")
        if len(report['isolates_removed_info']) > 10:
            lines.append(f"    ... and {len(report['isolates_removed_info']) - 10} more")
    else:
        lines.append("  No isolates excluded for missing data.")
    
    lines.extend([
        "",
        "MISSING DATA BY ANTIBIOTIC (After Cleaning)",
        "-" * 50,
    ])
    
    for ab, missing_pct in sorted(report.get('missing_data_stats', {}).items(), key=lambda x: x[1]):
        lines.append(f"  {ab}: {missing_pct:.1f}% missing ({100-missing_pct:.1f}% complete)")
    
    lines.extend([
        "",
        "CLEANING ACTIONS LOG",
        "-" * 50,
    ])
    
    cleaning_log = report.get('cleaning_log', {})
    lines.append(f"  Values standardized: {cleaning_log.get('values_standardized', 0)}")
    lines.append(f"  Invalid values nullified: {cleaning_log.get('invalid_values_nullified', 0)}")
    
    lines.extend([
        "",
        "=" * 70,
        "END OF REPORT",
        "=" * 70,
    ])
    
    report_text = '\n'.join(lines)
    
    if output_path:
        with open(output_path, 'w') as f:
            f.write(report_text)
    
    return report_text


if __name__ == "__main__":
    import os
    from pathlib import Path
    from .data_ingestion import create_unified_dataset
    
    project_root = Path(__file__).parent.parent.parent
    
    # Load unified dataset
    input_dir = project_root
    unified_path = project_root / "data" / "processed" / "unified_raw_dataset.csv"
    
    if not unified_path.exists():
        df = create_unified_dataset(str(input_dir), str(unified_path))
    else:
        df = pd.read_csv(unified_path)
    
    # Clean dataset with formal missing data strategy (70% antibiotic coverage, 30% isolate threshold)
    df_clean, report = clean_dataset(df, min_antibiotic_coverage=70.0, max_isolate_missing=30.0)
    
    # Save cleaned dataset
    clean_path = project_root / "data" / "processed" / "cleaned_dataset.csv"
    df_clean.to_csv(clean_path, index=False)
    print(f"\nCleaned dataset saved to: {clean_path}")
    
    # Generate and save report
    report_path = project_root / "data" / "processed" / "cleaning_report.txt"
    report_text = generate_cleaning_report(report, str(report_path))
    print(f"Cleaning report saved to: {report_path}")
