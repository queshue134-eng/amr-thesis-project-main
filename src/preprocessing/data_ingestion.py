"""
Data Ingestion and Consolidation Module for AMR Thesis Project
Phase 2.1 - Load, merge, and add metadata to CSV files from different regions and sites

This module enforces explicit metadata standardization to ensure traceability,
reproducibility, and stratified analysis capability.

Required metadata columns at ingestion:
- REGION: Geographic region (e.g., BARMM, Region VIII)
- SITE: Specific sampling site within region
- ENVIRONMENT: Environmental category (Water, Fish, Hospital)
- SAMPLING_SOURCE: Detailed sampling source (e.g., Drinking Water, Fish Tilapia)
"""

import pandas as pd
import numpy as np
import os
import re
import sys
from pathlib import Path
from typing import Tuple, Dict, List, Optional
import warnings

# Console utilities
sys.path.insert(0, str(Path(__file__).parent.parent))
try:
    from utils.console import console, Colors
except ImportError:
    class FallbackConsole:
        def header(self, t, s=None): print(f"\n{'='*50}\n{t}\n{'='*50}")
        def info(self, m): print(f"  → {m}")
        def success(self, m): print(f"  ✓ {m}")
        def warning(self, m): print(f"  ⚠ {m}")
        def kv(self, k, v, i=2): print(f"{'  '*i}{k}: {v}")
    console = FallbackConsole()


# Standard antibiotic abbreviations we expect
STANDARD_ANTIBIOTICS = [
    'AM', 'AMC', 'CPT', 'CN', 'CF', 'CPD', 'CTX', 'CFO', 'CFT', 'CZA',
    'IPM', 'AN', 'GM', 'N', 'NAL', 'ENR', 'MRB', 'PRA', 'DO', 'TE',
    'FT', 'C', 'SXT',
    # Alternate names
    'AMP', 'AMO', 'CFX', 'CFP', 'CFA', 'CFV', 'CTF', 'CFZ', 'IME',
    'AMI', 'GEN', 'NEO', 'NLA', 'MAR', 'DOX', 'TET', 'NIT', 'CHL'
]

# Required metadata columns for valid ingestion
REQUIRED_METADATA_COLUMNS = ['REGION', 'SITE', 'ENVIRONMENT', 'SAMPLING_SOURCE']

# Environment categorization mapping (sampling source -> environment type)
ENVIRONMENT_MAPPING = {
    'Drinking Water': 'Water',
    'Lake Water': 'Water',
    'River Water': 'Water',
    'Fish Banak': 'Fish',
    'Fish Gusaw': 'Fish',
    'Fish Tilapia': 'Fish',
    'Fish Kaolang': 'Fish',
    'Effluent Water Untreated': 'Hospital',
    'Effluent Water Treated': 'Hospital',
}

# =============================================================================
# SPECIES STANDARDIZATION MAPPING
# =============================================================================
# Purpose: Consolidate subspecies and standardize species names for robust analysis
# 
# Rationale:
# 1. Klebsiella pneumoniae subspecies (ssp/spp ozaenae) → K. pneumoniae
#    - Subspecies not relevant for general AMR analysis
#    - "spp" is a typo for "ssp" (subspecies)
# 
# 2. Salmonella subspecies → Salmonella species
#    - S. enterica ssp diarizonae is a subspecies
#    - Consolidation increases sample size for Salmonella category
#
# 3. Hybrid/ambiguous entries → EXCLUDE
#    - "Enterobacter cloacae complex/E COLI" suggests mixed culture
#    - Including in either group introduces contamination bias
#
# 4. Vibrio vulnificus → KEEP (single sample, but valid waterborne pathogen)
# =============================================================================

SPECIES_STANDARDIZATION = {
    # Klebsiella pneumoniae consolidation
    'Klebsiella pneumoniae ssp ozaenae': 'Klebsiella pneumoniae',
    'Klebsiella pneumoniae spp ozaenae': 'Klebsiella pneumoniae',  # typo: spp -> ssp
    
    # Salmonella consolidation
    'Salmonella enterica spp diarizonae': 'Salmonella species',
    'Salmonella enterica ssp diarizonae': 'Salmonella species',
    'Salmonella group': 'Salmonella species',
    
    # Keep these as-is (explicit mapping for clarity)
    'Escherichia coli': 'Escherichia coli',
    'Klebsiella pneumoniae': 'Klebsiella pneumoniae',
    'Enterobacter cloacae': 'Enterobacter cloacae',
    'Enterobacter aerogenes': 'Enterobacter aerogenes',
    'Vibrio vulnificus': 'Vibrio vulnificus',
}

# Entries to EXCLUDE from analysis (ambiguous/contaminated samples)
SPECIES_EXCLUSIONS = [
    'Enterobacter cloacae complex/E COLI',  # Mixed culture - laboratory ambiguity
]


def parse_isolate_code(code: str, strict: bool = False) -> Dict[str, str]:
    """
    Parse isolate code to extract metadata based on naming convention.
    
    Isolate code format: [Species Prefix]_[National Site][Local Site][Sample Source][Replicate][Colony]
    Example: EC_OADWR1C3
    - EC = Escherichia coli (Species Prefix)
    - O = Ormoc (National Site)
    - A = Alegria (Local Site)  
    - DW = Drinking Water (Sample Source)
    - R1 = Replicate 1
    - C3 = Colony 3
    
    Parameters
    ----------
    code : str
        Isolate code string to parse
    strict : bool, optional
        If True, raise ValueError on malformed or incomplete codes.
        If False (default), return dict with None values for unparsed fields.
        
        RECOMMENDATION: Use strict=True during initial data ingestion to catch
        data quality issues early. Use strict=False for exploratory analysis
        where some malformed codes are acceptable.
    
    Returns
    -------
    dict
        Dictionary containing parsed metadata fields:
        - national_site: Mapped national site name
        - local_site: Mapped local site name
        - sample_source: Detailed sampling source (standardized name)
        - environment: Environment category (Water, Fish, Hospital)
        - replicate: Replicate number
        - colony: Colony number
    
    Raises
    ------
    ValueError
        If strict=True and code cannot be fully parsed
    """
    metadata = {
        'national_site': None,
        'local_site': None,
        'sample_source': None,
        'environment': None,
        'replicate': None,
        'colony': None
    }
    
    if not code or not isinstance(code, str):
        if strict:
            raise ValueError(
                f"Invalid isolate code: '{code}'. "
                f"Code must be a non-empty string. "
                f"Expected format: [SpeciesPrefix]_[NationalSite][LocalSite][SampleSource]R[Rep]C[Col]"
            )
        return metadata
    
    # Remove prefix like EC_, VC_, SAL if present
    code_clean = re.sub(r'^[A-Z]+_', '', code.strip())
    
    # National site mapping
    national_site_map = {
        'O': 'Ormoc',
        'P': 'Pampanga',
        'M': 'Marawi'
    }
    
    # Local site mapping (second letter) - base mapping
    # Note: For Marawi (BARMM), 'A' maps to 'APMC', not 'Alegria'
    local_site_map_default = {
        'A': 'Alegria',
        'L': 'Larrazabal',
        'G': 'Gabriel',
        'R': 'Roque',
        'D': 'Dayawan',
        'T': 'Tuca Kialdan',
        'P': 'APMC',
        'K': 'K'  # Marawi K site
    }
    
    # Region-specific local site overrides (Marawi/BARMM)
    local_site_map_marawi = {
        'A': 'APMC',       # BARMM 'A' = APMC, not Alegria
        'E': 'E',          # Marawi sites
        'R': 'Roque',
        'D': 'Dayawan',
        'T': 'Tuca Kialdan',
        'K': 'K'
    }
    
    # Sample source mapping
    sample_source_map = {
        'DW': 'Drinking Water',
        'LW': 'Lake Water',
        'FB': 'Fish Banak',
        'FG': 'Fish Gusaw',
        'RW': 'River Water',
        'FT': 'Fish Tilapia',
        'EWU': 'Effluent Water Untreated',
        'EWT': 'Effluent Water Treated',
        'FK': 'Fish Kaolang'
    }
    
    # Parse national site (first character)
    national_char = None
    if len(code_clean) > 0:
        national_char = code_clean[0].upper()
        metadata['national_site'] = national_site_map.get(national_char, national_char)
    
    # Parse local site (second character) - use region-aware mapping
    if len(code_clean) > 1:
        local_char = code_clean[1].upper()
        # Use Marawi-specific mapping if national site is 'M' (Marawi/BARMM)
        if national_char == 'M':
            metadata['local_site'] = local_site_map_marawi.get(local_char, local_char)
        else:
            metadata['local_site'] = local_site_map_default.get(local_char, local_char)
    
    # Parse sample source (two letters after local site, before replicate R or colony C)
    # Pattern: [National][Local][SampleSource][Replicate][Colony]
    # Example: OLDWR3C1 -> OL = national+local, DW = sample source, R3 = replicate, C1 = colony
    # The sample source can be 2-3 letters (DW, LW, FB, FG, RW, FT, EWU, EWT, FK)
    sample_match = re.search(r'^[A-Z]{2}([A-Z]{2,3})(?:R\d|C\d)', code_clean.upper())
    if sample_match:
        sample_code = sample_match.group(1)
        metadata['sample_source'] = sample_source_map.get(sample_code, sample_code)
        # Derive environment from sample source
        metadata['environment'] = ENVIRONMENT_MAPPING.get(metadata['sample_source'], 'Unknown')
    else:
        # Fallback: try to find known sample source codes anywhere in the string
        for code in sample_source_map.keys():
            if code in code_clean.upper():
                metadata['sample_source'] = sample_source_map[code]
                metadata['environment'] = ENVIRONMENT_MAPPING.get(metadata['sample_source'], 'Unknown')
                break
    
    # Parse replicate number (R followed by digit)
    replicate_match = re.search(r'R(\d)', code_clean.upper())
    if replicate_match:
        metadata['replicate'] = int(replicate_match.group(1))
    
    # Parse colony number (C followed by digits)
    colony_match = re.search(r'C(\d+)', code_clean.upper())
    if colony_match:
        metadata['colony'] = int(colony_match.group(1))
    
    # =========================================================================
    # STRICT MODE VALIDATION
    # If strict=True, validate that critical fields were parsed successfully
    # =========================================================================
    if strict:
        missing_fields = []
        critical_fields = ['national_site', 'sample_source', 'environment']
        
        for field in critical_fields:
            if metadata[field] is None:
                missing_fields.append(field)
        
        # Check for 'Unknown' environment (indicates unmapped sample source)
        if metadata['environment'] == 'Unknown' and metadata['sample_source'] is not None:
            missing_fields.append(f"environment (sample_source '{metadata['sample_source']}' is not mapped)")
        
        if missing_fields:
            raise ValueError(
                f"Failed to fully parse isolate code '{code}'. "
                f"Missing/invalid fields: {missing_fields}. "
                f"Parsed so far: {metadata}. "
                f"Expected format: [SpeciesPrefix]_[NationalSite][LocalSite][SampleSource]R[Rep]C[Col]"
            )
    
    return metadata


def extract_region_from_filename(filename: str) -> Tuple[str, str]:
    """
    Extract region and site information from filename.
    """
    region = None
    site = None
    
    # Extract region
    if 'BARMM' in filename:
        region = 'BARMM'
    elif 'Region VIII' in filename:
        region = 'Region VIII - Eastern Visayas'
    elif 'Region III' in filename:
        region = 'Region III - Central Luzon'
    else:
        region_match = re.search(r'Region\s+([IVX]+(?:-[A-Za-z\s]+)?)', filename)
        if region_match:
            region = f'Region {region_match.group(1)}'
    
    # Extract site (after LOR-)
    site_match = re.search(r'LOR-([A-Z\s]+)\.csv', filename, re.IGNORECASE)
    if site_match:
        site = site_match.group(1).strip()
    
    return region, site


def process_csv_file(filepath: str) -> pd.DataFrame:
    """
    Process a single CSV file and extract structured data.
    
    CSV Structure:
    - Row 3: CODE, ISOLATE ID headers + summary columns at end
    - Row 4: ESBL + antibiotic names (AM, AMC, etc.)
    - Row 5: MIC/INT labels
    - Row 6+: Data rows
    """
    filename = os.path.basename(filepath)
    region, site = extract_region_from_filename(filename)
    
    try:
        # Read CSV with all columns
        df = pd.read_csv(filepath, header=None)
        
        # Find the row with CODE header (usually row 3)
        code_row = None
        for idx in range(min(10, len(df))):
            row_str = df.iloc[idx].astype(str).str.upper()
            if 'CODE' in row_str.values:
                code_row = idx
                break
        
        if code_row is None:
            print(f"Warning: Could not find CODE header in {filename}")
            return pd.DataFrame()
        
        # The antibiotic row is the next row (row 4 in typical files)
        antibiotic_row = code_row + 1
        data_start = code_row + 3  # Skip header, antibiotic names, and MIC/INT labels
        
        # Get column indices from code row
        code_header = df.iloc[code_row]
        code_col_idx = None
        isolate_col_idx = None
        scored_res_idx = None
        num_ab_idx = None
        mar_idx = None
        
        for idx, val in enumerate(code_header):
            val_str = str(val).strip().upper()
            if val_str == 'CODE':
                code_col_idx = idx
            elif val_str == 'ISOLATE ID':
                isolate_col_idx = idx
            elif 'SCORED' in val_str and 'RESISTANCE' in val_str:
                scored_res_idx = idx
            elif 'NO.' in val_str and 'ANTIBIOTIC' in val_str:
                num_ab_idx = idx
            elif 'MAR' in val_str and 'INDEX' in val_str:
                mar_idx = idx
        
        # Get antibiotic names from antibiotic row
        ab_header = df.iloc[antibiotic_row]
        antibiotics = []  # List of (column_idx, antibiotic_name)
        esbl_col_idx = None
        
        for idx, val in enumerate(ab_header):
            val_str = str(val).strip().upper()
            if val_str == 'ESBL':
                esbl_col_idx = idx
            elif val_str in STANDARD_ANTIBIOTICS:
                antibiotics.append((idx, val_str))
        
        # Process data rows
        processed_rows = []
        
        for row_idx in range(data_start, len(df)):
            row = df.iloc[row_idx]
            row_data = {}
            
            # Get CODE
            if code_col_idx is not None:
                row_data['CODE'] = row.iloc[code_col_idx]
            
            # Skip rows without valid CODE
            if pd.isna(row_data.get('CODE')) or str(row_data.get('CODE')).strip() == '':
                continue
            
            # Get ISOLATE_ID (species)
            if isolate_col_idx is not None:
                row_data['ISOLATE_ID'] = row.iloc[isolate_col_idx]
            
            # Get ESBL
            if esbl_col_idx is not None:
                row_data['ESBL'] = row.iloc[esbl_col_idx]
            
            # Get antibiotic INT values (interpretation: S, I, R)
            # For each antibiotic at column idx, the INT value is at idx+1
            for ab_col_idx, ab_name in antibiotics:
                int_idx = ab_col_idx + 1  # INT column is right after the antibiotic name/MIC column
                if int_idx < len(row):
                    value = row.iloc[int_idx]
                    # Only store if it's a valid resistance value
                    if pd.notna(value):
                        val_str = str(value).strip().upper()
                        if val_str in ['S', 'I', 'R', '*R']:
                            row_data[ab_name] = val_str.replace('*', '')
            
            # Get summary columns
            if scored_res_idx is not None:
                row_data['SCORED_RESISTANCE'] = row.iloc[scored_res_idx]
            if num_ab_idx is not None:
                row_data['NUM_ANTIBIOTICS_TESTED'] = row.iloc[num_ab_idx]
            if mar_idx is not None:
                row_data['MAR_INDEX'] = row.iloc[mar_idx]
            
            # Add metadata from filename
            row_data['REGION'] = region
            row_data['SITE'] = site
            row_data['SOURCE_FILE'] = filename
            
            # Parse isolate code for additional metadata
            code_metadata = parse_isolate_code(str(row_data.get('CODE', '')))
            row_data.update({
                'NATIONAL_SITE': code_metadata['national_site'],
                'LOCAL_SITE': code_metadata['local_site'],
                'SAMPLING_SOURCE': code_metadata['sample_source'],  # Explicit standardized name
                'ENVIRONMENT': code_metadata['environment'],  # Environment category
                'REPLICATE': code_metadata['replicate'],
                'COLONY': code_metadata['colony']
            })
            
            processed_rows.append(row_data)
        
        return pd.DataFrame(processed_rows)
        
    except Exception as e:
        print(f"Error processing {filename}: {e}")
        import traceback
        traceback.print_exc()
        return pd.DataFrame()


def load_all_csv_files(data_dir: str) -> pd.DataFrame:
    """
    Load and consolidate all CSV files from the data directory.
    
    Parameters:
    -----------
    data_dir : str
        Path to directory containing CSV files
    
    Returns:
    --------
    pd.DataFrame
        Unified raw dataset with metadata columns
    """
    all_data = []
    csv_files = list(Path(data_dir).glob('*.csv'))
    
    print(f"Found {len(csv_files)} CSV files")
    
    for csv_file in csv_files:
        print(f"Processing: {csv_file.name}")
        df = process_csv_file(str(csv_file))
        if not df.empty:
            all_data.append(df)
            print(f"  -> Loaded {len(df)} isolates")
    
    if not all_data:
        print("No data loaded!")
        return pd.DataFrame()
    
    # Concatenate all dataframes
    master_df = pd.concat(all_data, ignore_index=True)
    
    print(f"\nTotal isolates loaded: {len(master_df)}")
    print(f"Columns: {list(master_df.columns)}")
    
    return master_df


def validate_required_metadata(df: pd.DataFrame) -> Dict[str, any]:
    """
    Validate that required metadata columns exist and report coverage.
    
    This ensures traceability, reproducibility, and stratified analysis capability.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe
    
    Returns:
    --------
    dict
        Validation report with coverage statistics
    """
    validation_report = {
        'columns_present': [],
        'columns_missing': [],
        'coverage': {},
        'warnings': [],
        'is_valid': True
    }
    
    for col in REQUIRED_METADATA_COLUMNS:
        if col in df.columns:
            validation_report['columns_present'].append(col)
            # Calculate coverage (non-null percentage)
            non_null_count = df[col].notna().sum()
            coverage = (non_null_count / len(df)) * 100 if len(df) > 0 else 0
            validation_report['coverage'][col] = coverage
            
            if coverage < 80:
                validation_report['warnings'].append(
                    f"Warning: {col} has only {coverage:.1f}% coverage (below 80% threshold)"
                )
        else:
            validation_report['columns_missing'].append(col)
            validation_report['is_valid'] = False
    
    return validation_report


def standardize_species_names(df: pd.DataFrame, species_column: str = 'ISOLATE_ID') -> pd.DataFrame:
    """
    Standardize species names by consolidating subspecies and excluding ambiguous entries.
    
    This function addresses data quality issues in species identification:
    1. Consolidates Klebsiella pneumoniae subspecies (ssp/spp ozaenae) → K. pneumoniae
    2. Consolidates Salmonella subspecies → Salmonella species
    3. Excludes ambiguous entries (e.g., "Enterobacter cloacae complex/E COLI")
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe with species column
    species_column : str
        Name of column containing species names (default: 'ISOLATE_ID')
    
    Returns:
    --------
    pd.DataFrame
        Dataframe with standardized species names and excluded entries removed
    
    Notes:
    ------
    - Excluded entries are logged but not saved to output
    - Original species names are preserved in 'ORIGINAL_SPECIES' column
    """
    if species_column not in df.columns:
        warnings.warn(f"Species column '{species_column}' not found. Skipping standardization.")
        return df
    
    df = df.copy()
    
    # Store original species names for traceability
    df['ORIGINAL_SPECIES'] = df[species_column].copy()
    
    # Count entries before exclusion
    initial_count = len(df)
    
    # Apply exclusions first
    excluded_mask = df[species_column].isin(SPECIES_EXCLUSIONS)
    excluded_count = excluded_mask.sum()
    
    if excluded_count > 0:
        excluded_entries = df.loc[excluded_mask, [species_column, 'CODE']].to_dict('records')
        print(f"\n  Species exclusions applied ({excluded_count} rows removed):")
        for entry in excluded_entries:
            print(f"    - '{entry[species_column]}' (CODE: {entry.get('CODE', 'N/A')})")
        df = df[~excluded_mask].copy()
    
    # Apply standardization mapping
    standardized_count = 0
    for original, standardized in SPECIES_STANDARDIZATION.items():
        mask = df[species_column] == original
        count = mask.sum()
        if count > 0 and original != standardized:
            df.loc[mask, species_column] = standardized
            standardized_count += count
            print(f"    Mapped: '{original}' -> '{standardized}' ({count} rows)")
    
    # Report summary
    final_count = len(df)
    final_species = df[species_column].nunique()
    print(f"\n  Species standardization complete:")
    print(f"    - Initial rows: {initial_count}")
    print(f"    - Excluded rows: {excluded_count}")
    print(f"    - Final rows: {final_count}")
    print(f"    - Standardized mappings applied: {standardized_count}")
    print(f"    - Unique species (after): {final_species}")
    
    return df


def create_unified_dataset(input_dir: str, output_path: Optional[str] = None) -> pd.DataFrame:
    """
    Main function to create unified raw dataset with explicit metadata standardization.
    
    Ensures traceability, reproducibility, and stratified analysis capability
    by enforcing required metadata columns at ingestion.
    
    Parameters:
    -----------
    input_dir : str
        Directory containing raw CSV files
    output_path : str, optional
        Path to save the unified dataset
    
    Returns:
    --------
    pd.DataFrame
        Unified raw dataset with validated metadata columns
    """
    print("=" * 50)
    print("PHASE 2.1: Data Ingestion and Consolidation")
    print("=" * 50)
    
    # Load all CSV files
    master_df = load_all_csv_files(input_dir)
    
    if master_df.empty:
        return master_df
    
    # Validate required metadata columns
    print("\nValidating required metadata columns...")
    validation_report = validate_required_metadata(master_df)
    
    if validation_report['columns_present']:
        print(f"  Present: {', '.join(validation_report['columns_present'])}")
    if validation_report['columns_missing']:
        warnings.warn(f"Missing required metadata columns: {validation_report['columns_missing']}")
        print(f"  Missing: {', '.join(validation_report['columns_missing'])}")
    
    print("\nMetadata coverage:")
    for col, coverage in validation_report['coverage'].items():
        print(f"  {col}: {coverage:.1f}%")
    
    for warning in validation_report['warnings']:
        print(f"  {warning}")
    
    # Apply species standardization
    print("\nApplying species standardization...")
    master_df = standardize_species_names(master_df, species_column='ISOLATE_ID')
    
    # Reorder columns for clarity (with new metadata columns)
    metadata_cols = ['CODE', 'ISOLATE_ID', 'ORIGINAL_SPECIES', 'REGION', 'SITE', 'ENVIRONMENT', 
                     'SAMPLING_SOURCE', 'NATIONAL_SITE', 'LOCAL_SITE', 
                     'REPLICATE', 'COLONY', 'ESBL', 'SOURCE_FILE']
    
    summary_cols = ['SCORED_RESISTANCE', 'NUM_ANTIBIOTICS_TESTED', 'MAR_INDEX']
    
    # Get antibiotic columns (everything else)
    all_cols = set(master_df.columns)
    existing_metadata = [c for c in metadata_cols if c in all_cols]
    existing_summary = [c for c in summary_cols if c in all_cols]
    antibiotic_cols = sorted([c for c in all_cols if c not in metadata_cols + summary_cols])
    
    # Reorder
    new_order = existing_metadata + antibiotic_cols + existing_summary
    master_df = master_df[[c for c in new_order if c in master_df.columns]]
    
    # Save if output path provided
    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        master_df.to_csv(output_path, index=False)
        print(f"\nUnified dataset saved to: {output_path}")
    
    return master_df


if __name__ == "__main__":
    # Example usage
    project_root = Path(__file__).parent.parent.parent
    input_dir = project_root
    output_path = project_root / "data" / "processed" / "unified_raw_dataset.csv"
    
    df = create_unified_dataset(str(input_dir), str(output_path))
    if not df.empty:
        print("\nDataset Summary:")
        print(f"Shape: {df.shape}")
        print(f"Regions: {df['REGION'].unique()}")
        print(f"Sites: {df['SITE'].unique()}")
        if 'ENVIRONMENT' in df.columns:
            print(f"Environments: {df['ENVIRONMENT'].unique()}")
        if 'SAMPLING_SOURCE' in df.columns:
            print(f"Sampling Sources: {df['SAMPLING_SOURCE'].unique()}")
