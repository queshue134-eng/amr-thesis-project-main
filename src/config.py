"""
Centralized Configuration for AMR Thesis Project
=================================================

This module provides a single source of truth for all configurable parameters,
constants, and mappings used across the pipeline. Centralizing configuration:

1. Eliminates hard-coded values scattered across modules
2. Ensures consistency in parameter values
3. Makes the pipeline easily configurable for different datasets
4. Improves reproducibility by documenting all parameters in one place

Usage:
    from config import CONFIG, ANTIBIOTIC_CLASSES, RANDOM_STATE
    
    # Access nested config
    clustering_method = CONFIG['clustering']['linkage_method']
    
    # Use constants directly
    seed = RANDOM_STATE
"""

from pathlib import Path
from typing import Dict, List, Any


# =============================================================================
# PROJECT PATHS (pathlib for OS-agnostic compatibility)
# =============================================================================
PROJECT_ROOT = Path(__file__).parent.parent
DATA_DIR = PROJECT_ROOT / "data"
RAW_DATA_DIR = DATA_DIR / "raw"  # Raw CSVs moved to data/raw
PROCESSED_DATA_DIR = DATA_DIR / "processed"
FIGURES_DIR = PROCESSED_DATA_DIR / "figures"
MODELS_DIR = PROJECT_ROOT / "models"
ARTIFACTS_DIR = PROCESSED_DATA_DIR / "clustering_artifacts"

# =============================================================================
# REPRODUCIBILITY
# =============================================================================
RANDOM_STATE = 42  # Fixed seed for all random operations

# =============================================================================
# CONSOLE OUTPUT CONTROL
# =============================================================================
QUIET_MODE = False  # Set True to suppress verbose module output during pipeline

# =============================================================================
# DATA CLEANING THRESHOLDS
# =============================================================================
MIN_ANTIBIOTIC_COVERAGE = 70.0  # Minimum % of isolates an antibiotic must be tested on
MAX_ISOLATE_MISSING = 30.0  # Maximum % of missing antibiotic data per isolate

# =============================================================================
# RESISTANCE ENCODING
# =============================================================================
RESISTANCE_ENCODING = {
    'S': 0,  # Susceptible
    'I': 1,  # Intermediate
    'R': 2   # Resistant
}
RESISTANCE_DECODING = {v: k for k, v in RESISTANCE_ENCODING.items()}
RESISTANCE_THRESHOLD = 2  # Value >= this is considered "Resistant"

# =============================================================================
# MDR CLASSIFICATION (Magiorakos et al., 2012)
# =============================================================================
MDR_MIN_CLASSES = 3  # Minimum antibiotic classes for MDR classification
MDR_REFERENCE = "Magiorakos AP, et al. (2012). Clin Microbiol Infect, 18(3):268-281."

# -----------------------------------------------------------------------------
# IMPORTANT LIMITATION: SPECIES-AGNOSTIC MDR CLASSIFICATION
# -----------------------------------------------------------------------------
# The Magiorakos et al. (2012) framework specifies MDR definitions for each
# organism SEPARATELY, accounting for intrinsic resistances:
#
# - Pseudomonas aeruginosa: Has intrinsic penicillin resistance (chromosomal 
#   AmpC β-lactamase). Penicillin resistance should NOT count toward MDR.
#
# - Acinetobacter spp.: Has intrinsic resistance to aminopenicillins, 1st/2nd
#   generation cephalosporins. These should be excluded from MDR calculation.
#
# - Enterobacteriaceae (E. coli, Klebsiella, etc.): No significant intrinsic
#   resistances. All acquired resistances count toward MDR.
#
# CURRENT IMPLEMENTATION: This pipeline uses a UNIVERSAL antibiotic class
# mapping (ANTIBIOTIC_CLASSES below) applied to ALL species equally. This
# approach may:
# - OVERESTIMATE MDR rates for intrinsically resistant species (Pseudomonas)
# - Provide COMPARABLE rates for non-intrinsically resistant species (E. coli)
#
# SCIENTIFIC JUSTIFICATION: For environmental surveillance (vs. clinical), a
# universal classification simplifies cross-species comparison and pattern
# recognition. Future versions should implement species-specific class mappings
# for clinical applications.
#
# Reference for intrinsic resistances:
# - EUCAST Expert Rules 3.2 (2020): Intrinsic Resistance and Unusual Phenotypes
# - Poirel et al. (2018). Clin Microbiol Rev, 31(2):e00088-17
# -----------------------------------------------------------------------------

# =============================================================================
# MAR INDEX (Krumperman, 1983)
# =============================================================================
MAR_REFERENCE = "Krumperman PH. (1983). Appl Environ Microbiol, 46(1):165-170."

# =============================================================================
# ANTIBIOTIC CLASS DEFINITIONS
# Reference: CDC/CLSI Guidelines for Antimicrobial Susceptibility Testing
# This is the SINGLE SOURCE OF TRUTH for antibiotic class mappings
# =============================================================================
ANTIBIOTIC_CLASSES: Dict[str, str] = {
    # Penicillins
    'AM': 'Penicillins',
    'AMP': 'Penicillins',
    
    # Beta-lactam/Beta-lactamase inhibitor combinations
    'AMC': 'BL/BLI combinations',
    'AMO': 'BL/BLI combinations',
    'PRA': 'BL/BLI combinations',
    
    # Cephalosporins (1st gen)
    'CN': 'Cephalosporins-1st',
    'CF': 'Cephalosporins-1st',
    'CFX': 'Cephalosporins-1st',
    
    # Cephalosporins (3rd/4th gen)
    'CPD': 'Cephalosporins-3rd/4th',
    'CTX': 'Cephalosporins-3rd/4th',
    'CFT': 'Cephalosporins-3rd/4th',
    'CPT': 'Cephalosporins-3rd/4th',
    'CFA': 'Cephalosporins-3rd/4th',
    'CFV': 'Cephalosporins-3rd/4th',
    'CTF': 'Cephalosporins-3rd/4th',
    
    # Cephamycins
    'CFO': 'Cephamycins',
    
    # Cephalosporin/BLI combinations
    'CZA': 'Cephalosporin/BLI',
    'CFZ': 'Cephalosporin/BLI',
    
    # Carbapenems
    'IPM': 'Carbapenems',
    'MRB': 'Carbapenems',
    'IME': 'Carbapenems',
    'MAR': 'Carbapenems',
    
    # Aminoglycosides
    'AN': 'Aminoglycosides',
    'GM': 'Aminoglycosides',
    'N': 'Aminoglycosides',
    'AMI': 'Aminoglycosides',
    'GEN': 'Aminoglycosides',
    'NEO': 'Aminoglycosides',
    
    # Fluoroquinolones
    'NAL': 'Quinolones',
    'ENR': 'Fluoroquinolones',
    'NLA': 'Quinolones',
    
    # Tetracyclines
    'DO': 'Tetracyclines',
    'TE': 'Tetracyclines',
    'DOX': 'Tetracyclines',
    'TET': 'Tetracyclines',
    
    # Nitrofurans
    'FT': 'Nitrofurans',
    'NIT': 'Nitrofurans',
    
    # Phenicols
    'C': 'Phenicols',
    'CHL': 'Phenicols',
    
    # Folate pathway inhibitors
    'SXT': 'Folate pathway inhibitors',
}

# =============================================================================
# SPECIES-SPECIFIC MDR CLASS DEFINITIONS (Action Item 7)
# =============================================================================
# Magiorakos et al. (2012) specifies that MDR should be calculated using
# species-appropriate antibiotic categories, excluding intrinsic resistances.
#
# This dictionary maps species to their RELEVANT antibiotic classes for MDR.
# Classes NOT listed for a species are excluded from that species' MDR calculation.
# =============================================================================

MDR_CLASSES_BY_SPECIES: Dict[str, Dict[str, str]] = {
    # -------------------------------------------------------------------------
    # ENTEROBACTERIACEAE - No significant intrinsic resistances
    # All commonly tested classes are relevant
    # -------------------------------------------------------------------------
    'Escherichia coli': {
        'AM': 'Penicillins',
        'AMP': 'Penicillins',
        'AMC': 'BL/BLI combinations',
        'AMO': 'BL/BLI combinations',
        'CN': 'Cephalosporins-1st',
        'CF': 'Cephalosporins-1st',
        'CPD': 'Cephalosporins-3rd/4th',
        'CTX': 'Cephalosporins-3rd/4th',
        'CPT': 'Cephalosporins-3rd/4th',
        'CFO': 'Cephamycins',
        'IPM': 'Carbapenems',
        'MRB': 'Carbapenems',
        'AN': 'Aminoglycosides',
        'GM': 'Aminoglycosides',
        'NAL': 'Quinolones',
        'ENR': 'Fluoroquinolones',
        'TE': 'Tetracyclines',
        'DO': 'Tetracyclines',
        'FT': 'Nitrofurans',
        'C': 'Phenicols',
        'SXT': 'Folate pathway inhibitors',
    },
    
    'Klebsiella pneumoniae': {
        # Same as E. coli - no intrinsic resistance to tested agents
        'AM': 'Penicillins',
        'AMP': 'Penicillins',
        'AMC': 'BL/BLI combinations',
        'CPD': 'Cephalosporins-3rd/4th',
        'CTX': 'Cephalosporins-3rd/4th',
        'IPM': 'Carbapenems',
        'AN': 'Aminoglycosides',
        'GM': 'Aminoglycosides',
        'NAL': 'Quinolones',
        'ENR': 'Fluoroquinolones',
        'TE': 'Tetracyclines',
        'FT': 'Nitrofurans',
        'C': 'Phenicols',
        'SXT': 'Folate pathway inhibitors',
    },
    
    # -------------------------------------------------------------------------
    # PSEUDOMONAS AERUGINOSA - INTRINSIC RESISTANCES
    # Exclude: Penicillins (ampicillin), 1st/2nd gen cephalosporins, cephamycins
    # tetracyclines, trimethoprim-sulfamethoxazole, chloramphenicol
    # Reference: EUCAST Expert Rules 3.2 (2020)
    # -------------------------------------------------------------------------
    'Pseudomonas aeruginosa': {
        # Only anti-pseudomonal agents are relevant
        'CZA': 'Anti-pseudomonal cephalosporins',  # Ceftazidime-avibactam
        'CFZ': 'Anti-pseudomonal cephalosporins',  # Ceftazidime
        'PRA': 'Anti-pseudomonal penicillins',     # Piperacillin-tazobactam
        'IPM': 'Carbapenems',
        'MRB': 'Carbapenems',
        'AN': 'Aminoglycosides',
        'GM': 'Aminoglycosides',
        'ENR': 'Fluoroquinolones',                 # Enrofloxacin/Ciprofloxacin
        # NOTE: Excluded due to intrinsic resistance:
        # - Ampicillin, amoxicillin (chromosomal AmpC)
        # - 1st/2nd gen cephalosporins
        # - Tetracyclines (MexAB-OprM efflux)
        # - Trimethoprim-sulfamethoxazole
        # - Chloramphenicol
    },
    
    # -------------------------------------------------------------------------
    # ACINETOBACTER SPP. - INTRINSIC RESISTANCES
    # Exclude: Aminopenicillins, 1st/2nd gen cephalosporins, cephamycins
    # Reference: EUCAST Expert Rules 3.2 (2020)
    # -------------------------------------------------------------------------
    'Acinetobacter baumannii': {
        # Only certain classes are reliable for MDR assessment
        'CZA': 'Anti-Acinetobacter cephalosporins',
        'IPM': 'Carbapenems',
        'MRB': 'Carbapenems',
        'AN': 'Aminoglycosides',
        'GM': 'Aminoglycosides',
        'ENR': 'Fluoroquinolones',
        'TE': 'Tetracyclines',  # Variable susceptibility
        'DO': 'Tetracyclines',
        'SXT': 'Folate pathway inhibitors',  # Variable
        # NOTE: Excluded due to intrinsic resistance:
        # - Aminopenicillins
        # - 1st/2nd gen cephalosporins
        # - Cephamycins
    },
    
    # -------------------------------------------------------------------------
    # ENTEROBACTER SPP. - INDUCIBLE AmpC
    # Exclude: Ampicillin, 1st/2nd gen cephalosporins, cephamycins
    # -------------------------------------------------------------------------
    'Enterobacter cloacae': {
        'AMC': 'BL/BLI combinations',  # May still be useful
        'CPD': 'Cephalosporins-3rd/4th',
        'CTX': 'Cephalosporins-3rd/4th',
        'IPM': 'Carbapenems',
        'AN': 'Aminoglycosides',
        'GM': 'Aminoglycosides',
        'NAL': 'Quinolones',
        'ENR': 'Fluoroquinolones',
        'TE': 'Tetracyclines',
        'C': 'Phenicols',
        'SXT': 'Folate pathway inhibitors',
        # NOTE: Excluded - inducible AmpC makes these unreliable:
        # - Ampicillin
        # - 1st/2nd gen cephalosporins
        # - Cephamycins
    },
    
    # -------------------------------------------------------------------------
    # SALMONELLA SPP. - Environmental/enteric pathogen
    # All typically tested classes are relevant
    # -------------------------------------------------------------------------
    'Salmonella enterica': {
        'AM': 'Penicillins',
        'AMP': 'Penicillins',
        'CPD': 'Cephalosporins-3rd/4th',
        'CTX': 'Cephalosporins-3rd/4th',
        'AN': 'Aminoglycosides',
        'GM': 'Aminoglycosides',
        'NAL': 'Quinolones',
        'ENR': 'Fluoroquinolones',
        'TE': 'Tetracyclines',
        'C': 'Phenicols',
        'SXT': 'Folate pathway inhibitors',
    },
    
    # -------------------------------------------------------------------------
    # VIBRIO SPP. - Environmental aquatic pathogen
    # -------------------------------------------------------------------------
    'Vibrio': {
        'AM': 'Penicillins',
        'TE': 'Tetracyclines',
        'DO': 'Tetracyclines',
        'NAL': 'Quinolones',
        'ENR': 'Fluoroquinolones',
        'AN': 'Aminoglycosides',
        'C': 'Phenicols',
        'SXT': 'Folate pathway inhibitors',
    },
}

# Function to get appropriate class mapping for a species
def get_mdr_classes_for_species(species: str) -> Dict[str, str]:
    """
    Get the appropriate antibiotic class mapping for MDR calculation.
    
    Parameters
    ----------
    species : str
        Species name (will be matched flexibly)
    
    Returns
    -------
    dict
        Antibiotic-to-class mapping for the species.
        Falls back to universal ANTIBIOTIC_CLASSES if species not defined.
    """
    if not species or not isinstance(species, str):
        return ANTIBIOTIC_CLASSES
    
    species_lower = species.lower().strip()
    
    # Direct match
    for defined_species, classes in MDR_CLASSES_BY_SPECIES.items():
        if species_lower == defined_species.lower():
            return classes
    
    # Partial match (e.g., "E. coli" matches "Escherichia coli")
    for defined_species, classes in MDR_CLASSES_BY_SPECIES.items():
        if defined_species.lower() in species_lower or species_lower in defined_species.lower():
            return classes
    
    # Genus-level match
    species_genus = species_lower.split()[0] if ' ' in species_lower else species_lower
    for defined_species, classes in MDR_CLASSES_BY_SPECIES.items():
        defined_genus = defined_species.lower().split()[0]
        if species_genus == defined_genus:
            return classes
    
    # Fallback to universal mapping (with warning in calling code if needed)
    return ANTIBIOTIC_CLASSES

# =============================================================================
# CLUSTERING CONFIGURATION
# =============================================================================
CLUSTERING_CONFIG = {
    'linkage_method': 'ward',      # Ward linkage minimizes within-cluster variance - standard for phenotype clustering
    'distance_metric': 'euclidean', # Required by Ward linkage
    
    # K SELECTION: Determined AUTOMATICALLY by elbow + silhouette analysis
    # The pipeline evaluates all k values in [min_k, max_k] range and selects
    # the optimal k based on statistical criteria. NO preferred/default k is used.
    # See: validate_clustering.py:find_optimal_k() for implementation.
    'min_k': 3,                    # Minimum clusters to evaluate (avoid trivial solutions)
    'max_k': 8,                    # Maximum clusters to evaluate (algorithm selects optimal within range)
    'k_selection_method': 'combined',  # 'silhouette', 'elbow', or 'combined' (recommended)
    
    'silhouette_threshold': 0.40,  # Score >= this indicates "strong" cluster structure (Rousseeuw, 1987)
    'imputation_strategy': 'median', # Median is robust to outliers for missing value imputation
}

# =============================================================================
# SUPERVISED LEARNING CONFIGURATION
# =============================================================================
SUPERVISED_CONFIG = {
    'test_size': 0.20,  # 80/20 train-test split
    'cv_folds': 5,  # Cross-validation folds
    'min_samples_per_class': 2,  # Minimum samples required per class
    'random_state': RANDOM_STATE,
    'models': {
        'logistic_regression': {
            'max_iter': 1000,
            'solver': 'lbfgs',
        },
        'random_forest': {
            'n_estimators': 100,
            'n_jobs': -1,
        },
        'knn': {
            'n_neighbors': 5,
        }
    }
}

# =============================================================================
# ENVIRONMENT MAPPING
# =============================================================================
ENVIRONMENT_MAPPING: Dict[str, str] = {
    # Water sources
    'DW': 'Water',
    'DRINKING WATER': 'Water',
    'RIV': 'Water',
    'RIVER': 'Water',
    'STREAM': 'Water',
    'EW': 'Water',
    'EFFLUENT WATER': 'Hospital',  # Hospital effluent
    
    # Animal/Food sources
    'TL': 'Fish',
    'TILAPIA': 'Fish',
    'ML': 'Fish',
    'MILKFISH': 'Fish',
    'FISH': 'Fish',
    'CHICKEN': 'Animal',
    'POULTRY': 'Animal',
    
    # Other
    'SOIL': 'Soil',
    'SEDIMENT': 'Sediment',
}

# =============================================================================
# SPECIES ABBREVIATIONS TO FULL NAMES
# =============================================================================
SPECIES_ABBREVIATIONS: Dict[str, str] = {
    'EC': 'Escherichia coli',
    'KP': 'Klebsiella pneumoniae',
    'ECLA': 'Enterobacter cloacae',
    'EAER': 'Enterobacter aerogenes',
    'SAL': 'Salmonella group',
    'SEN': 'Salmonella enterica spp diarizonae',
    'KS': 'Klebsiella pneumoniae ssp ozaenae',
    'PS': 'Pseudomonas aeruginosa',
    'PRO': 'Proteus species',
    'SER': 'Serratia species',
}

# =============================================================================
# METADATA COLUMN DEFINITIONS
# =============================================================================
METADATA_COLUMNS = [
    'CODE', 'ISOLATE_ID', 'REGION', 'SITE', 'ENVIRONMENT',
    'SAMPLING_SOURCE', 'NATIONAL_SITE', 'LOCAL_SITE',
    'REPLICATE', 'COLONY', 'ESBL', 'SOURCE_FILE',
    'resistance_fingerprint', 'MAR_INDEX_COMPUTED',
    'RESISTANCE_COUNT', 'RESISTANT_CLASSES_COUNT',
    'MDR_FLAG', 'MDR_CATEGORY'
]

# =============================================================================
# OUTPUT FILE NAMING
# =============================================================================
OUTPUT_FILES = {
    'unified_raw': 'unified_raw_dataset.csv',
    'cleaned': 'cleaned_dataset.csv',
    'cleaning_report': 'cleaning_report.txt',
    'encoded': 'encoded_dataset.csv',
    'analysis_ready': 'analysis_ready_dataset.csv',
    'clustered': 'clustered_dataset.csv',
    'feature_matrix': 'feature_matrix_X.csv',
    'metadata': 'metadata.csv',
    'cluster_summary': 'cluster_summary_table.csv',
    'cluster_validation': 'cluster_validation.png',
    'dendrogram': 'dendrogram_highres.png',
    'heatmap': 'dendrogram_linked_heatmap.png',
}

# =============================================================================
# FULL CONFIG DICTIONARY (for passing to functions)
# =============================================================================
CONFIG: Dict[str, Any] = {
    'paths': {
        'project_root': PROJECT_ROOT,
        'data_dir': DATA_DIR,
        'raw_data_dir': RAW_DATA_DIR,
        'processed_dir': PROCESSED_DATA_DIR,
        'figures_dir': FIGURES_DIR,
        'models_dir': MODELS_DIR,
    },
    'random_state': RANDOM_STATE,
    'cleaning': {
        'min_antibiotic_coverage': MIN_ANTIBIOTIC_COVERAGE,
        'max_isolate_missing': MAX_ISOLATE_MISSING,
    },
    'encoding': {
        'resistance_encoding': RESISTANCE_ENCODING,
        'resistance_threshold': RESISTANCE_THRESHOLD,
    },
    'mdr': {
        'min_classes': MDR_MIN_CLASSES,
        'reference': MDR_REFERENCE,
    },
    'clustering': CLUSTERING_CONFIG,
    'supervised': SUPERVISED_CONFIG,
    'antibiotic_classes': ANTIBIOTIC_CLASSES,
    'environment_mapping': ENVIRONMENT_MAPPING,
    'species_abbreviations': SPECIES_ABBREVIATIONS,
    'output_files': OUTPUT_FILES,
}


def get_antibiotic_class(antibiotic_code: str) -> str:
    """
    Get the antibiotic class for a given antibiotic code.
    
    Parameters:
    -----------
    antibiotic_code : str
        Antibiotic abbreviation (e.g., 'AM', 'TE')
    
    Returns:
    --------
    str
        Antibiotic class name or 'Unknown' if not found
    """
    return ANTIBIOTIC_CLASSES.get(antibiotic_code.upper(), 'Unknown')


def detect_antibiotic_columns(df) -> List[str]:
    """
    Dynamically detect antibiotic columns from dataframe.
    
    This function identifies antibiotic columns using multiple strategies:
    1. Columns ending with '_encoded' (processed data)
    2. Columns matching known antibiotic codes
    3. Columns with S/I/R values (raw data)
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe
    
    Returns:
    --------
    List[str]
        List of detected antibiotic column names
    """
    import pandas as pd
    
    # Strategy 1: Look for _encoded suffix
    encoded_cols = [c for c in df.columns if c.endswith('_encoded')]
    if encoded_cols:
        return encoded_cols
    
    # Strategy 2: Match known antibiotic codes
    known_abs = set(ANTIBIOTIC_CLASSES.keys())
    matched_cols = [c for c in df.columns if c.upper() in known_abs]
    if matched_cols:
        return matched_cols
    
    # Strategy 3: Detect S/I/R columns
    sir_cols = []
    for col in df.columns:
        if df[col].dtype == 'object':
            unique_vals = set(df[col].dropna().astype(str).str.upper().unique())
            if unique_vals.issubset({'S', 'I', 'R', 'NaN', '', 'NA'}):
                if len(unique_vals.intersection({'S', 'I', 'R'})) > 0:
                    sir_cols.append(col)
    
    return sir_cols


if __name__ == "__main__":
    print("AMR Pipeline Configuration")
    print("=" * 50)
    print(f"Project Root: {PROJECT_ROOT}")
    print(f"Random State: {RANDOM_STATE}")
    print(f"Antibiotic Classes Defined: {len(ANTIBIOTIC_CLASSES)}")
    print(f"Clustering Method: {CLUSTERING_CONFIG['linkage_method']}")
    print(f"K Selection Range: {CLUSTERING_CONFIG['min_k']}-{CLUSTERING_CONFIG['max_k']} (auto-determined)")
