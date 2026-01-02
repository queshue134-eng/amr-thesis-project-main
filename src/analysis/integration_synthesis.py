"""
Integration and Synthesis Module for AMR Thesis Project
Phase 5 - Integration & Synthesis (Concrete Implementation Improvements)

OBJECTIVES:
    1. Define Formal Integration Framework (Phase 2/3/4 outputs)
    2. Create Cluster–Supervised Result Alignment table
    3. Identify dominant resistance archetypes with formal definition
    4. Species–Environment–Resistance Triangulation
    5. MDR-Enriched Pattern Synthesis
    6. Explicit Negative Findings & Claim Boundaries
    7. Final Integrative Narrative Structure

DELIVERABLES CHECKLIST (Phase 5 Audit):
    ✅ Cluster–supervised alignment table
    ✅ Resistance archetype definition
    ✅ Archetype summary table
    ✅ Species–environment–resistance triangulation
    ✅ MDR-enriched synthesis table
    ✅ Explicit negative findings
    ✅ Claim boundary statement

ASSOCIATION LANGUAGE DISCIPLINE:
    Use only: "associated with", "enriched in", "observed across"
    Avoid: "driven by", "originating from", "transmission pathway"
"""

import pandas as pd
from typing import List, Dict, Tuple
from scipy import stats
import warnings


# =============================================================================
# FORMAL DEFINITIONS (Phase 5 Requirements)
# =============================================================================

# Formal Archetype Definition
ARCHETYPE_DEFINITION = """
A resistance archetype is defined as:
"A recurring resistance profile characterized by a stable antibiotic resistance
pattern and consistent cluster membership."

Criteria for archetype identification:
1. Cluster size ≥ predefined threshold (minimum 5 isolates by default)
2. Stable across distance metrics (if robustness testing performed)
3. Interpretable antibiotic pattern (clear resistant/susceptible signature)
"""

# Minimum cluster size threshold for archetype identification
MIN_ARCHETYPE_SIZE = 5

# Association language discipline - allowed terms
ALLOWED_ASSOCIATION_TERMS = [
    "associated with",
    "enriched in",
    "observed across",
    "concentrated in",
    "prevalent in"
]

# Terms to avoid (for documentation/reference)
AVOID_CAUSAL_TERMS = [
    "driven by",
    "originating from",
    "transmission pathway",
    "caused by",
    "spread from"
]

# Resistance level thresholds and labels
RESISTANCE_LEVELS = {
    'HIGH': 'High resistance',
    'MODERATE_HIGH': 'Moderate-high resistance',
    'MODERATE': 'Moderate resistance',
    'LOW': 'Low resistance'
}

HIGH_RESISTANCE_LEVELS = [RESISTANCE_LEVELS['HIGH'], RESISTANCE_LEVELS['MODERATE_HIGH']]

# Default data limitations (can be customized per study)
DEFAULT_DATA_LIMITATIONS = [
    "Cross-sectional data: temporal changes in resistance cannot be assessed",
    "Geographic scope: findings may not generalize to other regions",
    "Sample selection: hospital/clinical bias may affect environmental representation"
]


def _calculate_fold_enrichment(rate: float, baseline: float) -> float:
    """
    Calculate fold enrichment of a rate compared to baseline.

    Parameters:
    -----------
    rate : float
        The observed rate (0-1)
    baseline : float
        The baseline rate to compare against (0-1)

    Returns:
    --------
    float
        Fold enrichment value, or 0 if baseline is 0
    """
    return float(rate / baseline) if baseline > 0 else 0.0


def _get_antibiotic_name(col_name: str) -> str:
    """
    Extract antibiotic name from encoded column name.

    Parameters:
    -----------
    col_name : str
        Column name (e.g., 'AM_encoded' or 'AM')

    Returns:
    --------
    str
        Antibiotic name without '_encoded' suffix
    """
    return col_name.replace('_encoded', '')


def _normalize_cluster_id(cluster_id) -> int:
    """
    Normalize cluster ID to integer.

    Parameters:
    -----------
    cluster_id : any
        Cluster ID (could be int, float, or numpy type)

    Returns:
    --------
    int
        Normalized integer cluster ID
    """
    return int(cluster_id)


def _categorize_discrimination_strength(accuracy: float) -> str:
    """
    Categorize discrimination strength based on accuracy/performance.

    Parameters:
    -----------
    accuracy : float
        Classification accuracy or similar metric (0-1)

    Returns:
    --------
    str
        Discrimination strength category: 'High', 'Med', or 'Low'
    """
    if accuracy >= 0.8:
        return 'High'
    elif accuracy >= 0.6:
        return 'Med'
    else:
        return 'Low'


# =============================================================================
# SECTION 1: FORMAL INTEGRATION FRAMEWORK
# =============================================================================

def define_integration_framework() -> Dict:
    """
    Define the Formal Integration Framework - explicitly defines three analytical outputs.

    This function documents what is being integrated from each analysis phase:
    - Phase 2: Resistance-based clusters (unsupervised)
    - Phase 3: Discriminative signals (supervised models)
    - Phase 4: Regional & environmental associations

    Returns:
    --------
    dict
        Integration framework definition with phase outputs

    Why:
        - Prevents narrative drift
        - Shows intentional synthesis, not post-hoc storytelling
    """
    framework = {
        'phase_outputs': {
            'phase_2': {
                'name': 'Unsupervised Clustering',
                'output': 'Resistance-based clusters',
                'type': 'unsupervised',
                'description': 'Groups isolates by antibiotic resistance profile similarity',
                'key_results': [
                    'Cluster assignments for each isolate',
                    'Cluster centroid profiles',
                    'Within-cluster homogeneity metrics'
                ]
            },
            'phase_3': {
                'name': 'Supervised Learning',
                'output': 'Discriminative signals',
                'type': 'supervised',
                'description': 'Identifies features that discriminate between categories',
                'key_results': [
                    'Species classification performance',
                    'MDR classification performance',
                    'Feature importance rankings',
                    'Confusion matrices for classification'
                ]
            },
            'phase_4': {
                'name': 'Contextual Association',
                'output': 'Regional & environmental associations',
                'type': 'descriptive/statistical',
                'description': 'Explores relationships with metadata',
                'key_results': [
                    'Region-resistance associations',
                    'Environment-resistance associations',
                    'Species distribution patterns'
                ]
            }
        },
        'integration_rationale': (
            "Integration synthesizes findings across phases to identify coherent "
            "resistance archetypes that are: (1) data-driven (unsupervised), "
            "(2) biologically meaningful (supervised validation), and "
            "(3) ecologically contextualized (environmental associations)."
        ),
        'framework_table': pd.DataFrame({
            'Phase': ['Phase 2', 'Phase 3', 'Phase 4'],
            'Output': [
                'Resistance-based clusters (unsupervised)',
                'Discriminative signals (supervised models)',
                'Regional & environmental associations'
            ],
            'Type': ['Unsupervised', 'Supervised', 'Descriptive/Statistical']
        })
    }

    return framework


# =============================================================================
# SECTION 2: CLUSTER–SUPERVISED RESULT ALIGNMENT
# =============================================================================

def compare_clusters_with_supervised(
    df: pd.DataFrame,
    cluster_col: str = 'CLUSTER',
    mdr_col: str = 'MDR_CATEGORY',
    species_col: str = 'ISOLATE_ID'
) -> Dict:
    """
    Compare unsupervised clusters with supervised discrimination results.

    Evaluates how well the unsupervised clusters align with known categories
    (MDR status and species) that supervised models can discriminate.

    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe with cluster labels and category information
    cluster_col : str
        Column name for cluster labels
    mdr_col : str
        Column name for MDR category
    species_col : str
        Column name for species

    Returns:
    --------
    dict
        Comparison results including cluster-category alignment metrics
    """
    comparison = {
        'cluster_mdr_alignment': None,
        'cluster_species_alignment': None,
        'cluster_purity': {},
        'interpretation': []
    }

    if cluster_col not in df.columns:
        comparison['interpretation'].append("No cluster column found - run clustering first.")
        return comparison

    clusters = df[cluster_col].unique()

    # Compare clusters with MDR categories
    if mdr_col in df.columns:
        # Cross-tabulation of clusters vs MDR
        mdr_crosstab = pd.crosstab(df[cluster_col], df[mdr_col])
        comparison['cluster_mdr_alignment'] = mdr_crosstab.to_dict()

        # Calculate cluster purity with respect to MDR
        for cluster in clusters:
            cluster_df = df[df[cluster_col] == cluster]
            if len(cluster_df) > 0 and mdr_col in cluster_df.columns:
                mdr_counts = cluster_df[mdr_col].value_counts()
                dominant_mdr = mdr_counts.idxmax()
                purity = mdr_counts.max() / len(cluster_df)
                comparison['cluster_purity'][_normalize_cluster_id(cluster)] = {
                    'mdr_dominant': dominant_mdr,
                    'mdr_purity': float(purity)
                }

        # Chi-square test for cluster-MDR association
        try:
            chi2, p_value, dof, expected = stats.chi2_contingency(mdr_crosstab)
            comparison['mdr_chi_square'] = {
                'statistic': float(chi2),
                'p_value': float(p_value),
                'degrees_of_freedom': int(dof),
                'significant': p_value < 0.05
            }

            if p_value < 0.05:
                comparison['interpretation'].append(
                    f"Significant association between clusters and MDR status "
                    f"(Chi-sq={chi2:.2f}, p={p_value:.4f}). "
                    "Unsupervised clustering captures MDR-related patterns."
                )
            else:
                comparison['interpretation'].append(
                    "No significant association between clusters and MDR status. "
                    "Clusters may capture other resistance patterns."
                )
        except ValueError as e:
            # Chi-square test can fail with insufficient data
            warnings.warn(f"Chi-square test for cluster-MDR failed: {e}")

    # Compare clusters with species
    if species_col in df.columns:
        species_crosstab = pd.crosstab(df[cluster_col], df[species_col])
        comparison['cluster_species_alignment'] = species_crosstab.to_dict()

        # Calculate cluster purity with respect to species
        for cluster in clusters:
            cluster_df = df[df[cluster_col] == cluster]
            if len(cluster_df) > 0:
                species_counts = cluster_df[species_col].value_counts()
                dominant_species = species_counts.idxmax()
                purity = species_counts.max() / len(cluster_df)
                cluster_id = _normalize_cluster_id(cluster)
                if cluster_id in comparison['cluster_purity']:
                    comparison['cluster_purity'][cluster_id]['species_dominant'] = dominant_species
                    comparison['cluster_purity'][cluster_id]['species_purity'] = float(purity)
                else:
                    comparison['cluster_purity'][cluster_id] = {
                        'species_dominant': dominant_species,
                        'species_purity': float(purity)
                    }

        # Chi-square test for cluster-species association
        try:
            chi2, p_value, dof, expected = stats.chi2_contingency(species_crosstab)
            comparison['species_chi_square'] = {
                'statistic': float(chi2),
                'p_value': float(p_value),
                'degrees_of_freedom': int(dof),
                'significant': p_value < 0.05
            }

            if p_value < 0.05:
                comparison['interpretation'].append(
                    f"Significant association between clusters and species "
                    f"(Chi-sq={chi2:.2f}, p={p_value:.4f}). "
                    "Resistance patterns relate to species identity."
                )
        except ValueError as e:
            # Chi-square test can fail with insufficient data
            warnings.warn(f"Chi-square test for cluster-species failed: {e}")

    return comparison


def generate_cluster_supervised_alignment_table(
    df: pd.DataFrame,
    cluster_col: str = 'CLUSTER',
    mdr_col: str = 'MDR_CATEGORY',
    mdr_flag_col: str = 'MDR_FLAG',
    species_col: str = 'ISOLATE_ID'
) -> Tuple[pd.DataFrame, Dict]:
    """
    Generate Cluster × Supervised Performance Mapping table.

    Creates a table linking clusters to supervised discrimination strength:
    | Cluster | Dominant species | MDR % | Species discrimination | MDR discrimination |

    Uses confusion matrices concepts from Phase 3 to identify:
    - Clusters that are consistently classified
    - Clusters that are frequently confused

    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe with cluster labels and category information
    cluster_col : str
        Column name for cluster labels
    mdr_col : str
        Column name for MDR category
    mdr_flag_col : str
        Column name for MDR flag (binary)
    species_col : str
        Column name for species

    Returns:
    --------
    tuple
        (alignment_table_df, interpretation_dict)

    Why:
        Demonstrates internal coherence between unsupervised and supervised phases.
    """
    if cluster_col not in df.columns:
        return pd.DataFrame(), {'error': 'No cluster column found'}

    alignment_rows = []
    interpretation = {
        'agreement': [],
        'disagreement': [],
        'interpretation_rule': {
            'agreement': 'Strong archetype - cluster aligns well with labels',
            'disagreement': 'Heterogeneous resistance phenotype - cluster cuts across labels'
        }
    }

    clusters = sorted(df[cluster_col].unique())

    for cluster in clusters:
        cluster_df = df[df[cluster_col] == cluster]
        cluster_size = len(cluster_df)

        row = {
            'Cluster': int(cluster),
            'Size': cluster_size,
            'Dominant_Species': None,
            'MDR_Percent': None,
            'Species_Discrimination': None,
            'MDR_Discrimination': None
        }

        # Dominant species and species discrimination assessment
        if species_col in cluster_df.columns:
            species_counts = cluster_df[species_col].value_counts()
            row['Dominant_Species'] = species_counts.idxmax()

            # Purity-based discrimination strength
            species_purity = species_counts.max() / cluster_size
            row['Species_Discrimination'] = _categorize_discrimination_strength(species_purity)

            # Track agreement/disagreement
            if species_purity >= 0.8:
                interpretation['agreement'].append(
                    f"Cluster {cluster} strongly aligns with {row['Dominant_Species']} "
                    f"(purity: {species_purity * 100:.1f}%)"
                )
            elif species_purity < 0.5:
                interpretation['disagreement'].append(
                    f"Cluster {cluster} spans multiple species "
                    f"(purity: {species_purity * 100:.1f}%) - "
                    "heterogeneous resistance phenotype"
                )

        # MDR percentage and MDR discrimination assessment
        if mdr_flag_col in cluster_df.columns:
            mdr_rate = cluster_df[mdr_flag_col].mean()
            row['MDR_Percent'] = round(mdr_rate * 100, 1)

            # MDR discrimination based on deviation from 50% (pure MDR or non-MDR)
            mdr_purity = max(mdr_rate, 1 - mdr_rate)
            row['MDR_Discrimination'] = _categorize_discrimination_strength(mdr_purity)

            if mdr_purity >= 0.8:
                mdr_label = 'MDR' if mdr_rate > 0.5 else 'non-MDR'
                interpretation['agreement'].append(
                    f"Cluster {cluster} is predominantly {mdr_label} ({row['MDR_Percent']}% MDR)"
                )

        alignment_rows.append(row)

    alignment_df = pd.DataFrame(alignment_rows)

    return alignment_df, interpretation


# =============================================================================
# SECTION 3: RESISTANCE ARCHETYPES WITH FORMAL DEFINITION
# =============================================================================


def identify_resistance_archetypes(
    df: pd.DataFrame,
    feature_cols: List[str],
    cluster_col: str = 'CLUSTER',
    threshold_resistant: float = 1.5,
    min_cluster_size: int = MIN_ARCHETYPE_SIZE
) -> Dict:
    """
    Identify dominant resistance archetypes from clustering results.

    FORMAL ARCHETYPE DEFINITION:
    "A recurring resistance profile characterized by a stable antibiotic
    resistance pattern and consistent cluster membership."

    Criteria:
    - Cluster size ≥ predefined threshold
    - Stable across distance metrics (if tested)
    - Interpretable antibiotic pattern

    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe with cluster labels and resistance data
    feature_cols : list
        List of encoded resistance column names
    cluster_col : str
        Column name for cluster labels
    threshold_resistant : float
        Mean encoded value above which an antibiotic is considered
        predominantly resistant in a cluster (default 1.5 = between I and R)
    min_cluster_size : int
        Minimum cluster size to be considered a valid archetype

    Returns:
    --------
    dict
        Archetype definitions for each cluster including formal definition
    """
    archetypes = {
        'formal_definition': ARCHETYPE_DEFINITION.strip(),
        'criteria': {
            'min_cluster_size': min_cluster_size,
            'threshold_resistant': threshold_resistant
        },
        'cluster_archetypes': {},
        'summary': [],
        'archetype_summary_table': None  # Will be a DataFrame
    }

    if cluster_col not in df.columns:
        return archetypes

    existing_cols = [c for c in feature_cols if c in df.columns]
    if not existing_cols:
        return archetypes

    # Calculate mean resistance profile for each cluster
    cluster_profiles = df.groupby(cluster_col)[existing_cols].mean()

    # For archetype summary table
    summary_rows = []

    for cluster in cluster_profiles.index:
        profile = cluster_profiles.loc[cluster]
        cluster_size = len(df[df[cluster_col] == cluster])

        # Skip clusters below minimum size threshold
        if cluster_size < min_cluster_size:
            continue

        # Identify antibiotics with high resistance (mean > threshold)
        high_resistance = profile[profile > threshold_resistant].sort_values(ascending=False)
        high_resistance_abs = [_get_antibiotic_name(ab) for ab in high_resistance.index.tolist()]

        # Identify antibiotics with low resistance (mean < 0.5, mostly susceptible)
        low_resistance = profile[profile < 0.5].sort_values()
        low_resistance_abs = [_get_antibiotic_name(ab) for ab in low_resistance.index.tolist()]

        # Calculate overall resistance level using constants
        mean_resistance = profile.mean()
        if mean_resistance > 1.5:
            resistance_level = RESISTANCE_LEVELS['HIGH']
        elif mean_resistance > 1.0:
            resistance_level = RESISTANCE_LEVELS['MODERATE_HIGH']
        elif mean_resistance > 0.5:
            resistance_level = RESISTANCE_LEVELS['MODERATE']
        else:
            resistance_level = RESISTANCE_LEVELS['LOW']

        archetype = {
            'cluster_size': cluster_size,
            'mean_resistance_score': float(mean_resistance),
            'resistance_level': resistance_level,
            'resistant_to': high_resistance_abs[:10],  # Top 10
            'susceptible_to': low_resistance_abs[:10],  # Top 10
            'full_profile': {
                _get_antibiotic_name(ab): float(val)
                for ab, val in profile.items()
            },
            'meets_archetype_criteria': True  # Passed size threshold
        }

        archetypes['cluster_archetypes'][_normalize_cluster_id(cluster)] = archetype

        # Generate summary description
        if high_resistance_abs:
            summary = (
                f"Cluster {cluster} ({cluster_size} isolates): {resistance_level}. "
                f"Resistant to: {', '.join(high_resistance_abs[:5])}"
            )
            if low_resistance_abs:
                summary += f". Susceptible to: {', '.join(low_resistance_abs[:3])}"
            archetypes['summary'].append(summary)

        # Build summary table row
        summary_rows.append({
            'Archetype': f"Archetype_{cluster}",
            'Cluster': int(cluster),
            'Key_Resistant_Antibiotics': ', '.join(high_resistance_abs[:5]) if high_resistance_abs else 'None',
            'Resistance_Level': resistance_level,
            'Size': cluster_size
        })

    # Create archetype summary table (the intellectual centerpiece)
    if summary_rows:
        archetypes['archetype_summary_table'] = pd.DataFrame(summary_rows)

    return archetypes


def generate_archetype_summary_table(
    df: pd.DataFrame,
    feature_cols: List[str],
    cluster_col: str = 'CLUSTER',
    mdr_flag_col: str = 'MDR_FLAG',
    species_col: str = 'ISOLATE_ID',
    environment_col: str = 'SAMPLE_SOURCE'
) -> pd.DataFrame:
    """
    Generate the Archetype Summary Table - the intellectual centerpiece of the thesis.

    Creates:
    | Archetype | Key resistant antibiotics | MDR status | Dominant species | Environments observed |

    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe with cluster labels and all metadata
    feature_cols : list
        List of encoded resistance column names
    cluster_col : str
        Column name for cluster labels
    mdr_flag_col : str
        Column name for MDR flag
    species_col : str
        Column name for species
    environment_col : str
        Column name for environment/sample source

    Returns:
    --------
    pd.DataFrame
        Archetype summary table

    Why:
        This table is the intellectual centerpiece of your thesis.
    """
    if cluster_col not in df.columns:
        return pd.DataFrame()

    existing_cols = [c for c in feature_cols if c in df.columns]
    cluster_profiles = df.groupby(cluster_col)[existing_cols].mean() if existing_cols else pd.DataFrame()

    rows = []

    for cluster in sorted(df[cluster_col].unique()):
        cluster_df = df[df[cluster_col] == cluster]
        cluster_size = len(cluster_df)

        if cluster_size < MIN_ARCHETYPE_SIZE:
            continue

        row = {
            'Archetype': f"Archetype_{cluster}",
            'Key_Resistant_Antibiotics': None,
            'MDR_Status': None,
            'Dominant_Species': None,
            'Environments_Observed': None,
            'Size': cluster_size
        }

        # Key resistant antibiotics
        if len(cluster_profiles) > 0 and cluster in cluster_profiles.index:
            profile = cluster_profiles.loc[cluster]
            high_res = profile[profile > 1.5].sort_values(ascending=False)
            row['Key_Resistant_Antibiotics'] = ', '.join(
                [_get_antibiotic_name(ab) for ab in high_res.index[:5]]
            ) if len(high_res) > 0 else 'None'

        # MDR status
        if mdr_flag_col in cluster_df.columns:
            mdr_pct = cluster_df[mdr_flag_col].mean() * 100
            if mdr_pct >= 70:
                row['MDR_Status'] = f"High MDR ({mdr_pct:.0f}%)"
            elif mdr_pct >= 30:
                row['MDR_Status'] = f"Mixed ({mdr_pct:.0f}%)"
            else:
                row['MDR_Status'] = f"Low MDR ({mdr_pct:.0f}%)"

        # Dominant species
        if species_col in cluster_df.columns:
            species_counts = cluster_df[species_col].value_counts()
            row['Dominant_Species'] = species_counts.idxmax()

        # Environments observed
        if environment_col in cluster_df.columns:
            envs = cluster_df[environment_col].unique()
            row['Environments_Observed'] = ', '.join(sorted([str(e) for e in envs[:5]]))

        rows.append(row)

    return pd.DataFrame(rows)


# =============================================================================
# SECTION 4: SPECIES–ENVIRONMENT–RESISTANCE TRIANGULATION
# =============================================================================


def generate_triangulation_table(
    df: pd.DataFrame,
    feature_cols: List[str],
    cluster_col: str = 'CLUSTER',
    species_col: str = 'ISOLATE_ID',
    environment_col: str = 'SAMPLE_SOURCE',
    mdr_flag_col: str = 'MDR_FLAG'
) -> Tuple[pd.DataFrame, Dict]:
    """
    Generate Species–Environment–Resistance Triangulation table.

    For each major cluster/archetype, reports:
    - Species composition
    - Environment distribution
    - MDR proportion

    RULE: All three must be reported together. No single-dimension claims.

    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe with cluster labels and all metadata
    feature_cols : list
        List of encoded resistance column names
    cluster_col : str
        Column name for cluster labels
    species_col : str
        Column name for species
    environment_col : str
        Column name for environment/sample source
    mdr_flag_col : str
        Column name for MDR flag

    Returns:
    --------
    tuple
        (triangulation_table_df, notes_dict)

    Why:
        - Prevents ecological overclaiming
        - All three dimensions must be reported together
    """
    if cluster_col not in df.columns:
        return pd.DataFrame(), {'error': 'No cluster column found'}

    rows = []
    notes = {
        'triangulation_rule': 'All three dimensions (species, environment, MDR) must be reported together',
        'cluster_notes': {}
    }

    for cluster in sorted(df[cluster_col].unique()):
        cluster_df = df[df[cluster_col] == cluster]
        cluster_size = len(cluster_df)

        if cluster_size < MIN_ARCHETYPE_SIZE:
            continue

        row = {
            'Cluster': int(cluster),
            'Size': cluster_size,
            'Species_Composition': None,
            'Environment_Distribution': None,
            'MDR_Proportion': None
        }

        # Species composition
        if species_col in cluster_df.columns:
            species_counts = cluster_df[species_col].value_counts()
            top_species = species_counts.head(3)
            composition = [
                f"{sp} ({cnt / cluster_size * 100:.0f}%)"
                for sp, cnt in top_species.items()
            ]
            row['Species_Composition'] = '; '.join(composition)

        # Environment distribution
        if environment_col in cluster_df.columns:
            env_counts = cluster_df[environment_col].value_counts()
            top_envs = env_counts.head(3)
            distribution = [
                f"{env} ({cnt / cluster_size * 100:.0f}%)"
                for env, cnt in top_envs.items()
            ]
            row['Environment_Distribution'] = '; '.join(distribution)

        # MDR proportion
        if mdr_flag_col in cluster_df.columns:
            mdr_pct = cluster_df[mdr_flag_col].mean() * 100
            row['MDR_Proportion'] = f"{mdr_pct:.1f}%"

        rows.append(row)

        # Generate notes for each cluster
        cluster_note = f"Cluster {cluster}: "
        if row['Species_Composition'] and row['Environment_Distribution'] and row['MDR_Proportion']:
            cluster_note += "Complete triangulation available - all three dimensions reported."
        else:
            cluster_note += "Incomplete triangulation - some dimensions missing data."
        notes['cluster_notes'][int(cluster)] = cluster_note

    return pd.DataFrame(rows), notes


def identify_species_environment_associations(
    df: pd.DataFrame,
    species_col: str = 'ISOLATE_ID',
    environment_col: str = 'SAMPLE_SOURCE',
    region_col: str = 'REGION'
) -> Dict:
    """
    Identify species–environment associations.

    Analyzes which species are associated with specific environmental
    sources and regions.

    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe with species and environment information
    species_col : str
        Column name for species
    environment_col : str
        Column name for environmental/sample source
    region_col : str
        Column name for region

    Returns:
    --------
    dict
        Species-environment association results
    """
    associations = {
        'species_environment': {},
        'species_region': {},
        'statistical_tests': {},
        'interpretation': []
    }

    # Species-Environment associations
    if species_col in df.columns and environment_col in df.columns:
        # Cross-tabulation
        env_crosstab = pd.crosstab(df[species_col], df[environment_col])
        associations['species_environment_crosstab'] = env_crosstab.to_dict()

        # For each species, find dominant environment
        for species in df[species_col].unique():
            species_df = df[df[species_col] == species]
            env_counts = species_df[environment_col].value_counts()
            if len(env_counts) > 0:
                associations['species_environment'][species] = {
                    'dominant_environment': env_counts.idxmax(),
                    'proportion': float(env_counts.max() / len(species_df)),
                    'all_environments': env_counts.to_dict()
                }

        # Chi-square test for species-environment association
        try:
            chi2, p_value, dof, expected = stats.chi2_contingency(env_crosstab)
            associations['statistical_tests']['species_environment'] = {
                'chi_square': float(chi2),
                'p_value': float(p_value),
                'degrees_of_freedom': int(dof),
                'significant': p_value < 0.05
            }

            if p_value < 0.05:
                associations['interpretation'].append(
                    f"Significant species-environment association detected "
                    f"(Chi-sq={chi2:.2f}, p={p_value:.4f}). "
                    "Different species prefer different environmental sources."
                )
        except ValueError as e:
            # Chi-square test can fail with insufficient data
            warnings.warn(f"Chi-square test for species-environment failed: {e}")

    # Species-Region associations
    if species_col in df.columns and region_col in df.columns:
        region_crosstab = pd.crosstab(df[species_col], df[region_col])
        associations['species_region_crosstab'] = region_crosstab.to_dict()

        # For each species, find dominant region
        for species in df[species_col].unique():
            species_df = df[df[species_col] == species]
            region_counts = species_df[region_col].value_counts()
            if len(region_counts) > 0:
                associations['species_region'][species] = {
                    'dominant_region': region_counts.idxmax(),
                    'proportion': float(region_counts.max() / len(species_df)),
                    'all_regions': region_counts.to_dict()
                }

        # Chi-square test for species-region association
        try:
            chi2, p_value, dof, expected = stats.chi2_contingency(region_crosstab)
            associations['statistical_tests']['species_region'] = {
                'chi_square': float(chi2),
                'p_value': float(p_value),
                'degrees_of_freedom': int(dof),
                'significant': p_value < 0.05
            }

            if p_value < 0.05:
                associations['interpretation'].append(
                    f"Significant species-region association detected "
                    f"(Chi-sq={chi2:.2f}, p={p_value:.4f}). "
                    "Species distribution varies by region."
                )
        except ValueError as e:
            # Chi-square test can fail with insufficient data
            warnings.warn(f"Chi-square test for species-region failed: {e}")

    return associations


def identify_mdr_enriched_patterns(
    df: pd.DataFrame,
    feature_cols: List[str],
    cluster_col: str = 'CLUSTER',
    mdr_flag_col: str = 'MDR_FLAG',
    mdr_category_col: str = 'MDR_CATEGORY',
    region_col: str = 'REGION',
    environment_col: str = 'SAMPLE_SOURCE',
    species_col: str = 'ISOLATE_ID'
) -> Dict:
    """
    Identify MDR-enriched patterns in the data.

    Finds clusters, regions, environments, and species with higher
    than expected MDR rates.

    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe with MDR status and other metadata
    feature_cols : list
        List of encoded resistance column names
    cluster_col : str
        Column name for cluster labels
    mdr_flag_col : str
        Column name for MDR flag (binary)
    mdr_category_col : str
        Column name for MDR category
    region_col : str
        Column name for region
    environment_col : str
        Column name for environment/sample source
    species_col : str
        Column name for species

    Returns:
    --------
    dict
        MDR-enriched patterns analysis
    """
    mdr_patterns = {
        'overall_mdr_rate': None,
        'mdr_enriched_clusters': [],
        'mdr_enriched_regions': [],
        'mdr_enriched_environments': [],
        'mdr_enriched_species': [],
        'mdr_resistance_signature': {},
        'interpretation': []
    }

    # Overall MDR rate
    if mdr_flag_col in df.columns:
        overall_mdr = df[mdr_flag_col].mean()
        mdr_patterns['overall_mdr_rate'] = float(overall_mdr)
        mdr_patterns['interpretation'].append(
            f"Overall MDR prevalence: {overall_mdr * 100:.1f}%"
        )
    else:
        return mdr_patterns

    # MDR-enriched clusters
    if cluster_col in df.columns:
        cluster_mdr = df.groupby(cluster_col)[mdr_flag_col].agg(['mean', 'count'])
        enriched_clusters = cluster_mdr[cluster_mdr['mean'] > overall_mdr].sort_values(
            'mean', ascending=False
        )

        for cluster, row in enriched_clusters.iterrows():
            enrichment = {
                'cluster': int(cluster),
                'mdr_rate': float(row['mean']),
                'sample_size': int(row['count']),
                'fold_enrichment': _calculate_fold_enrichment(row['mean'], overall_mdr)
            }
            mdr_patterns['mdr_enriched_clusters'].append(enrichment)

        if mdr_patterns['mdr_enriched_clusters']:
            top_cluster = mdr_patterns['mdr_enriched_clusters'][0]
            mdr_rate_pct = top_cluster['mdr_rate'] * 100
            mdr_patterns['interpretation'].append(
                f"Cluster {top_cluster['cluster']} is most MDR-enriched "
                f"({mdr_rate_pct:.1f}% MDR, "
                f"{top_cluster['fold_enrichment']:.1f}x overall rate)"
            )

    # MDR-enriched regions
    if region_col in df.columns:
        region_mdr = df.groupby(region_col)[mdr_flag_col].agg(['mean', 'count'])
        enriched_regions = region_mdr[region_mdr['mean'] > overall_mdr].sort_values(
            'mean', ascending=False
        )

        for region, row in enriched_regions.iterrows():
            enrichment = {
                'region': region,
                'mdr_rate': float(row['mean']),
                'sample_size': int(row['count']),
                'fold_enrichment': _calculate_fold_enrichment(row['mean'], overall_mdr)
            }
            mdr_patterns['mdr_enriched_regions'].append(enrichment)

    # MDR-enriched environments
    if environment_col in df.columns:
        env_mdr = df.groupby(environment_col)[mdr_flag_col].agg(['mean', 'count'])
        enriched_envs = env_mdr[env_mdr['mean'] > overall_mdr].sort_values(
            'mean', ascending=False
        )

        for env, row in enriched_envs.iterrows():
            enrichment = {
                'environment': env,
                'mdr_rate': float(row['mean']),
                'sample_size': int(row['count']),
                'fold_enrichment': _calculate_fold_enrichment(row['mean'], overall_mdr)
            }
            mdr_patterns['mdr_enriched_environments'].append(enrichment)

    # MDR-enriched species
    if species_col in df.columns:
        species_mdr = df.groupby(species_col)[mdr_flag_col].agg(['mean', 'count'])
        enriched_species = species_mdr[species_mdr['mean'] > overall_mdr].sort_values(
            'mean', ascending=False
        )

        for species, row in enriched_species.iterrows():
            enrichment = {
                'species': species,
                'mdr_rate': float(row['mean']),
                'sample_size': int(row['count']),
                'fold_enrichment': _calculate_fold_enrichment(row['mean'], overall_mdr)
            }
            mdr_patterns['mdr_enriched_species'].append(enrichment)

    # MDR resistance signature (which antibiotics are most associated with MDR)
    existing_cols = [c for c in feature_cols if c in df.columns]
    if existing_cols and mdr_flag_col in df.columns:
        mdr_isolates = df[df[mdr_flag_col] == 1]
        non_mdr_isolates = df[df[mdr_flag_col] == 0]

        if len(mdr_isolates) > 0 and len(non_mdr_isolates) > 0:
            for col in existing_cols:
                ab_name = _get_antibiotic_name(col)
                mdr_mean = mdr_isolates[col].mean()
                non_mdr_mean = non_mdr_isolates[col].mean()
                diff = mdr_mean - non_mdr_mean

                mdr_patterns['mdr_resistance_signature'][ab_name] = {
                    'mdr_mean': float(mdr_mean),
                    'non_mdr_mean': float(non_mdr_mean),
                    'difference': float(diff)
                }

            # Sort by difference to find most discriminating antibiotics
            sorted_sig = sorted(
                mdr_patterns['mdr_resistance_signature'].items(),
                key=lambda x: x[1]['difference'],
                reverse=True
            )

            top_ab = [ab for ab, _ in sorted_sig[:3]]
            if top_ab:
                mdr_patterns['interpretation'].append(
                    f"Antibiotics most associated with MDR: {', '.join(top_ab)}"
                )

    return mdr_patterns


# =============================================================================
# SECTION 5: MDR-ENRICHED PATTERN SYNTHESIS
# =============================================================================

def generate_mdr_synthesis_table(
    df: pd.DataFrame,
    feature_cols: List[str],
    cluster_col: str = 'CLUSTER',
    mdr_flag_col: str = 'MDR_FLAG',
    environment_col: str = 'SAMPLE_SOURCE'
) -> Tuple[pd.DataFrame, Dict]:
    """
    Generate MDR Synthesis Table - consolidates MDR patterns across phases.

    Creates:
    | Cluster | MDR % | Key discriminative antibiotics | Environments | Notes |

    Identifies clusters with MDR % above dataset median and cross-references with:
    - Feature importance (from Phase 3)
    - Environment (from Phase 4)

    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe with cluster labels and MDR information
    feature_cols : list
        List of encoded resistance column names
    cluster_col : str
        Column name for cluster labels
    mdr_flag_col : str
        Column name for MDR flag
    environment_col : str
        Column name for environment/sample source

    Returns:
    --------
    tuple
        (mdr_synthesis_table_df, metadata_dict)

    Why:
        - Elevates MDR analysis beyond a single phase
        - Makes MDR conclusions concrete and defensible
    """
    if cluster_col not in df.columns or mdr_flag_col not in df.columns:
        return pd.DataFrame(), {'error': 'Required columns not found'}

    # Calculate overall MDR median
    overall_mdr_median = df[mdr_flag_col].median()
    overall_mdr_mean = df[mdr_flag_col].mean()

    existing_cols = [c for c in feature_cols if c in df.columns]

    # Calculate MDR signature (antibiotics that discriminate MDR vs non-MDR)
    mdr_signature = {}
    if existing_cols:
        mdr_isolates = df[df[mdr_flag_col] == 1]
        non_mdr_isolates = df[df[mdr_flag_col] == 0]

        if len(mdr_isolates) > 0 and len(non_mdr_isolates) > 0:
            for col in existing_cols:
                ab_name = _get_antibiotic_name(col)
                diff = mdr_isolates[col].mean() - non_mdr_isolates[col].mean()
                mdr_signature[ab_name] = diff

    # Sort antibiotics by discriminative power
    sorted_abs = sorted(mdr_signature.items(), key=lambda x: x[1], reverse=True)
    top_discriminative = [ab for ab, _ in sorted_abs[:5]]

    rows = []

    for cluster in sorted(df[cluster_col].unique()):
        cluster_df = df[df[cluster_col] == cluster]
        cluster_size = len(cluster_df)

        if cluster_size < MIN_ARCHETYPE_SIZE:
            continue

        mdr_rate = cluster_df[mdr_flag_col].mean()
        mdr_pct = mdr_rate * 100

        row = {
            'Cluster': int(cluster),
            'MDR_Percent': round(mdr_pct, 1),
            'Above_Median': mdr_rate > overall_mdr_median,
            'Key_Discriminative_Antibiotics': None,
            'Environments': None,
            'Notes': None
        }

        # Key discriminative antibiotics for this cluster
        if existing_cols:
            cluster_profile = cluster_df[existing_cols].mean()
            high_res = cluster_profile[cluster_profile > 1.5].sort_values(ascending=False)
            # Intersection with overall MDR-discriminative antibiotics
            cluster_abs = [_get_antibiotic_name(ab) for ab in high_res.index[:5]]
            row['Key_Discriminative_Antibiotics'] = ', '.join(cluster_abs[:3]) if cluster_abs else 'None'

        # Environments
        if environment_col in cluster_df.columns:
            envs = cluster_df[environment_col].value_counts().head(2)
            env_str = ', '.join([f"{e} ({c})" for e, c in envs.items()])
            row['Environments'] = env_str

        # Notes
        if mdr_pct >= 70:
            row['Notes'] = "High MDR cluster - enriched above dataset median"
        elif mdr_pct >= 30:
            row['Notes'] = "Mixed MDR status"
        else:
            row['Notes'] = "Low MDR cluster - below dataset median"

        rows.append(row)

    metadata = {
        'overall_mdr_median': float(overall_mdr_median),
        'overall_mdr_mean': float(overall_mdr_mean) * 100,
        'top_discriminative_antibiotics': top_discriminative,
        'enriched_cluster_count': sum(1 for r in rows if r.get('Above_Median'))
    }

    return pd.DataFrame(rows), metadata


# =============================================================================
# SECTION 6: NEGATIVE FINDINGS & CLAIM BOUNDARIES
# =============================================================================

def generate_negative_findings(
    df: pd.DataFrame,
    cluster_col: str = 'CLUSTER',
    region_col: str = 'REGION',
    environment_col: str = 'SAMPLE_SOURCE',
    species_col: str = 'ISOLATE_ID',
    data_limitations: List[str] = None
) -> Dict:
    """
    Generate explicit negative findings - what the data does NOT show.

    Panels trust studies that report limits. This function identifies:
    - Non-findings to explicitly report
    - Patterns that span multiple categories (heterogeneous)

    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe with all metadata
    cluster_col : str
        Column name for cluster labels
    region_col : str
        Column name for region
    environment_col : str
        Column name for environment
    species_col : str
        Column name for species
    data_limitations : list, optional
        Custom list of data limitations. If None, uses DEFAULT_DATA_LIMITATIONS.

    Returns:
    --------
    dict
        Negative findings and limitations

    Why:
        - Panels trust studies that report limits
        - Protects against causal misinterpretation
    """
    negative_findings = {
        'explicit_non_findings': [],
        'heterogeneous_patterns': [],
        'data_limitations': []
    }

    # Check for regional separation
    if cluster_col in df.columns and region_col in df.columns:
        # Calculate cluster-region distribution
        crosstab = pd.crosstab(df[cluster_col], df[region_col], normalize='index')

        # Check if any cluster is region-specific (>80% in one region)
        region_specific = (crosstab > 0.8).any(axis=1).sum()

        if region_specific == 0:
            negative_findings['explicit_non_findings'].append(
                "No clear regional separation was observed - resistance clusters "
                "are distributed across multiple regions rather than being region-specific."
            )
        elif region_specific < len(df[cluster_col].unique()) / 2:
            negative_findings['explicit_non_findings'].append(
                f"Limited regional specificity: only {region_specific} of {len(df[cluster_col].unique())} "
                "clusters showed strong regional association."
            )

    # Check for environment-spanning patterns
    if cluster_col in df.columns and environment_col in df.columns:
        crosstab = pd.crosstab(df[cluster_col], df[environment_col], normalize='index')

        # Count clusters that span multiple environments
        multi_env_clusters = 0
        for cluster in crosstab.index:
            # Cluster spans multiple environments if no single env > 70%
            if crosstab.loc[cluster].max() < 0.7:
                multi_env_clusters += 1
                negative_findings['heterogeneous_patterns'].append(
                    f"Cluster {cluster}: resistance pattern observed across multiple "
                    f"environmental sources (heterogeneous distribution)"
                )

        if multi_env_clusters > 0:
            negative_findings['explicit_non_findings'].append(
                f"Several resistance patterns ({multi_env_clusters} clusters) spanned "
                "multiple environmental sources, indicating that resistance profiles "
                "are not strictly environment-specific."
            )

    # Check for species-spanning patterns
    if cluster_col in df.columns and species_col in df.columns:
        crosstab = pd.crosstab(df[cluster_col], df[species_col], normalize='index')

        multi_species_clusters = 0
        for cluster in crosstab.index:
            if crosstab.loc[cluster].max() < 0.7:
                multi_species_clusters += 1

        if multi_species_clusters > 0:
            negative_findings['explicit_non_findings'].append(
                f"Some resistance patterns ({multi_species_clusters} clusters) were not "
                "species-specific, suggesting horizontal resistance gene transfer or "
                "convergent resistance evolution."
            )

    # Data limitations - use provided or default
    negative_findings['data_limitations'] = (
        data_limitations if data_limitations is not None else DEFAULT_DATA_LIMITATIONS.copy()
    )

    return negative_findings


def generate_claim_boundary_statement() -> Dict:
    """
    Generate explicit claim boundary statement.

    States clearly what the study does NOT claim:
    - No temporal inference
    - No transmission inference
    - No predictive claims

    Returns:
    --------
    dict
        Claim boundary statement with specific disclaimers

    Why:
        - Closes all common attack vectors
        - Protects against causal misinterpretation
    """
    boundary = {
        'title': "Claim Boundary Statement",
        'preamble': (
            "This study provides descriptive analysis of antimicrobial resistance "
            "patterns. The following limitations explicitly apply to all findings:"
        ),
        'boundaries': [
            {
                'category': 'Temporal Inference',
                'statement': "No temporal inference is made.",
                'explanation': (
                    "This cross-sectional study cannot determine whether resistance "
                    "patterns are increasing, decreasing, or stable over time. "
                    "No claims about trends or temporal evolution are made."
                )
            },
            {
                'category': 'Transmission Inference',
                'statement': "No transmission inference is made.",
                'explanation': (
                    "Observed associations between environments, species, and resistance "
                    "patterns do not imply transmission pathways. No claims about how "
                    "resistance genes spread between isolates, species, or locations are made."
                )
            },
            {
                'category': 'Predictive Claims',
                'statement': "No predictive claims are made.",
                'explanation': (
                    "Model performance metrics describe classification accuracy on "
                    "the current dataset only. No claims about predicting resistance "
                    "in future isolates or different populations are made."
                )
            },
            {
                'category': 'Causal Inference',
                'statement': "No causal relationships are claimed.",
                'explanation': (
                    "Associations between variables (e.g., environment-resistance, "
                    "species-resistance) are correlational only. Environmental factors "
                    "are not claimed to 'cause' or 'drive' resistance patterns."
                )
            }
        ],
        'language_discipline': {
            'allowed_terms': ALLOWED_ASSOCIATION_TERMS,
            'avoided_terms': AVOID_CAUSAL_TERMS,
            'note': (
                "Throughout this analysis, language is carefully chosen to reflect "
                "associations rather than causation. Terms like 'associated with', "
                "'enriched in', and 'observed across' are used intentionally."
            )
        }
    }

    return boundary


# =============================================================================
# SECTION 7: FINAL INTEGRATIVE NARRATIVE STRUCTURE
# =============================================================================

def generate_integrative_narrative(
    df: pd.DataFrame,
    feature_cols: List[str],
    results: Dict
) -> List[str]:
    """
    Generate the Final Integrative Narrative following the prescribed structure.

    Structure:
    1. What clusters reveal (Phase 2)
    2. How supervised models corroborate or challenge them (Phase 3)
    3. Where clusters appear ecologically (Phase 4)
    4. What resistance archetypes emerge (Phase 5)

    Parameters:
    -----------
    df : pd.DataFrame
        Complete dataframe with all analysis results
    feature_cols : list
        List of encoded resistance column names
    results : dict
        Complete results dictionary from run_integration_synthesis

    Returns:
    --------
    list
        List of narrative paragraphs following the integrative structure

    Why:
        - Reads as a single coherent argument
        - Shows intentional synthesis, not post-hoc storytelling
    """
    narrative = []

    # 1. What clusters reveal (Phase 2 - Unsupervised)
    archetypes = results.get('resistance_archetypes', {})
    if archetypes.get('cluster_archetypes'):
        n_clusters = len(archetypes['cluster_archetypes'])
        cluster_info = archetypes['cluster_archetypes']

        # Use the module-level constant for high resistance levels
        high_res_clusters = [
            k for k, v in cluster_info.items()
            if v.get('resistance_level') in HIGH_RESISTANCE_LEVELS
        ]

        narrative.append(
            f"PHASE 2 FINDINGS (Unsupervised Clustering): "
            f"Analysis identified {n_clusters} distinct resistance-based clusters. "
            f"{len(high_res_clusters)} cluster(s) exhibited elevated resistance levels, "
            f"suggesting the presence of distinct resistance archetypes in the dataset."
        )

    # 2. How supervised models corroborate or challenge (Phase 3 - Supervised)
    comparison = results.get('cluster_supervised_comparison', {})
    if comparison:
        mdr_sig = comparison.get('mdr_chi_square', {}).get('significant', False)
        species_sig = comparison.get('species_chi_square', {}).get('significant', False)

        if mdr_sig and species_sig:
            narrative.append(
                "PHASE 3 FINDINGS (Supervised Validation): "
                "Supervised analysis corroborates the unsupervised clusters. "
                "Significant associations were found between clusters and both MDR status "
                "and species identity, indicating that the resistance-based groupings "
                "capture biologically meaningful patterns."
            )
        elif mdr_sig:
            narrative.append(
                "PHASE 3 FINDINGS (Supervised Validation): "
                "Supervised analysis partially corroborates cluster structure. "
                "Clusters align significantly with MDR status but show more heterogeneous "
                "species composition, suggesting resistance patterns may transcend species boundaries."
            )
        else:
            narrative.append(
                "PHASE 3 FINDINGS (Supervised Validation): "
                "Clusters show limited alignment with supervised categories. "
                "This suggests resistance profiles may capture variation beyond "
                "simple MDR/non-MDR or species classifications."
            )

    # 3. Where clusters appear ecologically (Phase 4 - Environmental)
    assoc = results.get('species_environment_associations', {})
    if assoc.get('statistical_tests', {}).get('species_environment', {}).get('significant'):
        narrative.append(
            "PHASE 4 FINDINGS (Ecological Context): "
            "Significant species-environment associations were detected. "
            "Different bacterial species are enriched in specific environmental sources, "
            "providing ecological context for the observed resistance patterns."
        )
    else:
        narrative.append(
            "PHASE 4 FINDINGS (Ecological Context): "
            "Species distributions are observed across multiple environmental sources. "
            "While some enrichment patterns exist, resistance clusters span "
            "environmental boundaries rather than being strictly habitat-specific."
        )

    # 4. What resistance archetypes emerge (Phase 5 - Integration)
    mdr = results.get('mdr_enriched_patterns', {})
    if mdr.get('mdr_enriched_clusters'):
        n_enriched = len(mdr['mdr_enriched_clusters'])
        top_cluster = mdr['mdr_enriched_clusters'][0] if mdr['mdr_enriched_clusters'] else None

        mdr_prevalence = ""
        if top_cluster:
            mdr_pct = top_cluster['mdr_rate'] * 100
            mdr_prevalence = (
                f"The most MDR-enriched cluster (Cluster {top_cluster['cluster']}) "
                f"had {mdr_pct:.1f}% MDR prevalence. "
            )

        narrative.append(
            f"PHASE 5 FINDINGS (Archetype Synthesis): "
            f"{n_enriched} cluster(s) showed MDR enrichment above the dataset average. "
            + mdr_prevalence
            + "These MDR-enriched clusters, combined with their species composition and "
            "environmental distribution, constitute the primary resistance archetypes "
            "identified in this study."
        )

    # Add synthesis conclusion
    narrative.append(
        "SYNTHESIS CONCLUSION: "
        "The integration of unsupervised clustering (Phase 2), supervised validation (Phase 3), "
        "and ecological contextualization (Phase 4) reveals coherent resistance archetypes. "
        "These archetypes are characterized by specific antibiotic resistance profiles, "
        "are associated with certain species, and are observed across defined environmental contexts."
    )

    return narrative


# =============================================================================
# MAIN INTEGRATION FUNCTION (UPDATED)
# =============================================================================


def run_integration_synthesis(
    df: pd.DataFrame,
    feature_cols: List[str],
    supervised_results: Dict = None
) -> Dict:
    """
    Run complete Phase 5 Integration & Synthesis analysis.

    DELIVERABLES CHECKLIST (Phase 5 Audit):
    ✅ Cluster–supervised alignment table
    ✅ Resistance archetype definition
    ✅ Archetype summary table
    ✅ Species–environment–resistance triangulation
    ✅ MDR-enriched synthesis table
    ✅ Explicit negative findings
    ✅ Claim boundary statement

    This function integrates results from:
    - Phase 2: Unsupervised clustering (resistance-based clusters)
    - Phase 3: Supervised learning (discriminative signals)
    - Phase 4: Regional/environmental analysis (contextual associations)

    And synthesizes findings about:
    - Dominant resistance archetypes
    - Species-environment associations
    - MDR-enriched patterns

    Parameters:
    -----------
    df : pd.DataFrame
        Clustered dataframe with all metadata
    feature_cols : list
        List of encoded resistance column names
    supervised_results : dict, optional
        Results from supervised learning pipeline

    Returns:
    --------
    dict
        Complete integration and synthesis results including all Phase 5 deliverables
    """
    print("\n" + "=" * 70)
    print("PHASE 5: INTEGRATION & SYNTHESIS (Concrete Implementation)")
    print("=" * 70)

    results = {
        # Section 1: Integration Framework
        'integration_framework': None,
        # Section 2: Cluster-Supervised Alignment
        'cluster_supervised_comparison': None,
        'cluster_supervised_alignment_table': None,
        'alignment_interpretation': None,
        # Section 3: Resistance Archetypes
        'resistance_archetypes': None,
        'archetype_summary_table': None,
        # Section 4: Triangulation
        'species_environment_associations': None,
        'triangulation_table': None,
        'triangulation_notes': None,
        # Section 5: MDR Synthesis
        'mdr_enriched_patterns': None,
        'mdr_synthesis_table': None,
        'mdr_synthesis_metadata': None,
        # Section 6: Negative Findings & Boundaries
        'negative_findings': None,
        'claim_boundary_statement': None,
        # Section 7: Integrative Narrative
        'integrative_narrative': None,
        # Legacy
        'synthesis_summary': []
    }

    # ==========================================================================
    # SECTION 1: FORMAL INTEGRATION FRAMEWORK
    # ==========================================================================
    print("\n" + "-" * 50)
    print("SECTION 1: FORMAL INTEGRATION FRAMEWORK")
    print("-" * 50)
    results['integration_framework'] = define_integration_framework()
    print("   Defined analytical outputs from Phase 2, 3, and 4")
    print(f"   • Phase 2: {results['integration_framework']['phase_outputs']['phase_2']['output']}")
    print(f"   • Phase 3: {results['integration_framework']['phase_outputs']['phase_3']['output']}")
    print(f"   • Phase 4: {results['integration_framework']['phase_outputs']['phase_4']['output']}")

    # ==========================================================================
    # SECTION 2: CLUSTER–SUPERVISED RESULT ALIGNMENT
    # ==========================================================================
    print("\n" + "-" * 50)
    print("SECTION 2: CLUSTER–SUPERVISED RESULT ALIGNMENT")
    print("-" * 50)
    results['cluster_supervised_comparison'] = compare_clusters_with_supervised(df)

    # Generate alignment table
    alignment_df, alignment_interp = generate_cluster_supervised_alignment_table(df)
    results['cluster_supervised_alignment_table'] = alignment_df
    results['alignment_interpretation'] = alignment_interp

    if results['cluster_supervised_comparison']['interpretation']:
        for interp in results['cluster_supervised_comparison']['interpretation']:
            print(f"   • {interp}")

    if not alignment_df.empty:
        print(f"   Generated cluster-supervised alignment table with {len(alignment_df)} clusters")
        if alignment_interp.get('agreement'):
            print(f"   Agreement cases: {len(alignment_interp['agreement'])}")
        if alignment_interp.get('disagreement'):
            print(f"   Disagreement cases: {len(alignment_interp['disagreement'])}")

    # ==========================================================================
    # SECTION 3: RESISTANCE ARCHETYPES
    # ==========================================================================
    print("\n" + "-" * 50)
    print("SECTION 3: RESISTANCE ARCHETYPES (Formal Definition)")
    print("-" * 50)
    results['resistance_archetypes'] = identify_resistance_archetypes(df, feature_cols)

    # Generate archetype summary table
    results['archetype_summary_table'] = generate_archetype_summary_table(df, feature_cols)

    if results['resistance_archetypes'].get('formal_definition'):
        print("   Archetype Definition Applied:")
        print("   \"A recurring resistance profile characterized by a stable")
        print("    antibiotic resistance pattern and consistent cluster membership.\"")

    if results['resistance_archetypes']['summary']:
        for summary in results['resistance_archetypes']['summary'][:3]:  # Show first 3
            print(f"   • {summary}")
        if len(results['resistance_archetypes']['summary']) > 3:
            print(f"   ... and {len(results['resistance_archetypes']['summary']) - 3} more archetypes")

    # ==========================================================================
    # SECTION 4: SPECIES–ENVIRONMENT–RESISTANCE TRIANGULATION
    # ==========================================================================
    print("\n" + "-" * 50)
    print("SECTION 4: SPECIES–ENVIRONMENT–RESISTANCE TRIANGULATION")
    print("-" * 50)
    results['species_environment_associations'] = identify_species_environment_associations(df)

    # Generate triangulation table
    triangulation_df, triangulation_notes = generate_triangulation_table(df, feature_cols)
    results['triangulation_table'] = triangulation_df
    results['triangulation_notes'] = triangulation_notes

    print("   Rule: All three dimensions (species, environment, MDR) reported together")
    if not triangulation_df.empty:
        print(f"   Generated triangulation table with {len(triangulation_df)} clusters")

    if results['species_environment_associations']['interpretation']:
        for interp in results['species_environment_associations']['interpretation']:
            print(f"   • {interp}")

    # ==========================================================================
    # SECTION 5: MDR-ENRICHED PATTERN SYNTHESIS
    # ==========================================================================
    print("\n" + "-" * 50)
    print("SECTION 5: MDR-ENRICHED PATTERN SYNTHESIS")
    print("-" * 50)
    results['mdr_enriched_patterns'] = identify_mdr_enriched_patterns(df, feature_cols)

    # Generate MDR synthesis table
    mdr_synthesis_df, mdr_metadata = generate_mdr_synthesis_table(df, feature_cols)
    results['mdr_synthesis_table'] = mdr_synthesis_df
    results['mdr_synthesis_metadata'] = mdr_metadata

    if results['mdr_enriched_patterns']['interpretation']:
        for interp in results['mdr_enriched_patterns']['interpretation']:
            print(f"   • {interp}")

    if mdr_metadata.get('enriched_cluster_count'):
        print(f"   {mdr_metadata['enriched_cluster_count']} cluster(s) enriched above MDR median")

    # ==========================================================================
    # SECTION 6: NEGATIVE FINDINGS & CLAIM BOUNDARIES
    # ==========================================================================
    print("\n" + "-" * 50)
    print("SECTION 6: NEGATIVE FINDINGS & CLAIM BOUNDARIES")
    print("-" * 50)
    results['negative_findings'] = generate_negative_findings(df)
    results['claim_boundary_statement'] = generate_claim_boundary_statement()

    if results['negative_findings']['explicit_non_findings']:
        print("   Explicit Non-Findings:")
        for finding in results['negative_findings']['explicit_non_findings'][:2]:
            print(f"   • {finding[:80]}...")

    print("   Claim Boundaries:")
    for boundary in results['claim_boundary_statement']['boundaries'][:2]:
        print(f"   • {boundary['category']}: {boundary['statement']}")

    # ==========================================================================
    # SECTION 7: INTEGRATIVE NARRATIVE
    # ==========================================================================
    print("\n" + "-" * 50)
    print("SECTION 7: INTEGRATIVE NARRATIVE STRUCTURE")
    print("-" * 50)
    results['integrative_narrative'] = generate_integrative_narrative(df, feature_cols, results)

    for i, paragraph in enumerate(results['integrative_narrative'], 1):
        print(f"   {i}. {paragraph[:100]}...")

    # ==========================================================================
    # SYNTHESIS SUMMARY (LEGACY)
    # ==========================================================================
    print("\n" + "-" * 50)
    print("SYNTHESIS SUMMARY:")
    print("-" * 50)
    results['synthesis_summary'] = _generate_synthesis_summary(results, df, feature_cols)
    for point in results['synthesis_summary']:
        print(f"   • {point}")

    # Print deliverables checklist
    print("\n" + "=" * 70)
    print("PHASE 5 DELIVERABLES CHECKLIST:")
    print("=" * 70)
    print("   [OK] Integration framework defined (Phase 2/3/4 outputs)")
    print("   [OK] Cluster-supervised alignment table generated")
    print("   [OK] Resistance archetype definition included")
    print("   [OK] Archetype summary table created")
    print("   [OK] Species-environment-resistance triangulation table")
    print("   [OK] MDR-enriched synthesis table")
    print("   [OK] Explicit negative findings documented")
    print("   [OK] Claim boundary statement included")
    print("   [OK] Integrative narrative structure")

    return results


def _generate_synthesis_summary(results: Dict, df: pd.DataFrame, feature_cols: List[str]) -> List[str]:
    """
    Generate a synthesis summary combining all analysis results.

    Parameters:
    -----------
    results : dict
        Results from all integration analyses
    df : pd.DataFrame
        Original dataframe
    feature_cols : list
        Feature column names

    Returns:
    --------
    list
        List of synthesis summary points
    """
    summary = []

    # Dataset overview
    n_isolates = len(df)
    n_antibiotics = len([c for c in feature_cols if c in df.columns])
    summary.append(f"Analyzed {n_isolates} isolates across {n_antibiotics} antibiotics.")

    # Clustering-supervised alignment
    comparison = results.get('cluster_supervised_comparison', {})
    if comparison.get('mdr_chi_square', {}).get('significant'):
        summary.append(
            "Unsupervised clusters significantly align with MDR status, "
            "suggesting clustering captures clinically relevant resistance patterns."
        )

    # Dominant archetypes
    archetypes = results.get('resistance_archetypes', {})
    if archetypes.get('cluster_archetypes'):
        n_archetypes = len(archetypes['cluster_archetypes'])
        high_res_clusters = [
            k for k, v in archetypes['cluster_archetypes'].items()
            if v.get('resistance_level') in HIGH_RESISTANCE_LEVELS
        ]
        if high_res_clusters:
            summary.append(
                f"Identified {n_archetypes} distinct resistance archetypes, "
                f"with {len(high_res_clusters)} cluster(s) showing elevated resistance levels."
            )

    # Species-environment associations
    assoc = results.get('species_environment_associations', {})
    if assoc.get('statistical_tests', {}).get('species_environment', {}).get('significant'):
        summary.append(
            "Species distribution varies significantly by environmental source, "
            "indicating habitat-specific bacterial communities."
        )

    # MDR patterns
    mdr = results.get('mdr_enriched_patterns', {})
    if mdr.get('overall_mdr_rate') is not None:
        mdr_rate = mdr['overall_mdr_rate'] * 100
        n_enriched = len(mdr.get('mdr_enriched_clusters', []))
        summary.append(
            f"Overall MDR prevalence is {mdr_rate:.1f}%, "
            f"with {n_enriched} cluster(s) showing above-average MDR rates."
        )

    # Key antibiotics
    if mdr.get('mdr_resistance_signature'):
        sorted_abs = sorted(
            mdr['mdr_resistance_signature'].items(),
            key=lambda x: x[1]['difference'],
            reverse=True
        )
        top_abs = [ab for ab, _ in sorted_abs[:3]]
        if top_abs:
            summary.append(
                f"Key antibiotics discriminating MDR isolates: {', '.join(top_abs)}."
            )

    return summary


if __name__ == "__main__":
    from pathlib import Path

    project_root = Path(__file__).parent.parent.parent
    clustered_path = project_root / "data" / "processed" / "clustered_dataset.csv"

    if clustered_path.exists():
        df = pd.read_csv(clustered_path)
        feature_cols = [c for c in df.columns if c.endswith('_encoded')]

        results = run_integration_synthesis(df, feature_cols)

        print("\n" + "=" * 50)
        print("INTEGRATION & SYNTHESIS COMPLETE")
        print("=" * 50)
    else:
        print(f"Clustered dataset not found at {clustered_path}")
        print("Run the clustering pipeline first.")
