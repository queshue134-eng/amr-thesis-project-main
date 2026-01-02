# Study Limitations

> **AMR Thesis Project — Study Limitations**  
> **Last Updated:** December 28, 2025  
> **Data Source:** INOHAC AMR Project Two

---

## Table of Contents

1. [Scope Boundaries](#1-scope-boundaries)
2. [Statistical Caveats](#2-statistical-caveats)
3. [Missing Data Handling](#3-missing-data-handling)
4. [Generalizability Constraints](#4-generalizability-constraints)
5. [Methodological Considerations](#5-methodological-considerations)

---

## 1. Scope Boundaries

| Boundary             | Description                                                                      |
| -------------------- | -------------------------------------------------------------------------------- |
| **Research Purpose** | Exploratory pattern recognition only — NOT clinical decision support             |
| **Temporal Scope**   | Single snapshot — no longitudinal/temporal data available                        |
| **Geographic Scope** | Limited to 3 Philippine regions (Ormoc, Pampanga, Marawi/BARMM)                  |
| **Sample Size**      | 491 isolates, 6 bacterial species (after cleaning)                               |
| **Data Type**        | Environmental and hospital isolates — not representative of clinical populations |

### Species Coverage

| Species                | Count | Percentage |
| ---------------------- | ----- | ---------- |
| Escherichia coli       | 227   | 46.2%      |
| Klebsiella pneumoniae  | 149   | 30.3%      |
| Enterobacter cloacae   | 68    | 13.8%      |
| Enterobacter aerogenes | 23    | 4.7%       |
| Salmonella species     | 23    | 4.7%       |
| Vibrio vulnificus      | 1     | 0.2%       |

> **Note:** _Vibrio vulnificus_ is represented by only 1 isolate, making species-specific conclusions unreliable.

---

## 2. Statistical Caveats

### Clustering Quality

| Metric                 | Value         | Interpretation                  |
| ---------------------- | ------------- | ------------------------------- |
| Silhouette Score (k=4) | 0.466         | Moderate structure (not strong) |
| Rousseeuw Range        | 0.26 – 0.50   | Weak-to-moderate clustering     |
| Cluster Boundaries     | Probabilistic | Not deterministic separations   |

### Silhouette Interpretation Reference (Rousseeuw, 1987)

| Range           | Interpretation                           |
| --------------- | ---------------------------------------- |
| 0.71 – 1.00     | Strong structure                         |
| 0.51 – 0.70     | Reasonable structure                     |
| **0.26 – 0.50** | **Weak-moderate structure** ← This study |
| ≤ 0.25          | No meaningful structure                  |

### Sample Size Per Cluster (k=4)

| Cluster      | n   | %     | Statistical Adequacy           |
| ------------ | --- | ----- | ------------------------------ |
| 1            | 23  | 4.7%  | Small — interpret with caution |
| 2            | 93  | 18.9% | Adequate                       |
| 3            | 123 | 25.1% | Adequate                       |
| 4            | 252 | 51.3% | Adequate (dominant cluster)    |
| **Smallest** | 23  | —     | Below ideal minimum of 30      |

---

## 3. Missing Data Handling

### Single Imputation Caveat

> **IMPORTANT:** This pipeline uses **median imputation** for missing resistance values after threshold-based exclusion. Single imputation has inherent limitations:

| Limitation                      | Consequence                                                                                                                                             |
| ------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Uncertainty Underestimation** | Imputed values treated as observed; true values unknown. Standard errors and confidence intervals are artificially deflated.                            |
| **Variance Reduction**          | All missing values for an antibiotic receive identical median value, reducing variability and potentially making clusters appear more cohesive.         |
| **MAR Assumption**              | Median imputation assumes Missing At Random (MAR). AST missingness may be informative (selective testing based on suspected resistance), violating MAR. |

### Why Acceptable for This Study

1. **Conservative thresholds** (70% antibiotic coverage, 30% max isolate missing) exclude severely sparse data before imputation
2. **Exploratory purpose** — pattern recognition, not clinical decision-making
3. **Sensitivity robustness** — alternative imputation methods yielded consistent cluster structures (not shown)

### Recommended Improvement

For confirmatory analyses, consider **Multiple Imputation (MICE)** to properly propagate uncertainty.

---

## 4. Generalizability Constraints

| Constraint                         | Implication                                                                      |
| ---------------------------------- | -------------------------------------------------------------------------------- |
| **Water-Fish-Human Nexus Context** | Findings specific to environmental surveillance of aquaculture and water sources |
| **Regional Specificity**           | Patterns may not transfer to other Philippine regions or countries               |
| **Non-Clinical Isolates**          | Hospital effluent samples differ from clinical patient isolates                  |
| **Temporal Snapshot**              | No ability to assess trends, seasonality, or epidemic dynamics                   |

### Transfer Limitations

These findings **should NOT** be used to:

- Make clinical treatment decisions
- Predict resistance patterns in clinical settings
- Establish national antimicrobial policy
- Replace laboratory susceptibility testing

---

## 5. Methodological Considerations

### Encoding Simplification

| Concern                    | Details                                                                                            |
| -------------------------- | -------------------------------------------------------------------------------------------------- |
| **Ordinal Encoding**       | S=0, I=1, R=2 treats intermediate as exactly halfway between susceptible and resistant             |
| **Biological Reality**     | Clinical breakpoints are antibiotic-specific; some I interpretations may be clinically susceptible |
| **Alternative Approaches** | Binary (S vs. R+I), one-hot encoding — not explored in this study                                  |

### MDR Definition

- Uses Magiorakos et al. (2012) threshold: Resistant to ≥3 antibiotic classes
- Classification is based on phenotype only, not genotype
- Intermediate results (I) are counted as non-resistant for MDR calculation

### Clustering Algorithm Selection

- Ward's linkage chosen for variance minimization
- Alternative methods (complete, average, single linkage) may yield different groupings
- Distance metric (Euclidean) treats ordinal data as continuous

### Supervised Learning Caveats

| Concern                   | Note                                                                               |
| ------------------------- | ---------------------------------------------------------------------------------- |
| **Small Species Classes** | _Vibrio vulnificus_ (n=1) and _Salmonella_ (n=23) have limited training data       |
| **Cross-Validation**      | 5-fold stratified CV used, but may still overestimate performance for rare species |
| **Feature Importance**    | Random Forest importance reflects correlation, not causation                       |

### Species-Agnostic MDR Classification

> **IMPORTANT LIMITATION:** The Magiorakos et al. (2012) framework specifies MDR definitions for each organism SEPARATELY, accounting for intrinsic resistances.

| Species                  | Intrinsic Resistance                         | Impact on MDR        |
| ------------------------ | -------------------------------------------- | -------------------- |
| _Pseudomonas aeruginosa_ | Penicillins (AmpC β-lactamase)               | May OVERESTIMATE MDR |
| _Acinetobacter spp._     | Aminopenicillins, 1st/2nd gen cephalosporins | May OVERESTIMATE MDR |
| _E. coli_, _Klebsiella_  | No significant intrinsic resistances         | Comparable rates     |

**This pipeline uses a UNIVERSAL antibiotic class mapping applied to ALL species equally.** For clinical applications, species-specific class mappings should be implemented.

**References:**

- EUCAST Expert Rules 3.2 (2020): Intrinsic Resistance and Unusual Phenotypes
- Poirel et al. (2018). Clin Microbiol Rev, 31(2):e00088-17

---

## Disclaimer

> ⚠️ **This tool is intended for exploratory pattern recognition and surveillance analysis only.**
>
> - It should NOT be used for clinical decision support
> - No patient-level identifiers are processed
> - Findings are hypothesis-generating, not confirmatory

---

## Document Version

| Version | Date              | Changes                                                                           |
| ------- | ----------------- | --------------------------------------------------------------------------------- |
| 1.1     | December 28, 2025 | Updated cluster sizes to k=4 actual values; added species-agnostic MDR limitation |
| 1.0     | December 20, 2025 | Initial limitations documentation                                                 |

---

_This document is part of Phase 8 — Documentation & Reporting for the AMR Thesis Project._

**See also**: [TECHNICAL_REFERENCE.md](TECHNICAL_REFERENCE.md) | [README.md](../README.md)
