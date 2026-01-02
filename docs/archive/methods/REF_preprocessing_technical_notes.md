# Preprocessing Methods Documentation

## Phase 2: Data Preprocessing

This document provides comprehensive documentation of all preprocessing decisions, rules, and justifications for the AMR thesis project.

---

## Table of Contents

1. [Decision Log](#1-decision-log)
2. [Data Ingestion (Phase 2.1)](#2-data-ingestion-phase-21)
3. [Data Cleaning (Phase 2.2-2.3)](#3-data-cleaning-phase-22-23)
4. [Resistance Encoding (Phase 2.4)](#4-resistance-encoding-phase-24)
5. [Feature Engineering (Phase 2.5)](#5-feature-engineering-phase-25)
6. [Exclusion Summary](#6-exclusion-summary)

---

## 1. Decision Log

### 1.1 Preprocessing Decision Summary Table

| Decision | Rule Applied | Justification |
|----------|--------------|---------------|
| **Antibiotic inclusion** | ≥70% coverage | Ensures comparability across isolates; antibiotics tested in fewer samples lack statistical power |
| **MDR definition** | ≥3 antibiotic classes | Standard AMR definition per Magiorakos et al. (2012) |
| **MAR Index formula** | Resistant / Tested | Krumperman (1983) standard formula |
| **Missing data per isolate** | ≤30% threshold | Maintains data quality while preserving sample size |
| **Duplicate handling** | Remove by CODE | Keep first occurrence; log all removals |
| **Resistance encoding** | S=0, I=1, R=2 | Ordinal encoding preserves biological meaning |
| **Species standardization** | Controlled vocabulary | Ensures consistent naming across sources |
| **Invalid values** | Set to NULL | Only {S, I, R} allowed; multi-label entries rejected |

### 1.2 Threshold Justification

| Threshold | Value | Rationale | Alternative Considered |
|-----------|-------|-----------|------------------------|
| Antibiotic coverage | 70% | Balance between data completeness and sample retention | 50% (too permissive), 90% (too restrictive) |
| Isolate missing data | 30% | Ensures sufficient data per isolate for pattern discrimination | 50% (lower quality), 20% (excessive data loss) |

---

## 2. Data Ingestion (Phase 2.1)

### 2.1 Objectives

- Load and merge AST data from multiple CSV files
- Extract metadata from filenames and isolate codes
- Standardize antibiotic abbreviations across sources
- Create a unified raw dataset with validated metadata

### 2.2 Required Metadata Columns

| Column | Description | Source | Validation Rule |
|--------|-------------|--------|-----------------|
| **REGION** | Geographic region | Extracted from filename | Must be non-empty |
| **SITE** | Sampling site | Extracted from filename | Must be non-empty |
| **ENVIRONMENT** | Environmental category | Derived from sampling source | Must be {Water, Fish, Hospital} |
| **SAMPLING_SOURCE** | Detailed sampling source | Parsed from isolate code | Must match controlled vocabulary |

### 2.3 Isolate Code Parsing Rules

```
Format: [Species Prefix]_[National Site][Local Site][Sample Source][Replicate][Colony]

Example: EC_OADWR1C3
├── EC    = Escherichia coli (Species Prefix)
├── O     = Ormoc (National Site)
├── A     = Alegria (Local Site)
├── DW    = Drinking Water (Sample Source)
├── R1    = Replicate 1
└── C3    = Colony 3
```

### 2.4 Environment Categorization Rules

| Sampling Source | Environment Category |
|-----------------|---------------------|
| Drinking Water (DW), Lake Water (LW), River Water (RW) | Water |
| Fish Banak (FB), Fish Gusaw (FG), Fish Tilapia (FT), Fish Kaolang (FK) | Fish |
| Effluent Water Untreated (EWU), Effluent Water Treated (EWT) | Hospital |

### 2.5 Output

- **File**: `unified_raw_dataset.csv`
- **Contents**: Consolidated dataset with all isolates and validated metadata

---

## 3. Data Cleaning (Phase 2.2-2.3)

### 3.1 Validation Rules

| Rule | Description | Action | Logged |
|------|-------------|--------|--------|
| Valid resistance values | Only {S, I, R} allowed | Invalid values → NULL | Yes |
| No multi-label entries | No "S/R", "S,I", etc. | Invalid values → NULL | Yes |
| Species name variants | Standardize to controlled vocabulary | Map to standard name | Yes |
| Duplicate isolates | Same CODE value | Remove, keep first | Yes |

### 3.2 Species Name Standardization (Controlled Vocabulary)

> **Note**: After species standardization, isolates are reduced from 492 to 491 (1 excluded), and unique species from 10 to 6.

#### Consolidation Rules

| Original Species Name | Standardized Name | Action |
|----------------------|-------------------|--------|
| Klebsiella pneumoniae ssp ozaenae | Klebsiella pneumoniae | Subspecies consolidated |
| Klebsiella pneumoniae spp ozaenae | Klebsiella pneumoniae | Typo variant consolidated |
| Salmonella enterica ssp diarizonae | Salmonella species | Subspecies consolidated |
| Salmonella group | Salmonella species | Generic consolidated |
| Escherichia coli | Escherichia coli | Retained |
| Enterobacter cloacae | Enterobacter cloacae | Retained |
| Enterobacter aerogenes | Enterobacter aerogenes | Retained |
| Vibrio vulnificus | Vibrio vulnificus | Retained |

#### Exclusion Rules

| Excluded Species | Reason | Count |
|------------------|--------|-------|
| Enterobacter cloacae complex/E COLI | Ambiguous hybrid identification | 1 |

#### Final Species List (6 species)

| Species | Count | Percentage |
|---------|-------|------------|
| Escherichia coli | 227 | 46.2% |
| Klebsiella pneumoniae | 149 | 30.3% |
| Enterobacter cloacae | 68 | 13.8% |
| Enterobacter aerogenes | 23 | 4.7% |
| Salmonella species | 23 | 4.7% |
| Vibrio vulnificus | 1 | 0.2% |

### 3.3 Missing Data Strategy

| Parameter | Threshold | Rationale |
|-----------|-----------|-----------|
| **Minimum antibiotic coverage** | ≥70% | Antibiotics tested in fewer than 70% of isolates excluded |
| **Maximum missing per isolate** | ≤30% | Isolates with >30% missing AST values excluded |

**Procedure**:
1. Compute antibiotic test coverage (% isolates tested per antibiotic)
2. Retain antibiotics meeting ≥70% threshold
3. Remove isolates exceeding 30% missing-value threshold
4. Document all exclusions in cleaning report

### 3.4 Duplicate Detection and Resolution

**Strategy**: Exact match on CODE column
**Action**: Remove duplicates, retain first occurrence
**Documentation**: Log all duplicate removals with index and CODE

### 3.5 Output

- **File**: `cleaned_dataset.csv`
- **Report**: `cleaning_report.txt` containing:
  - Thresholds applied
  - Data retention summary
  - Validation summary
  - Antibiotic test coverage table
  - Exclusion summary table
  - Cleaning actions log

---

## 4. Resistance Encoding (Phase 2.4)

### 4.1 Encoding Scheme

| Original Value | Encoded Value | Biological Interpretation |
|----------------|---------------|---------------------------|
| S (Susceptible) | 0 | No resistance detected |
| I (Intermediate) | 1 | Intermediate resistance |
| R (Resistant) | 2 | Full resistance |

### 4.2 Rationale for Ordinal Encoding

Ordinal encoding preserves the biological meaning of resistance levels:
- **Enables meaningful distance calculations** for hierarchical clustering
- **Supports numerical operations** for supervised learning
- **Maintains interpretability** in resistance indices

### 4.3 Output

- **File**: `encoded_dataset.csv`
- **New columns**: `{ANTIBIOTIC}_encoded` for each antibiotic

---

## 5. Feature Engineering (Phase 2.5)

### 5.1 MAR Index (Multiple Antibiotic Resistance Index)

**Formula**:
```
MAR = a / b

Where:
  a = Number of antibiotics to which the isolate is resistant (R)
  b = Total number of antibiotics tested on the isolate
```

**Reference**: Krumperman PH. (1983). Multiple antibiotic resistance indexing of *Escherichia coli* to identify high-risk sources of fecal contamination of foods. *Applied and Environmental Microbiology*, 46(1), 165-170.

**Interpretation**:
| MAR Value | Interpretation |
|-----------|----------------|
| > 0.2 | High-risk contamination source |
| = 0 | Fully susceptible isolate |
| = 1 | Pan-resistant isolate |

### 5.2 MDR Classification

**Definition**: Multi-Drug Resistant (MDR) if resistant to ≥3 antimicrobial categories.

**Reference**: Magiorakos AP, et al. (2012). Multidrug-resistant, extensively drug-resistant and pandrug-resistant bacteria. *Clinical Microbiology and Infection*, 18(3), 268-281.

**Antimicrobial Categories for MDR Calculation**:

| Category | Antibiotics |
|----------|-------------|
| Penicillins | AM, AMP |
| β-lactam/BLI combinations | AMC, PRA |
| Cephalosporins (1st generation) | CN, CF |
| Cephalosporins (3rd/4th generation) | CPD, CTX, CFT, CPT |
| Cephamycins | CFO |
| Cephalosporin/BLI combinations | CZA |
| Carbapenems | IPM, MRB |
| Aminoglycosides | AN, GM, N |
| Quinolones/Fluoroquinolones | NAL, ENR |
| Tetracyclines | DO, TE |
| Nitrofurans | FT |
| Phenicols | C |
| Folate pathway inhibitors | SXT |

### 5.3 Derived Features Summary

| Feature | Formula | Description |
|---------|---------|-------------|
| **MAR_INDEX_COMPUTED** | a / b (Krumperman, 1983) | Multiple Antibiotic Resistance index (0-1) |
| **RESISTANCE_COUNT** | Count where encoded = 2 | Total number of resistant antibiotics |
| **RESISTANT_CLASSES_COUNT** | Count of unique resistant classes | Number of antimicrobial categories with resistance |
| **MDR_FLAG** | Boolean: Classes ≥ 3 (Magiorakos, 2012) | Multi-Drug Resistant indicator |
| **MDR_CATEGORY** | "MDR" or "Non-MDR" | Categorical MDR status |
| **{AB}_RESISTANT** | Binary: 1 if R, 0 if S/I | Per-antibiotic binary resistance indicator |

### 5.4 Output

- **File**: `analysis_ready_dataset.csv`
- **Additional outputs**:
  - `feature_matrix_X.csv`: Encoded resistance values only (for clustering/ML)
  - `metadata.csv`: Sample identification and derived features

---

## 6. Exclusion Summary

### 6.1 Exclusion Summary Table Template

| Exclusion Type | Count | Percentage | Reason |
|----------------|-------|------------|--------|
| Low-coverage antibiotics | [N] | [%] | <70% coverage |
| High-missing isolates | [N] | [%] | >30% missing data |
| Invalid resistance values | [N] | [%] | Not in {S, I, R} |
| Duplicate isolates | [N] | [%] | Same CODE value |

### 6.2 Data Retention Summary Template

| Stage | Records In | Records Out | Retention Rate |
|-------|------------|-------------|----------------|
| Raw ingestion | [N] | [N] | 100% |
| After cleaning | [N] | [N] | [%] |
| After encoding | [N] | [N] | [%] |
| Analysis-ready | [N] | [N] | [%] |

---

## Document Version

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2024 | Initial preprocessing documentation |

---

*This document is part of Phase 8 — Documentation & Reporting for the AMR Thesis Project.*

**See also**: [METHODOLOGY.md](../METHODOLOGY.md) | [DOCUMENTATION.md](../DOCUMENTATION.md)
