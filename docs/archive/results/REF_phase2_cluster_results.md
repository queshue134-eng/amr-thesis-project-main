# Phase 2 Results: Data Preprocessing

## Preprocessing and Data Preparation Results

This document provides the actual Phase 2 preprocessing results for the AMR thesis project.

---

## Table of Contents

1. [Data Ingestion Results](#1-data-ingestion-results)
2. [Data Cleaning Results](#2-data-cleaning-results)
3. [Encoding Results](#3-encoding-results)
4. [Feature Engineering Results](#4-feature-engineering-results)
5. [Final Dataset Summary](#5-final-dataset-summary)

---

## 1. Data Ingestion Results

### 1.1 Source Files Processed

| Region | Isolates Loaded | Percentage |
|--------|-----------------|------------|
| BARMM | 249 | 50.7% |
| Region III - Central Luzon | 140 | 28.5% |
| Region VIII - Eastern Visayas | 102 | 20.8% |
| **Total** | **491** | **100%** |

> **Note**: 1 isolate excluded during species standardization (ambiguous hybrid entry: *Enterobacter cloacae complex/E COLI*).

### 1.2 Metadata Extraction Summary

| Metadata Field | Coverage | Notes |
|----------------|----------|-------|
| REGION | 100% | All isolates have region |
| SITE | 100% | Extracted from isolate codes |
| ENVIRONMENT | 100% | Mapped from sampling source |
| SAMPLING_SOURCE | 100% | 9 distinct sources |
| SPECIES (ISOLATE_ID) | 100% | 6 species identified (after standardization) |

### 1.3 Antibiotic Columns Identified

| Category | Antibiotics | Count |
|----------|-------------|-------|
| Retained | AM, AMC, AN, C, CF, CFO, CFT, CN, CPD, CTX, CZA, DO, ENR, FT, GM, IPM, MRB, N, PRA, SXT, TE | 21 |
| Excluded (low coverage) | AMI, CFA, CFV, CPT, CTF, GEN, IME, MAR | 8 |
| **Total tested** | - | **29** |

---

## 2. Data Cleaning Results

### 2.1 Cleaning Parameters Applied

| Parameter | Value Applied | Justification |
|-----------|---------------|---------------|
| Minimum antibiotic coverage | 70% | Balances data retention with completeness |
| Maximum missing per isolate | 30% | Ensures reliable resistance profiles |
| Duplicate handling | Remove by CODE | Prevents replicate bias |

### 2.2 Data Retention Summary

| Stage | Records In | Records Out | Retention Rate |
|-------|------------|-------------|----------------|
| Raw ingestion | 584 | 584 | 100% |
| After duplicate removal | 584 | 582 | 99.7% |
| After antibiotic filtering | 582 | 582 | 100% |
| After isolate filtering | 582 | 492 | 84.5% |
| After species standardization | 492 | **491** | 99.8% |
| **Final cleaned dataset** | - | **491** | **84.0%** |

### 2.3 Exclusion Summary

| Exclusion Type | Count | Percentage | Details |
|----------------|-------|------------|---------|
| Low-coverage antibiotics | 8 | 27.6% of tested | AMI, CFA, CFV, CPT, CTF, GEN, IME, MAR |
| High-missing isolates | 90 | 15.4% | Exceeded 30% missing threshold |
| Invalid resistance values | 0 | 0% | All values valid S/I/R |
| Duplicate isolates | 2 | 0.3% | Removed by CODE |

### 2.4 Antibiotic Coverage Table

| Antibiotic | Coverage (%) | Retained? |
|------------|--------------|-----------|
| GM | 89.9% | ✅ Yes |
| IPM | 89.9% | ✅ Yes |
| MRB | 89.7% | ✅ Yes |
| AN | 89.5% | ✅ Yes |
| CFO | 89.3% | ✅ Yes |
| C | 89.0% | ✅ Yes |
| TE | 89.0% | ✅ Yes |
| SXT | 89.0% | ✅ Yes |
| AMC | 86.9% | ✅ Yes |
| DO | 85.4% | ✅ Yes |
| AM | 73.2% | ✅ Yes |
| CPT | 69.1% | ❌ No |
| AMI, CFA, CFV, CTF, GEN, IME, MAR | 0.2% | ❌ No |

### 2.5 Species Distribution (After Standardization)

| Species | Count | Percentage |
|---------|-------|------------|
| Escherichia coli | 227 | 46.2% |
| Klebsiella pneumoniae | 149 | 30.3% |
| Enterobacter cloacae | 68 | 13.8% |
| Enterobacter aerogenes | 23 | 4.7% |
| Salmonella species | 23 | 4.7% |
| Vibrio vulnificus | 1 | 0.2% |
| **Total** | **491** | **100%** |

> **Species Standardization Applied**:
> - *Klebsiella pneumoniae ssp/spp ozaenae* → *Klebsiella pneumoniae* (3 merged)
> - *Salmonella enterica ssp diarizonae* + *Salmonella group* → *Salmonella species* (2 merged)
> - *Enterobacter cloacae complex/E COLI* → Excluded (1 ambiguous hybrid)

---

## 3. Encoding Results

### 3.1 Encoding Scheme Applied

| Original | Encoded | Count | Percentage |
|----------|---------|-------|------------|
| S (Susceptible) | 0 | ~7,200 | ~66% |
| I (Intermediate) | 1 | ~300 | ~3% |
| R (Resistant) | 2 | ~3,300 | ~31% |

*Note: Approximate counts across all 22 encoded antibiotic columns (491 × 22 = 10,802 cells, minus missing)*

### 3.2 Example Encoded Columns

| Column Name | S (0) | I (1) | R (2) |
|-------------|-------|-------|-------|
| AM_encoded | 167 | 2 | 231 |
| TE_encoded | 366 | 6 | 119 |
| DO_encoded | 363 | 13 | 106 |

**Total encoded columns**: 22

---

## 4. Feature Engineering Results

### 4.1 MAR Index Distribution

| Statistic | Value |
|-----------|-------|
| Mean | 0.1020 |
| Median | 0.0476 |
| Std Dev | 0.1153 |
| Min | 0.0000 |
| Max | 0.6190 |

**Interpretation**: Mean MAR = 0.102 indicates average 10.2% resistant across tested antibiotics. Values > 0.2 indicate high resistance.

### 4.2 MDR Classification Results

| MDR Status | Count | Percentage |
|------------|-------|------------|
| MDR (≥3 classes) | 94 | 19.1% |
| Non-MDR (<3 classes) | 397 | 80.9% |
| **Total** | **491** | **100%** |

### 4.3 Resistant Classes Distribution

| Classes Resistant | Isolates | Percentage | MDR? |
|-------------------|----------|------------|------|
| 0 | ~200 | ~40% | No |
| 1-2 | ~200 | ~40% | No |
| 3+ | 94 | 19.1% | Yes |

---

## 5. Final Dataset Summary

### 5.1 Analysis-Ready Dataset Characteristics

| Characteristic | Value |
|----------------|-------|
| Total isolates | 491 |
| Total species | 6 |
| Total antibiotics (retained) | 22 |
| Encoded columns | 22 |
| Metadata columns | 8 |
| MDR prevalence | 19.1% |
| Mean MAR index | 0.1019 |

### 5.2 Regional Distribution

| Region | Isolates | Percentage |
|--------|----------|------------|
| BARMM | 249 | 50.7% |
| Region III - Central Luzon | 140 | 28.5% |
| Region VIII - Eastern Visayas | 102 | 20.8% |
| **Total** | **491** | **100%** |

### 5.3 Environment Distribution

| Environment | Isolates | Percentage |
|-------------|----------|------------|
| Fish | 274 | 55.8% |
| Water | 176 | 35.8% |
| Hospital | 41 | 8.4% |
| **Total** | **491** | **100%** |

### 5.4 Sampling Source Details

| Sampling Source | Count | Percentage |
|-----------------|-------|------------|
| Fish Tilapia | 115 | 23.4% |
| Drinking Water | 82 | 16.7% |
| Fish Gusaw | 70 | 14.2% |
| River Water | 56 | 11.4% |
| Fish Kaolang | 52 | 10.6% |
| Lake Water | 39 | 7.9% |
| Fish Banak | 37 | 7.5% |
| Effluent Water (Untreated) | 36 | 7.3% |
| Effluent Water (Treated) | 5 | 1.0% |

### 5.5 Output Files Generated

| File | Description | Location |
|------|-------------|----------|
| unified_raw_dataset.csv | Raw consolidated data | data/processed/ |
| cleaned_dataset.csv | Cleaned and standardized | data/processed/ |
| cleaning_report.txt | Cleaning documentation | data/processed/ |
| encoded_dataset.csv | Numerically encoded | data/processed/ |
| analysis_ready_dataset.csv | Final with features | data/processed/ |
| feature_matrix_X.csv | Features only | data/processed/ |
| metadata.csv | Metadata only | data/processed/ |

---

## Interpretation Notes

> **Scope**: These preprocessing results describe the data preparation steps performed on 491 environmental bacterial isolates (6 species, after standardization) from three Philippine regions. Decisions regarding exclusion thresholds (70%/30%), encoding schemes (S=0, I=1, R=2), and species consolidation follow established standards.

---

## Document Version

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2024 | Initial results template |
| 2.0 | 2025-12-18 | **Populated with actual pipeline results** |

---

*This document is part of Phase 8 — Documentation & Reporting for the AMR Thesis Project.*

**See also**: [preprocessing.md](../methods/preprocessing.md) | [phase3_discrimination.md](phase3_discrimination.md)
