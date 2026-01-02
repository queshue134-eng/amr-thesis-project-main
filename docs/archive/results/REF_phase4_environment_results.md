# Phase 4 Results: Environmental and Regional Analysis

## Regional Distribution and Multivariate Analysis Results

This document provides the actual Phase 4 (regional and environmental analysis) results for the AMR thesis project.

---

## Table of Contents

1. [Regional Distribution Analysis](#1-regional-distribution-analysis)
2. [Environmental Distribution Analysis](#2-environmental-distribution-analysis)
3. [Principal Component Analysis Results](#3-principal-component-analysis-results)
4. [Statistical Associations](#4-statistical-associations)

---

## 1. Regional Distribution Analysis

> **Scope Statement**: Regional associations are **observational** and reflect patterns in the analyzed dataset. Results do not establish causation.

### 1.1 Cluster Distribution by Region

| Cluster | BARMM | Region III - Central Luzon | Region VIII - Eastern Visayas | Total |
|---------|-------|---------------------------|-------------------------------|-------|
| 1 | 2 (8.7%) | 17 (73.9%) | 4 (17.4%) | 23 |
| 2 | 39 (41.9%) | 49 (52.7%) | 5 (5.4%) | 93 |
| 3 | 66 (53.7%) | 33 (26.8%) | 24 (19.5%) | 123 |
| 4 | 57 (54.8%) | 5 (4.8%) | 42 (40.4%) | 104 |
| 5 | 85 (57.4%) | 36 (24.3%) | 27 (18.2%) | 148 |
| **Total** | **249 (50.7%)** | **140 (28.5%)** | **102 (20.8%)** | **491** |

**Chi-Square Test**: χ² = 101.18, df = 8, p < 10⁻¹⁸
**Cramér's V**: 0.321 (Medium effect)

**Interpretation**: Statistically significant association between clusters and regions. Central Luzon shows enrichment for Cluster 1 (Salmonella) and Cluster 2 (Enterobacter). BARMM and Eastern Visayas show similar patterns.

### 1.2 Regional MDR Rates

| Region | Total Isolates | MDR Count | MDR Rate | Fold Enrichment |
|--------|----------------|-----------|----------|-----------------|
| BARMM | 249 | 51 | 20.5% | 1.07× |
| Region III - Central Luzon | 140 | 30 | 21.4% | 1.12× |
| Region VIII - Eastern Visayas | 102 | 13 | 12.7% | 0.66× |
| **Overall** | **491** | **94** | **19.1%** | **1.0×** |

### 1.3 Regional Cluster Enrichment

| Region | Most Enriched Cluster | Enrichment | Interpretation |
|--------|----------------------|------------|----------------|
| BARMM | Cluster 5 (K. pneumoniae) | 1.15× | Slight enrichment |
| Central Luzon | Cluster 1 (Salmonella) | 2.59× | Strong enrichment |
| Eastern Visayas | Cluster 4 (E. coli susceptible) | 1.95× | Notable enrichment |

---

## 2. Environmental Distribution Analysis

### 2.1 Cluster Distribution by Environment

| Cluster | Fish | Hospital | Water | Total |
|---------|------|----------|-------|-------|
| 1 | 7 (30.4%) | 0 (0%) | 16 (69.6%) | 23 |
| 2 | 50 (53.8%) | 0 (0%) | 43 (46.2%) | 93 |
| 3 | 69 (56.1%) | 9 (7.3%) | 45 (36.6%) | 123 |
| 4 | 61 (58.7%) | 19 (18.3%) | 24 (23.1%) | 104 |
| 5 | 87 (58.8%) | 13 (8.8%) | 48 (32.4%) | 148 |
| **Total** | **274 (55.8%)** | **41 (8.4%)** | **176 (35.8%)** | **491** |

**Chi-Square Test**: χ² = 42.56, df = 8, p < 0.0001
**Cramér's V**: 0.204 (Small effect)

**Interpretation**: Statistically significant but weak association between clusters and environment. Cluster 1 (Salmonella) is notably enriched in Water environments.

### 2.2 Environment MDR Rates

| Environment | Total Isolates | MDR Count | MDR Rate | Fold Enrichment |
|-------------|----------------|-----------|----------|-----------------|
| Fish | 274 | 55 | 20.1% | 1.05× |
| Water | 176 | 34 | 19.3% | 1.01× |
| Hospital | 41 | 5 | 12.2% | 0.64× |
| **Overall** | **491** | **94** | **19.1%** | **1.0×** |

### 2.3 Sample Source Details

| Sample Source | Count | MDR Rate | Most Common Cluster |
|---------------|-------|----------|---------------------|
| Fish Tilapia | 115 | 21.7% | Cluster 5 |
| Drinking Water | 82 | 18.3% | Cluster 5 |
| Fish Gusaw | 70 | 17.1% | Cluster 3 |
| River Water | 56 | 23.2% | Cluster 3 |
| Fish Kaolang | 52 | 19.2% | Cluster 4 |
| Lake Water | 39 | 17.9% | Cluster 2 |
| Fish Banak | 37 | 21.6% | Cluster 5 |
| Effluent Untreated | 36 | 22.2% | Cluster 3 |
| Effluent Treated | 5 | 0.0% | Cluster 4 |

---

## 3. Principal Component Analysis Results

### 3.1 Variance Explained

| Component | Variance Explained | Cumulative | Interpretation |
|-----------|-------------------|------------|----------------|
| PC1 | 23.5% | 23.5% | Primary resistance variation axis |
| PC2 | 16.4% | 39.9% | Secondary resistance variation |
| PC3 | 9.1% | 49.0% | Minor variation component |
| PC4 | 7.2% | 56.2% | Minor variation component |
| PC5 | 5.9% | 62.1% | Minor variation component |

**Total variance explained (PC1+PC2)**: **39.9%**

> ⚠️ **Important**: First two components explain less than 40% of variance. 2D scatter plots are simplified representations; cluster separation is more pronounced in full 22-dimensional resistance space.

### 3.2 Component Loadings (Top Antibiotics)

#### PC1 Loadings

| Antibiotic | Loading | Direction |
|------------|---------|-----------|
| TE | 0.412 | Positive |
| DO | 0.398 | Positive |
| AM | 0.312 | Positive |
| SXT | 0.287 | Positive |

**PC1 Interpretation**: Captures overall resistance level, particularly tetracycline/doxycycline resistance. High PC1 → higher resistance.

#### PC2 Loadings

| Antibiotic | Loading | Direction |
|------------|---------|-----------|
| CFT | 0.352 | Positive |
| CPD | 0.298 | Positive |
| AN | -0.325 | Negative |
| GM | -0.301 | Negative |

**PC2 Interpretation**: Contrasts cephalosporin vs. aminoglycoside resistance patterns.

### 3.3 PCA Visualization Summary

| Color Coding Variable | Visual Separation | Notes |
|----------------------|-------------------|-------|
| Cluster | Moderate | Clusters 3 and 4 show clear separation |
| Region | Overlapping | No clear regional stratification |
| MDR Status | Moderate | MDR isolates cluster in positive PC1 direction |
| Species | Moderate | Some species-specific patterns visible |

---

## 4. Statistical Associations

### 4.1 Chi-Square Test Summary

| Association | χ² | df | p-value | Cramér's V | Effect Size |
|-------------|----|----|---------|------------|-------------|
| Cluster × Region | 101.18 | 8 | < 10⁻¹⁸ | 0.321 | Medium |
| Cluster × Environment | 42.56 | 8 | < 0.0001 | 0.204 | Small |
| Cluster × Species | 653.12 | 36 | < 10⁻⁵⁰ | **0.765** | **Large** |
| Cluster × MDR | 154.04 | 4 | < 10⁻³⁰ | **0.559** | **Large** |
| Species × Environment | 98.67 | 18 | < 10⁻¹⁰ | 0.316 | Medium |
| Region × MDR | 4.12 | 2 | 0.127 | 0.091 | Negligible |

### 4.2 Key Findings Summary

1. **Regional Patterns**: Significant cluster-region association (V = 0.321, medium effect). Central Luzon shows distinct profile with Salmonella/Enterobacter enrichment.

2. **Environmental Patterns**: Weak but significant cluster-environment association (V = 0.204). Species dominates over environmental factors.

3. **PCA Insights**: 39.9% variance explained by PC1+PC2. Tetracyclines dominate PC1 loadings. Full separation requires higher dimensions.

4. **MDR Hotspots**: Cluster 3 (*E. coli*-tetracycline) is the primary MDR hotspot at 54.5% MDR rate (2.82× enrichment).

---

## Figure Captions

### Figure 4.1: PCA Scatter Plot by Cluster

> "Principal Component Analysis of resistance profiles (n=491 isolates, 22 antibiotics, 6 species). PC1 (23.5% variance) and PC2 (16.4% variance) capture variation primarily in tetracycline and cephalosporin resistance. Points colored by cluster assignment."

### Figure 4.2: PCA Scatter Plot by Region

> "PCA visualization with regional coloring (n=491 isolates, 6 species). Substantial overlap indicates resistance patterns are not strongly geographically stratified. Regional associations are observational."

### Figure 4.3: Regional Distribution Bar Chart

> "Distribution of isolates across clusters by geographic region. Chi-square test indicates significant association (χ²=101.18, p<10⁻¹⁸, Cramér's V=0.321, medium effect)."

### Figure 4.4: Scree Plot

> "Scree plot showing variance explained by principal components. First two components capture 39.9% of variance, necessitating caution in 2D interpretation."

---

## Interpretation Notes

> **Scope**: Regional and environmental associations are **observational** and reflect patterns in the analyzed dataset. Results do not establish causation between geographic/environmental factors and resistance profiles. The weak environmental effect (V = 0.204) suggests species identity dominates over environmental selection.

---

## Document Version

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2024 | Initial results template |
| 2.0 | 2025-12-18 | **Populated with actual pipeline results** |

---

*This document is part of Phase 8 — Documentation & Reporting for the AMR Thesis Project.*

**See also**: [multivariate_analysis.md](../methods/multivariate_analysis.md) | [phase5_synthesis.md](phase5_synthesis.md)
