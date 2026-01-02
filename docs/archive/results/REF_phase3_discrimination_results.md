# Phase 3 Results: Clustering and Pattern Discrimination

> **⚠️ DEPRECATED — HISTORICAL DOCUMENT**  
> **This document reflects the original k=5 clustering analysis.**  
> **Current analysis uses k=4** (determined via comprehensive multi-metric validation).  
> **For current cluster selection rationale, see:** [REF_clustering_technical_notes.md § 4.3](../methods/REF_clustering_technical_notes.md#43-comprehensive-validation-metrics-comparison)  
> **This file is preserved for historical reference only.**

---

## Hierarchical Clustering and Supervised Learning Results

This document provides the actual Phase 3 (clustering) and Phase 4 (supervised discrimination) results for the AMR thesis project.

---

## Table of Contents

1. [Clustering Results (Phase 3)](#1-clustering-results-phase-3)
2. [Supervised Discrimination Results (Phase 4)](#2-supervised-discrimination-results-phase-4)
3. [Feature Importance Analysis](#3-feature-importance-analysis)
4. [Cluster-Supervised Comparison](#4-cluster-supervised-comparison)

---

## 1. Clustering Results (Phase 3)

### 1.1 Clustering Parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Linkage method | Ward | Minimizes within-cluster variance |
| Distance metric | Euclidean | Captures ordinal resistance differences |
| Number of clusters | **5** | Selected via combined elbow + silhouette analysis |
| Imputation | Median | Robust to outliers |
| Silhouette score | **0.488** | Good clustering structure (> 0.4 threshold) |

### 1.1.1 Cluster Number Selection (k=5 Justification)

The choice of k=5 clusters was determined using a **combined methodology** balancing multiple criteria:

| k | Silhouette | WCSS | Interpretation |
|---|------------|------|----------------|
| 2 | 0.378 | 2395 | Moderate (below elbow) |
| 3 | 0.417 | 1769 | Strong (below elbow) |
| 4 | 0.465 | 1486 | Strong (elbow point) |
| **5** | **0.488** | **1238** | **\* SELECTED (elbow + sample size) \*** |
| 6 | 0.517 | 1013 | Higher sil. (small clusters) |
| 7 | 0.527 | 895 | Higher sil. (small clusters) |
| 8 | 0.551 | 796 | Higher sil. (small clusters) |
| 9 | 0.573 | 727 | Higher sil. (small clusters) |
| 10 | 0.585 | 660 | Higher sil. (small clusters) |

**Interpretation Key:**

| Label | Meaning |
|-------|---------|
| **(below elbow)** | k values before WCSS diminishing returns point; too few clusters to capture resistance heterogeneity |
| **(elbow point)** | Where WCSS reduction begins to slow significantly |
| **\* SELECTED \*** | Chosen for optimal balance of elbow position, sample size, and silhouette |
| **(small clusters)** | Higher silhouette but creates clusters with <20 isolates, reducing statistical reliability |

**Why k=5 instead of maximizing silhouette (k=10)?**

1. **Elbow Method:** WCSS reduction slows significantly around k=4-5 (diminishing returns)
2. **Sample Size Constraint:** At k=5, smallest cluster has 23 isolates (adequate for chi-square tests). At k=10, smallest clusters have <15 isolates, reducing statistical power.
3. **Silhouette Trade-off:** k=5 (0.488) exceeds the 0.4 "good clustering" threshold. The marginal improvement to 0.585 at k=10 (+0.097) does not justify the loss of interpretability and statistical reliability.
4. **Literature Alignment:** AMR phenotype studies typically identify 3-7 major resistance patterns; 5 clusters is biologically interpretable.

> **Note:** Negative per-sample silhouette values (visible in detailed silhouette plots) indicate isolates near cluster boundaries. This is expected and does not invalidate the overall clustering quality as long as the average silhouette remains above threshold.

### 1.2 Cluster Summary

| Cluster | N | Percentage | MDR Rate | Mean MAR | Dominant Species | Major Environment |
|---------|---|------------|----------|----------|------------------|-------------------|
| 1 | 23 | 4.7% | 26.1% | 0.196 | *Salmonella* (100.0%) | Water (69.6%) |
| 2 | 93 | 18.9% | 20.4% | 0.164 | *E. cloacae* (71.0%) | Fish (53.2%) |
| 3 | 123 | 25.1% | **54.5%** | 0.181 | *E. coli* (77.2%) | Fish (56.1%) |
| 4 | 104 | 21.2% | **0.0%** | 0.002 | *E. coli* (98.1%) | Fish (58.7%) |
| 5 | 148 | 30.1% | 1.4% | 0.053 | *K. pneumoniae* (77.0%) | Fish (58.8%) |
| **Total** | **491** | **100%** | **19.1%** | **0.102** | - | - |

### 1.3 Cluster Profiles

#### Cluster 1: Salmonella-Aminoglycoside

| Attribute | Value |
|-----------|-------|
| **Size** | 23 isolates (4.7% of dataset) |
| **MDR Rate** | 26.1% |
| **Mean MAR** | 0.196 |
| **Resistance Level** | Moderate |

**High-Resistance Antibiotics**: AN, CN, GM (aminoglycosides)
**Low-Resistance Antibiotics**: AMC, CTX, CPT

**Characterization**: Small cluster dominated by *Salmonella* from water environments with aminoglycoside resistance pattern.

---

#### Cluster 2: Enterobacter-Beta-lactam

| Attribute | Value |
|-----------|-------|
| **Size** | 93 isolates (18.9% of dataset) |
| **MDR Rate** | 20.4% |
| **Mean MAR** | 0.164 |
| **Resistance Level** | Moderate |

**High-Resistance Antibiotics**: AM, CF, CN (ampicillin/cephalosporins)
**Low-Resistance Antibiotics**: AN, CTX, CPT

**Characterization**: *Enterobacter cloacae*-dominated cluster with beta-lactam resistance, mixed fish/water sources.

---

#### Cluster 3: E. coli-Tetracycline-MDR

| Attribute | Value |
|-----------|-------|
| **Size** | 123 isolates (25.0% of dataset) |
| **MDR Rate** | **54.5%** (highest) |
| **Mean MAR** | 0.181 |
| **Resistance Level** | **High** |

**High-Resistance Antibiotics**: TE, DO, AM (tetracyclines/ampicillin)
**Low-Resistance Antibiotics**: AN, CTX, CZA

**Characterization**: MDR-enriched *E. coli* cluster with tetracycline dominance. Aquaculture-associated.

---

#### Cluster 4: E. coli-Susceptible

| Attribute | Value |
|-----------|-------|
| **Size** | 104 isolates (21.1% of dataset) |
| **MDR Rate** | **0.0%** (lowest) |
| **Mean MAR** | 0.002 |
| **Resistance Level** | **Susceptible** |

**High-Resistance Antibiotics**: None (broadly susceptible)
**Low-Resistance Antibiotics**: CF, C, CFT

**Characterization**: Nearly pan-susceptible *E. coli* cluster. May represent commensal strains unexposed to antibiotic selection.

---

#### Cluster 5: Klebsiella-Intermediate

| Attribute | Value |
|-----------|-------|
| **Size** | 148 isolates (30.1% of dataset) |
| **MDR Rate** | 1.4% |
| **Mean MAR** | 0.053 |
| **Resistance Level** | Low |

**High-Resistance Antibiotics**: AM, FT, CN (low levels)
**Low-Resistance Antibiotics**: Most antibiotics

**Characterization**: *K. pneumoniae*-dominated cluster with minimal resistance. Largest cluster.

---

### 1.4 Cluster Quality Metrics

| Metric | Value | Interpretation |
|--------|-------|----------------|
| Silhouette score | 0.488 | Strong clustering (> 0.4) |
| Average cluster size | 98.4 | Well-distributed |
| Size standard deviation | 45.3 | Moderate variability |
| Robustness ARI | > 0.8 | Stable across methods |

### 1.5 k Selection Justification

| k | Silhouette Score | Notes |
|---|------------------|-------|
| 3 | 0.421 | Acceptable |
| 4 | 0.465 | Strong |
| **5** | **0.488** | **Selected** |
| 6 | 0.517 | Higher but fragmented |

**Rationale**: k=5 selected based on:
1. Strong silhouette (0.488)
2. Biological interpretability (distinct species/resistance patterns)
3. Adequate cluster sizes (min 23, sufficient for statistics)

---

## 2. Supervised Discrimination Results (Phase 4)

> **Scope Statement**: Results reflect pattern discrimination within the analyzed dataset and are **not** predictive of performance on external datasets.

### 2.1 Task A: Species Discrimination

#### Data Split

| Set | Samples | Percentage |
|-----|---------|------------|
| Training | 392 | 80% |
| Test | 99 | 20% |
| **Total** | **491** | **100%** |

#### Model Evaluation Metrics

| Model | Accuracy | Precision (Macro) | Recall (Macro) | F1-Score (Macro) |
|-------|----------|-------------------|----------------|------------------|
| Logistic Regression | 0.727 | 0.698 | 0.712 | 0.705 |
| Random Forest | **0.736** | **0.721** | **0.728** | **0.724** |
| k-Nearest Neighbors | 0.689 | 0.665 | 0.674 | 0.669 |

**Best Model**: Random Forest (F1 = 0.724)

#### Interpretation

> "The Random Forest model demonstrates **moderate discriminative capacity** for species identification based on resistance fingerprints (F1 = 0.724). *E. coli* and *K. pneumoniae* are well-distinguished; *Enterobacter* species show some confusion due to similar resistance profiles."

### 2.2 Task B: Cluster Discrimination

#### Model Evaluation Metrics

| Model | Accuracy | F1-Score (Macro) |
|-------|----------|------------------|
| Logistic Regression | 0.891 | 0.876 |
| **Random Forest** | **0.923** | **0.912** |
| k-Nearest Neighbors | 0.867 | 0.854 |

**Best Model**: Random Forest (F1 = 0.912)

#### Interpretation

> "Supervised models effectively discriminate clusters based on resistance features (F1 = 0.912), validating that clusters represent distinct, learnable patterns."

---

## 3. Feature Importance Analysis

### 3.1 Species Discrimination - Top Antibiotics

| Rank | Antibiotic | Importance Score | Biological Context |
|------|------------|------------------|-------------------|
| 1 | TE | 0.142 | Tetracycline - aquaculture associated |
| 2 | DO | 0.128 | Doxycycline - tetracycline class |
| 3 | AM | 0.115 | Ampicillin - intrinsic resistance marker |
| 4 | CFT | 0.098 | Cefixime - cephalosporin marker |
| 5 | SXT | 0.087 | Trimethoprim-sulfamethoxazole |

### 3.2 Feature Importance Summary

| Antibiotic Class | Combined Importance | Notes |
|------------------|---------------------|-------|
| Tetracyclines (TE, DO) | 0.270 | Highest discriminative power |
| Penicillins (AM) | 0.115 | Species-specific intrinsic patterns |
| Folate inhibitors (SXT) | 0.087 | Co-resistance marker |

> **Note**: Importance scores indicate associative contribution to group separation, **not** causal relationships.

---

## 4. Cluster-Supervised Comparison

### 4.1 Clusters vs. MDR Status

| Cluster | Non-MDR | MDR | Total | MDR Rate | Enrichment |
|---------|---------|-----|-------|----------|------------|
| 1 | 17 | 6 | 23 | 26.1% | 1.37× |
| 2 | 74 | 19 | 93 | 20.4% | 1.07× |
| 3 | 56 | **67** | 123 | **54.5%** | **2.86×** |
| 4 | 104 | 0 | 104 | 0.0% | 0.00× |
| 5 | 146 | 2 | 148 | 1.4% | 0.07× |
| **Total** | **397** | **94** | **491** | **19.1%** | **1.0×** |

**Chi-Square Test**: χ² = 154.04, df = 4, p < 0.0001
**Cramér's V**: 0.559 (Large effect)

**Interpretation**: Highly significant association between clusters and MDR status. Cluster 3 is 2.82× enriched for MDR compared to overall rate.

### 4.2 Clusters vs. Species

| Cluster | E. coli | K. pneum. | E. cloacae | Other | Dominant |
|---------|---------|-----------|------------|-------|----------|
| 1 | 0 | 0 | 2 | 21 | Salmonella (100.0%) |
| 2 | 6 | 22 | 66 | 0 | E. cloacae (70.2%) |
| 3 | 95 | 15 | 0 | 13 | E. coli (77.2%) |
| 4 | 102 | 0 | 0 | 2 | E. coli (98.1%) |
| 5 | 24 | 109 | 0 | 15 | K. pneumoniae (77.0%) |

**Chi-Square Test**: χ² = 653.12, df = 36, p < 0.0001
**Cramér's V**: 0.765 (Large effect)

**Interpretation**: Very strong species-cluster association. Species identity is the dominant factor in cluster membership.

---

## Figure Captions

### Figure 3.1: Hierarchical Clustering Dendrogram

> "Dendrogram showing hierarchical clustering of resistance profiles (n=491 isolates, 6 species) using Ward linkage with Euclidean distance. Horizontal line indicates cut point for k=5 clusters. Silhouette score = 0.488."

### Figure 3.2: Resistance Heatmap with Clustering

> "Heatmap illustrating hierarchical clustering of resistance profiles. Clusters were derived **solely from resistance features** (S=0, I=1, R=2); metadata variables are shown for **post-hoc interpretation only**."

### Figure 4.1: Confusion Matrix (Species Discrimination)

> "Confusion matrix for species discrimination using Random Forest classifier (n=99 test samples). F1-macro = 0.724. Results reflect pattern consistency, not predictive performance on external data."

---

## Document Version

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2024 | Initial results template |
| 2.0 | 2025-12-18 | **Populated with actual pipeline results** |

---

*This document is part of Phase 8 — Documentation & Reporting for the AMR Thesis Project.*

**See also**: [clustering.md](../methods/clustering.md) | [supervised_models.md](../methods/supervised_models.md)
