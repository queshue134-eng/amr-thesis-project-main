# Integration Methods Documentation

## Phase 6: Integration and Synthesis

This document provides comprehensive documentation of integration methods that synthesize results from clustering, supervised learning, and multivariate analysis for the AMR thesis project.

---

## Table of Contents

1. [Integration Objectives](#1-integration-objectives)
2. [Cluster-Supervised Comparison](#2-cluster-supervised-comparison)
3. [Resistance Archetype Identification](#3-resistance-archetype-identification)
4. [Species-Environment Association Analysis](#4-species-environment-association-analysis)
5. [MDR-Enriched Pattern Identification](#5-mdr-enriched-pattern-identification)
6. [Phase Linking Guidelines](#6-phase-linking-guidelines)

---

## 1. Integration Objectives

### 1.1 Purpose

The integration phase synthesizes findings from all analysis phases to:

1. Evaluate alignment between unsupervised (clustering) and supervised (discrimination) results
2. Identify characteristic resistance archetypes for each cluster
3. Characterize species-environment associations
4. Identify MDR-enriched patterns across clusters, regions, and environments

### 1.2 Integration Framework

```
┌─────────────────────────────────────────────────────────────────┐
│                    INTEGRATION FRAMEWORK                         │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  INPUTS                                                          │
│  ├── Phase 3: Clustered dataset with cluster labels              │
│  ├── Phase 4: Supervised model results and feature importance    │
│  └── Phase 5: PCA results and regional distributions             │
│                                                                 │
│  INTEGRATION METHODS                                             │
│  ├── Cluster-supervised comparison                               │
│  ├── Resistance archetype identification                         │
│  ├── Species-environment association analysis                    │
│  └── MDR-enriched pattern identification                         │
│                                                                 │
│  OUTPUTS                                                         │
│  ├── Cluster purity analysis                                     │
│  ├── Archetype characterization tables                           │
│  ├── Association matrices                                        │
│  └── MDR enrichment summaries                                    │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

---

## 2. Cluster-Supervised Comparison

### 2.1 Objective

Evaluate how well unsupervised clusters align with supervised discrimination results (species identity, MDR status).

### 2.2 Methods

#### 2.2.1 Cross-Tabulation

| Cross-Tab | Purpose |
|-----------|---------|
| Clusters × Species | Do clusters correspond to species groups? |
| Clusters × MDR Status | Do clusters capture MDR patterns? |

#### 2.2.2 Cluster Purity Calculation

**Definition**:
```
Purity(cluster) = max(category_count) / cluster_size

Where:
- category_count = count of most common category in cluster
- cluster_size = total isolates in cluster
```

**Interpretation**:
| Purity | Interpretation |
|--------|----------------|
| > 0.9 | Very high; cluster nearly homogeneous |
| 0.7 - 0.9 | High; cluster dominated by one category |
| 0.5 - 0.7 | Moderate; cluster has mixed composition |
| < 0.5 | Low; cluster is heterogeneous |

### 2.3 Interpretation Guidelines

| Finding | Interpretation |
|---------|----------------|
| High cluster-MDR association | Clustering captures MDR-related resistance patterns |
| High cluster-species association | Clustering reflects species-specific resistance profiles |
| High cluster purity | Clusters represent homogeneous groups |
| Low purity but significant association | Clusters capture gradient rather than discrete categories |

### 2.4 Reporting Template

```markdown
### Cluster-Supervised Comparison

**Objective**: Assess alignment between unsupervised clusters and supervised categories.

#### Cluster × MDR Cross-Tabulation

| Cluster | Non-MDR | MDR | Total | MDR Rate | Purity |
|---------|---------|-----|-------|----------|--------|
| 1 | [n] | [n] | [N] | [%] | [value] |
| 2 | [n] | [n] | [N] | [%] | [value] |

**Chi-square**: χ² = [value], df = [df], p = [p-value]

**Interpretation**: [Significant/Not significant] association between clusters and MDR status. [Additional observations about purity and patterns.]
```

---

## 3. Resistance Archetype Identification

### 3.1 Definition

An **archetype** is a characteristic resistance profile that defines each cluster, determined by the mean resistance pattern across cluster members.

### 3.2 Archetype Characterization Method

For each cluster, compute:

| Metric | Calculation | Purpose |
|--------|-------------|---------|
| Mean resistance profile | Average encoded value (0-2) per antibiotic | Define cluster archetype |
| High-resistance antibiotics | Antibiotics with mean > 1.5 | Identify dominant resistance |
| Low-resistance antibiotics | Antibiotics with mean < 0.5 | Identify susceptibility |
| Resistance breadth | Count of high-resistance antibiotics | Characterize resistance extent |

### 3.3 Resistance Level Classification

| Mean Encoded Value | Classification | Description |
|--------------------|----------------|-------------|
| > 1.5 | High resistance | Approaching "R" (resistant) |
| 1.0 - 1.5 | Moderate-high resistance | Mixed I/R |
| 0.5 - 1.0 | Moderate resistance | Mixed S/I |
| < 0.5 | Low resistance | Approaching "S" (susceptible) |

### 3.4 Archetype Table Template

| Cluster | N | MDR Rate | MAR Index | High-Resistance ABs | Low-Resistance ABs | Archetype Label |
|---------|---|----------|-----------|---------------------|-------------------|-----------------|
| 1 | [n] | [%] | [mean] | [AB list] | [AB list] | [descriptive label] |
| 2 | [n] | [%] | [mean] | [AB list] | [AB list] | [descriptive label] |

### 3.5 Archetype Labeling Convention

| Pattern Characteristics | Suggested Label |
|------------------------|-----------------|
| High β-lactam, high aminoglycoside | "Multi-class resistant" |
| High β-lactam, low others | "β-lactam focused resistance" |
| Low resistance across all | "Susceptible profile" |
| High fluoroquinolone focus | "Quinolone-resistant profile" |

---

## 4. Species-Environment Association Analysis

### 4.1 Objective

Identify associations between bacterial species and environmental sources to characterize ecological patterns.

### 4.2 Methods

#### 4.2.1 Cross-Tabulation Analysis

Cross-tabulate:
- Species × Sample Source
- Species × Environment Category
- Species × Region

#### 4.2.2 Dominant Environment Identification

For each species, identify:
- Primary environment (highest count)
- Secondary environment (second highest)
- Breadth (number of environments with ≥10% of species isolates)

### 4.3 Association Table Template

| Species | Primary Env | Primary % | Secondary Env | Secondary % | Breadth |
|---------|-------------|-----------|---------------|-------------|---------|
| E. coli | Water | 45% | Fish | 30% | 3 |
| K. pneumoniae | Hospital | 60% | Water | 25% | 2 |

### 4.4 Statistical Testing

**Chi-square test** for species-environment independence:
- H₀: Species and environment are independent
- H₁: Species distribution differs by environment
- Significance: α = 0.05

---

## 5. MDR-Enriched Pattern Identification

### 5.1 Objective

Identify clusters, regions, environments, and species with higher-than-expected MDR rates.

### 5.2 Enrichment Calculation

**Fold Enrichment**:
```
Fold_Enrichment = Group_MDR_Rate / Overall_MDR_Rate

Example:
- Overall MDR rate: 40%
- Cluster 3 MDR rate: 80%
- Fold enrichment: 80% / 40% = 2.0×
```

### 5.3 MDR Enrichment Table Template

| Grouping Variable | Category | N | MDR Count | MDR Rate | Fold Enrichment |
|-------------------|----------|---|-----------|----------|-----------------|
| Cluster | Cluster 1 | [n] | [n] | [%] | [×] |
| Cluster | Cluster 2 | [n] | [n] | [%] | [×] |
| Region | BARMM | [n] | [n] | [%] | [×] |
| Environment | Hospital | [n] | [n] | [%] | [×] |
| Species | E. coli | [n] | [n] | [%] | [×] |

**Overall MDR rate**: [%]

### 5.4 MDR Resistance Signature

Identify antibiotics most associated with MDR by comparing mean resistance:

| Antibiotic | MDR Mean | Non-MDR Mean | Difference | MDR Associated? |
|------------|----------|--------------|------------|-----------------|
| [AB] | [value] | [value] | [diff] | Yes/No |

**MDR-associated**: Difference > 0.5 and MDR Mean > 1.5

---

## 6. Phase Linking Guidelines

### 6.1 Phase Linking Sentences (MANDATORY)

**At the start of each Results section**, include a linking sentence that connects to previous phases:

| Current Phase | Linking Sentence Template |
|---------------|---------------------------|
| **Phase 3 (Clustering)** | "Following data preprocessing (Phase 2), hierarchical clustering was performed to identify natural groupings in resistance profiles..." |
| **Phase 4 (Supervised)** | "Building on the resistance-based clusters identified in Phase 3, supervised learning evaluated how resistance fingerprints discriminate known categories..." |
| **Phase 5 (Regional)** | "With cluster assignments from Phase 3, regional and environmental distribution analysis characterized the geographic and ecological patterns..." |
| **Phase 6 (Integration)** | "Integrating results from clustering (Phase 3), supervised discrimination (Phase 4), and regional analysis (Phase 5), this synthesis identifies dominant resistance archetypes..." |

### 6.2 Cross-Reference Format

When referring to results from another phase:

```markdown
As shown in Phase 3 (Section 3.2, Table 3.1), Cluster 2 exhibited the highest 
MDR proportion (72%). Supervised discrimination (Phase 4) confirmed that 
resistance fingerprints effectively distinguish MDR status (F1=0.85, Table 4.2).
```

### 6.3 Integration Summary Template

```markdown
## Phase 6: Integration and Synthesis Summary

### Key Integrated Findings

1. **Cluster-MDR Alignment**: [Finding about how clusters relate to MDR]

2. **Resistance Archetypes**: [Summary of identified archetypes]

3. **Species-Environment Patterns**: [Key species-environment associations]

4. **MDR Enrichment Hotspots**: [Clusters/regions/environments with elevated MDR]

### Coherence Statement

> The integration of unsupervised clustering (Phase 3), supervised discrimination 
> (Phase 4), and multivariate analysis (Phase 5) reveals consistent patterns: 
> [summary of cross-phase consistency]. These findings support [main conclusion] 
> while acknowledging [key limitations].
```

---

## Document Version

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2024 | Initial integration documentation |

---

*This document is part of Phase 8 — Documentation & Reporting for the AMR Thesis Project.*

**See also**: [clustering.md](clustering.md) | [supervised_models.md](supervised_models.md) | [multivariate_analysis.md](multivariate_analysis.md)
