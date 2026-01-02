# Multivariate Analysis Documentation

## Phase 5: Regional and Environmental Analysis

This document provides comprehensive documentation of multivariate analysis methods, including Principal Component Analysis (PCA), cross-tabulation, and statistical testing for the AMR thesis project.

---

## Table of Contents

1. [Principal Component Analysis (PCA)](#1-principal-component-analysis-pca)
2. [Cross-Tabulation Analysis](#2-cross-tabulation-analysis)
3. [Statistical Testing](#3-statistical-testing)
4. [Reporting Standards](#4-reporting-standards)

---

## 1. Principal Component Analysis (PCA)

### 1.1 Objective

Reduce dimensionality of resistance profiles and visualize patterns in a lower-dimensional space, enabling identification of dominant sources of variation in the dataset.

### 1.2 PCA Parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| **Number of components** | 2 (default) | Sufficient for visualization; can be extended |
| **Preprocessing** | StandardScaler | Zero mean, unit variance for each feature |
| **Missing value handling** | Median imputation | Applied before scaling |

### 1.3 PCA Procedure

```
1. DATA PREPARATION
   ├── Extract encoded resistance columns ({AB}_encoded)
   ├── Impute missing values (median strategy)
   └── Standardize features (StandardScaler)

2. PCA COMPUTATION
   ├── Fit PCA model to standardized data
   ├── Extract principal components
   ├── Compute explained variance ratios
   └── Calculate component loadings

3. INTERPRETATION
   ├── Identify antibiotics with highest absolute loadings
   ├── Interpret components in terms of resistance patterns
   └── Visualize in 2D space with metadata coloring
```

### 1.4 Loading Interpretation

| Loading Value | Interpretation |
|---------------|----------------|
| > 0.3 | Strong positive contribution to component |
| < -0.3 | Strong negative contribution to component |
| -0.3 to 0.3 | Weak contribution |

### 1.5 Component Interpretation Template

> "**PC1** (XX% variance explained) showed high loadings for [antibiotics], suggesting this axis captures variation in [resistance pattern description]. **PC2** (YY% variance explained) showed high loadings for [antibiotics], representing [pattern description]."

### 1.6 Visualization Types

| Plot Type | Description | Color Coding |
|-----------|-------------|--------------|
| **PCA Scatter Plot** | 2D projection of isolates | By cluster, region, or MDR status |
| **PCA Biplot** | Scatter plot with loading vectors | Shows antibiotic contributions |

---

## 2. Cross-Tabulation Analysis

### 2.1 Objective

Analyze distributions and associations between clusters and metadata variables (region, environment, species).

### 2.2 Cross-Tabulations Computed

| Cross-Tabulation | Purpose |
|------------------|---------|
| Clusters × Regions | Geographic distribution of resistance patterns |
| Clusters × Sample Sources | Environmental associations |
| Clusters × Species | Species composition of clusters |
| Clusters × MDR Status | MDR enrichment by cluster |

### 2.3 Cross-Tab Table Template

**Example: Cluster × Region Cross-Tabulation**

| Cluster | BARMM | Region III | Region VIII | Total |
|---------|-------|------------|-------------|-------|
| Cluster 1 | [n] | [n] | [n] | [N] |
| Cluster 2 | [n] | [n] | [n] | [N] |
| ... | ... | ... | ... | ... |
| **Total** | [N] | [N] | [N] | **[Total]** |

### 2.4 Enrichment Analysis

**Fold Enrichment Calculation**:
```
Fold_Enrichment = Group_Rate / Overall_Rate

Example: MDR enrichment in Cluster 3
- Cluster 3 MDR rate: 75%
- Overall MDR rate: 45%
- Fold enrichment: 75% / 45% = 1.67×
```

**Interpretation Guide**:
| Fold Enrichment | Interpretation |
|-----------------|----------------|
| > 2.0 | Strong enrichment |
| 1.5 - 2.0 | Moderate enrichment |
| 1.0 - 1.5 | Slight enrichment |
| < 1.0 | Depletion |

---

## 3. Statistical Testing

### 3.1 Chi-Square Test of Independence

**Purpose**: Test for statistical associations between categorical variables.

**Hypotheses**:
| Test | H₀ (Null) | H₁ (Alternative) |
|------|-----------|------------------|
| Cluster-Region | Clusters and regions are independent | Clusters differ by region |
| Cluster-Environment | Clusters and sample sources are independent | Clusters differ by environment |
| Cluster-Species | Clusters and species are independent | Clusters differ by species |
| Cluster-MDR | Clusters and MDR status are independent | Clusters differ by MDR rate |

**Decision Rule**:
- **Significance level**: α = 0.05
- If p < 0.05: Reject H₀, association is statistically significant
- If p ≥ 0.05: Fail to reject H₀, no significant association

### 3.2 Chi-Square Assumptions and Limitations

| Assumption | Description | Mitigation |
|------------|-------------|------------|
| Expected counts ≥ 5 | Each cell should have expected count ≥ 5 | Report if violated; consider Fisher's exact test |
| Independence | Observations are independent | Ensured by study design |
| Random sampling | Samples represent target population | Acknowledge as limitation if not met |

### 3.3 Statistical Results Table Template

| Association | Chi-Square (χ²) | df | p-value | Significant? |
|-------------|-----------------|----|---------|--------------| 
| Cluster × Region | [value] | [df] | [p] | Yes/No |
| Cluster × Environment | [value] | [df] | [p] | Yes/No |
| Cluster × Species | [value] | [df] | [p] | Yes/No |
| Cluster × MDR | [value] | [df] | [p] | Yes/No |

---

## 4. Reporting Standards

### 4.1 PCA Results Reporting Template

```markdown
### Principal Component Analysis Results

**Objective**: Reduce dimensionality and visualize resistance pattern variation.

**Method**: PCA performed on standardized encoded resistance features (n=[N] antibiotics).

#### Variance Explained

| Component | Variance Explained | Cumulative |
|-----------|-------------------|------------|
| PC1 | [X]% | [X]% |
| PC2 | [Y]% | [X+Y]% |

#### Component Loadings (Top Antibiotics)

**PC1**:
- [AB1]: [loading]
- [AB2]: [loading]

**PC2**:
- [AB1]: [loading]
- [AB2]: [loading]

**Interpretation**: [Component interpretation in terms of resistance patterns]
```

### 4.2 Distribution Analysis Reporting Template

```markdown
### Cluster Distribution by [Variable]

**Objective**: Characterize cluster composition by [variable].

**Method**: Cross-tabulation with chi-square test of independence.

#### Cross-Tabulation Results

[Cross-tab table]

#### Statistical Test

- Chi-square: [χ²]
- Degrees of freedom: [df]
- p-value: [p]
- **Interpretation**: [Significant/Not significant] association between clusters and [variable].

#### Key Findings

- Cluster [X] showed enrichment in [category] ([fold]× enrichment)
- [Additional observations]
```

### 4.3 Figure Caption Requirements

**PCA Plot Caption Elements**:
1. Sample size and feature count
2. Variance explained by displayed components
3. Color coding variable
4. Interpretation scope

**Example Caption**:
> "**Figure 5.1**: PCA visualization of resistance profiles (n=187 isolates, 18 antibiotics). PC1 and PC2 explain 32% and 18% of total variance, respectively. Points colored by cluster assignment. High loadings on PC1 reflect variation in β-lactam resistance; PC2 captures aminoglycoside resistance variation."

**Cross-Tab Figure Caption Elements**:
1. Variables cross-tabulated
2. Sample size
3. Statistical test result
4. Interpretation boundary

**Example Caption**:
> "**Figure 5.2**: Cluster distribution by geographic region (n=187 isolates). Chi-square test indicates significant association between clusters and regions (χ²=45.2, df=8, p<0.001). Distribution reflects observed patterns and does not imply regional causation of resistance profiles."

### 4.4 Scope Statement for Regional Analysis

> "Regional and environmental associations are **observational** and reflect patterns in the analyzed dataset. Results do not establish causation between geographic/environmental factors and resistance profiles. Temporal and confounding factors were not controlled."

---

## Document Version

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2024 | Initial multivariate analysis documentation |

---

*This document is part of Phase 8 — Documentation & Reporting for the AMR Thesis Project.*

**See also**: [clustering.md](clustering.md) | [supervised_models.md](supervised_models.md) | [integration.md](integration.md)
