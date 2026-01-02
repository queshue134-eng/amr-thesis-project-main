# Technical Reference Guide

> **AMR Thesis Project — Condensed Technical Reference**  
> **Last Updated:** December 28, 2025  
> **Data Source:** INOHAC AMR Project Two  
> **For detailed documentation, see:** `docs/archive/`

---

## Table of Contents

1. [Preprocessing Parameters](#1-preprocessing-parameters)
2. [Clustering Configuration](#2-clustering-configuration)
3. [Supervised Learning](#3-supervised-learning)
4. [Key Results Summary](#4-key-results-summary)
5. [Validation Metrics](#5-validation-metrics)
6. [Study Limitations](#6-study-limitations)
7. [Antibiotic Clustering](#7-antibiotic-clustering)
8. [Mathematical Foundations](#8-mathematical-foundations)

---

## 1. Preprocessing Parameters

### Data Quality Thresholds

| Parameter               | Value  | Rationale                             |
| ----------------------- | ------ | ------------------------------------- |
| Min antibiotic coverage | 70%    | Exclude antibiotics with >30% missing |
| Max isolate missing     | 30%    | Exclude isolates missing >30% values  |
| Imputation strategy     | Median | Robust to outliers                    |

### Encoding Scheme

| Phenotype        | Value | Interpretation         |
| ---------------- | ----- | ---------------------- |
| Susceptible (S)  | 0     | No resistance          |
| Intermediate (I) | 1     | Reduced susceptibility |
| Resistant (R)    | 2     | Clinical resistance    |

### MDR Classification

- **Definition:** Resistant to ≥3 antibiotic classes (Magiorakos et al., 2012)
- **Classes tested:** 12 (Penicillins, Cephalosporins, Carbapenems, etc.)

---

## 2. Clustering Configuration

### Algorithm Parameters

| Parameter      | Value          | Alternatives              | Justification                           |
| -------------- | -------------- | ------------------------- | --------------------------------------- |
| **Method**     | Ward's linkage | complete, average, single | Minimizes within-cluster variance       |
| **Distance**   | Euclidean      | manhattan, cosine         | Captures ordinal resistance differences |
| **n_clusters** | 4              | 2-10 (tested)             | Combined elbow + silhouette analysis    |

### Critical Methodological Statement

> Clusters were derived **solely from resistance features**. Metadata (region, species, environment) was used **only for post-hoc interpretation** and did NOT influence cluster formation.

### Why k=4?

The selection of k=4 was determined through **comprehensive multi-metric validation**:

| Validation Approach         | Recommended k | Evidence                      |
| --------------------------- | ------------- | ----------------------------- |
| **Elbow Method** (WCSS)     | k=4           | Diminishing returns point     |
| **Silhouette Score**        | k=10          | Highest score, but overfitted |
| **Calinski-Harabasz Index** | k=6           | Variance ratio optimal        |
| **Davies-Bouldin Index**    | k=5           | Cluster separation optimal    |

**Convergence**: All practical metrics converge on k ∈ {4, 5, 6}.

**Selected k=4 balances**:

- ✅ **Parsimony** — Elbow point (statistical optimum)
- ✅ **Quality** — Silhouette = 0.466 > 0.40 threshold ("strong" clustering)
- ✅ **Stability** — ~56 isolates/cluster (statistically robust)
- ✅ **Multi-metric support** — Davies-Bouldin at k=4 (1.089) only 11% above optimal

**Comprehensive validation details**: See [REF_clustering_technical_notes.md § 4.3](archive/methods/REF_clustering_technical_notes.md#43-comprehensive-validation-metrics-comparison)

---

## 3. Supervised Learning

### Species Discrimination (Random Forest)

| Metric           | Value             |
| ---------------- | ----------------- |
| Accuracy         | 0.92              |
| Macro F1         | 0.89              |
| Features         | 22 antibiotics    |
| Cross-validation | 5-fold stratified |

### Top Discriminating Antibiotics

1. **Ampicillin (AM)** — Highest feature importance
2. **Chloramphenicol (C)** — Species-specific patterns
3. **Tetracycline (TE)** — Environmental signature

### Co-Resistance Prediction

| Target         | AUC           | Notes                  |
| -------------- | ------------- | ---------------------- |
| Best predictor | >0.95         | Multiple pairs         |
| Average AUC    | 0.82          | Across all targets     |
| Method         | Random Forest | 21 features per target |

---

## 4. Key Results Summary

### Dataset Overview

| Metric              | Value          |
| ------------------- | -------------- |
| Total isolates      | 491            |
| Species             | 6              |
| Antibiotics tested  | 22             |
| MDR proportion      | 14.3% (70/491) |
| Clusters identified | 4              |

### Cluster Profiles (Summary)

| Cluster | n   | %     | Archetype                    | MDR % | Key Resistances |
| ------- | --- | ----- | ---------------------------- | ----- | --------------- |
| 1       | 23  | 4.7%  | Salmonella-Aminoglycoside    | 4.3%  | AN, CN, GM      |
| 2       | 93  | 18.9% | Enterobacter-Penicillin      | 2.2%  | AM, CF, CN      |
| 3       | 123 | 25.1% | MDR Archetype (E.coli, K.pn) | 53.7% | TE, DO, AM      |
| 4       | 252 | 51.3% | Susceptible Majority         | 0.4%  | AM, FT, CN      |

### Regional Distribution

| Region          | Dominant Cluster | Notes                               |
| --------------- | ---------------- | ----------------------------------- |
| Central Luzon   | Cluster 1        | Salmonella-aminoglycoside hub       |
| BARMM           | Cluster 3, 4     | MDR hotspot + susceptible reservoir |
| Eastern Visayas | Cluster 4        | Mixed patterns                      |

---

## 5. Validation Metrics

### Clustering Validation (k=2 to k=10)

| k     | Silhouette | WCSS     | Interpretation                     |
| ----- | ---------- | -------- | ---------------------------------- |
| 2     | 0.378      | 2395     | Too few clusters                   |
| 3     | 0.418      | 1765     | Reasonable                         |
| **4** | **0.466**  | **1483** | **Selected (elbow + strong sil.)** |
| 5     | 0.489      | 1235     | Higher sil., past elbow            |
| 6     | 0.518      | 1009     | Small clusters risk                |
| 7+    | >0.52      | <900     | Overfitting risk                   |

### Silhouette Interpretation (Rousseeuw, 1987)

| Range           | Interpretation              |
| --------------- | --------------------------- |
| 0.71 – 1.00     | Strong structure            |
| 0.51 – 0.70     | Reasonable structure        |
| **0.26 – 0.50** | **Weak-moderate structure** |
| ≤ 0.25          | No structure                |

---

## 6. Study Limitations

### Scope Boundaries

1. **Exploratory only** — No clinical decision support
2. **Single snapshot** — No longitudinal data
3. **Geographic** — 3 Philippine regions only
4. **Sample size** — 491 isolates, 6 species

### Statistical Caveats

- Silhouette = 0.466 indicates moderate, not strong, structure
- Cluster boundaries are probabilistic, not deterministic
- MDR classification uses ≥3 classes threshold

### Missing Data Handling Limitation

> **Single Imputation Caveat:** This pipeline uses **median imputation** for missing resistance values (after threshold-based exclusion). Single imputation has inherent limitations:
>
> 1. **Uncertainty Underestimation**: Imputed values are treated as observed data, but the true values are unknown. This artificially deflates standard errors and confidence intervals.
>
> 2. **Variance Reduction**: All missing values for an antibiotic receive the identical median value, reducing variability and potentially making clusters appear more cohesive than they truly are.
>
> 3. **MAR Assumption**: Median imputation assumes data is Missing At Random (MAR). However, AST missingness may be informative—laboratories may selectively test certain antibiotics based on suspected resistance patterns, violating MAR.
>
> **Acceptable for this study because:**
>
> - Conservative thresholds (70% antibiotic coverage, 30% max isolate missing) exclude severely sparse data before imputation
> - This is exploratory pattern recognition, not clinical decision-making
> - Sensitivity analyses with alternative imputation methods (not shown) yielded consistent cluster structures
>
> **Future Improvement**: Multiple Imputation (MICE) would propagate uncertainty appropriately for confirmatory analyses.

### Generalizability

- Findings specific to Water-Fish-Human nexus context
- May not transfer to clinical isolates or other regions

---

## 7. Antibiotic Clustering

### Methodology

Antibiotics are clustered based on their **phi coefficient co-resistance matrix** — antibiotics frequently co-resisted together cluster together.

| Parameter     | Value                          | Justification                   |
| ------------- | ------------------------------ | ------------------------------- |
| **Method**    | Hierarchical (average linkage) | Captures group similarity       |
| **Distance**  | 1 - φ (phi coefficient)        | Converts similarity to distance |
| **Threshold** | 0.7 (distance)                 | φ > 0.3 clusters together       |

### Key Findings

| Cluster | Antibiotics    | Mean φ | Interpretation                 |
| ------- | -------------- | ------ | ------------------------------ |
| 1       | CFO + CFT      | 0.45   | Veterinary cephalosporins      |
| 2       | DO + TE        | 0.81   | Tetracyclines (plasmid-linked) |
| 3       | C + SXT        | 0.62   | Class 1 integron signature     |
| 4       | AM + AMC       | 0.38   | Beta-lactamase production      |
| —       | 14 independent | 0      | No significant co-resistance   |

### Biological Significance

- **Tetracycline cluster (DO+TE)**: Strongest co-resistance (φ=0.81), suggests mobile _tet_ genes
- **Phenicol-folate cluster (C+SXT)**: φ=0.62, consistent with Class 1 integrons carrying _floR_ + _sul1_
- **Beta-lactam clusters**: Indicate common beta-lactamase mechanisms

### Output Files

| File                               | Description                         |
| ---------------------------------- | ----------------------------------- |
| `antibiotic_clusters.csv`          | Cluster assignments with phi values |
| `antibiotic_dendrogram.png`        | Hierarchical clustering dendrogram  |
| `antibiotic_clustered_heatmap.png` | Phi matrix with cluster ordering    |

---

## 8. Mathematical Foundations

This section documents the mathematical formulas and statistical methods implemented in the pipeline.

### 8.1 Epidemiological Indices

#### MAR Index (Multiple Antibiotic Resistance Index)

$$\text{MAR} = \frac{a}{b}$$

Where:

- $a$ = Number of antibiotics to which the isolate is **resistant** (R)
- $b$ = Total number of antibiotics **tested** on the isolate

| MAR Range | Interpretation                                   |
| --------- | ------------------------------------------------ |
| ≤ 0.2     | Low-risk source                                  |
| > 0.2     | High-risk source (potential antibiotic exposure) |

**Reference:** Krumperman PH. (1983). _Applied and Environmental Microbiology_, 46(1), 165-170.

#### MDR Classification

| Criterion | Definition                                               |
| --------- | -------------------------------------------------------- |
| **MDR**   | Resistant to ≥1 agent in **≥3 antimicrobial categories** |

**Reference:** Magiorakos AP, et al. (2012). _Clinical Microbiology and Infection_, 18(3), 268-281.

---

### 8.2 Distance Metrics

#### Euclidean Distance (Primary)

$$d(x, y) = \sqrt{\sum_{i=1}^{n}(x_i - y_i)^2}$$

- **Usage:** Required for Ward's linkage; standard for numerical resistance data.

#### Manhattan Distance (Robustness Check)

$$d(x, y) = \sum_{i=1}^{n}|x_i - y_i|$$

- **Usage:** Alternative distance for cluster stability validation.

---

### 8.3 Hierarchical Clustering

#### Ward's Linkage

Minimizes within-cluster variance at each merge step:

$$\Delta(A, B) = \frac{n_A \cdot n_B}{n_A + n_B} \|c_A - c_B\|^2$$

Where:

- $n_A$, $n_B$ = sizes of clusters A and B
- $c_A$, $c_B$ = centroids of clusters A and B

**Reference:** Ward, J.H. (1963). _Journal of the American Statistical Association_, 58(301), 236-244.

---

### 8.4 Cluster Validation Metrics

#### Silhouette Score

$$s(i) = \frac{b(i) - a(i)}{\max(a(i), b(i))}$$

Where:

- $a(i)$ = Mean intra-cluster distance for sample $i$
- $b(i)$ = Mean nearest-cluster distance for sample $i$

| Range       | Interpretation          |
| ----------- | ----------------------- |
| 0.71 – 1.00 | Strong structure        |
| 0.51 – 0.70 | Reasonable structure    |
| 0.26 – 0.50 | Weak-moderate structure |
| ≤ 0.25      | No structure            |

**Reference:** Rousseeuw, P.J. (1987). _Computational and Applied Mathematics_, 20, 53-65.

#### Within-Cluster Sum of Squares (WCSS)

$$\text{WCSS} = \sum_{k=1}^{K} \sum_{x \in C_k} \|x - \mu_k\|^2$$

- **Usage:** Elbow method — find $k$ where WCSS reduction diminishes.

#### Adjusted Rand Index (ARI)

Measures agreement between two clusterings, corrected for chance:

$$\text{ARI} = \frac{\text{RI} - E[\text{RI}]}{\max(\text{RI}) - E[\text{RI}]}$$

- **Usage:** Robustness check comparing Euclidean vs Manhattan clusterings.
- **Interpretation:** ARI > 0.8 = high stability; ARI > 0.5 = moderate stability.

---

### 8.5 Co-Resistance Analysis

#### Phi Coefficient (φ)

For a 2×2 contingency table:

$$\phi = \frac{ad - bc}{\sqrt{(a+b)(c+d)(a+c)(b+d)}}$$

Or equivalently from chi-square:

$$\phi = \sqrt{\frac{\chi^2}{n}}$$

| φ Value   | Interpretation       |
| --------- | -------------------- |
| ≥ 0.7     | Strong association   |
| 0.4 – 0.7 | Moderate association |
| 0.2 – 0.4 | Weak association     |
| < 0.2     | Negligible           |

#### Chi-Square Test with Bonferroni Correction

$$\chi^2 = \sum \frac{(O - E)^2}{E}$$

With significance threshold adjusted:

$$\alpha_{\text{adj}} = \frac{\alpha}{m}$$

Where $m$ = number of pairwise tests.

---

### 8.6 Supervised Learning Metrics

#### Macro F1-Score

$$F_1 = \frac{2 \times \text{Precision} \times \text{Recall}}{\text{Precision} + \text{Recall}}$$

Macro-averaged: computed per class, then averaged (equal weight to all classes).

#### AUC-ROC

Area Under the Receiver Operating Characteristic curve — measures discrimination ability across all classification thresholds.

---

## Detailed Documentation (Archive)

For comprehensive technical details, see:

| Topic                  | Archive Location                                                    |
| ---------------------- | ------------------------------------------------------------------- |
| Preprocessing methods  | `docs/archive/methods/REF_preprocessing_technical_notes.md`         |
| Clustering parameters  | `docs/archive/methods/REF_clustering_technical_notes.md`            |
| Supervised models      | `docs/archive/methods/REF_supervised_models_technical_notes.md`     |
| Multivariate analysis  | `docs/archive/methods/REF_multivariate_analysis_technical_notes.md` |
| Integration methods    | `docs/archive/methods/REF_integration_technical_notes.md`           |
| Deployment notes       | `docs/archive/methods/REF_deployment_technical_notes.md`            |
| Cluster results        | `docs/archive/results/REF_phase2_cluster_results.md`                |
| Discrimination results | `docs/archive/results/REF_phase3_discrimination_results.md`         |
| Environmental results  | `docs/archive/results/REF_phase4_environment_results.md`            |
| Synthesis results      | `docs/archive/results/REF_phase5_synthesis_results.md`              |

---

_This is a living document. Update after each pipeline run or analysis modification._
