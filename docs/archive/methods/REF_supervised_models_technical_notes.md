# Supervised Models Documentation

## Phase 4: Supervised Learning for Pattern Discrimination

This document provides comprehensive documentation of supervised learning methods, terminology enforcement, and metric interpretation guidelines for the AMR thesis project.

---

## Table of Contents

1. [Terminology Enforcement](#1-terminology-enforcement)
2. [Model Selection and Rationale](#2-model-selection-and-rationale)
3. [Data Splitting and Leakage Prevention](#3-data-splitting-and-leakage-prevention)
4. [Model Evaluation Metrics](#4-model-evaluation-metrics)
5. [Feature Importance Analysis](#5-feature-importance-analysis)
6. [Reporting Standards](#6-reporting-standards)

---

## 1. Terminology Enforcement

### 1.1 Language Control Table (CRITICAL)

| ✅ ALLOWED | ❌ NOT ALLOWED | Rationale |
|------------|----------------|-----------|
| **Discrimination** | ~~Prediction~~ | We assess pattern consistency, not future outcomes |
| **Consistency** | ~~Accuracy on unseen data~~ | Test set evaluates generalization of patterns |
| **Alignment** | ~~Forecast~~ | Resistance fingerprints align with categories |
| **Pattern discrimination** | ~~Predictive power~~ | Descriptive, not prognostic |
| **Discriminative capacity** | ~~Predictive accuracy~~ | Within-dataset evaluation |
| **Model shows consistent alignment** | ~~Model performs well~~ | Neutral, precise language |
| **Demonstrates discriminative capacity** | ~~Predicts accurately~~ | Pattern-focused phrasing |
| **Strong pattern consistency** | ~~High predictive accuracy~~ | Avoids forecasting implications |
| **Associative importance** | ~~Predictive importance~~ | Features help discriminate, not predict |

### 1.2 Metric Interpretation Template (USE EVERY TIME)

> **Standard Statement (use verbatim)**:
> "Evaluation metrics quantify the consistency with which resistance fingerprints align with predefined categories. These metrics assess **pattern discrimination** within the analyzed dataset and are **not** predictive of outcomes in external datasets."

### 1.3 Forbidden Terms Checklist

Before submitting any results section, search for and replace:

| Search For | Replace With |
|------------|--------------|
| `predict*` | `discriminate`, `align`, `separate` |
| `forecast*` | `characterize`, `describe` |
| `classify (implying deployment)` | `discriminate`, `distinguish` |
| `accuracy (implying prediction)` | `alignment consistency`, `discrimination rate` |
| `model performs well` | `model shows consistent alignment` |

---

## 2. Model Selection and Rationale

### 2.1 Rationalized Model Set

| Model | Category | Key Hyperparameters | Purpose |
|-------|----------|---------------------|---------|
| **Logistic Regression** | Linear | max_iter=1000, random_state=42 | Linear baseline with coefficient interpretation |
| **Random Forest** | Tree-based | n_estimators=100, random_state=42 | Nonlinear model with Gini feature importance |
| **k-Nearest Neighbors** | Distance-based | n_neighbors=5 | Distance-based consistency check |

### 2.2 Model Category Descriptions

| Category | Description | Interpretability |
|----------|-------------|------------------|
| **Linear** | Learns linear decision boundaries in feature space | High - coefficient magnitudes indicate feature contributions |
| **Tree-based** | Partitions feature space using decision rules | Medium - Gini importance shows feature split contribution |
| **Distance-based** | Classifies based on similarity to training instances | Low - no native feature importance |

### 2.3 Why This Model Set

| Design Choice | Rationale |
|---------------|-----------|
| **Three models only** | Sufficient diversity without excessive complexity |
| **One per category** | Represents different learning paradigms |
| **Interpretable options** | Supports biological interpretation of results |
| **Standard implementations** | Reproducible with scikit-learn defaults |

---

## 3. Data Splitting and Leakage Prevention

### 3.1 Leakage-Safe Pipeline

**CRITICAL**: Split **BEFORE** any preprocessing to prevent data leakage.

```
┌─────────────────────────────────────────────────────────────┐
│              LEAKAGE-SAFE PREPROCESSING ORDER               │
├─────────────────────────────────────────────────────────────┤
│                                                             │
│  Step 1: SPLIT FIRST                                        │
│  ├── 80% Training Set (stratified by target)                │
│  └── 20% Test Set (stratified by target)                    │
│                                                             │
│  Step 2: FIT ON TRAIN ONLY                                  │
│  ├── Imputer.fit(X_train)                                   │
│  └── Scaler.fit(X_train)                                    │
│                                                             │
│  Step 3: TRANSFORM BOTH                                     │
│  ├── X_train = Imputer.transform(X_train)                   │
│  ├── X_train = Scaler.transform(X_train)                    │
│  ├── X_test = Imputer.transform(X_test)                     │
│  └── X_test = Scaler.transform(X_test)                      │
│                                                             │
└─────────────────────────────────────────────────────────────┘
```

### 3.2 Split Parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| **Train-Test Ratio** | 80%-20% | Standard split for model evaluation |
| **Stratification** | By target variable | Preserves class distribution |
| **Random State** | 42 (fixed) | Reproducibility |

### 3.3 Strict Feature-Label Separation

**Input Features** (ALLOWED):
- Encoded resistance fingerprints: `{ANTIBIOTIC}_encoded` columns only

**Excluded Variables** (NEVER USE AS FEATURES):
| Variable | Reason |
|----------|--------|
| REGION | Metadata - would cause contextual leakage |
| SITE | Metadata - would cause contextual leakage |
| SPECIES | Target variable in species discrimination task |
| ENVIRONMENT | Metadata - would cause contextual leakage |
| SAMPLE_SOURCE | Metadata - would cause contextual leakage |
| CODE | Identifier - not a feature |

### 3.4 MDR Target Transparency Statement

> **IMPORTANT**: The MDR label is derived from the SAME resistance features used as input. MDR discrimination is treated as **self-consistency discrimination**—evaluating how consistently resistance fingerprints align with MDR status. This explicitly acknowledges the "predicted MDR from MDR features" relationship and is **not** predictive in nature.

---

## 4. Model Evaluation Metrics

### 4.1 Metrics Definitions

| Metric | Formula | Interpretation |
|--------|---------|----------------|
| **Accuracy** | (TP + TN) / Total | Overall correct discriminations |
| **Precision** | TP / (TP + FP) | Of positive discriminations, how many correct |
| **Recall** | TP / (TP + FN) | Of actual positives, how many correctly discriminated |
| **F1-Score** | 2 × (P × R) / (P + R) | Harmonic mean of precision and recall |

### 4.2 Averaging Method: MACRO

**Definition**:
```
Macro_Metric = (1/n_classes) × Σ(class_metric)
```

**Rationale**: All classes treated equally regardless of sample size, preventing class imbalance bias.

### 4.3 Confusion Matrix Interpretation

| Actual \ Discriminated | Positive | Negative |
|------------------------|----------|----------|
| **Positive** | True Positive (TP) | False Negative (FN) |
| **Negative** | False Positive (FP) | True Negative (TN) |

### 4.4 Results Scope Statement (MANDATORY)

**Add to every results table**:
> "**Scope**: Results reflect pattern discrimination within the analyzed dataset and are **not** predictive of performance on external datasets."

---

## 5. Feature Importance Analysis

### 5.1 Importance Extraction by Model Type

| Model | Method | Interpretation |
|-------|--------|----------------|
| **Random Forest** | Gini importance (mean decrease in impurity) | Contribution to splits |
| **Logistic Regression** | Absolute coefficient magnitude | Weight in linear combination |
| **k-NN** | Not available | Distance-based; no feature weights |

### 5.2 Interpretation Guidelines

| Statement Type | ✅ Allowed | ❌ Not Allowed |
|----------------|------------|----------------|
| Importance claim | "AM shows high **associative importance** for MDR discrimination" | ~~"AM **predicts** MDR status"~~ |
| Ranking | "Top discriminative antibiotics were AM, AMC, GM" | ~~"Most predictive antibiotics"~~ |
| Biological link | "High importance may relate to common resistance mechanisms" | ~~"High importance indicates AM causes MDR"~~ |

### 5.3 Biological Restraint Statement

> "Feature importance scores indicate **associative importance**—the degree to which each antibiotic contributes to group separation. High importance does **NOT** imply causation or mechanism. Biological interpretation should consider known resistance mechanisms and clinical context."

### 5.4 Feature Importance Table Template

| Rank | Antibiotic | Importance Score | Task | Interpretation |
|------|------------|------------------|------|----------------|
| 1 | [AB] | [Score] | MDR/Species | [Biological context] |
| 2 | [AB] | [Score] | MDR/Species | [Biological context] |
| ... | ... | ... | ... | ... |

---

## 6. Reporting Standards

### 6.1 Task Separation

Each discrimination task runs as an **independent experiment**:

| Task | Target | Type | Separate Report |
|------|--------|------|-----------------|
| Task A: Species discrimination | SPECIES | Multi-class | Yes |
| Task B: MDR discrimination | MDR_FLAG | Binary | Yes |

### 6.2 Results Section Template

```markdown
### [Task Name] Discrimination Results

**Objective**: Evaluate how resistance fingerprints discriminate [target variable].

**Method**: [Model name] trained on 80% stratified sample, evaluated on held-out 20%.

**Scope**: Results reflect pattern discrimination within analyzed dataset.

#### Evaluation Metrics

| Model | Accuracy | Precision | Recall | F1-Score |
|-------|----------|-----------|--------|----------|
| [Model 1] | [value] | [value] | [value] | [value] |
| [Model 2] | [value] | [value] | [value] | [value] |

**Interpretation**: Evaluation metrics quantify the consistency with which resistance fingerprints align with predefined categories.

#### Feature Importance

[Feature importance table]

**Note**: Importance scores indicate associative contribution to group separation, not predictive or causal relationships.
```

### 6.3 Figure Caption Requirements

Every supervised learning figure must include:

1. **Data source**: "Based on N isolates with [target] labels"
2. **Method**: Model type and key parameters
3. **Interpretation boundary**: Discrimination vs. prediction distinction

**Example Caption**:
> "**Figure 4.1**: Confusion matrix for MDR discrimination using Random Forest classifier. Results reflect pattern consistency within the analyzed dataset (n=187 isolates, 80/20 train-test split). Metrics quantify alignment between resistance fingerprints and MDR categories, not predictive performance on external data."

---

## Document Version

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2024 | Initial supervised models documentation |

---

*This document is part of Phase 8 — Documentation & Reporting for the AMR Thesis Project.*

**See also**: [clustering.md](clustering.md) | [multivariate_analysis.md](multivariate_analysis.md) | [METHODOLOGY.md](../METHODOLOGY.md)
