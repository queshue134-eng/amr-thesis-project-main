# Sensitivity Analysis Report

**Generated:** 2025-12-25 21:17:05

## Overview

This report compares different train-test split ratios and cross-validation 
configurations to assess the stability of supervised validation results.

---

## Split Ratio Comparison

| Split | Model | F1 Score | Accuracy | Stability (std) |
|-------|-------|----------|----------|-----------------|
| 70/30 | Logistic Regression | 0.923 +/- 0.026 | 0.965 | 0.0257 |
| 70/30 | Random Forest | 0.944 +/- 0.016 | 0.974 | 0.0159 |
| 70/30 | KNN | 0.917 +/- 0.034 | 0.964 | 0.0336 |
| 80/20 | Logistic Regression | 0.952 +/- 0.031 | 0.978 | 0.0314 |
| 80/20 | Random Forest | 0.973 +/- 0.022 | 0.988 | 0.0217 |
| 80/20 | KNN | 0.956 +/- 0.028 | 0.980 | 0.0281 |
| 90/10 | Logistic Regression | 0.983 +/- 0.020 | 0.992 | 0.0205 |
| 90/10 | Random Forest | 0.991 +/- 0.018 | 0.996 | 0.0177 |
| 90/10 | KNN | 0.958 +/- 0.026 | 0.980 | 0.0263 |

### Best Split per Model

- **Logistic Regression**: 90/10 (F1=0.983)
- **Random Forest**: 90/10 (F1=0.991)
- **KNN**: 90/10 (F1=0.958)

---

## Cross-Validation Comparison

| CV Folds | Model | F1 Score | Accuracy | Stability (std) |
|----------|-------|----------|----------|-----------------|
| 5-fold | Logistic Regression | 0.953 +/- 0.009 | 0.978 | 0.0090 |
| 5-fold | Random Forest | 0.947 +/- 0.025 | 0.976 | 0.0248 |
| 5-fold | KNN | 0.933 +/- 0.038 | 0.970 | 0.0383 |
| 10-fold | Logistic Regression | 0.932 +/- 0.045 | 0.970 | 0.0451 |
| 10-fold | Random Forest | 0.955 +/- 0.035 | 0.980 | 0.0354 |
| 10-fold | KNN | 0.952 +/- 0.037 | 0.978 | 0.0373 |

---

## Interpretation Guidelines

### Split Ratio Selection
- **70/30**: More test data -> more reliable estimate of generalization, but less training data
- **80/20**: Standard balance (used in main pipeline)
- **90/10**: More training data -> potentially better model, but less reliable test estimate

### Cross-Validation
- **5-fold**: Standard choice, good bias-variance tradeoff
- **10-fold**: Lower bias but higher variance, may be unstable with small datasets

### Stability (std)
- **std < 0.02**: Very stable results
- **std 0.02-0.05**: Acceptable stability  
- **std > 0.05**: Results may be sensitive to random initialization

---

## Recommendation

Based on the analysis, the original **80/20 split** with **5-fold cross-validation** 
provides a good balance of:
1. Sufficient training data for model learning
2. Reliable test set for evaluation
3. Stable results across random seeds

If results are highly consistent across configurations (low std), the choice of 
split ratio has minimal impact on conclusions.
