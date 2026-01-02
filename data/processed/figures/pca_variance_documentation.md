# PCA Variance Analysis Documentation

> **AMR Thesis Project — PCA Variance Explained Report**  
> **Last Updated:** December 20, 2025  
> **Data Source:** INOHAC AMR Project Two

---

## Summary Statistics

| Metric | Value |
|--------|-------|
| PC1 Explained Variance | 23.5% |
| PC2 Explained Variance | 16.4% |
| PC1 + PC2 Cumulative | 39.9% |
| PCs for 50% variance | 3 |
| PCs for 75% variance | 7 |
| PCs for 90% variance | 12 |

---

## Interpretation

The first two principal components capture only **39.9%** of total resistance variation. This limited proportion suggests 2D plots should be interpreted cautiously as they represent only a fraction of the full 22-dimensional resistance space.

> **Limitation Acknowledgment:** Full resistance space is multi-dimensional; these projections emphasize the two dominant axes of variation but may not fully represent cluster separation in higher dimensions.

---

## PC Loadings Interpretation

### PC1 Axis

| Characteristic | Value |
|----------------|-------|
| Primary driver | ENR (Enrofloxacin) |
| Axis interpretation | Dominant resistance variation |

### PC2 Axis

| Characteristic | Value |
|----------------|-------|
| Primary driver | CF (Cephalothin) |
| Axis interpretation | Secondary patterns orthogonal to PC1 |

---

## Figure Caption Template

> Principal Component Analysis of resistance profiles (n=491 isolates, 22 antibiotics, 6 species). **PC1 explains 23.5% of variance, PC2 explains 16.4% (cumulative 39.9%)**. Points colored by [cluster/region/environment].

---

## Variance Breakdown by Component

| Component | Individual (%) | Cumulative (%) |
|-----------|----------------|----------------|
| PC1 | 23.5 | 23.5 |
| PC2 | 16.4 | 39.9 |
| PC3 | 9.1 | 49.0 |
| PC4 | 7.2 | 56.2 |
| PC5 | 5.9 | 62.1 |

---

*This document is auto-generated as part of the AMR Thesis Project pipeline.*

**See also**: [TECHNICAL_REFERENCE.md](../../docs/TECHNICAL_REFERENCE.md) | [REF_multivariate_analysis_technical_notes.md](../../docs/archive/methods/REF_multivariate_analysis_technical_notes.md)
