# Threshold Sensitivity Analysis Report

**Generated:** 2026-01-02 20:13:36

## Executive Summary

**Stability Level:** LOW

**Mean Pairwise ARI:** 0.8313
**ARI Range:** [0.5782, 1.0000]

**Interpretation:** ⚠️ Clustering results show LOW STABILITY - they are THRESHOLD-DEPENDENT. Consider conducting further investigation to determine optimal thresholds.

---

## Threshold Combinations Tested

| AB Coverage (%) | Isolate Missing (%) | Isolates Retained | Antibiotics Retained | Silhouette |
|-----------------|---------------------|-------------------|---------------------|------------|
| 50 | 40 | 491 (100.0%) | 22 (100.0%) | 0.4889 |
| 60 | 35 | 491 (100.0%) | 22 (100.0%) | 0.4889 |
| 70 | 30 | 490 (99.8%) | 22 (100.0%) | 0.4879 |
| 80 | 25 | 488 (99.4%) | 22 (100.0%) | 0.4862 |
| 90 | 20 | 488 (99.4%) | 20 (90.9%) | 0.5607 |

---

## Pairwise ARI Matrix

| Threshold | (50.0, 40.0) | (60.0, 35.0) | (70.0, 30.0) | (80.0, 25.0) | (90.0, 20.0) |
|---|---|---|---|---|---|
| (50.0, 40.0)| 1.000 | 1.000 | 1.000 | 1.000 | 0.578 ||
| (60.0, 35.0)| 1.000 | 1.000 | 1.000 | 1.000 | 0.578 ||
| (70.0, 30.0)| 1.000 | 1.000 | 1.000 | 1.000 | 0.578 ||
| (80.0, 25.0)| 1.000 | 1.000 | 1.000 | 1.000 | 0.578 ||
| (90.0, 20.0)| 0.578 | 0.578 | 0.578 | 0.578 | 1.000 ||

---

## Recommendation

**Optimal Threshold:** (90.0, 20.0)

**Rationale:** Based on combined analysis of silhouette score (clustering quality), data retention, and stability vs. baseline, the optimal threshold is (90.0, 20.0).

### Ranked Thresholds

| Rank | Threshold | Silhouette | Retention (%) | Combined Score |
|------|-----------|------------|---------------|----------------|
| 1 | (90.0, 20.0) | 0.5607 | 99.4 | 0.6156 |
| 2 | (50.0, 40.0) | 0.4889 | 100.0 | 0.5179 |
| 3 | (60.0, 35.0) | 0.4889 | 100.0 | 0.5179 |
| 4 | (70.0, 30.0) | 0.4879 | 99.8 | 0.4114 |
| 5 | (80.0, 25.0) | 0.4862 | 99.4 | 0.2000 |

---

## Conclusion

This sensitivity analysis validates the choice of data cleaning thresholds by demonstrating 
that clustering results are robust (or sensitive) to threshold variations. 

For thesis defense, this analysis provides evidence-based justification for the selected 
thresholds and documents the stability of the identified resistance phenotypes.