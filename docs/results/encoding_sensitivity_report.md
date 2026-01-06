# Encoding Sensitivity Analysis Report

**Generated:** 2026-01-07 06:30:45

## Executive Summary

**Stability Level:** MODERATE

**Mean Pairwise ARI:** 0.9043
**ARI Range:** [0.7690, 1.0000]

**Interpretation:** Clustering results show MODERATE STABILITY. Most encoding schemes produce similar groupings, though binary encoding may differ.

---

## Encoding Schemes Tested

| Encoding | S | I | R | Silhouette | Description |
|----------|---|---|---|------------|-------------|
| linear | 0 | 1 | 2 | 0.4889 | Default linear encoding (equal intervals) |
| clinical_weighted | 0.0 | 1.5 | 3.0 | 0.4889 | Clinical-weighted encoding (larger I→R gap) |
| binary | 0 | 1 | 1 | 0.3878 | Binary encoding (I=R, non-susceptible) |
| threshold_normalized | 0.0 | 0.5 | 1.0 | 0.4889 | Normalized threshold encoding [0-1] |
| exponential | 0 | 1 | 4 | 0.5372 | Exponential encoding (emphasizes R) |

---

## ARI vs Baseline (Linear)

| Encoding | ARI vs Baseline |
|----------|-----------------|
| linear | 1.0000 |
| clinical_weighted | 1.0000 |
| binary | 0.7923 |
| threshold_normalized | 1.0000 |
| exponential | 0.9656 |

---

## Recommendation

**Recommended Encoding:** linear

**Rationale:** Linear encoding is recommended (default) as results are stable. Best silhouette achieved by 'exponential' (0.5372).

---

## Conclusion

This analysis validates the robustness of clustering results to different 
resistance encoding schemes. For thesis defense, this provides evidence that 
the identified resistance phenotypes are not artifacts of the encoding choice.