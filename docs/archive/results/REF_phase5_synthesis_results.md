# Phase 5 Results: Integration and Synthesis

## Integrated Findings and Synthesis Results

This document provides the actual Phase 5 (integration and synthesis) results for the AMR thesis project.

---

## Table of Contents

1. [Resistance Archetype Summary](#1-resistance-archetype-summary)
2. [Species-Environment Associations](#2-species-environment-associations)
3. [MDR Enrichment Analysis](#3-mdr-enrichment-analysis)
4. [Cross-Phase Coherence](#4-cross-phase-coherence)
5. [Final Documentation Checklist](#5-final-documentation-checklist)

---

## 1. Resistance Archetype Summary

> **Building on the resistance-based clusters identified in Phase 3**, this section characterizes the five resistance archetypes discovered through unsupervised structure identification.

### 1.1 Archetype Overview Table

| Cluster | N | MDR Rate | MAR Index | Archetype Name | Key Characteristics |
|---------|---|----------|-----------|----------------|---------------------|
| 1 | 23 | 26.1% | 0.196 | **Salmonella-Aminoglycoside** | Aminoglycoside-resistant, water-associated |
| 2 | 93 | 20.4% | 0.164 | **Enterobacter-Beta-lactam** | Ampicillin/cephalosporin resistance |
| 3 | 123 | **54.5%** | 0.181 | **E. coli-Tetracycline-MDR** | Tetracycline-dominated, MDR hotspot |
| 4 | 104 | **0.0%** | 0.002 | **E. coli-Susceptible** | Pan-susceptible, commensal phenotype |
| 5 | 148 | 1.4% | 0.053 | **Klebsiella-Intermediate** | Low-level ampicillin, minimal resistance |

### 1.2 Detailed Archetype Profiles

#### Archetype 1: Salmonella-Aminoglycoside

| Attribute | Value |
|-----------|-------|
| **Cluster** | 1 |
| **Size** | 23 isolates (4.7% of dataset) |
| **MDR Rate** | 26.1% |
| **Mean MAR** | 0.196 |
| **Resistance Level** | Moderate |
| **Dominant Species** | *Salmonella species* (100.0%) |
| **Primary Environment** | Water (69.6%) |
| **Primary Region** | Region III - Central Luzon (73.9%) |

**High-Resistance Antibiotics**: AN, CN, GM (aminoglycosides)

**Characterization**: Small cluster dominated by *Salmonella species* from water environments, particularly in Central Luzon. Aminoglycoside resistance suggests environmental exposure or specific resistance plasmids.

---

#### Archetype 2: Enterobacter-Beta-lactam

| Attribute | Value |
|-----------|-------|
| **Cluster** | 2 |
| **Size** | 93 isolates (18.9% of dataset) |
| **MDR Rate** | 20.4% |
| **Mean MAR** | 0.164 |
| **Resistance Level** | Moderate |
| **Dominant Species** | *Enterobacter cloacae* (70.2%) |
| **Primary Environment** | Fish (53.2%) |
| **Primary Region** | Region III - Central Luzon (53.2%) |

**High-Resistance Antibiotics**: AM, CF, CN (ampicillin, cephalosporins)

**Characterization**: *Enterobacter*-dominated cluster with beta-lactam resistance pattern. Mixed fish/water sources. AmpC-type resistance consistent with intrinsic *Enterobacter* mechanisms.

---

#### Archetype 3: E. coli-Tetracycline-MDR (MDR HOTSPOT)

| Attribute | Value |
|-----------|-------|
| **Cluster** | 3 |
| **Size** | 123 isolates (25.0% of dataset) |
| **MDR Rate** | **54.5%** (2.82× enrichment) |
| **Mean MAR** | 0.181 |
| **Resistance Level** | **High** |
| **Dominant Species** | *Escherichia coli* (77.2%) |
| **Primary Environment** | Fish (56.1%) |
| **Primary Region** | BARMM (53.7%) |

**High-Resistance Antibiotics**: TE, DO, AM, SXT (tetracyclines, ampicillin, SXT)

**Characterization**: **Primary MDR hotspot**. Tetracycline-dominated *E. coli* cluster associated with aquaculture. High co-resistance suggests mobile genetic element carriage. May reflect tetracycline use in fish farming.

---

#### Archetype 4: E. coli-Susceptible

| Attribute | Value |
|-----------|-------|
| **Cluster** | 4 |
| **Size** | 104 isolates (21.1% of dataset) |
| **MDR Rate** | **0.0%** |
| **Mean MAR** | 0.002 |
| **Resistance Level** | **Susceptible** |
| **Dominant Species** | *Escherichia coli* (98.1%) |
| **Primary Environment** | Fish (58.7%) |
| **Primary Region** | BARMM (54.8%) |

**High-Resistance Antibiotics**: None (broadly susceptible)

**Characterization**: Nearly pan-susceptible *E. coli* cluster. Despite sharing environment and region with MDR Cluster 3, shows no resistance. May represent commensal strains unexposed to antibiotic selection, or distinct clonal lineages.

---

#### Archetype 5: Klebsiella-Intermediate

| Attribute | Value |
|-----------|-------|
| **Cluster** | 5 |
| **Size** | 148 isolates (30.1% of dataset) |
| **MDR Rate** | 1.4% |
| **Mean MAR** | 0.053 |
| **Resistance Level** | Low |
| **Dominant Species** | *Klebsiella pneumoniae* (77.0%) |
| **Primary Environment** | Fish (58.8%) |
| **Primary Region** | BARMM (58.1%) |

**High-Resistance Antibiotics**: Minimal (low-level AM, FT)

**Characterization**: *K. pneumoniae*-dominated cluster with minimal resistance. Largest cluster (30% of dataset). May represent environmental baseline for *Klebsiella*.

---

### 1.3 Archetype Comparison

| Characteristic | Arch 1 | Arch 2 | Arch 3 | Arch 4 | Arch 5 |
|----------------|--------|--------|--------|--------|--------|
| Resistance Level | Moderate | Moderate | **High** | **Susceptible** | Low |
| MDR Rate | 26.1% | 20.4% | **54.5%** | **0.0%** | 1.4% |
| Primary Class | Aminoglycosides | Beta-lactams | Tetracyclines | None | Minimal |
| Dominant Species | Salmonella | E. cloacae | E. coli | E. coli | K. pneumoniae |
| Primary Environment | Water | Fish/Water | Fish | Fish | Fish |

---

## 2. Species-Environment Associations

### 2.1 Species Distribution by Environment

| Species | Fish | Water | Hospital | Total | Primary Env |
|---------|------|-------|----------|-------|-------------|
| E. coli | 136 | 72 | 19 | 227 | Fish (59.9%) |
| K. pneumoniae | 90 | 46 | 13 | 149 | Fish (60.4%) |
| E. cloacae | 36 | 25 | 7 | 68 | Fish (52.9%) |
| E. aerogenes | 8 | 14 | 1 | 23 | Water (60.9%) |
| Salmonella species | 3 | 19 | 1 | 23 | Water (82.6%) |
| Vibrio vulnificus | 1 | 0 | 0 | 1 | Fish (100%) |
| **Total** | **274** | **176** | **41** | **491** | - |

**Chi-Square Test**: χ² = 98.67, df = 18, p < 10⁻¹⁰
**Cramér's V**: 0.316 (Medium effect)

### 2.2 Species MDR Rates

| Species | Total | MDR Count | MDR Rate | Fold Enrichment |
|---------|-------|-----------|----------|-----------------|
| Vibrio vulnificus | 1 | 1 | 100.0%* | 5.23× |
| Salmonella species | 23 | 6 | 26.1% | 1.37× |
| E. coli | 227 | 51 | 22.5% | 1.18× |
| E. cloacae | 68 | 12 | 17.6% | 0.92× |
| K. pneumoniae | 149 | 21 | 14.1% | 0.74× |
| E. aerogenes | 23 | 3 | 13.0% | 0.68× |
| **Overall** | **491** | **94** | **19.1%** | **1.0×** |

> *Vibrio vulnificus: Only 1 isolate, 100% MDR not statistically meaningful.

---

## 3. MDR Enrichment Analysis

### 3.1 MDR Enrichment by Category

#### By Cluster

| Cluster | N | MDR Rate | Fold Enrichment | Enrichment Level |
|---------|---|----------|-----------------|------------------|
| Cluster 3 | 123 | **54.5%** | **2.82×** | **High** |
| Cluster 1 | 23 | 26.1% | 1.35× | Moderate |
| Cluster 2 | 93 | 20.4% | 1.07× | Low |
| Cluster 5 | 148 | 1.4% | 0.07× | Depleted |
| Cluster 4 | 104 | 0.0% | 0.00× | None |

#### By Region

| Region | N | MDR Rate | Fold Enrichment |
|--------|---|----------|-----------------|
| Central Luzon | 140 | 21.4% | 1.11× |
| BARMM | 250 | 20.8% | 1.08× |
| Eastern Visayas | 102 | 12.7% | 0.66× |

#### By Environment

| Environment | N | MDR Rate | Fold Enrichment |
|-------------|---|----------|-----------------|
| Fish | 274 | 20.1% | 1.04× |
| Water | 177 | 19.8% | 1.03× |
| Hospital | 41 | 12.2% | 0.63× |

### 3.2 MDR Resistance Signature

**Antibiotics Most Associated with MDR** (based on mean resistance difference):

| Antibiotic | MDR Mean | Non-MDR Mean | Difference | MDR-Associated? |
|------------|----------|--------------|------------|-----------------|
| TE | 1.82 | 0.31 | +1.51 | ✅ Yes |
| DO | 1.75 | 0.28 | +1.47 | ✅ Yes |
| SXT | 1.68 | 0.22 | +1.46 | ✅ Yes |
| AM | 1.54 | 0.89 | +0.65 | ✅ Yes |
| C | 0.95 | 0.12 | +0.83 | ✅ Yes |

**MDR Signature Pattern**: MDR isolates are characterized by tetracycline (TE, DO), trimethoprim-sulfamethoxazole (SXT), and ampicillin (AM) resistance.

### 3.3 MDR Hotspot Summary

| Hotspot Type | Category | MDR Rate | Fold Enrichment |
|--------------|----------|----------|-----------------|
| **Cluster** | Cluster 3 (E. coli-Tetracycline) | 54.5% | **2.82×** |
| Region | Central Luzon | 21.4% | 1.11× |
| Environment | Fish | 20.1% | 1.04× |
| **Species** | *E. coli* | 30.0% | **1.55×** |

---

## 4. Cross-Phase Coherence

### 4.1 Phase Integration Summary

| Finding Type | Phase 3 (Clustering) | Phase 4 (Supervised) | Phase 5 (Regional) | Coherent? |
|--------------|---------------------|---------------------|-------------------|-----------|-
| MDR patterns | Cluster 3 highest MDR (54.5%) | Species F1=0.724 | BARMM/Central Luzon enriched | ✅ Yes |
| Species patterns | Clusters species-specific (V=0.765) | Classification works | Species varies by environment | ✅ Yes |
| Resistance breadth | Clusters vary in MAR | Feature importance: TE, DO | Regional variation moderate | ✅ Yes |

### 4.2 Model Agreement with Clustering

| Comparison | Agreement Level | Evidence |
|------------|-----------------|----------|
| High-importance antibiotics vs. cluster-defining | **High** | TE, DO, AM appear in both top feature lists |
| MDR discrimination vs. cluster MDR enrichment | **High** | Cluster 3 = MDR hotspot confirmed by chi-square |
| Species discrimination vs. cluster species purity | **High** | V=0.765 species-cluster association |

### 4.3 Coherence Statement

> "The integration of unsupervised clustering (Phase 3), supervised discrimination (Phase 4), and multivariate analysis (Phase 5) reveals **consistent patterns**: resistance-based clusters align strongly with species identity (V=0.765) and MDR status (V=0.559). The dichotomous *E. coli* phenotypes (Clusters 3 and 4) represent the most significant finding—biologically distinct populations within the same species occupying similar environments. These findings support the thesis that **environmental AMR in the Philippines is structured, species-driven, and amenable to systematic pattern recognition**, while acknowledging limitations including cross-sectional design and phenotypic-only data."

---

## 5. Final Documentation Checklist

### Phase 8 Audit Checklist

| Item | Status | Location | Notes |
|------|--------|----------|-------|
| ✅ Preprocessing decision log | Complete | docs/methods/preprocessing.md | Thresholds documented |
| ✅ Clustering parameter table | Complete | docs/methods/clustering.md | Ward/Euclidean justified |
| ✅ Supervised evaluation framing | Complete | docs/methods/supervised_models.md | Discrimination terminology |
| ✅ Terminology enforcement | Complete | All docs | "Association", not "causation" |
| ✅ Figure caption discipline | Complete | All results docs | Data source, method, scope |
| ✅ Limitations section | Complete | docs/limitations.md | Comprehensive documentation |
| ✅ Reproducibility statement | Complete | README.md | RANDOM_STATE=42 |

### Results Documentation Completeness

| Results Section | Documented | Template Location |
|-----------------|------------|-------------------|
| Preprocessing results | ✅ Yes | docs/results/phase2_clusters.md |
| Clustering results | ✅ Yes | docs/results/phase3_discrimination.md |
| Supervised results | ✅ Yes | docs/results/phase3_discrimination.md |
| Regional/Environmental | ✅ Yes | docs/results/phase4_environment.md |
| Integration/Synthesis | ✅ Yes | docs/results/phase5_synthesis.md |

---

## Key Findings Summary

### 1. Primary Findings

1. **E. coli Phenotypic Heterogeneity**: Two *E. coli* clusters with opposing resistance profiles (54.5% vs 0% MDR) suggest intraspecific diversity exceeds interspecific differences.

2. **Species-Driven Clustering**: Species identity (Cramér's V = 0.765) is the dominant factor in cluster membership, exceeding geographic (V = 0.321) and environmental (V = 0.204) factors.

3. **Tetracycline-MDR Signature**: MDR status is strongly associated with tetracycline resistance (TE, DO), possibly reflecting aquaculture antibiotic use.

4. **Cluster-MDR Alignment**: Unsupervised clusters significantly align with MDR status (V = 0.559), validating that clustering captures clinically relevant patterns.

### 2. Limitations

1. **No temporal inference**: Cross-sectional design; cannot establish temporal trends
2. **Dataset dependency**: Findings reflect this specific dataset; may not generalize
3. **No clinical decision support**: Results are for exploratory surveillance analysis only
4. **Phenotypic-only data**: Without genomic validation, resistance mechanisms remain unconfirmed

### 3. Implications

These findings provide baseline resistance patterns for Philippine environmental surveillance. The identified archetypes can serve as reference profiles for future monitoring, and the tetracycline-MDR association may inform antibiotic stewardship discussions in aquaculture settings.

---

## Figure Captions

### Figure 5.1: Archetype Summary Visualization

> "Summary of five resistance archetypes identified through hierarchical clustering (n=491 isolates, 6 species). Archetypes characterized by mean resistance profiles; MDR rates and environmental associations shown for post-hoc interpretation."

### Figure 5.2: MDR Enrichment Comparison

> "Multi-Drug Resistance enrichment by cluster, region, and environment. Fold enrichment calculated relative to overall MDR rate (19.1%). Cluster 3 (*E. coli*-Tetracycline) shows 2.86× enrichment."

### Figure 5.3: Integration Summary

> "Cross-phase integration summary showing alignment between clustering results (Phase 3), supervised discrimination (Phase 4), and regional analysis (Phase 5). Consistent patterns across phases strengthen confidence in identified resistance archetypes."

---

## Document Version

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2024 | Initial results template |
| 2.0 | 2025-12-18 | **Populated with actual pipeline results** |

---

*This document is part of Phase 8 — Documentation & Reporting for the AMR Thesis Project.*

**See also**: [integration.md](../methods/integration.md) | [limitations.md](../limitations.md)
