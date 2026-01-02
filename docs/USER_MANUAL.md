# AMR Thesis Project - User Manual

> **Version:** 2.1  
> **Last Updated:** December 28, 2025  
> **Author:** Queshue  
> **Data Source:** INOHAC AMR Project Two

---

## Table of Contents

1. [Quick Start](#1-quick-start)
2. [Prerequisites](#2-prerequisites)
3. [CLI Reference](#3-cli-reference)
4. [Pipeline Details](#4-pipeline-details)
5. [Input Data Format](#5-input-data-format)
6. [Output Files](#6-output-files)
7. [Dashboard Guide](#7-dashboard-guide)
8. [Configuration](#8-configuration)
9. [Troubleshooting](#9-troubleshooting)
10. [Advanced Usage](#10-advanced-usage)

---

## 1. Quick Start

```bash
# 1. Navigate to project directory
cd amr_thesis_project_code

# 2. Activate virtual environment
.\venv\Scripts\Activate.ps1     # Windows PowerShell
source venv/bin/activate         # Linux/Mac

# 3. Run full pipeline
python main.py --all

# 4. Launch dashboard
python main.py --app
```

**That's it!** The pipeline will process all data and launch the interactive dashboard.

---

## 2. Prerequisites

### 2.1 System Requirements

| Requirement | Minimum | Recommended |
| ----------- | ------- | ----------- |
| Python      | 3.9+    | 3.11+       |
| RAM         | 4 GB    | 8 GB        |
| Storage     | 500 MB  | 1 GB        |

### 2.2 Python Dependencies

Install all dependencies:

```bash
pip install -r requirements.txt
```

Key packages:

- `pandas`, `numpy` — Data processing
- `scipy`, `scikit-learn` — Clustering & ML
- `matplotlib`, `seaborn` — Visualization
- `streamlit` — Interactive dashboard
- `joblib` — Model serialization

---

## 3. CLI Reference

**All commands go through `main.py`** — never run scripts directly.

### 3.1 Available Commands

| Command                        | Description                                                  |
| ------------------------------ | ------------------------------------------------------------ |
| `python main.py --pipeline`    | Core data flow: Ingestion → Cleaning → Encoding → Clustering |
| `python main.py --validate`    | Run validation scripts (cluster quality, co-resistance)      |
| `python main.py --sensitivity` | Run sensitivity analysis (encoding, thresholds)              |
| `python main.py --analyze`     | Run analysis modules (regional, synthesis, supervised)       |
| `python main.py --viz`         | Regenerate all visualizations                                |
| `python main.py --app`         | Launch Streamlit dashboard                                   |
| `python main.py --all`         | Run everything in sequence                                   |
| `python main.py --help`        | Show help                                                    |

### 3.2 Optional Flags

| Flag    | Description            | Example                           |
| ------- | ---------------------- | --------------------------------- |
| `--k N` | Override cluster count | `python main.py --pipeline --k 4` |

### 3.3 Command Examples

```bash
# Full pipeline with default settings (auto-detect k)
python main.py --all

# Re-run with 6 clusters
python main.py --pipeline --k 6 --validate --analyze --viz

# Only regenerate figures
python main.py --viz

# Launch dashboard only (data must exist)
python main.py --app
```

---

## 4. Pipeline Details

### 4.1 Pipeline Stages

```
┌─────────────────────────────────────────────────────────────────┐
│                         --pipeline                               │
├─────────────────────────────────────────────────────────────────┤
│  Step 1: Data Ingestion                                          │
│          - Reads all *.csv from project root                     │
│          - Extracts metadata from isolate codes                  │
│          - Output: unified_raw_dataset.csv                       │
├─────────────────────────────────────────────────────────────────┤
│  Step 2: Data Cleaning                                           │
│          - Applies coverage thresholds (70% antibiotic, 30% max) │
│          - Handles missing values                                │
│          - Output: cleaned_dataset.csv, cleaning_report.txt      │
├─────────────────────────────────────────────────────────────────┤
│  Step 3: Resistance Encoding                                     │
│          - S=0, I=1, R=2                                         │
│          - Output: encoded_dataset.csv                           │
├─────────────────────────────────────────────────────────────────┤
│  Step 4: Feature Engineering                                     │
│          - Computes MAR index, MDR status                        │
│          - Separates features from metadata                      │
│          - Output: analysis_ready_dataset.csv                    │
├─────────────────────────────────────────────────────────────────┤
│  Step 5: Cluster Validation                                      │
│          - Elbow + Silhouette analysis                           │
│          - Selects optimal k (constrained 3-6)                   │
│          - Output: cluster_validation.png                        │
├─────────────────────────────────────────────────────────────────┤
│  Step 6: Hierarchical Clustering                                 │
│          - Ward linkage, Euclidean distance                      │
│          - Output: clustered_dataset.csv, linkage_matrix.pkl     │
└─────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────┐
│                         --validate                               │
├─────────────────────────────────────────────────────────────────┤
│  - validate_clustering.py: Cluster quality metrics               │
│  - coresistance_analysis.py: Co-resistance network analysis      │
│  - antibiotic_clustering.py: Antibiotic clustering by φ matrix   │
└─────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────┐
│                        --sensitivity                             │
├─────────────────────────────────────────────────────────────────┤
│  - encoding_sensitivity.py: Tests impact of S/I/R values         │
│  - threshold_sensitivity.py: Tests cleaning threshold impact     │
└─────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────┐
│                         --analyze                                │
├─────────────────────────────────────────────────────────────────┤
│  - Regional & Environmental distribution analysis                │
│  - Integration & Synthesis (archetypes, patterns)                │
│  - Species Classifier (Random Forest)                            │
└─────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────┐
│                           --viz                                  │
├─────────────────────────────────────────────────────────────────┤
│  - Resistance heatmaps                                           │
│  - Dendrograms                                                   │
│  - PCA biplots                                                   │
│  - Distribution charts                                           │
│  - All saved to data/processed/figures/                          │
└─────────────────────────────────────────────────────────────────┘
```

### 4.2 Execution Order

When running `--all`, stages execute in this order:

1. `--pipeline` (data processing + clustering)
2. `--validate` (validation scripts)
3. `--sensitivity` (robustness checks)
4. `--analyze` (analysis modules)
5. `--viz` (visualization regeneration)

---

## 5. Input Data Format

### 5.1 Required CSV Structure

Place raw CSV files in the `data/raw/` directory.

| Column           | Description               | Example                 |
| ---------------- | ------------------------- | ----------------------- |
| `CODE`           | Unique isolate identifier | `EC_OADWR1C3`           |
| Antibiotic codes | S/I/R values              | `AM`, `AMC`, `TE`, etc. |

### 5.2 Isolate Code Convention

Format: `{SPECIES}_{SITE}{SOURCE}{REPLICATE}C{COLONY}`

| Component     | Values                                                |
| ------------- | ----------------------------------------------------- |
| Species       | EC (E. coli), KP (Klebsiella), PS (Pseudomonas), etc. |
| National Site | O (Ormoc), P (Pampanga), M (Marawi)                   |
| Local Site    | A (Alegria), L (Larrazabal), G (Gabriel), etc.        |
| Source        | DW (Drinking Water), RW (River), FB (Fish), etc.      |
| Replicate     | 1, 2, 3                                               |
| Colony        | C1, C2, C3, etc.                                      |

Example: `EC_OADWR1C3` = _E. coli_, Ormoc-Alegria, Drinking Water, Replicate 1, Colony 3

### 5.3 Resistance Values

| Value | Meaning      | Encoded As |
| ----- | ------------ | ---------- |
| S     | Susceptible  | 0          |
| I     | Intermediate | 1          |
| R     | Resistant    | 2          |

---

## 6. Output Files

All outputs are saved to `data/processed/`:

### 6.1 Datasets

| File                         | Description              |
| ---------------------------- | ------------------------ |
| `unified_raw_dataset.csv`    | Consolidated raw data    |
| `cleaned_dataset.csv`        | After quality filtering  |
| `encoded_dataset.csv`        | Numeric encoding applied |
| `analysis_ready_dataset.csv` | With computed features   |
| `clustered_dataset.csv`      | With CLUSTER column      |

### 6.2 Figures (`data/processed/figures/`)

| File                               | Description                        |
| ---------------------------------- | ---------------------------------- |
| `cluster_validation.png`           | Elbow + Silhouette plot            |
| `dendrogram_highres.png`           | Hierarchical clustering tree       |
| `dendrogram_linked_heatmap.png`    | Dendrogram with resistance heatmap |
| `cluster_profiles.png`             | Mean resistance by cluster         |
| `pca_biplot.png`                   | PCA with antibiotic loadings       |
| `pca_by_cluster.png`               | PCA colored by cluster             |
| `coresistance_network.png`         | Co-resistance network graph        |
| `antibiotic_clusters.csv`          | Antibiotic cluster assignments     |
| `antibiotic_dendrogram.png`        | Antibiotic hierarchical clustering |
| `antibiotic_clustered_heatmap.png` | Phi coefficient matrix (clustered) |
| `mdr_distribution.png`             | MDR by cluster                     |
| `mar_distribution.png`             | MAR index by cluster               |

### 6.3 Models (`models/`)

| File                        | Description                      |
| --------------------------- | -------------------------------- |
| `species_classifier.joblib` | Trained Random Forest classifier |

### 6.4 Artifacts (`data/processed/clustering_artifacts/`)

| File                   | Description          |
| ---------------------- | -------------------- |
| `linkage_matrix.pkl`   | Scipy linkage matrix |
| `clustering_info.pkl`  | Cluster metadata     |
| `dendrogram_order.pkl` | Leaf ordering        |

---

## 7. Dashboard Guide

### 7.1 Launching

```bash
python main.py --app
```

The dashboard opens at `http://localhost:8501`.

### 7.2 Available Analyses

| Tab                       | Description                                     |
| ------------------------- | ----------------------------------------------- |
| **Overview**              | Dataset summary, isolate counts                 |
| **Enhanced Summary**      | Cluster distribution table                      |
| **Resistance Heatmap**    | Visual resistance profiles                      |
| **Cluster Analysis**      | Cluster statistics, MDR by cluster              |
| **Cluster Validation**    | Dynamic k selection justification with metrics  |
| **Species Classifier**    | Confusion matrix, feature importance            |
| **Co-Resistance Network** | Significant antibiotic associations             |
| **Antibiotic Clustering** | Antibiotics clustered by co-resistance patterns |
| **PCA Analysis**          | Dimensionality reduction plots                  |
| **Regional Distribution** | Cluster × Region cross-tabulation               |
| **Model Evaluation**      | Performance metrics summary                     |
| **Methodology**           | Technical methods documentation                 |
| **Limitations**           | Study limitations                               |

### 7.3 Dashboard Screenshots

#### Overview Tab

![Dashboard Overview](docs/images/dashboard_overview.png)
_Dataset summary with species distribution and MDR statistics_

#### Cluster Validation

![Cluster Validation](docs/images/dashboard_cluster_validation.png)
_Dynamic k selection showing elbow + silhouette analysis with k=4 selected_

#### Resistance Heatmap

![Resistance Heatmap](docs/images/dashboard_heatmap.png)
_Clustered resistance heatmap showing antibiotic susceptibility patterns_

#### PCA Analysis

![PCA Analysis](docs/images/dashboard_pca.png)
_Principal Component Analysis colored by cluster membership_

#### Model Evaluation

![Model Evaluation](docs/images/dashboard_model_eval.png)
_Species classifier confusion matrix and performance metrics_

#### Antibiotic Clustering

![Antibiotic Clustering](docs/images/dashboard_antibiotic_clustering.png)
_Antibiotic co-resistance patterns and phi coefficient matrix_

### 7.4 Refreshing Data

The dashboard reads from `data/processed/`. To update:

1. Stop the dashboard (Ctrl+C)
2. Run `python main.py --pipeline --validate --analyze --viz`
3. Restart `python main.py --app`

---

## 8. Configuration

### 8.1 Central Config File

All parameters are in `src/config.py`:

| Parameter                             | Default | Description                    |
| ------------------------------------- | ------- | ------------------------------ |
| `RANDOM_STATE`                        | 42      | Reproducibility seed           |
| `MIN_ANTIBIOTIC_COVERAGE`             | 70.0    | Min % coverage per antibiotic  |
| `MAX_ISOLATE_MISSING`                 | 30.0    | Max % missing per isolate      |
| `MDR_MIN_CLASSES`                     | 3       | Classes for MDR classification |
| `CLUSTERING_CONFIG['min_k']`          | 3       | Minimum k for auto-selection   |
| `CLUSTERING_CONFIG['max_k']`          | 6       | Maximum k for auto-selection   |
| `CLUSTERING_CONFIG['linkage_method']` | ward    | Clustering method              |

> **Dynamic k Selection:** The optimal number of clusters (k) is **automatically determined** by the pipeline using combined elbow + silhouette analysis. The algorithm evaluates all k values in the [min_k, max_k] range and selects the optimal value based on statistical criteria.

### 8.2 Overriding Cluster Count

Use the `--k` flag only if you need to override automatic selection:

```bash
python main.py --pipeline --k 4  # Force k=4 (not recommended unless justified)
```

---

## 9. Troubleshooting

### 9.1 Common Issues

| Issue                          | Solution                                     |
| ------------------------------ | -------------------------------------------- |
| `ModuleNotFoundError`          | Activate venv: `.\venv\Scripts\Activate.ps1` |
| `streamlit: command not found` | Install: `pip install streamlit` or use venv |
| `No data loaded`               | Ensure CSV files are in project root         |
| `Matplotlib thread error`      | Set `MPLBACKEND=Agg` before running          |
| Dashboard shows old data       | Re-run pipeline, restart dashboard           |

### 9.2 Matplotlib Backend Issue

If you see "main thread is not in main loop":

```bash
$env:MPLBACKEND='Agg'; python main.py --all
```

### 9.3 Check Pipeline Output

Verify data exists:

```bash
ls data/processed/*.csv
ls data/processed/figures/*.png
```

---

## 10. Advanced Usage

### 10.1 Importing Modules Directly

```python
# In Python script or Jupyter
import sys
sys.path.insert(0, 'src')

from preprocessing.data_ingestion import create_unified_dataset
from clustering.hierarchical_clustering import run_clustering_pipeline
from config import ANTIBIOTIC_CLASSES, RANDOM_STATE
```

### 10.2 Custom k Analysis

Run with different k values and compare:

```bash
python main.py --pipeline --k 3 --viz
# Inspect figures, then:
python main.py --pipeline --k 5 --viz
```

### 10.3 Adding New Data

1. Add new CSV files to project root
2. Ensure they follow the CODE naming convention
3. Run `python main.py --all`

### 10.4 Re-training Models Only

```bash
python main.py --analyze  # Runs species classifier
```

---

## File Structure Reference

```
amr_thesis_project_code/
├── main.py                    # 🔴 CENTRAL ENTRY POINT
├── requirements.txt           # Python dependencies
├── README.md                  # Project overview
├── USER_MANUAL.md             # This file
├── REFACTOR_LOG.md            # Change history
│
├── src/                       # Source code modules
│   ├── config.py              # 🔧 Central configuration (single source of truth)
│   ├── preprocessing/         # Data processing (ingestion, cleaning, encoding)
│   ├── clustering/            # Hierarchical clustering
│   ├── supervised/            # ML models (species classifier)
│   ├── analysis/              # Regional/environmental analysis
│   ├── visualization/         # Plot generation
│   └── utils/                 # Shared utilities (console styling)
│
├── scripts/                   # Standalone analysis scripts
│   ├── validate_clustering.py    # Silhouette + elbow analysis
│   ├── coresistance_analysis.py  # Co-resistance network
│   ├── antibiotic_clustering.py  # Antibiotic grouping
│   ├── encoding_sensitivity.py   # Encoding scheme comparison
│   └── threshold_sensitivity.py  # Data threshold analysis
│
├── app/                       # Streamlit dashboard
│   └── streamlit_app.py
│
├── models/                    # Trained ML models
│   └── species_classifier.joblib
│
├── data/
│   └── processed/             # Pipeline outputs
│       ├── figures/           # All visualizations (39+ files)
│       └── clustering_artifacts/  # Reproducibility artifacts
│
└── docs/                      # Reference documentation
    ├── TECHNICAL_REFERENCE.md # Condensed technical summary
    ├── limitations.md         # Study scope & limitations
    └── archive/               # Historical reference notes
```

---

## Support

For questions or issues:

1. Check this manual first
2. Review `REFACTOR_LOG.md` for recent changes
3. Check `docs/` reference files for technical details

---

**Remember:** Always use `python main.py [flags]` — never run scripts directly!
