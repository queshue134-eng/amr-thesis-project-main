# Deployment Documentation

## Phase 7: Interactive Dashboard

This document provides comprehensive documentation of the deployment architecture, dashboard components, and usage guidelines for the AMR thesis project interactive dashboard.

---

## Table of Contents

1. [Dashboard Architecture](#1-dashboard-architecture)
2. [System Requirements](#2-system-requirements)
3. [Installation and Setup](#3-installation-and-setup)
4. [Dashboard Components](#4-dashboard-components)
5. [Usage Guidelines](#5-usage-guidelines)
6. [Reproducibility Statement](#6-reproducibility-statement)

---

## 1. Dashboard Architecture

### 1.1 Technology Stack

| Component | Technology | Version | Purpose |
|-----------|------------|---------|---------|
| **Framework** | Streamlit | ≥1.24.0 | Interactive web dashboard |
| **Data Processing** | pandas | ≥1.5.0 | Data manipulation |
| **Visualization** | matplotlib, seaborn | ≥3.6.0, ≥0.12.0 | Static plots |
| **Machine Learning** | scikit-learn | ≥1.1.0 | ML algorithms |
| **Scientific Computing** | scipy | ≥1.9.0 | Clustering, statistics |
| **Model Persistence** | joblib | ≥1.2.0 | Save/load models |

### 1.2 Architecture Diagram

```
┌─────────────────────────────────────────────────────────────────┐
│                    DASHBOARD ARCHITECTURE                        │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  USER INTERFACE (Streamlit)                                      │
│  ├── Sidebar: Navigation, file upload, parameters                │
│  └── Main area: Visualizations, tables, analysis results         │
│                                                                 │
│  DATA LAYER                                                      │
│  ├── File upload (CSV)                                           │
│  ├── Default datasets (if available)                             │
│  └── Session state caching                                       │
│                                                                 │
│  ANALYSIS MODULES                                                │
│  ├── Preprocessing (data loading, encoding)                      │
│  ├── Clustering (hierarchical clustering)                        │
│  ├── Visualization (heatmaps, dendrograms)                       │
│  ├── Supervised (model evaluation)                               │
│  └── PCA (dimensionality reduction)                              │
│                                                                 │
│  OUTPUT                                                          │
│  ├── Interactive plots                                           │
│  ├── Summary tables                                              │
│  └── Downloadable results                                        │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

---

## 2. System Requirements

### 2.1 Hardware Requirements

| Requirement | Minimum | Recommended |
|-------------|---------|-------------|
| **RAM** | 4 GB | 8 GB |
| **Disk Space** | 500 MB | 1 GB |
| **CPU** | Dual-core | Quad-core |

### 2.2 Software Requirements

| Software | Version | Notes |
|----------|---------|-------|
| **Python** | 3.8+ | Required |
| **pip** | Latest | Package manager |
| **Web Browser** | Modern | Chrome, Firefox, Safari, Edge |

### 2.3 Operating System Compatibility

| OS | Supported |
|----|-----------|
| Windows 10/11 | ✅ |
| macOS 10.14+ | ✅ |
| Linux (Ubuntu 18.04+) | ✅ |

---

## 3. Installation and Setup

### 3.1 Installation Steps

```bash
# 1. Clone or download repository
cd amr-thesis-project

# 2. Create virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # Linux/macOS
# OR
venv\Scripts\activate     # Windows

# 3. Install dependencies
pip install -r requirements.txt

# 4. Verify installation
python -c "import streamlit; print('Streamlit version:', streamlit.__version__)"
```

### 3.2 Running the Dashboard

```bash
# Start dashboard
streamlit run app/streamlit_app.py

# Dashboard will open at: http://localhost:8501
```

### 3.3 Alternative Port

```bash
# If port 8501 is in use
streamlit run app/streamlit_app.py --server.port 8502
```

---

## 4. Dashboard Components

### 4.1 Available Pages/Views

| Page | Description | Key Features |
|------|-------------|--------------|
| **Overview** | Dataset summary and preview | Row/column counts, data types |
| **Resistance Heatmap** | Visual resistance profiles | Color-coded S/I/R patterns |
| **Cluster Analysis** | Cluster statistics and distribution | Cluster sizes, MDR rates |
| **PCA Analysis** | Dimensionality reduction | 2D scatter, loadings |
| **Regional Distribution** | Geographic analysis | Region-cluster cross-tabs |
| **Model Evaluation** | Supervised learning results | Metrics, confusion matrices |
| **Integration & Synthesis** | Combined findings | Archetype summaries |

### 4.2 Interactive Features

| Feature | Description |
|---------|-------------|
| **File Upload** | Upload custom CSV datasets |
| **Parameter Selection** | Adjust clustering parameters |
| **Color Coding** | Select variables for plot coloring |
| **Filtering** | Filter data by species, region, etc. |
| **Export** | Download visualizations and tables |

---

## 5. Usage Guidelines

### 5.1 Data Requirements

| Requirement | Description |
|-------------|-------------|
| **Format** | CSV file |
| **Required columns** | CODE, resistance columns (S/I/R values) |
| **Recommended** | Metadata columns (REGION, SPECIES, etc.) |

### 5.2 Workflow

```
1. Upload Data
   └── Use sidebar file uploader or default dataset

2. View Overview
   └── Check data summary and quality

3. Explore Patterns
   ├── Resistance heatmap
   └── Cluster analysis

4. Analyze Dimensions
   └── PCA visualization

5. Examine Distributions
   └── Regional/environmental analysis

6. Review Model Results
   └── Supervised learning metrics
```

### 5.3 Interpretation Caveats

> **IMPORTANT**: Dashboard results are for **exploratory analysis only**. They should not be used for clinical decision support. All findings reflect patterns within the uploaded dataset and are not predictive of external data.

---

## 6. Reproducibility Statement

### 6.1 Software Versions

| Package | Version Used | Purpose |
|---------|--------------|---------|
| Python | 3.8+ | Runtime |
| pandas | ≥1.5.0 | Data manipulation |
| numpy | ≥1.23.0 | Numerical operations |
| scipy | ≥1.9.0 | Clustering |
| scikit-learn | ≥1.1.0 | Machine learning |
| matplotlib | ≥3.6.0 | Visualization |
| seaborn | ≥0.12.0 | Statistical plots |
| streamlit | ≥1.24.0 | Dashboard |
| joblib | ≥1.2.0 | Model persistence |

### 6.2 Random Seeds

| Component | Random State | Purpose |
|-----------|--------------|---------|
| Train-test split | 42 | Reproducible data partitioning |
| Random Forest | 42 | Reproducible model training |
| Logistic Regression | 42 | Reproducible model training |
| SVM | 42 | Reproducible model training |
| Decision Tree | 42 | Reproducible model training |

### 6.3 Environment Isolation

**Recommended**: Use virtual environments to isolate dependencies.

```bash
# Create isolated environment
python -m venv amr_env

# Activate and install
source amr_env/bin/activate  # Linux/macOS
pip install -r requirements.txt

# Document environment
pip freeze > environment_snapshot.txt
```

### 6.4 Reproducibility Checklist

| Item | Status | Notes |
|------|--------|-------|
| ✅ Fixed random seeds | 42 | Consistent across all stochastic operations |
| ✅ Versioned dependencies | requirements.txt | Minimum versions specified |
| ✅ Documented parameters | docs/methods/ | All parameters documented |
| ✅ Data provenance | cleaning_report.txt | Data transformations logged |
| ✅ Environment isolation | Virtual environment | Recommended practice |

### 6.5 Full Reproducibility Statement

> "All analyses were performed using Python 3.8+ with scikit-learn (≥1.1.0), scipy (≥1.9.0), and pandas (≥1.5.0). Random states were fixed at 42 for all stochastic operations to ensure reproducibility. The complete software environment is documented in requirements.txt. Data preprocessing decisions are logged in cleaning_report.txt, and all methodological parameters are documented in the docs/methods/ directory."

---

## Document Version

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2024 | Initial deployment documentation |

---

*This document is part of Phase 8 — Documentation & Reporting for the AMR Thesis Project.*

**See also**: [RUNNING_THE_SYSTEM.md](../RUNNING_THE_SYSTEM.md) | [ARCHITECTURE.md](../ARCHITECTURE.md)
