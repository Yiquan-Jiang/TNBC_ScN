# Chemotherapy-Induced Senescent Macrophages Drive Stromal Fibrosis and Immune Exclusion via the IL-1β-matCAF Axis in TNBC

## Chemotherapy-Induced Cellular Senescence and Immune Microenvironment Analysis

[![R Version](https://img.shields.io/badge/R-≥4.0-blue.svg)](https://www.r-project.org/)
[![Python Version](https://img.shields.io/badge/Python-≥3.8-green.svg)](https://www.python.org/)
[![Seurat](https://img.shields.io/badge/Seurat-v5.0+-orange.svg)](https://satijalab.org/seurat/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

---

## 📖 Overview

This repository contains comprehensive analysis code for investigating **chemotherapy-induced cellular senescence** and its effects on the **tumor immune microenvironment** in breast cancer. The study integrates multiple data modalities across different model systems:

### Key Scientific Questions

1. **Does chemotherapy induce cellular senescence in the tumor microenvironment?**
2. **How do senescent cells reshape immune cell composition and function?**
3. **What are the spatial relationships between IL1R1+ CAFs, CXCL13+ T cells, and CDKN1A+ macrophages?**
4. **Can senolytic therapy (ABT-263) reverse chemotherapy-induced immunosuppression?**

### Key Findings

- 🔬 **Myeloid Cell Senescence**: Chemotherapy induces p21/CDKN1A-mediated senescence in tumor-associated macrophages
- 🔥 **SASP Activation**: Senescent macrophages secrete IL1B and other SASP factors
- 🧬 **CAF Remodeling**: Enhanced ECM production and IL1R1+ CAF expansion post-chemotherapy
- 🛡️ **T Cell Dysfunction**: Altered CXCL13+ T cell proportions and increased exhaustion markers
- 🎯 **Spatial Colocalization**: IL1R1+ CAFs colocalize with CXCL13+ exhausted T cells in human tumors

---

## 📁 Repository Structure

This repository is organized into **multiple branches**, each containing analysis code for different datasets:

```
Main Branch (Current)
├── README.md                              # This file
├── Senescence_NAC_vs_naive_scRNA_GSE169246.R  # Human TNBC scRNA-seq (GSE169246)
│
├── mouse_data_github_release/             # Branch: mouse-scrnaseq
│   ├── scripts/                           # R analysis scripts (01-07)
│   ├── renv.lock                          # R environment lockfile
│   └── README.md                          # Mouse data documentation
│
├── Patient_data_R_analysis/               # Branch: patient-scrnaseq
│   ├── scripts/FINAL_ANALYSIS_CODE.R      # Complete patient analysis
│   ├── renv.lock                          # R environment lockfile
│   └── README.md                          # Patient data documentation
│
└── publication_code_brca_spatial_transcriptomes/  # Branch: spatial-transcriptomics
    ├── 01_spatial_colocalization_base.py  # Base spatial analysis module
    ├── 02_il1r1_cxcl13_extended_analysis.py
    ├── 03_comprehensive_analysis_with_tumor_and_cdkn1a.py
    ├── 04_cellchat_analysis.R             # Cell-cell communication
    └── README.md                          # Spatial data documentation
```

---

## 🔬 Datasets & Analysis Branches

### 1. 🐭 Mouse scRNA-seq Data (`mouse_data_github_release/`)

**Description**: Single-cell RNA sequencing analysis of mouse tumor models treated with chemotherapy and senolytic combinations.

**Treatment Groups**:
| Group | Description |
|-------|-------------|
| Untreated | Control group |
| Chemo | Chemotherapy treatment |
| Chemo_Abt263 | Chemotherapy + ABT-263 (senolytic) |
| Chemo_Abt263_IL1B | Chemotherapy + ABT-263 + IL1B treatment |

**Analysis Pipeline**:
1. `01_data_processing.R` - Data loading, QC, and CAF subtyping
2. `02_umap_visualization.R` - UMAP plots using SCP package
3. `03_caf_analysis.R` - CAF subtype classification (iCAF, myCAF, apCAF, matCAF)
4. `04_differential_enrichment.R` - DE analysis and pathway enrichment
5. `05_comprehensive_visualization.R` - Publication figures
6. `06_merged_data_analysis.R` - Multi-dataset integration
7. `07_nebulosa_expression.R` - Expression density plots

**CAF Subtype Markers**:
| Subtype | Key Markers |
|---------|-------------|
| iCAF (Inflammatory) | Il6, Il1b, Cxcl1, Cxcl2, Ccl2, Tnf |
| myCAF (Myofibroblastic) | Acta2, Tagln, Myh11, Postn, Lox |
| apCAF (Antigen-presenting) | H2-Aa, H2-Ab1, Cd74, Ciita |
| matCAF (Matrix) | Col1a1, Col1a2, Col3a1, Fn1, Sparc |

---

### 2. 👥 Patient scRNA-seq Data (`Patient_data_R_analysis/`)

**Description**: Single-cell RNA sequencing analysis of human TNBC samples comparing treatment-naive (CTRL) vs neoadjuvant chemotherapy (NAC) patients.

**Data Source**: 8 samples (4 CTRL, 4 NAC) from clinical TNBC specimens

**Key Analyses**:
- Cellular senescence scoring (CDKN1A, CDKN2A)
- SASP factor expression analysis
- T cell exhaustion and cytotoxicity profiling
- Fibroblast ECM remodeling assessment
- Senescence-IL1B correlation analysis

**Output Panels**:
| Panel | Description | Size |
|-------|-------------|------|
| A | UMAP by Cell Type | 10×8 in |
| B | UMAP by Treatment Group | 10×8 in |
| C | CXCL13+ T Cell Proportion | 4×6 in |
| D-E | CD8+ Cytotoxicity/Exhaustion Scores | 5×6 in |
| F-H | Senescence Marker DotPlots | 6×10 in |
| I | Senescence-IL1B Correlation | 8×8 in |
| J-L | ECM and T Cell Marker Expression | varies |

---

### 3. 🗺️ Spatial Transcriptomics (`publication_code_brca_spatial_transcriptomes/`)

**Description**: Spatial analysis of breast cancer tissue sections to investigate cell-cell colocalization and spatial relationships.

**Analysis Focus**:
1. **IL1R1+ CAF vs CXCL13+ T cell colocalization**
2. **CDKN1A+ Macrophage vs CAF spatial relationships**
3. **Tumor microenvironment cell-cell interactions**

**Code Structure**:
```
SpatialColocalizationAnalyzer (base)
    └── IL1R1_CXCL13_ColocalizationAnalyzer (extended)
        └── ComprehensiveSpatialAnalyzer (full analysis)
```

**Cell Type Markers**:
| Cell Type | Marker Genes |
|-----------|--------------|
| CAF | COL1A1, COL1A2, COL3A1, DCN, LUM, FAP, PDGFRA |
| T cells | CD3D, CD3E, CD8A, TRAC, LCK, ZAP70 |
| CXCL13+ T | CXCL13, PDCD1, CTLA4, TIGIT, LAG3, TOX |
| Macrophages | CD68, CD163, CD14, MARCO, CSF1R |

---

### 4. 📊 Public Dataset Analysis (`Senescence_NAC_vs_naive_scRNA_GSE169246.R`)

**Description**: Validation analysis using publicly available TNBC scRNA-seq data (GSE169246) comparing treatment-naive vs chemotherapy-treated samples.

**Key Features**:
- RPCA-based data integration
- scATOMIC cell type annotation
- Senescence cell identification (CDKN1A+/CDKN2A+ cells in G1 phase)
- pySCENIC transcription factor analysis
- MDSC annotation and proportion analysis

---

## 🛠️ Installation

### R Environment Setup

```r
# Install renv for reproducible environment
install.packages("renv")

# Navigate to desired analysis directory and restore
setwd("mouse_data_github_release")  # or Patient_data_R_analysis
renv::restore()
```

### Core R Packages

```r
# BiocManager installation
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Core packages
install.packages(c("Seurat", "dplyr", "ggplot2", "tidyr", 
                   "patchwork", "ggpubr", "ggsci", "harmony"))

# SCP for high-quality visualization
remotes::install_github("zhanghao-njmu/SCP")

# Bioconductor packages
BiocManager::install(c("ComplexHeatmap", "clusterProfiler", 
                       "enrichplot", "Nebulosa"))
```

### Python Environment Setup

```bash
# Create conda environment
conda create -n spatial_analysis python=3.9
conda activate spatial_analysis

# Install dependencies
pip install scanpy>=1.9.0 squidpy>=1.2.0 pandas numpy 
pip install matplotlib seaborn scikit-learn scipy anndata
```

---

## 🚀 Quick Start

### Mouse Data Analysis

```r
# Load and process mouse scRNA-seq data
setwd("mouse_data_github_release/scripts")
source("01_data_processing.R")
source("02_umap_visualization.R")
source("03_caf_analysis.R")
```

### Patient Data Analysis

```r
# Run complete patient analysis
setwd("Patient_data_R_analysis/scripts")
source("FINAL_ANALYSIS_CODE.R")
```

### Spatial Transcriptomics

```bash
cd publication_code_brca_spatial_transcriptomes

# Run comprehensive Python analysis
python 03_comprehensive_analysis_with_tumor_and_cdkn1a.py

# Run CellChat analysis (R)
Rscript 04_cellchat_analysis.R
```

---

## 📊 Figure Specifications

All figures are generated with publication-quality settings:

| Parameter | Value |
|-----------|-------|
| Resolution | 300 DPI |
| Format | PDF (vector, editable) |
| Font Type | Type 42 (editable text) |
| Font Family | Arial / Helvetica |
| Color Scheme | Colorblind-friendly |

### R Settings

```r
options(pdf.fonttype = 42, ps.fonttype = 42)
```

### Python Settings

```python
plt.rcParams.update({
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'font.family': 'Arial',
    'figure.dpi': 300
})
```

---

## 📈 Statistical Methods

### Single-Cell Analysis
- **Normalization**: SCTransform / LogNormalize
- **Batch Correction**: Harmony / RPCA integration
- **Differential Expression**: Wilcoxon rank-sum test
- **Cell Type Scoring**: AddModuleScore

### Spatial Analysis
- **Colocalization**: Pearson/Spearman correlation
- **Spatial Proximity**: Permutation test (n=1000)
- **Neighbor Counting**: Mann-Whitney U test

### Multiple Testing
- **Adjustment Method**: Benjamini-Hochberg FDR
- **Significance Threshold**: p.adj < 0.05

---

## 📦 Data Availability

| Dataset | Accession | Description |
|---------|-----------|-------------|
| Human TNBC scRNA-seq | GSE169246 | NAC vs treatment-naive TNBC |
| Spatial Transcriptomics | VISDS000449-451 | 10X Visium breast cancer |
| Mouse scRNA-seq | [Contact authors] | ScN treatment study |

---

## 🔧 Troubleshooting

### Common Issues

1. **Memory Error**
   ```r
   options(future.globals.maxSize = 8000 * 1024^2)
   ```

2. **SCP Installation Failed**
   ```r
   BiocManager::install(c("SingleCellExperiment", "scater"))
   remotes::install_github("zhanghao-njmu/SCP", dependencies = TRUE)
   ```

3. **DotPlot Color Mapping Issues**
   - Use `scale = FALSE` parameter
   - Create manual grouping columns instead of `split.by`

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 📧 Contact

- **Author**: Yi-Quan Jiang
- **Email**: jiangyq@sysucc.org.cn
- **Institution**: Sun Yat-Sen University Cancer Center

For questions about the analysis code, please open an issue in this repository.

---

## Acknowledgments

- [Seurat](https://satijalab.org/seurat/) team for the scRNA-seq analysis framework
- [SCP](https://github.com/zhanghao-njmu/SCP) for high-quality visualization tools
- [Harmony](https://github.com/immunogenomics/harmony) for batch correction
- [Squidpy](https://squidpy.readthedocs.io/) for spatial transcriptomics analysis
- [CellChat](https://github.com/sqjin/CellChat) for cell-cell communication analysis

---

**Last Updated**: December 2025

