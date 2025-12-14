# Single-Cell RNA-seq Analysis Pipeline

## Chemotherapy-Induced Senescence in Triple-Negative Breast Cancer (TNBC)

[![R Version](https://img.shields.io/badge/R-≥4.0-blue.svg)](https://www.r-project.org/)
[![Seurat](https://img.shields.io/badge/Seurat-v5.0+-green.svg)](https://satijalab.org/seurat/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

---

## Overview

This repository contains the R analysis pipeline for investigating the effects of neoadjuvant chemotherapy (NAC) on the tumor immune microenvironment in triple-negative breast cancer (TNBC) using single-cell RNA sequencing (scRNA-seq).

### Key Findings

1. **Myeloid Cell Senescence**: NAC induces cellular senescence in tumor-associated macrophages via the p21/CDKN1A pathway
2. **SASP Activation**: Senescent macrophages secrete IL1B and other SASP factors
3. **Fibroblast Remodeling**: Enhanced ECM production and CAF activation post-chemotherapy
4. **T Cell Exhaustion**: Increased exhaustion markers and altered CXCL13+ T cell proportions

---

## Repository Structure

```
R_analysis/
├── README.md                    # This file
├── renv.lock                    # R environment lockfile (reproducibility)
├── environment.yaml             # Conda environment (alternative)
├── scripts/
│   └── FINAL_ANALYSIS_CODE.R    # Main analysis pipeline
├── data/
│   └── .gitkeep                 # Placeholder (data not included)
└── output/
    ├── figures/                 # Publication-quality figures
    └── statistics/              # Statistical results (CSV)
```

---

## Requirements

### System Requirements

- **OS**: macOS, Linux, or Windows (WSL2 recommended)
- **RAM**: ≥32 GB (64 GB recommended for large datasets)
- **Disk**: ≥50 GB free space

### Software Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| R | ≥4.0 | Base environment |
| Seurat | ≥5.0 | scRNA-seq analysis |
| SCP | latest | High-quality visualization |
| harmony | ≥1.0 | Batch correction |
| dplyr | ≥1.0 | Data manipulation |
| ggplot2 | ≥3.4 | Plotting |
| patchwork | ≥1.1 | Figure composition |
| ggpubr | ≥0.6 | Statistical annotations |
| ggsci | ≥3.0 | Journal color palettes |

---

## Installation

### Option 1: Using renv (Recommended)

```r
# Install renv if not already installed
install.packages("renv")

# Restore the project environment
renv::restore()
```

### Option 2: Manual Installation

```r
# Install BiocManager
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install CRAN packages
install.packages(c("Seurat", "dplyr", "ggplot2", "tidyr", 
                   "patchwork", "ggpubr", "ggsci", "harmony"))

# Install SCP from GitHub
if (!require("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("zhanghao-njmu/SCP")
```

### Option 3: Using Conda

```bash
conda env create -f environment.yaml
conda activate scrnaseq-tnbc
```

---

## Data Preparation

### Input Data Format

The pipeline expects 10X Genomics scRNA-seq data in the following structure:

```
data/
├── sample1/
│   └── output/
│       └── filter_matrix/
│           ├── barcodes.tsv.gz
│           ├── features.tsv.gz
│           └── matrix.mtx.gz
├── sample2/
│   └── ...
```

### Data Availability

Raw sequencing data have been deposited in the Gene Expression Omnibus (GEO) under accession number **GSE169246**.

---

## Usage

### Quick Start

```r
# Set working directory
setwd("/path/to/your/data")

# Source the analysis script
source("scripts/FINAL_ANALYSIS_CODE.R")
```

### Step-by-Step Execution

The script is organized into 9 sections that can be run sequentially:

| Section | Description | Output |
|---------|-------------|--------|
| 1 | Environment Setup | Loaded packages, directories |
| 2 | Data Loading & QC | Merged Seurat object |
| 3 | Normalization & Clustering | UMAP embeddings |
| 4 | Cell Type Annotation | Annotated cells |
| 5 | Hypothesis Testing | Statistical results |
| 6 | Supplementary Analyses | Module scores |
| 7 | Visualization | Publication figures (PDF) |
| 8 | Export Statistics | Summary tables (CSV) |
| 9 | Session Info | Reproducibility record |

---

## Output Files

### Figures (PDF, 300 DPI)

| Panel | Description | Dimensions |
|-------|-------------|------------|
| A | UMAP by Cell Type | 10×8 in |
| B | UMAP by Treatment Group | 10×8 in |
| C | CXCL13+ T Cell Proportion | 4×6 in |
| D | CD8+ Cytotoxicity Score | 5×6 in |
| E | CD8+ Exhaustion Score | 5×6 in |
| F | Senescence Markers DotPlot | 6×10 in |
| G | SASP (Senescent vs Non) | 6×10 in |
| H | SASP (CTRL vs NAC) | 6×10 in |
| I | Senescence-IL1B Correlation | 8×8 in |
| J | Fibroblast ECM Score | 5×6 in |
| K | CD8+ T Cell Markers | 6.5×3.5 in |
| L | Fibroblast ECM Markers | 5×3.5 in |

### Statistical Summary

All statistical tests and results are exported to `output/statistics/statistical_summary.csv`.

---

## Key Methods

### Batch Correction
- **Harmony** integration across samples

### Cell Type Annotation
- Marker-based scoring using canonical gene sets

### Statistical Tests
- **Wilcoxon rank-sum test**: Non-parametric comparisons
- **Fisher's exact test**: Proportion comparisons
- **Pearson correlation**: Expression correlations

### Visualization Notes

> **Important**: For DotPlot visualizations, we use `group.by` instead of `split.by` and set `scale = FALSE` to ensure proper color mapping and display absolute expression values.

---

## Troubleshooting

### Common Issues

1. **Memory Error**
   ```r
   # Increase memory limit (Linux/macOS)
   options(future.globals.maxSize = 8000 * 1024^2)
   ```

2. **SCP Installation Failed**
   ```r
   # Try installing dependencies first
   BiocManager::install(c("SingleCellExperiment", "scater"))
   remotes::install_github("zhanghao-njmu/SCP", dependencies = TRUE)
   ```

3. **DotPlot Color Mapping Issues**
   - Use `scale = FALSE` parameter
   - Create manual grouping columns instead of `split.by`

---

## Citation

If you use this code, please cite:

```bibtex
@article{YourLastName2024,
  title={Chemotherapy-induced senescence reshapes the tumor immune microenvironment in triple-negative breast cancer},
  author={Your Name and Collaborators},
  journal={Journal Name},
  year={2024},
  doi={10.xxxx/xxxxx}
}
```

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Contact

- **Author**: [Your Name]
- **Email**: [your.email@institution.edu]
- **Institution**: [Your Institution]

---

## Acknowledgments

- [Seurat](https://satijalab.org/seurat/) team for the analysis framework
- [SCP](https://github.com/zhanghao-njmu/SCP) for visualization tools
- [Harmony](https://github.com/immunogenomics/harmony) for batch correction

