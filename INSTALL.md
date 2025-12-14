# Installation Guide

## System Requirements

- **R version**: 4.3.0 or higher (tested on R 4.5)
- **Operating System**: macOS, Linux, or Windows
- **RAM**: Minimum 16GB recommended for single-cell analysis

## Quick Start

### Option 1: Using renv (Recommended)

The `renv.lock` file contains all package versions used in this analysis. This ensures exact reproducibility.

```r
# Install renv if not already installed
install.packages("renv")

# Restore all packages from renv.lock
renv::restore()
```

### Option 2: Manual Installation

If you prefer to install packages manually:

#### Core Packages

```r
# CRAN packages
install.packages(c(
  "Seurat",
  "ggplot2",
  "dplyr",
  "tidyr",
  "patchwork",
  "cowplot",
  "viridis",
  "RColorBrewer",
  "scales",
  "tibble",
  "ggpubr",
  "ggsignif",
  "ggbeeswarm",
  "ggridges",
  "ggdist",
  "wesanderson",
  "ggsci",
  "paletteer",
  "rstatix"
))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "clusterProfiler",
  "enrichplot",
  "DOSE",
  "org.Mm.eg.db",
  "ComplexHeatmap",
  "Nebulosa"
))
```

#### SCP Package (Single Cell Pipeline)

SCP package is available from GitHub:

```r
# Install devtools if not installed
install.packages("devtools")

# Install SCP from GitHub
devtools::install_github("zhanghao-njmu/SCP")
```

Alternatively, download the SCP package source from:
https://github.com/zhanghao-njmu/SCP

## Verification

To verify the installation:

```r
# Check if all packages are available
required_packages <- c(
  "Seurat", "SCP", "ggplot2", "dplyr", "patchwork",
  "ComplexHeatmap", "Nebulosa", "clusterProfiler"
)

for(pkg in required_packages) {
  if(requireNamespace(pkg, quietly = TRUE)) {
    cat("✓", pkg, "installed\n")
  } else {
    cat("✗", pkg, "NOT installed\n")
  }
}
```

## Troubleshooting

### Common Issues

1. **BiocManager version conflicts**
   ```r
   BiocManager::install(version = "3.18")  # or appropriate version
   ```

2. **SCP installation fails**
   - Try installing dependencies first
   - Check GitHub issues: https://github.com/zhanghao-njmu/SCP/issues

3. **Memory errors during analysis**
   - Increase R memory limit
   - Use raster mode for large datasets: `DimPlot(..., raster = TRUE)`

4. **Font issues in PDF output**
   - Install system fonts (Arial/Helvetica)
   - Use `pdf.options(useDingbats = FALSE)` before plotting

## Package Versions

Key package versions used in this analysis (check renv.lock for complete list):

| Package | Version |
|---------|---------|
| Seurat | 5.0+ |
| SCP | latest |
| ggplot2 | 3.4+ |
| ComplexHeatmap | 2.14+ |
| Nebulosa | 1.8+ |
| clusterProfiler | 4.6+ |

## Contact

For installation issues, please open an issue in this repository.
