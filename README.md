# ScN Mouse Single-Cell RNA-seq Analysis Pipeline

## Overview

This repository contains the analysis scripts for single-cell RNA sequencing (scRNA-seq) data analysis of the Senolytic Combination (ScN) mouse tumor microenvironment study. The analysis focuses on characterizing cellular heterogeneity and treatment response across different cell populations including Cancer-Associated Fibroblasts (CAFs), Myeloid cells, and T cells.

## Project Structure

```
github_release/
├── README.md                    # This file
├── renv.lock                    # R package version lock file
├── .gitignore                   # Git ignore configuration
├── scripts/
│   ├── 01_data_processing.R     # Data loading and initial processing
│   ├── 02_umap_visualization.R  # UMAP visualization using SCP package
│   ├── 03_caf_analysis.R        # CAF subtype classification and analysis
│   ├── 04_differential_enrichment.R  # Differential expression and enrichment
│   ├── 05_comprehensive_visualization.R  # Publication-quality figures
│   ├── 06_merged_data_analysis.R    # Merged dataset analysis
│   └── 07_nebulosa_expression.R     # Expression density plots
```

## Analysis Pipeline

### 1. Data Processing (`01_data_processing.R`)
- Loading Seurat objects
- Sample filtering and treatment group assignment
- UMAP recalculation for filtered data
- CAF subtype classification using marker gene module scores

### 2. UMAP Visualization (`02_umap_visualization.R`)
- High-quality UMAP plots using SCP package
- Cell type and treatment group visualizations
- CAF subtype distribution plots
- Cell statistics plots

### 3. CAF Analysis (`03_caf_analysis.R`)
- CAF subtype classification (iCAF, myCAF, apCAF, matCAF)
- IL1R1 expression analysis across CAF subtypes
- ECM remodeling and inflammatory capacity scoring
- Matrix CAF and CXCL13+ CAF proportion analysis

### 4. Differential Expression & Enrichment (`04_differential_enrichment.R`)
- Differential expression analysis between CAF subtypes
- GO and KEGG pathway enrichment
- Feature heatmaps with functional annotations
- Volcano plots and enrichment network visualizations

### 5. Comprehensive Visualization (`05_comprehensive_visualization.R`)
- Myeloid cell senescence scoring
- IL1B expression analysis
- T cell reactivity (PDCD1+/GZMK+ tumor-reactive T cells)
- CAF ratio and proportion analysis
- Publication-ready multi-panel figures

### 6. Merged Data Analysis (`06_merged_data_analysis.R`)
- Integration of multiple datasets
- Cross-dataset cell type comparison
- Treatment response analysis across datasets

### 7. Expression Density Plots (`07_nebulosa_expression.R`)
- Nebulosa-based expression density visualization
- IL1R1 and COL1A1 expression patterns
- Co-expression analysis in CAF subtypes

## CAF Subtype Markers

| Subtype | Key Markers |
|---------|-------------|
| iCAF (Inflammatory) | Il6, Il1b, Cxcl1, Cxcl2, Ccl2, Tnf |
| myCAF (Myofibroblastic) | Acta2, Tagln, Myh11, Postn, Lox |
| apCAF (Antigen-presenting) | H2-Aa, H2-Ab1, Cd74, Ciita |
| matCAF (Matrix) | Col1a1, Col1a2, Col3a1, Fn1, Sparc |

## Treatment Groups

- **Untreated**: Control group
- **Chemo**: Chemotherapy treatment
- **Chemo_Abt263**: Chemotherapy + ABT-263 (senolytic)
- **Chemo_Abt263_IL1B**: Chemotherapy + ABT-263 + IL1B treatment

## Key Findings Visualized

1. **Myeloid Senescence**: Cellular aging signature distribution across treatment groups
2. **IL1B Expression**: Pro-inflammatory cytokine expression in myeloid cells
3. **T Cell Reactivity**: PDCD1+GZMK+ exhausted/reactive T cell proportions
4. **CAF Heterogeneity**: Distribution and functional characterization of CAF subtypes
5. **IL1R1 Expression**: IL1B receptor expression patterns in CAF subtypes

## Dependencies

### R Packages (Core)
- Seurat (≥5.0)
- SeuratObject
- SCP (Single Cell Pipeline)
- ggplot2, dplyr, tidyr
- patchwork, cowplot
- ComplexHeatmap
- Nebulosa
- clusterProfiler, enrichplot

### Additional Visualization Packages
- viridis, RColorBrewer
- ggsci, wesanderson, paletteer
- ggpubr, ggsignif
- ggridges, ggbeeswarm

## Reproducibility

This project uses `renv` for R package management. To restore the exact package versions used:

```r
# Install renv if not already installed
install.packages("renv")

# Restore packages from renv.lock
renv::restore()
```

## Figure Specifications

All figures are generated following publication standards:
- **Resolution**: 300 DPI minimum
- **Format**: Editable PDF (Type 42 fonts embedded)
- **Font**: Helvetica/Arial, consistent sizing
- **Color**: Colorblind-friendly palettes
- **Aspect Ratio**: Golden ratio (1.618) when appropriate

## Citation

If you use this code, please cite the associated publication.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contact

For questions regarding the analysis pipeline, please open an issue in this repository.
