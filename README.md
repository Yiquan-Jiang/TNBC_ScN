# Breast Cancer Spatial Transcriptomics Analysis

## Publication Code Repository

This repository contains the analysis code and visualization scripts for the spatial transcriptomics study of breast cancer, focusing on:

1. **IL1R1+ CAF (Cancer-Associated Fibroblasts) and CXCL13+ T cell colocalization**
2. **CDKN1A+ Macrophage and CAF spatial relationships**
3. **Tumor microenvironment cell-cell interactions**

---

## Table of Contents

- [Requirements](#requirements)
- [Data Structure](#data-structure)
- [Code Structure](#code-structure)
- [Usage](#usage)
- [Output Files](#output-files)
- [Visualization Settings](#visualization-settings)
- [Citation](#citation)

---

## Requirements

### Python Dependencies

```bash
pip install -r 00_requirements.txt
```

**Main packages:**
- scanpy >= 1.9.0
- pandas >= 1.5.0
- numpy >= 1.21.0
- matplotlib >= 3.5.0
- seaborn >= 0.11.0
- squidpy >= 1.2.0
- scikit-learn >= 1.1.0
- scipy >= 1.9.0
- h5py >= 3.7.0
- anndata >= 0.8.0

### R Dependencies

```r
# Install via Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("CellChat", "ComplexHeatmap"))
install.packages(c("Seurat", "dplyr", "ggplot2", "patchwork", "circlize", "future"))
```

---

## Data Structure

```
BRCA_spatial_transcriptom/
├── download_page_VISDS000449_spatial/    # 10X Visium spatial data
├── download_page_VISDS000450_spatial/
├── download_page_VISDS000451_spatial/
├── ...                                    # Additional samples
├── BRCA_GSE176078_TS.h5ad                # Single-cell reference (AnnData)
├── BRCA_GSE176078_TS.rds                 # Single-cell reference (Seurat)
└── publication_code/                      # Analysis scripts
```

Each spatial data folder should contain:
- `*.h5` - Filtered feature barcode matrix
- `spatial/` - Spatial coordinates and images

---

## Code Structure

### Analysis Pipeline

| Script | Description |
|--------|-------------|
| `spatial_colocalization_base.py` | Base module for spatial data loading and cell type annotation |
| `il1r1_cxcl13_extended_analysis.py` | IL1R1+ CAF and CXCL13+ T cell identification and colocalization |
| `comprehensive_analysis_with_tumor_and_cdkn1a.py` | Complete analysis including tumor cells and CDKN1A+ macrophages |
| `cellchat_analysis.R` | CellChat cell-cell communication analysis |
| `run_analysis.py` | Main runner script for executing analysis pipelines |

### Inheritance Structure

```
SpatialColocalizationAnalyzer (spatial_colocalization_base.py)
    └── IL1R1_CXCL13_ColocalizationAnalyzer (il1r1_cxcl13_extended_analysis.py)
        └── ComprehensiveSpatialAnalyzer (comprehensive_analysis_with_tumor_and_cdkn1a.py)
```

---

## Usage

### Quick Start: Run Complete Analysis

```bash
# Navigate to publication_code directory
cd publication_code

# Run comprehensive Python analysis
python 03_comprehensive_analysis_with_tumor_and_cdkn1a.py

# Run CellChat analysis (R)
Rscript 04_cellchat_analysis.R
```

### Step-by-Step Analysis

#### 1. Basic Spatial Colocalization

```python
from spatial_colocalization_base import SpatialColocalizationAnalyzer

analyzer = SpatialColocalizationAnalyzer(data_dir="..")
analyzer.load_spatial_data()
analyzer.annotate_cell_types()

# Colocalization analysis
coloc_results = analyzer.calculate_colocalization('IL1R1_pos_Fibroblasts', 'T_cells')
analyzer.plot_colocalization_results(coloc_results, 'IL1R1+ Fibroblasts vs T cells')
```

#### 2. Extended IL1R1/CXCL13 Analysis

```python
from il1r1_cxcl13_extended_analysis import IL1R1_CXCL13_ColocalizationAnalyzer

analyzer = IL1R1_CXCL13_ColocalizationAnalyzer(data_dir="..")
results = analyzer.run_extended_analysis()
```

#### 3. Comprehensive Analysis

```python
from comprehensive_analysis_with_tumor_and_cdkn1a import ComprehensiveSpatialAnalyzer

analyzer = ComprehensiveSpatialAnalyzer(data_dir="..")
results = analyzer.run_comprehensive_analysis()
```

---

## Output Files

### Visualization Outputs (PDF/PNG, 300 DPI)

| File | Description |
|------|-------------|
| `enhanced_all_samples_with_tumor_and_HE.pdf` | IL1R1+ CAF vs CXCL13+ T cell spatial distribution with H&E |
| `enhanced_all_samples_CDKN1A_Macrophage_CAF_with_HE.pdf` | CDKN1A+ Macrophage vs CAF spatial distribution with H&E |
| `comparative_spatial_analysis_summary.pdf` | Comparative analysis summary |
| `IL1R1_CAF_vs_CXCL13_T_colocalization_analysis.pdf` | Detailed IL1R1-CXCL13 analysis |
| `CDKN1A_Macrophage_vs_CAF_colocalization_analysis.pdf` | Detailed CDKN1A-CAF analysis |
| `spatial_distribution_*.pdf` | Per-sample spatial distributions |

### Data Outputs (CSV)

| File | Description |
|------|-------------|
| `IL1R1_CAF_CXCL13_T_colocalization_results.csv` | IL1R1-CXCL13 correlation results |
| `IL1R1_CAF_CXCL13_T_spatial_proximity_results.csv` | Spatial proximity analysis |
| `CDKN1A_Macrophage_CAF_colocalization_results.csv` | CDKN1A-CAF correlation results |
| `analysis_summary.csv` | Overall analysis summary |

### CellChat Outputs

| File | Description |
|------|-------------|
| `cellchat_object.rds` | Complete CellChat object |
| `communication_summary.csv` | Cell communication summary |
| `signaling_pathways.csv` | Signaling pathway analysis |
| `ligand_receptor_pairs.csv` | Ligand-receptor pairs |
| `cellchat_*.pdf` | Various CellChat visualizations |

---

## Visualization Settings

All figures are generated with publication-quality settings:

```python
plt.rcParams.update({
    'pdf.fonttype': 42,      # Editable text in PDF
    'ps.fonttype': 42,       # PostScript compatibility
    'font.family': 'Arial',
    'font.size': 12,
    'axes.titlesize': 25,
    'axes.labelsize': 20,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 20,
    'figure.titlesize': 25
})
```

**Output specifications:**
- **DPI**: 300
- **Format**: PDF (vector, editable) + PNG (raster)
- **Font**: Arial
- **Color scheme**: Colorblind-friendly

---

## Gene Signatures

### Cell Type Markers

| Cell Type | Marker Genes |
|-----------|--------------|
| **CAF** | COL1A1, COL1A2, COL3A1, DCN, LUM, POSTN, FAP, PDGFRA, PDGFRB, ACTA2, VIM, THY1 |
| **T cells** | CD3D, CD3E, CD3G, CD4, CD8A, CD8B, TRAC, TRBC1, TRBC2, LCK, ZAP70 |
| **CXCL13+ T** | CXCL13, PDCD1, CTLA4, TIGIT, LAG3, HAVCR2, TOX, ICOS |
| **Macrophages** | CD68, CD163, CD14, MARCO, MSR1, MRC1, FCGR3A, CSF1R |
| **Epithelial** | EPCAM, KRT7, KRT8, KRT18, KRT19, CDH1, MUC1 |

### Target Cell Identification

| Cell Population | Criteria |
|-----------------|----------|
| **IL1R1+ CAF** | CAF signature > 75th percentile AND IL1R1 expression > 75th percentile |
| **CXCL13+ T cells** | T cell signature > 75th percentile AND CXCL13 expression > 75th percentile |
| **CDKN1A+ Macrophages** | Macrophage signature > 75th percentile AND CDKN1A expression > 75th percentile |
| **Tumor cells** | Epithelial signature > 75th percentile |

---

## Statistical Methods

### Colocalization Analysis

1. **Score Correlation**: Pearson and Spearman correlations between composite scores
2. **Spatial Proximity**: Permutation test (n=1000) for nearest neighbor distances
3. **Neighbor Counting**: Mann-Whitney U test for cells within radius

### Multiple Sample Analysis

- Per-sample correlations with significance testing
- Overall pooled correlation analysis
- Abundance correlation between cell types

---

## Key Findings (Template)

After running the analysis, key results will be summarized:

1. **IL1R1+ CAF vs CXCL13+ T cells**:
   - Overall Pearson r: [value]
   - p-value: [value]
   - Samples with significant colocalization: [X/N]

2. **CDKN1A+ Macrophage vs CAF**:
   - Overall Pearson r: [value]
   - p-value: [value]
   - Samples with significant colocalization: [X/N]

---

## Citation

If you use this code in your research, please cite:

```
[Your paper citation here]
```

---

## License

[Add license information]

---

## Contact

For questions about the analysis code, please contact:
- [Author name and email]

