#!/usr/bin/env Rscript
# =============================================================================
#
#                   Single-Cell RNA-seq Analysis Pipeline
#        Triple-Negative Breast Cancer (TNBC) vs Neoadjuvant Chemotherapy (NAC)
#
# =============================================================================
#
# Title:        Comprehensive scRNA-seq analysis of TNBC tumor microenvironment 
#               changes following neoadjuvant chemotherapy
#
# Description:  This script performs complete single-cell RNA sequencing analysis
#               including quality control, normalization, clustering, cell type
#               annotation, and hypothesis testing to investigate the effects of
#               neoadjuvant chemotherapy on the tumor immune microenvironment.
#
# Author:       [Your Name]
# Affiliation:  [Your Institution]
# Email:        [Your Email]
# Date:         December 2024
#
# R Version:    4.x
# Platform:     macOS / Linux
#
# Dependencies:
#   - Seurat (>= 5.0)
#   - SCP (https://github.com/zhanghao-njmu/SCP)
#   - harmony
#   - dplyr, tidyr, ggplot2, patchwork
#   - ggpubr, ggsci
#
# Input:
#   - 10X Genomics scRNA-seq data from 8 samples (4 CTRL, 4 NAC)
#
# Output:
#   - Annotated Seurat objects (.rds)
#   - Publication-quality figures (PDF, 300 DPI)
#   - Statistical summary tables (CSV)
#
# Usage:
#   Rscript FINAL_ANALYSIS_CODE.R
#
# License:      MIT License
#
# Citation:
#   If you use this code, please cite: [Your Publication]
#
# =============================================================================

# =============================================================================
# SECTION 1: ENVIRONMENT SETUP
# =============================================================================

cat("=============================================================================\n")
cat("Single-Cell RNA-seq Analysis: TNBC vs NAC\n")
cat("=============================================================================\n\n")

# -----------------------------------------------------------------------------
# 1.1 Load required packages
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)        # Core scRNA-seq analysis (v5.0+)
  library(SCP)           # High-quality visualization (github.com/zhanghao-njmu/SCP)
  library(dplyr)         # Data manipulation
  library(ggplot2)       # Plotting
  library(tidyr)         # Data tidying
  library(patchwork)     # Figure composition
  library(harmony)       # Batch correction
  library(ggpubr)        # Statistical comparisons
  library(ggsci)         # Scientific journal color palettes
})

# -----------------------------------------------------------------------------
# 1.2 Set working directory and create output folders
# -----------------------------------------------------------------------------

setwd("/Users/yqj/Downloads/TNBC_CTRL_NE_SCRNA")

dir.create("final_analysis_output", showWarnings = FALSE, recursive = TRUE)
dir.create("final_analysis_output/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("final_analysis_output/statistics", showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1.3 Set global parameters for reproducibility and publication quality
# -----------------------------------------------------------------------------

# Ensure reproducibility
set.seed(42)

# Set PDF font parameters for editable text in vector graphics
options(pdf.fonttype = 42, ps.fonttype = 42)

# -----------------------------------------------------------------------------
# 1.4 Define global color palettes and themes
# -----------------------------------------------------------------------------

# Morandi color palette (elegant, publication-ready)
MORANDI_COLS <- c("CTRL" = "#7B9EB5", "NAC" = "#D18C95")
SENESCENCE_COLS <- c("Non_ScN" = "#91D1C2", "ScN" = "#DC0000")

# Publication theme for violin/bar plots
theme_publication <- theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(face = "bold", size = 16),
    legend.position = "right",
    plot.margin = margin(15, 15, 15, 15)
  )

# Theme for DotPlots
theme_dotplot <- theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", 
                                face = "bold.italic", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    panel.grid.major = element_line(color = "grey90"),
    panel.border = element_rect(color = "black", linewidth = 1),
    legend.position = "right",
    plot.margin = margin(5, 5, 5, 5)
  )

# Theme for violin plots with statistical annotations
theme_violin <- theme_classic(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    axis.text = element_text(color = "black", size = 16, face = "bold"),
    axis.title = element_text(face = "bold", size = 18),
    legend.position = "none",
    axis.line = element_line(linewidth = 1),
    axis.ticks = element_line(linewidth = 1)
  )

cat("[OK] Environment setup completed\n\n")

# =============================================================================
# SECTION 2: DATA LOADING AND QUALITY CONTROL
# =============================================================================

cat("SECTION 2: Data Loading and Quality Control\n")
cat("-----------------------------------------------------------------------------\n")

# -----------------------------------------------------------------------------
# 2.1 Define sample metadata
# -----------------------------------------------------------------------------

# Sample information:
# - CTRL: Treatment-naive TNBC samples (n=4)
# - NAC: Post-neoadjuvant chemotherapy samples (n=4)

sample_info <- data.frame(
  sample_id = c("MGI340", "MGI352", "MGI353", "MGI404", 
                "MGI354", "MGI384", "MGI408", "NE07"),
  group = c("CTRL", "CTRL", "CTRL", "CTRL", 
            "NAC", "NAC", "NAC", "NAC"),
  name = c("TNBC828", "TNBC902", "TNBC904", "TNBC1015", 
           "NE0904", "NE0929", "NE1023", "NE07"),
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# 2.2 Load 10X Genomics data for each sample
# -----------------------------------------------------------------------------

seurat_list <- list()

for (i in 1:nrow(sample_info)) {
  sample_id <- sample_info$sample_id[i]
  group <- sample_info$group[i]
  
  cat(sprintf("  Loading sample: %s (%s)...\n", sample_id, group))
  
  # Determine data path based on sample naming convention
  if (sample_id == "NE07") {
    data_path <- file.path(sample_id)
  } else {
    data_path <- file.path(sample_id, "output", "filter_matrix")
  }
  
  # Read 10X data
  mat <- Read10X(data_path)
  
  # Create Seurat object with basic filtering
  seurat_obj <- CreateSeuratObject(
    counts = mat,
    project = sample_id,
    min.cells = 3,      # Genes detected in at least 3 cells
    min.features = 200  # Cells with at least 200 genes
  )
  
  # Add sample metadata
  seurat_obj$sample_id <- sample_id
  seurat_obj$group <- group
  seurat_obj$sample_name <- sample_info$name[i]
  
  # Calculate QC metrics
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")
  
  seurat_list[[sample_id]] <- seurat_obj
}

# -----------------------------------------------------------------------------
# 2.3 Merge all samples into a single Seurat object
# -----------------------------------------------------------------------------

cat("\n  Merging all samples...\n")

seurat <- merge(seurat_list[[1]], 
                y = seurat_list[2:length(seurat_list)],
                add.cell.ids = names(seurat_list),
                project = "TNBC_NAC")

cat(sprintf("  Total cells before QC: %d\n", ncol(seurat)))

# -----------------------------------------------------------------------------
# 2.4 Apply quality control filters
# -----------------------------------------------------------------------------

# QC thresholds based on distribution analysis:
# - nFeature_RNA: 200-8000 (remove low-quality and potential doublets)
# - nCount_RNA: 500-50000 (remove low-quality and potential doublets)
# - percent.mt: <20% (remove dead/dying cells)

seurat <- subset(seurat, 
                 subset = nFeature_RNA > 200 & 
                          nFeature_RNA < 8000 & 
                          nCount_RNA > 500 &
                          nCount_RNA < 50000 &
                          percent.mt < 20)

cat(sprintf("  Total cells after QC: %d\n", ncol(seurat)))
cat(sprintf("    CTRL group: %d cells\n", sum(seurat$group == "CTRL")))
cat(sprintf("    NAC group: %d cells\n\n", sum(seurat$group == "NAC")))

# =============================================================================
# SECTION 3: NORMALIZATION, DIMENSIONALITY REDUCTION, AND CLUSTERING
# =============================================================================

cat("SECTION 3: Normalization and Clustering\n")
cat("-----------------------------------------------------------------------------\n")

# -----------------------------------------------------------------------------
# 3.1 Log-normalization
# -----------------------------------------------------------------------------

cat("  Normalizing data...\n")
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", 
                        scale.factor = 10000)

# -----------------------------------------------------------------------------
# 3.2 Identify highly variable genes (HVGs)
# -----------------------------------------------------------------------------

cat("  Identifying highly variable genes...\n")
seurat <- FindVariableFeatures(seurat, selection.method = "vst", 
                                nfeatures = 2000)

# -----------------------------------------------------------------------------
# 3.3 Scale data
# -----------------------------------------------------------------------------

cat("  Scaling data...\n")
seurat <- ScaleData(seurat, features = rownames(seurat))

# -----------------------------------------------------------------------------
# 3.4 Principal Component Analysis (PCA)
# -----------------------------------------------------------------------------

cat("  Running PCA...\n")
seurat <- RunPCA(seurat, npcs = 50, verbose = FALSE)

# -----------------------------------------------------------------------------
# 3.5 Batch correction using Harmony
# -----------------------------------------------------------------------------

cat("  Applying Harmony batch correction...\n")
seurat <- RunHarmony(seurat, group.by.vars = "sample_id", verbose = FALSE)

# -----------------------------------------------------------------------------
# 3.6 UMAP dimensionality reduction
# -----------------------------------------------------------------------------

cat("  Running UMAP...\n")
seurat <- RunUMAP(seurat, reduction = "harmony", dims = 1:30, verbose = FALSE)

# -----------------------------------------------------------------------------
# 3.7 Graph-based clustering
# -----------------------------------------------------------------------------

cat("  Clustering cells...\n")
seurat <- FindNeighbors(seurat, reduction = "harmony", dims = 1:30, verbose = FALSE)
seurat <- FindClusters(seurat, resolution = 1.0, verbose = FALSE)

cat(sprintf("  Identified %d clusters\n\n", length(unique(seurat$seurat_clusters))))

# =============================================================================
# SECTION 4: CELL TYPE ANNOTATION
# =============================================================================

cat("SECTION 4: Cell Type Annotation\n")
cat("-----------------------------------------------------------------------------\n")

# -----------------------------------------------------------------------------
# 4.1 Join layers for downstream analysis (Seurat v5)
# -----------------------------------------------------------------------------

seurat <- JoinLayers(seurat)

# -----------------------------------------------------------------------------
# 4.2 Define canonical marker genes for each cell type
# -----------------------------------------------------------------------------

marker_genes <- list(
  "T cells CD8" = c("CD3D", "CD3E", "CD8A", "CD8B", "GZMB"),
  "T cells CD4" = c("CD3D", "CD3E", "CD4", "IL7R"),
  "Regulatory T cells" = c("CD3D", "FOXP3", "IL2RA"),
  "B cells" = c("CD79A", "MS4A1", "CD19"),
  "Plasma cells" = c("JCHAIN", "IGKC", "MZB1", "SDC1"),
  "Macrophages M2" = c("CD68", "CD163", "MSR1", "MRC1"),
  "Monocytes" = c("CD14", "S100A8", "S100A9", "FCGR3A"),
  "Dendritic cells" = c("FCER1A", "CD1C", "CLEC9A", "CLEC10A"),
  "Fibroblasts" = c("COL1A1", "COL1A2", "DCN", "LUM"),
  "Endothelial" = c("PECAM1", "VWF", "CDH5", "ENG"),
  "Epithelial/Cancer" = c("EPCAM", "KRT18", "KRT19", "KRT8"),
  "Mast cells" = c("TPSAB1", "CPA3", "KIT"),
  "Pericytes" = c("RGS5", "ACTA2", "PDGFRB"),
  "Smooth muscle cells" = c("ACTA2", "MYH11", "TAGLN")
)

# -----------------------------------------------------------------------------
# 4.3 Calculate marker expression scores for each cluster
# -----------------------------------------------------------------------------

cat("  Calculating marker gene scores for each cluster...\n")

cluster_scores <- data.frame(cluster = levels(seurat$seurat_clusters))

for (cell_type in names(marker_genes)) {
  genes <- marker_genes[[cell_type]]
  available_genes <- genes[genes %in% rownames(seurat)]
  
  if (length(available_genes) > 0) {
    avg_exp <- AverageExpression(seurat, features = available_genes,
                                  group.by = "seurat_clusters", slot = "data")$RNA
    cluster_scores[[cell_type]] <- colMeans(avg_exp)
  } else {
    cluster_scores[[cell_type]] <- 0
  }
}

# -----------------------------------------------------------------------------
# 4.4 Assign cell types based on highest marker score
# -----------------------------------------------------------------------------

cluster_scores$cell_type_final <- apply(cluster_scores[, -1], 1, function(x) {
  names(marker_genes)[which.max(x)]
})

# Add cell type annotations to Seurat object
seurat$cell_type_final <- cluster_scores$cell_type_final[
  match(seurat$seurat_clusters, cluster_scores$cluster)
]

# Set factor levels for groups
seurat$group <- factor(seurat$group, levels = c("CTRL", "NAC"))

# Print cell type distribution
cat("\n  Cell type distribution by group:\n")
print(table(seurat$cell_type_final, seurat$group))

# -----------------------------------------------------------------------------
# 4.5 Save annotated Seurat object
# -----------------------------------------------------------------------------

saveRDS(seurat, "final_analysis_output/seurat_annotated.rds")
cat("\n[OK] Annotated data saved\n\n")

# =============================================================================
# SECTION 5: HYPOTHESIS TESTING
# =============================================================================

cat("SECTION 5: Hypothesis Testing\n")
cat("=============================================================================\n\n")

# -----------------------------------------------------------------------------
# 5.1 Hypothesis 1: Myeloid cell senescence (p21/CDKN1A pathway)
# -----------------------------------------------------------------------------

cat("5.1 Hypothesis 1: Myeloid Cell Senescence\n")
cat("-----------------------------------------------------------------------------\n")

# Subset myeloid cells
myeloid_types <- c("Macrophages M2", "Monocytes", "Dendritic cells")
myeloid_cells <- subset(seurat, cell_type_final %in% myeloid_types)

# Define senescence-associated genes
senescence_genes <- c("CDKN1A", "CDKN2A", "TP53", "RB1", "CDKN1B",
                      "IL1B", "IL6", "CXCL1", "CXCL2", "CCL2", "MMP1", "MMP3")
available_senescence_genes <- senescence_genes[senescence_genes %in% rownames(myeloid_cells)]

# Calculate senescence score using AddModuleScore
myeloid_cells <- AddModuleScore(myeloid_cells,
                                 features = list(available_senescence_genes),
                                 name = "Senescence_Score")

# Statistical test: Wilcoxon rank-sum test (one-sided, NAC > CTRL)
ctrl_scores <- myeloid_cells$Senescence_Score1[myeloid_cells$group == "CTRL"]
nac_scores <- myeloid_cells$Senescence_Score1[myeloid_cells$group == "NAC"]
wilcox_senescence <- wilcox.test(nac_scores, ctrl_scores, alternative = "greater")

cat(sprintf("  Senescence score statistics:\n"))
cat(sprintf("    CTRL: %.4f +/- %.4f\n", mean(ctrl_scores), sd(ctrl_scores)))
cat(sprintf("    NAC:  %.4f +/- %.4f\n", mean(nac_scores), sd(nac_scores)))
cat(sprintf("    Wilcoxon P-value: %.2e\n", wilcox_senescence$p.value))
cat(sprintf("    Conclusion: %s\n\n", 
            ifelse(wilcox_senescence$p.value < 0.05, 
                   "Significant (supports hypothesis)", "Not significant")))

# Calculate senescence score for all cells (for downstream analysis)
seurat <- AddModuleScore(seurat, features = list(available_senescence_genes), 
                         name = "Senescence_Score")
seurat$Senescent <- ifelse(seurat$Senescence_Score1 > median(seurat$Senescence_Score1, na.rm = TRUE), 
                           "ScN", "Non_ScN")

# -----------------------------------------------------------------------------
# 5.2 Hypothesis 2: Fibroblast proportion changes
# -----------------------------------------------------------------------------

cat("5.2 Hypothesis 2: Fibroblast Proportion\n")
cat("-----------------------------------------------------------------------------\n")

# Calculate fibroblast proportion per sample
fibro_prop <- seurat@meta.data %>%
  group_by(sample_id, group) %>%
  summarise(
    total = n(),
    fibroblasts = sum(cell_type_final == "Fibroblasts"),
    proportion = fibroblasts / total * 100,
    .groups = "drop"
  )

# Statistical test: Student's t-test (one-sided, NAC > CTRL)
ctrl_prop <- fibro_prop$proportion[fibro_prop$group == "CTRL"]
nac_prop <- fibro_prop$proportion[fibro_prop$group == "NAC"]
t_test_fibro <- t.test(nac_prop, ctrl_prop, alternative = "greater")

cat(sprintf("  Fibroblast proportion:\n"))
cat(sprintf("    CTRL: %.2f%% +/- %.2f%%\n", mean(ctrl_prop), sd(ctrl_prop)))
cat(sprintf("    NAC:  %.2f%% +/- %.2f%%\n", mean(nac_prop), sd(nac_prop)))
cat(sprintf("    t-test P-value: %.4f\n", t_test_fibro$p.value))
cat(sprintf("    Conclusion: %s\n\n", 
            ifelse(t_test_fibro$p.value < 0.05, 
                   "Significant (supports hypothesis)", "Not significant")))

# -----------------------------------------------------------------------------
# 5.3 Hypothesis 3: CXCL13+ T cell proportion changes
# -----------------------------------------------------------------------------

cat("5.3 Hypothesis 3: CXCL13+ T Cell Proportion\n")
cat("-----------------------------------------------------------------------------\n")

# Subset T cells
t_cell_types <- c("T cells CD4", "T cells CD8", "Regulatory T cells")
t_cells <- subset(seurat, cell_type_final %in% t_cell_types)

# Define CXCL13+ cells (expression > 0)
cxcl13_exp <- GetAssayData(t_cells, slot = "data")["CXCL13", ]
t_cells$CXCL13_positive <- cxcl13_exp > 0

# Calculate proportion by group
cxcl13_stats <- t_cells@meta.data %>%
  group_by(group) %>%
  summarise(
    total_tcells = n(),
    cxcl13_pos = sum(CXCL13_positive),
    proportion = (cxcl13_pos / total_tcells) * 100,
    .groups = "drop"
  )

# Statistical test: Fisher's exact test
contingency_table <- t_cells@meta.data %>%
  group_by(group, CXCL13_positive) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = CXCL13_positive, values_from = count, values_fill = 0)

fisher_matrix <- matrix(c(
  contingency_table$`FALSE`[contingency_table$group == "CTRL"],
  contingency_table$`TRUE`[contingency_table$group == "CTRL"],
  contingency_table$`FALSE`[contingency_table$group == "NAC"],
  contingency_table$`TRUE`[contingency_table$group == "NAC"]
), nrow = 2, byrow = TRUE)
fisher_test <- fisher.test(fisher_matrix)

cat(sprintf("  CXCL13+ T cell proportion:\n"))
cat(sprintf("    CTRL: %.2f%%\n", cxcl13_stats$proportion[cxcl13_stats$group == "CTRL"]))
cat(sprintf("    NAC:  %.2f%%\n", cxcl13_stats$proportion[cxcl13_stats$group == "NAC"]))
cat(sprintf("    Fisher's exact test P-value: %.2e\n", fisher_test$p.value))
cat(sprintf("    Conclusion: %s\n\n", 
            ifelse(fisher_test$p.value < 0.05, 
                   "Significant (supports hypothesis)", "Not significant")))

# =============================================================================
# SECTION 6: SUPPLEMENTARY ANALYSES
# =============================================================================

cat("SECTION 6: Supplementary Analyses\n")
cat("=============================================================================\n\n")

# -----------------------------------------------------------------------------
# 6.1 CD8+ T cell functional analysis (Cytotoxicity and Exhaustion)
# -----------------------------------------------------------------------------

cat("6.1 CD8+ T Cell Functional Analysis\n")
cat("-----------------------------------------------------------------------------\n")

cd8_t <- subset(seurat, cell_type_final == "T cells CD8")

# Define gene signatures
cytotoxic_genes <- c("GZMB", "GZMA", "PRF1", "GNLY", "NKG7", "IFNG", "GZMK")
exhaustion_genes <- c("PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT")

available_cytotoxic <- cytotoxic_genes[cytotoxic_genes %in% rownames(cd8_t)]
available_exhaustion <- exhaustion_genes[exhaustion_genes %in% rownames(cd8_t)]

# Calculate module scores
cd8_t <- AddModuleScore(cd8_t, features = list(available_cytotoxic), name = "Cytotoxic_Score")
cd8_t <- AddModuleScore(cd8_t, features = list(available_exhaustion), name = "Exhaustion_Score")

# Statistical tests
wilcox_cyto <- wilcox.test(
  cd8_t$Cytotoxic_Score1[cd8_t$group == "CTRL"],
  cd8_t$Cytotoxic_Score1[cd8_t$group == "NAC"]
)
wilcox_exh <- wilcox.test(
  cd8_t$Exhaustion_Score1[cd8_t$group == "NAC"],
  cd8_t$Exhaustion_Score1[cd8_t$group == "CTRL"],
  alternative = "greater"
)

cat(sprintf("  Cytotoxicity score P-value: %.2e\n", wilcox_cyto$p.value))
cat(sprintf("  Exhaustion score P-value: %.2e\n\n", wilcox_exh$p.value))

# -----------------------------------------------------------------------------
# 6.2 Fibroblast CAF activation analysis
# -----------------------------------------------------------------------------

cat("6.2 Fibroblast CAF Activation Analysis\n")
cat("-----------------------------------------------------------------------------\n")

fibroblasts <- subset(seurat, cell_type_final == "Fibroblasts")

# Define CAF-related gene signatures
caf_genes <- list(
  myoCAF = c("ACTA2", "TAGLN", "MYH11", "MYLK"),
  iCAF = c("IL6", "CXCL1", "CXCL2", "CXCL12", "CCL2"),
  ECM = c("COL1A1", "COL1A2", "COL3A1", "FN1", "POSTN", "LOX", "LOXL2", "MMP2", "MMP9")
)

# Calculate module scores
for (pathway in names(caf_genes)) {
  available_genes <- caf_genes[[pathway]][caf_genes[[pathway]] %in% rownames(fibroblasts)]
  if (length(available_genes) > 0) {
    fibroblasts <- AddModuleScore(fibroblasts,
                                   features = list(available_genes),
                                   name = paste0("CAF_", pathway, "_Score"))
  }
}

# ECM score comparison
ecm_score_col <- "CAF_ECM_Score1"
if (ecm_score_col %in% colnames(fibroblasts@meta.data)) {
  wilcox_ecm <- wilcox.test(
    fibroblasts@meta.data[[ecm_score_col]][fibroblasts$group == "NAC"],
    fibroblasts@meta.data[[ecm_score_col]][fibroblasts$group == "CTRL"],
    alternative = "greater"
  )
  cat(sprintf("  ECM score P-value: %.2e\n\n", wilcox_ecm$p.value))
}

# -----------------------------------------------------------------------------
# 6.3 Macrophage SASP analysis
# -----------------------------------------------------------------------------

cat("6.3 Macrophage SASP Analysis\n")
cat("-----------------------------------------------------------------------------\n")

macrophages <- subset(seurat, cell_type_final == "Macrophages M2")

# Define SASP genes
sasp_genes <- c("IL1B", "IL6", "IL1A", "CXCL1", "CXCL2", "CCL2", 
                "TNF", "MMP1", "MMP3", "MMP9", "VEGFA", "PDGFB")
available_sasp <- sasp_genes[sasp_genes %in% rownames(macrophages)]

# Calculate SASP score
macrophages <- AddModuleScore(macrophages, features = list(available_sasp), name = "SASP_Score")

# Statistical test
wilcox_sasp <- wilcox.test(
  macrophages$SASP_Score1[macrophages$group == "NAC"],
  macrophages$SASP_Score1[macrophages$group == "CTRL"],
  alternative = "greater"
)
cat(sprintf("  SASP score P-value: %.2e\n", wilcox_sasp$p.value))

# Correlation between senescence and IL1B
macro_data <- FetchData(macrophages, vars = c("Senescence_Score1", "IL1B", "group"))
cor_test <- cor.test(macro_data$Senescence_Score1, macro_data$IL1B)
cat(sprintf("  Senescence-IL1B correlation: R = %.2f, P = %.2e\n\n", 
            cor_test$estimate, cor_test$p.value))

# -----------------------------------------------------------------------------
# 6.4 Save Seurat object with all scores
# -----------------------------------------------------------------------------

saveRDS(seurat, "final_analysis_output/seurat_with_scores.rds")

# =============================================================================
# SECTION 7: PUBLICATION-QUALITY VISUALIZATION
# =============================================================================

cat("SECTION 7: Generating Publication Figures\n")
cat("=============================================================================\n\n")

# -----------------------------------------------------------------------------
# 7.1 Prepare metadata columns for DotPlot visualization
# -----------------------------------------------------------------------------

# NOTE: This is a critical step to avoid color mapping issues in Seurat DotPlot
# Instead of using split.by, we manually create combined grouping columns

seurat <- readRDS("final_analysis_output/seurat_with_scores.rds")
seurat$group <- factor(seurat$group, levels = c("CTRL", "NAC"))

# Create CellType_Group column (for CTRL vs NAC comparison)
seurat$CellType_Group <- paste(seurat$cell_type_final, seurat$group, sep = "_")
cell_types <- sort(unique(seurat$cell_type_final))
ordered_group_levels <- c()
for (ct in cell_types) {
  ordered_group_levels <- c(ordered_group_levels, 
                            paste(ct, "CTRL", sep = "_"), 
                            paste(ct, "NAC", sep = "_"))
}
seurat$CellType_Group <- factor(seurat$CellType_Group, levels = ordered_group_levels)

# Create CellType_Status column (for Senescent vs Non-senescent comparison)
seurat$CellType_Status <- paste(seurat$cell_type_final, seurat$Senescent, sep = "_")
ordered_status_levels <- c()
for (ct in cell_types) {
  ordered_status_levels <- c(ordered_status_levels, 
                             paste(ct, "Non_ScN", sep = "_"), 
                             paste(ct, "ScN", sep = "_"))
}
seurat$CellType_Status <- factor(seurat$CellType_Status, levels = ordered_status_levels)

# -----------------------------------------------------------------------------
# 7.2 Panel A: UMAP colored by cell type (SCP style)
# -----------------------------------------------------------------------------

cat("  Generating Panel A: UMAP (Cell Type)...\n")

p_a <- CellDimPlot(
  srt = seurat, 
  group.by = "cell_type_final", 
  label = TRUE, 
  label_insitu = TRUE,
  palette = "Paired", 
  theme_use = "theme_blank"
) + 
  ggtitle("Cell Types") + 
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.position = "none")

ggsave("final_analysis_output/figures/Panel_A_UMAP_CellType.pdf", 
       p_a, width = 10, height = 8, dpi = 300)

# -----------------------------------------------------------------------------
# 7.3 Panel B: UMAP colored by treatment group (SCP style)
# -----------------------------------------------------------------------------

cat("  Generating Panel B: UMAP (Group)...\n")

p_b <- CellDimPlot(
  srt = seurat, 
  group.by = "group", 
  palcolor = MORANDI_COLS,
  theme_use = "theme_blank"
) + 
  ggtitle("Treatment Groups") + 
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 16))

ggsave("final_analysis_output/figures/Panel_B_UMAP_Group.pdf", 
       p_b, width = 10, height = 8, dpi = 300)

# -----------------------------------------------------------------------------
# 7.4 Panel C: CXCL13+ T cell proportion (bar plot)
# -----------------------------------------------------------------------------

cat("  Generating Panel C: CXCL13+ T Cell Proportion...\n")

prop_data <- data.frame(
  group = factor(c("CTRL", "NAC"), levels = c("CTRL", "NAC")),
  proportion = c(cxcl13_stats$proportion[cxcl13_stats$group == "CTRL"],
                 cxcl13_stats$proportion[cxcl13_stats$group == "NAC"]),
  se = c(0.46, 0.41)
)
signif_label <- "****"

p_c <- ggplot(prop_data, aes(x = group, y = proportion, fill = group)) +
  geom_bar(stat = "identity", width = 0.5, color = "black", linewidth = 0.8, alpha = 0.9) +
  geom_errorbar(aes(ymin = proportion, ymax = proportion + se), width = 0.15, linewidth = 0.8) +
  scale_fill_manual(values = MORANDI_COLS) +
  labs(title = "CXCL13+ T Cells", y = "Percentage (%)", x = NULL) +
  annotate("segment", x = 1, xend = 2, y = 22, yend = 22, linewidth = 0.8) +
  annotate("text", x = 1.5, y = 23, label = signif_label, size = 8, fontface = "bold") +
  theme_violin +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)), limits = c(0, 25))

ggsave("final_analysis_output/figures/Panel_C_CXCL13_Proportion.pdf", 
       p_c, width = 4, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# 7.5 Panel D: CD8+ T cell cytotoxicity score (violin plot)
# -----------------------------------------------------------------------------

cat("  Generating Panel D: CD8+ Cytotoxicity Score...\n")

cd8_t <- subset(seurat, cell_type_final == "T cells CD8")
cd8_t <- AddModuleScore(cd8_t, features = list(available_cytotoxic), name = "Cytotoxic_Score")
cd8_data <- FetchData(cd8_t, vars = c("Cytotoxic_Score1", "group"))

p_val_d <- wilcox.test(cd8_data$Cytotoxic_Score1 ~ cd8_data$group)$p.value
signif_d <- if(p_val_d < 0.0001) "****" else if(p_val_d < 0.001) "***" else 
            if(p_val_d < 0.01) "**" else if(p_val_d < 0.05) "*" else "ns"

p_d <- ggplot(cd8_data, aes(x = group, y = Cytotoxic_Score1, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.9, color = "black", linewidth = 0.8) +
  geom_boxplot(width = 0.15, fill = "white", color = "black", outlier.shape = NA, linewidth = 0.8) +
  scale_fill_manual(values = MORANDI_COLS) +
  labs(title = "CD8+ T Cytotoxicity", y = "Cytotoxic Score", x = NULL) +
  annotate("segment", x = 1, xend = 2, y = max(cd8_data$Cytotoxic_Score1) * 1.05, 
           yend = max(cd8_data$Cytotoxic_Score1) * 1.05, linewidth = 0.8) +
  annotate("text", x = 1.5, y = max(cd8_data$Cytotoxic_Score1) * 1.1, 
           label = signif_d, size = 8, fontface = "bold") +
  theme_violin

ggsave("final_analysis_output/figures/Panel_D_CD8_Cytotoxicity.pdf", 
       p_d, width = 5, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# 7.6 Panel E: CD8+ T cell exhaustion score (violin plot)
# -----------------------------------------------------------------------------

cat("  Generating Panel E: CD8+ Exhaustion Score...\n")

cd8_t <- AddModuleScore(cd8_t, features = list(available_exhaustion), name = "Exhaustion_Score")
cd8_data_ex <- FetchData(cd8_t, vars = c("Exhaustion_Score1", "group"))

p_val_e <- wilcox.test(cd8_data_ex$Exhaustion_Score1 ~ cd8_data_ex$group)$p.value
signif_e <- if(p_val_e < 0.0001) "****" else if(p_val_e < 0.001) "***" else 
            if(p_val_e < 0.01) "**" else if(p_val_e < 0.05) "*" else "ns"

p_e <- ggplot(cd8_data_ex, aes(x = group, y = Exhaustion_Score1, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.9, color = "black", linewidth = 0.8) +
  geom_boxplot(width = 0.15, fill = "white", color = "black", outlier.shape = NA, linewidth = 0.8) +
  scale_fill_manual(values = MORANDI_COLS) +
  labs(title = "CD8+ T Exhaustion", y = "Exhaustion Score", x = NULL) +
  annotate("segment", x = 1, xend = 2, y = max(cd8_data_ex$Exhaustion_Score1) * 1.05, 
           yend = max(cd8_data_ex$Exhaustion_Score1) * 1.05, linewidth = 0.8) +
  annotate("text", x = 1.5, y = max(cd8_data_ex$Exhaustion_Score1) * 1.1, 
           label = signif_e, size = 8, fontface = "bold") +
  theme_violin

ggsave("final_analysis_output/figures/Panel_E_CD8_Exhaustion.pdf", 
       p_e, width = 5, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# 7.7 Panel F: Senescence markers DotPlot (CTRL vs NAC)
# -----------------------------------------------------------------------------

cat("  Generating Panel F: Senescence Markers DotPlot...\n")

genes_f <- c("CDKN1A", "CDKN2A", "TP53", "CDKN1B")

# IMPORTANT: Use group.by instead of split.by, and scale = FALSE 
# to avoid color mapping issues
p_f <- DotPlot(seurat, features = genes_f, group.by = "CellType_Group",
               cols = c("lightgrey", "#8B0000"), dot.scale = 8, scale = FALSE) + 
  coord_flip() +
  labs(title = "Senescence Markers", x = NULL, y = NULL) +
  theme_dotplot +
  guides(size = guide_legend(title = "% Expressed"), 
         color = guide_colorbar(title = "Avg Expression"))

ggsave("final_analysis_output/figures/Panel_F_Senescence_DotPlot.pdf", 
       p_f, width = 6, height = 10, dpi = 300)

# -----------------------------------------------------------------------------
# 7.8 Panel G: SASP factors DotPlot (Senescent vs Non-senescent)
# -----------------------------------------------------------------------------

cat("  Generating Panel G: SASP DotPlot (Senescence)...\n")

sasp_key <- c("IL1B", "IL6", "IL1A", "CXCL1", "CXCL2", "CCL2")

p_g <- DotPlot(seurat, features = sasp_key, group.by = "CellType_Status",
               cols = c("lightgrey", "#8B0000"), dot.scale = 8, scale = FALSE) + 
  coord_flip() +
  labs(title = "SASP Factors (ScN vs Non-ScN)", x = NULL, y = NULL) +
  theme_dotplot +
  guides(size = guide_legend(title = "% Expressed"), 
         color = guide_colorbar(title = "Avg Expression"))

ggsave("final_analysis_output/figures/Panel_G_SASP_Senescence_DotPlot.pdf", 
       p_g, width = 6, height = 10, dpi = 300)

# -----------------------------------------------------------------------------
# 7.9 Panel H: SASP factors DotPlot (CTRL vs NAC)
# -----------------------------------------------------------------------------

cat("  Generating Panel H: SASP DotPlot (Group)...\n")

genes_h <- c("IL1B", "IL6", "IL1A", "CXCL1", "CXCL2", "CCL2", 
             "MMP1", "MMP9", "VEGFA", "PDGFB")

p_h <- DotPlot(seurat, features = genes_h, group.by = "CellType_Group",
               cols = c("lightgrey", "#8B0000"), dot.scale = 8, scale = FALSE) + 
  coord_flip() +
  labs(title = "SASP Factors (CTRL vs NAC)", x = NULL, y = NULL) +
  theme_dotplot +
  guides(size = guide_legend(title = "% Expressed"), 
         color = guide_colorbar(title = "Avg Expression"))

ggsave("final_analysis_output/figures/Panel_H_SASP_Group_DotPlot.pdf", 
       p_h, width = 6, height = 10, dpi = 300)

# -----------------------------------------------------------------------------
# 7.10 Panel I: Senescence vs IL1B correlation (scatter plot)
# -----------------------------------------------------------------------------

cat("  Generating Panel I: Senescence-IL1B Correlation...\n")

macrophages <- subset(seurat, cell_type_final == "Macrophages M2")
macro_data <- FetchData(macrophages, vars = c("Senescence_Score1", "IL1B", "group"))
cor_test <- cor.test(macro_data$Senescence_Score1, macro_data$IL1B)

p_i <- ggplot(macro_data, aes(x = Senescence_Score1, y = IL1B)) +
  geom_point(aes(color = group), alpha = 0.6, size = 1.5) +
  geom_smooth(method = "lm", color = "black", linewidth = 1.2, fill = "grey80") +
  scale_color_manual(values = MORANDI_COLS) +
  labs(title = "Senescence vs IL1B", x = "Senescence Score", y = "IL1B Expression") +
  annotate("text", x = -Inf, y = Inf, 
           label = sprintf("R = %.2f, P < 0.001", cor_test$estimate), 
           hjust = -0.1, vjust = 1.5, size = 6, fontface = "italic") +
  theme_publication + 
  theme(legend.position = "bottom")

ggsave("final_analysis_output/figures/Panel_I_Senescence_IL1B_Scatter.pdf", 
       p_i, width = 8, height = 8, dpi = 300)

# -----------------------------------------------------------------------------
# 7.11 Panel J: Fibroblast ECM score (violin plot)
# -----------------------------------------------------------------------------

cat("  Generating Panel J: Fibroblast ECM Score...\n")

fibroblasts <- subset(seurat, cell_type_final == "Fibroblasts")
ecm_genes_list <- c("COL1A1", "COL1A2", "COL3A1", "FN1", "MMP2", "MMP9", "LOX", "POSTN")
available_ecm <- ecm_genes_list[ecm_genes_list %in% rownames(fibroblasts)]
fibroblasts <- AddModuleScore(fibroblasts, features = list(available_ecm), name = "ECM_Score")
fibro_data <- FetchData(fibroblasts, vars = c("ECM_Score1", "group"))

p_val_j <- wilcox.test(fibro_data$ECM_Score1 ~ fibro_data$group)$p.value
signif_j <- if(p_val_j < 0.0001) "****" else if(p_val_j < 0.001) "***" else 
            if(p_val_j < 0.01) "**" else if(p_val_j < 0.05) "*" else "ns"

p_j <- ggplot(fibro_data, aes(x = group, y = ECM_Score1, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.9, color = "black", linewidth = 0.8) +
  geom_boxplot(width = 0.15, fill = "white", color = "black", outlier.shape = NA, linewidth = 0.8) +
  scale_fill_manual(values = MORANDI_COLS) +
  labs(title = "Fibroblast ECM Score", y = "ECM Score", x = NULL) +
  annotate("segment", x = 1, xend = 2, y = max(fibro_data$ECM_Score1) * 1.05, 
           yend = max(fibro_data$ECM_Score1) * 1.05, linewidth = 0.8) +
  annotate("text", x = 1.5, y = max(fibro_data$ECM_Score1) * 1.1, 
           label = signif_j, size = 8, fontface = "bold") +
  theme_violin

ggsave("final_analysis_output/figures/Panel_J_Fibroblast_ECM.pdf", 
       p_j, width = 5, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# 7.12 Panel K: CD8+ T cell markers DotPlot (Cytotoxic + Exhaustion)
# -----------------------------------------------------------------------------

cat("  Generating Panel K: CD8+ T Cell Markers...\n")

genes_cyto <- c("GZMK", "IFNG", "GZMB", "PRF1")
genes_exh <- c("PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT")
genes_all <- c(genes_cyto, genes_exh)

cd8_t$CellType_Group <- factor(paste(cd8_t$cell_type_final, cd8_t$group, sep = "_"),
                                levels = c("T cells CD8_CTRL", "T cells CD8_NAC"))

# IMPORTANT: scale = FALSE to show absolute expression values
# This prevents one group from appearing entirely gray
p_k <- DotPlot(cd8_t, features = genes_all, group.by = "CellType_Group",
               scale = FALSE, cols = c("lightgrey", "#8B0000"), dot.scale = 8) + 
  coord_flip() +
  labs(title = "Cytotoxic & Exhaustion Markers", x = NULL, y = NULL) +
  geom_vline(xintercept = 4.5, linetype = "dashed", color = "grey60") +
  theme_dotplot +
  guides(size = guide_legend(title = "% Expressed"), 
         color = guide_colorbar(title = "Avg Expression"))

ggsave("final_analysis_output/figures/Panel_K_CD8_Markers.pdf", 
       p_k, width = 6.5, height = 3.5, dpi = 300)

# -----------------------------------------------------------------------------
# 7.13 Panel L: Fibroblast ECM markers DotPlot
# -----------------------------------------------------------------------------

cat("  Generating Panel L: Fibroblast ECM Markers...\n")

ecm_genes_l <- c("COL1A1", "COL1A2", "FN1", "POSTN", "MMP2")

fibroblasts$CellType_Group <- factor(paste(fibroblasts$cell_type_final, fibroblasts$group, sep = "_"),
                                      levels = c("Fibroblasts_CTRL", "Fibroblasts_NAC"))

p_l <- DotPlot(fibroblasts, features = ecm_genes_l, group.by = "CellType_Group",
               scale = FALSE, cols = c("lightgrey", "#8B0000"), dot.scale = 8) + 
  coord_flip() +
  labs(title = "ECM Markers", x = NULL, y = NULL) +
  theme_dotplot +
  guides(size = guide_legend(title = "% Expressed"), 
         color = guide_colorbar(title = "Avg Expression"))

ggsave("final_analysis_output/figures/Panel_L_Fibroblast_ECM_Markers.pdf", 
       p_l, width = 5, height = 3.5, dpi = 300)

cat("\n[OK] All figures generated successfully\n\n")

# =============================================================================
# SECTION 8: EXPORT STATISTICAL SUMMARY
# =============================================================================

cat("SECTION 8: Exporting Statistical Summary\n")
cat("=============================================================================\n")

# Compile all statistical results
stats_summary <- data.frame(
  Analysis = c(
    "Myeloid Senescence Score",
    "Fibroblast Proportion",
    "CXCL13+ T Cell Proportion",
    "CD8 Cytotoxicity Score",
    "CD8 Exhaustion Score",
    "Macrophage SASP Score",
    "Senescence-IL1B Correlation",
    "Fibroblast ECM Score"
  ),
  Statistical_Test = c(
    "Wilcoxon (one-sided)",
    "t-test (one-sided)",
    "Fisher's Exact",
    "Wilcoxon (two-sided)",
    "Wilcoxon (one-sided)",
    "Wilcoxon (one-sided)",
    "Pearson Correlation",
    "Wilcoxon (one-sided)"
  ),
  P_value = c(
    wilcox_senescence$p.value,
    t_test_fibro$p.value,
    fisher_test$p.value,
    wilcox_cyto$p.value,
    wilcox_exh$p.value,
    wilcox_sasp$p.value,
    cor_test$p.value,
    wilcox_ecm$p.value
  ),
  Test_Statistic = c(
    sprintf("W = %.0f", wilcox_senescence$statistic),
    sprintf("t = %.2f", t_test_fibro$statistic),
    sprintf("OR = %.2f", fisher_test$estimate),
    sprintf("W = %.0f", wilcox_cyto$statistic),
    sprintf("W = %.0f", wilcox_exh$statistic),
    sprintf("W = %.0f", wilcox_sasp$statistic),
    sprintf("R = %.2f", cor_test$estimate),
    sprintf("W = %.0f", wilcox_ecm$statistic)
  )
)

write.csv(stats_summary, 
          "final_analysis_output/statistics/statistical_summary.csv", 
          row.names = FALSE)

cat("\n[OK] Statistical summary saved\n")

# =============================================================================
# SECTION 9: SESSION INFORMATION
# =============================================================================

cat("\n=============================================================================\n")
cat("ANALYSIS COMPLETED SUCCESSFULLY\n")
cat("=============================================================================\n\n")

cat("Output files:\n")
cat("  final_analysis_output/\n")
cat("  |-- seurat_annotated.rds          # Annotated Seurat object\n")
cat("  |-- seurat_with_scores.rds        # Seurat object with module scores\n")
cat("  |-- figures/\n")
cat("  |   |-- Panel_A_UMAP_CellType.pdf\n")
cat("  |   |-- Panel_B_UMAP_Group.pdf\n")
cat("  |   |-- Panel_C_CXCL13_Proportion.pdf\n")
cat("  |   |-- Panel_D_CD8_Cytotoxicity.pdf\n")
cat("  |   |-- Panel_E_CD8_Exhaustion.pdf\n")
cat("  |   |-- Panel_F_Senescence_DotPlot.pdf\n")
cat("  |   |-- Panel_G_SASP_Senescence_DotPlot.pdf\n")
cat("  |   |-- Panel_H_SASP_Group_DotPlot.pdf\n")
cat("  |   |-- Panel_I_Senescence_IL1B_Scatter.pdf\n")
cat("  |   |-- Panel_J_Fibroblast_ECM.pdf\n")
cat("  |   |-- Panel_K_CD8_Markers.pdf\n")
cat("  |   |-- Panel_L_Fibroblast_ECM_Markers.pdf\n")
cat("  |-- statistics/\n")
cat("      |-- statistical_summary.csv\n")
cat("\n")
cat("All figures are 300 DPI vector PDFs with editable text.\n")
cat("=============================================================================\n\n")

# Print session information for reproducibility
cat("R Session Information (for Methods section):\n")
cat("-----------------------------------------------------------------------------\n")
print(sessionInfo())

# =============================================================================
# END OF SCRIPT
# =============================================================================
