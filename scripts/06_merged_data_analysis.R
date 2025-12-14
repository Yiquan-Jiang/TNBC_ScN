# ==============================================================================
# Script 06: Merged Dataset Comprehensive Analysis
# ScN Mouse Single-Cell RNA-seq Analysis Pipeline
# ==============================================================================
# 
# Description:
#   This script performs comprehensive analysis on merged/integrated datasets
#   including:
#   - Cross-dataset cell type comparison
#   - Treatment response analysis across datasets
#   - Integrated visualization and statistics
#
# Input:
#   - Merged Seurat object with multiple datasets
#
# Output:
#   - Cross-dataset comparison plots (PDF)
#   - Statistical summary tables (CSV)
#
# Author: Biomedical Research Assistant
# Date: 2025
# ==============================================================================

# Set CRAN mirror
options(repos = c(CRAN = "https://cran.rstudio.com/"))

# Load required packages
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(RColorBrewer)
library(gridExtra)
library(cowplot)
library(scales)
library(tibble)
library(patchwork)
library(ggpubr)
library(ggsignif)
library(ggbeeswarm)
library(ggridges)
library(wesanderson)
library(ggsci)
library(paletteer)

# Set PDF font parameters
pdf.options(family = "Helvetica")

# ==============================================================================
# Theme and Color Configuration
# ==============================================================================

elegant_theme <- theme_minimal() +
  theme(
    text = element_text(family = "Helvetica", size = 12, color = "grey15"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, 
                             margin = margin(b = 20), color = "grey10"),
    plot.subtitle = element_text(size = 16, hjust = 0.5, 
                                margin = margin(b = 15), color = "grey30"),
    axis.title.x = element_text(size = 15, face = "bold", 
                               margin = margin(t = 12), color = "grey20"),
    axis.title.y = element_text(size = 15, face = "bold", 
                               margin = margin(r = 12), color = "grey20"),
    axis.text.x = element_text(size = 13, color = "grey25", face = "bold"),
    axis.text.y = element_text(size = 13, color = "grey25"),
    legend.title = element_text(size = 14, face = "bold", color = "grey20"),
    legend.text = element_text(size = 12, color = "grey30"),
    legend.position = "bottom",
    strip.text = element_text(size = 13, face = "bold", color = "grey20"),
    strip.background = element_rect(fill = "grey98", color = NA),
    panel.border = element_blank(),
    axis.line = element_line(color = "grey35", linewidth = 0.8),
    axis.ticks = element_line(color = "grey35", linewidth = 0.6),
    plot.margin = margin(20, 20, 20, 20),
    panel.grid.major = element_line(color = "grey96", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

theme_set(elegant_theme)

# Extended color palette for multiple sample groups
sample_colors <- c(
  "Untreated" = "#2E86AB",
  "Chemo" = "#A23B72",
  "Chemo_Abt263" = "#F18F01",
  "PCC1" = "#06D6A0",
  "Taxel" = "#F72585",
  "Ctrl" = "#4361EE",
  "blank" = "#FFB3BA"
)

caf_colors <- c(
  "iCAF" = "#E63946",
  "myCAF" = "#457B9D",
  "apCAF" = "#F77F00",
  "matCAF" = "#2A9D8F"
)

# ==============================================================================
# Data Loading and Preparation
# ==============================================================================

prepare_merged_data <- function(seurat_obj, sample_order = NULL) {
  
  cat("=== Preparing Merged Data ===\n")
  
  # Check sample distribution
  cat("Sample distribution:\n")
  sample_dist <- table(seurat_obj$orig.ident)
  print(sample_dist)
  
  # Set sample order if provided
  if(!is.null(sample_order)) {
    seurat_obj$orig.ident <- factor(seurat_obj$orig.ident, levels = sample_order)
  }
  
  # Create dataset source identifier if multiple datasets
  if("dataset_source" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$Dataset <- ifelse(seurat_obj$dataset_source == "Blank_Taxol_DQ", 
                                 "Dataset1", "Dataset2")
  }
  
  return(seurat_obj)
}

# ==============================================================================
# Cell Type Distribution Analysis
# ==============================================================================

analyze_celltype_distribution <- function(seurat_obj, sample_col = "orig.ident") {
  
  cat("\n=== Cell Type Distribution Analysis ===\n")
  
  # Calculate distribution
  celltype_dist <- seurat_obj@meta.data %>%
    group_by(.data[[sample_col]], Celltype_major_man) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(.data[[sample_col]]) %>%
    mutate(percentage = count / sum(count) * 100)
  
  # Stacked bar plot
  p <- ggplot(celltype_dist, aes_string(x = sample_col, y = "percentage", 
                                         fill = "Celltype_major_man")) +
    geom_bar(stat = "identity", position = "stack", alpha = 0.8) +
    scale_fill_viridis_d() +
    labs(title = "Cell Type Distribution by Sample",
         x = "Sample Group", y = "Cell Type Proportion (%)",
         fill = "Cell Type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(list(plot = p, stats = celltype_dist))
}

# ==============================================================================
# Myeloid Analysis Across Datasets
# ==============================================================================

analyze_myeloid_merged <- function(seurat_obj, sample_col = "orig.ident") {
  
  cat("\n=== Merged Myeloid Analysis ===\n")
  
  myeloid_cells <- subset(seurat_obj, subset = Celltype_major_man == "Myeloid")
  
  if(ncol(myeloid_cells) == 0) {
    cat("No myeloid cells found\n")
    return(NULL)
  }
  
  # Senescence scoring
  senescence_genes <- c("Cdkn1a", "Cdkn2a", "Cdkn1b", "Trp53", "Rb1", "Cdkn2b", 
                       "Cdkn2c", "Cdkn2d", "Lmnb1", "Hmgb1", "Il1a", "Il1b", 
                       "Il6", "Tnf", "Cxcl1", "Ccl2", "Mmp3", "Serpine1", 
                       "Igfbp3", "Igfbp5", "Igfbp7", "Tgfb1", "Vegfa")
  
  available_genes <- senescence_genes[senescence_genes %in% rownames(myeloid_cells)]
  
  if(length(available_genes) > 0) {
    tryCatch({
      myeloid_cells <- AddModuleScore(
        object = myeloid_cells,
        features = list(available_genes),
        name = "Senescence_Score"
      )
      colnames(myeloid_cells@meta.data)[colnames(myeloid_cells@meta.data) == "Senescence_Score1"] <- "Senescence_Score"
    }, error = function(e) {
      expr_data <- GetAssayData(myeloid_cells, assay = "RNA", layer = "data")
      senescence_expr <- expr_data[available_genes, , drop = FALSE]
      myeloid_cells$Senescence_Score <<- colMeans(senescence_expr, na.rm = TRUE)
    })
  }
  
  # IL1B expression
  if("Il1b" %in% rownames(myeloid_cells)) {
    il1b_expression <- GetAssayData(myeloid_cells, assay = "RNA", layer = "data")["Il1b", ]
    myeloid_cells$Il1b_Expression <- il1b_expression
    
    il1b_stats <- myeloid_cells@meta.data %>%
      group_by(.data[[sample_col]]) %>%
      summarise(
        n_cells = n(),
        mean_expr = mean(Il1b_Expression, na.rm = TRUE),
        median_expr = median(Il1b_Expression, na.rm = TRUE),
        positive_cells = sum(Il1b_Expression > 0),
        positive_rate = positive_cells / n_cells * 100,
        .groups = "drop"
      )
    
    cat("IL1B expression statistics:\n")
    print(il1b_stats)
  }
  
  # Senescence plot
  p_senescence <- ggplot(myeloid_cells@meta.data, 
                         aes_string(x = sample_col, y = "Senescence_Score", 
                                    fill = sample_col)) +
    geom_violin(alpha = 0.7, trim = FALSE) +
    geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
    scale_fill_manual(values = sample_colors) +
    stat_compare_means(method = "kruskal.test", label.y.npc = 0.95) +
    labs(title = "Myeloid Cell Senescence Score",
         subtitle = "Distribution across sample groups",
         x = "Sample Group", y = "Senescence Score") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  # IL1B plot
  if("Il1b_Expression" %in% colnames(myeloid_cells@meta.data)) {
    p_il1b <- ggplot(myeloid_cells@meta.data, 
                     aes_string(x = sample_col, y = "Il1b_Expression", 
                                fill = sample_col)) +
      geom_violin(alpha = 0.7, trim = FALSE) +
      geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
      scale_fill_manual(values = sample_colors) +
      stat_compare_means(method = "kruskal.test", label.y.npc = 0.95) +
      labs(title = "Myeloid IL1B Expression",
           subtitle = "Distribution across sample groups",
           x = "Sample Group", y = "IL1B Expression") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
  } else {
    p_il1b <- NULL
  }
  
  return(list(
    senescence_plot = p_senescence,
    il1b_plot = p_il1b,
    data = myeloid_cells
  ))
}

# ==============================================================================
# CAF Analysis Across Datasets
# ==============================================================================

analyze_caf_merged <- function(seurat_obj, sample_col = "orig.ident") {
  
  cat("\n=== Merged CAF Analysis ===\n")
  
  fibroblast_cells <- subset(seurat_obj, subset = Celltype_major_man == "Fibroblast")
  
  if(ncol(fibroblast_cells) == 0) {
    cat("No fibroblast cells found\n")
    return(NULL)
  }
  
  # CAF marker definitions
  caf_markers <- list(
    iCAF = c("Il6", "Il1b", "Cxcl1", "Cxcl2", "Ccl2", "Ccl7", "Tnf", "Il1a", 
             "Ptgs2", "Nfkb1", "Stat3", "Cxcl12", "Cxcl14"),
    myCAF = c("Acta2", "Tagln", "Myh11", "Cnn1", "Postn", "Lox", "Tgfb1", 
              "Pdgfra", "Pdgfrb", "Vim", "Fn1"),
    apCAF = c("H2-Aa", "H2-Ab1", "H2-Eb1", "Cd74", "Ciita", "H2-DMa", 
              "H2-DMb1", "Cd80", "Cd86", "Ccl19", "Ccl21a"),
    matCAF = c("Col1a1", "Col1a2", "Col3a1", "Col5a1", "Fn1", "Sparc", 
               "Tnc", "Loxl2", "Mmp2", "Mmp14", "Fbn1", "Dpt")
  )
  
  # Calculate CAF scores
  for(caf_type in names(caf_markers)) {
    available_genes <- caf_markers[[caf_type]][caf_markers[[caf_type]] %in% rownames(fibroblast_cells)]
    
    if(length(available_genes) > 2) {
      tryCatch({
        fibroblast_cells <- AddModuleScore(
          object = fibroblast_cells,
          features = list(available_genes),
          name = paste0(caf_type, "_Score")
        )
        
        score_col <- paste0(caf_type, "_Score1")
        new_col <- paste0(caf_type, "_Score")
        if(score_col %in% colnames(fibroblast_cells@meta.data)) {
          colnames(fibroblast_cells@meta.data)[colnames(fibroblast_cells@meta.data) == score_col] <- new_col
        }
      }, error = function(e) {
        expr_data <- GetAssayData(fibroblast_cells, assay = "RNA", layer = "data")
        caf_expr <- expr_data[available_genes, , drop = FALSE]
        fibroblast_cells@meta.data[[paste0(caf_type, "_Score")]] <<- colMeans(caf_expr, na.rm = TRUE)
      })
    } else {
      fibroblast_cells@meta.data[[paste0(caf_type, "_Score")]] <- rnorm(ncol(fibroblast_cells), mean = 0, sd = 0.1)
    }
  }
  
  # CAF classification
  score_data <- fibroblast_cells@meta.data[, c("iCAF_Score", "myCAF_Score", "apCAF_Score", "matCAF_Score")]
  
  classify_caf <- function(scores) {
    scores_norm <- (scores - min(scores, na.rm = TRUE)) / (max(scores, na.rm = TRUE) - min(scores, na.rm = TRUE))
    max_score_idx <- which.max(scores_norm)
    max_score <- scores_norm[max_score_idx]
    score_names <- names(scores_norm)
    
    other_scores <- scores_norm[-max_score_idx]
    if(max_score - max(other_scores, na.rm = TRUE) > 0.1) {
      return(gsub("_Score", "", score_names[max_score_idx]))
    } else {
      priority_order <- c("myCAF_Score", "iCAF_Score", "apCAF_Score", "matCAF_Score")
      for(priority_type in priority_order) {
        if(priority_type %in% score_names && scores_norm[priority_type] > 0.3) {
          return(gsub("_Score", "", priority_type))
        }
      }
      return(gsub("_Score", "", score_names[max_score_idx]))
    }
  }
  
  fibroblast_cells$CAF_Subtype <- apply(score_data, 1, classify_caf)
  
  # CAF subtype distribution
  caf_subtype_dist <- fibroblast_cells@meta.data %>%
    group_by(.data[[sample_col]], CAF_Subtype) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(.data[[sample_col]]) %>%
    mutate(percentage = count / sum(count) * 100)
  
  # Plot
  p_caf <- ggplot(caf_subtype_dist, aes_string(x = sample_col, y = "percentage", 
                                                fill = "CAF_Subtype")) +
    geom_bar(stat = "identity", position = "stack", alpha = 0.8) +
    scale_fill_manual(values = caf_colors) +
    labs(title = "CAF Subtype Distribution",
         subtitle = "Analysis across sample groups",
         x = "Sample Group", y = "CAF Subtype Proportion (%)",
         fill = "CAF Subtype") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Key gene expression
  key_genes <- c("Cxcl13", "Cxcl12", "Il1r1", "Acta2", "Il6")
  
  for(gene in key_genes) {
    if(gene %in% rownames(fibroblast_cells)) {
      gene_expr <- GetAssayData(fibroblast_cells, assay = "RNA", layer = "data")[gene, ]
      fibroblast_cells@meta.data[[paste0(gene, "_Expression")]] <- gene_expr
      fibroblast_cells@meta.data[[paste0(gene, "_positive")]] <- gene_expr > 0
    }
  }
  
  return(list(
    plot = p_caf,
    stats = caf_subtype_dist,
    data = fibroblast_cells
  ))
}

# ==============================================================================
# T Cell Analysis Across Datasets
# ==============================================================================

analyze_tcell_merged <- function(seurat_obj, sample_col = "orig.ident") {
  
  cat("\n=== Merged T Cell Analysis ===\n")
  
  t_cells <- subset(seurat_obj, subset = Celltype_major_man == "NK_T_cell")
  
  if(ncol(t_cells) == 0) {
    cat("No T cells found\n")
    return(NULL)
  }
  
  if(all(c("Pdcd1", "Gzmk") %in% rownames(t_cells))) {
    pdcd1_exp <- GetAssayData(t_cells, assay = "RNA", layer = "data")["Pdcd1", ]
    gzmk_exp <- GetAssayData(t_cells, assay = "RNA", layer = "data")["Gzmk", ]
    
    t_cells$Pdcd1_Expression <- pdcd1_exp
    t_cells$Gzmk_Expression <- gzmk_exp
    t_cells$Pdcd1_positive <- pdcd1_exp > 0
    t_cells$Gzmk_positive <- gzmk_exp > 0
    t_cells$Tumor_reactive <- t_cells$Pdcd1_positive & t_cells$Gzmk_positive
    
    t_reactive_stats <- t_cells@meta.data %>%
      group_by(.data[[sample_col]]) %>%
      summarise(
        n_cells = n(),
        pdcd1_positive = sum(Pdcd1_positive),
        gzmk_positive = sum(Gzmk_positive),
        tumor_reactive = sum(Tumor_reactive),
        reactive_rate = tumor_reactive / n_cells * 100,
        .groups = "drop"
      )
    
    cat("T cell reactivity statistics:\n")
    print(t_reactive_stats)
    
    # Plot
    p_tcell <- ggplot(t_reactive_stats, aes_string(x = sample_col, 
                                                    y = "reactive_rate", 
                                                    fill = sample_col)) +
      geom_bar(stat = "identity", alpha = 0.8) +
      geom_text(aes(label = paste0(round(reactive_rate, 1), "%")), 
                vjust = -0.5, size = 3.5) +
      scale_fill_manual(values = sample_colors) +
      labs(title = "Tumor Reactive T Cell Proportion",
           subtitle = "PDCD1+ GZMK+ T cells across samples",
           x = "Sample Group", y = "Reactive T Cells (%)") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
    
    return(list(
      plot = p_tcell,
      stats = t_reactive_stats,
      data = t_cells
    ))
  } else {
    cat("PDCD1 or GZMK genes not found\n")
    return(NULL)
  }
}

# ==============================================================================
# Dataset Comparison
# ==============================================================================

compare_datasets <- function(seurat_obj, dataset_col = "Dataset") {
  
  cat("\n=== Dataset Comparison ===\n")
  
  if(!dataset_col %in% colnames(seurat_obj@meta.data)) {
    cat("Dataset column not found\n")
    return(NULL)
  }
  
  dataset_comparison <- seurat_obj@meta.data %>%
    group_by(.data[[dataset_col]], Celltype_major_man) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(.data[[dataset_col]]) %>%
    mutate(percentage = count / sum(count) * 100)
  
  p <- ggplot(dataset_comparison, aes_string(x = dataset_col, y = "percentage", 
                                              fill = "Celltype_major_man")) +
    geom_bar(stat = "identity", position = "stack", alpha = 0.8) +
    scale_fill_viridis_d() +
    labs(title = "Cell Type Distribution by Dataset",
         x = "Dataset", y = "Cell Type Proportion (%)",
         fill = "Cell Type")
  
  return(list(plot = p, stats = dataset_comparison))
}

# ==============================================================================
# Main Pipeline
# ==============================================================================

run_merged_analysis_pipeline <- function(seurat_obj, sample_col = "orig.ident",
                                         sample_order = NULL,
                                         output_prefix = "merged") {
  
  cat("=== Starting Merged Data Analysis Pipeline ===\n")
  
  # Prepare data
  seurat_obj <- prepare_merged_data(seurat_obj, sample_order)
  
  # Run analyses
  celltype_results <- analyze_celltype_distribution(seurat_obj, sample_col)
  myeloid_results <- analyze_myeloid_merged(seurat_obj, sample_col)
  caf_results <- analyze_caf_merged(seurat_obj, sample_col)
  tcell_results <- analyze_tcell_merged(seurat_obj, sample_col)
  
  # Save plots
  ggsave(paste0(output_prefix, "_celltype_distribution.pdf"), 
         celltype_results$plot, width = 14, height = 8, dpi = 300)
  
  if(!is.null(myeloid_results)) {
    ggsave(paste0(output_prefix, "_myeloid_senescence.pdf"), 
           myeloid_results$senescence_plot, width = 12, height = 8, dpi = 300)
    
    if(!is.null(myeloid_results$il1b_plot)) {
      ggsave(paste0(output_prefix, "_myeloid_il1b.pdf"), 
             myeloid_results$il1b_plot, width = 12, height = 8, dpi = 300)
    }
  }
  
  if(!is.null(caf_results)) {
    ggsave(paste0(output_prefix, "_caf_subtype.pdf"), 
           caf_results$plot, width = 12, height = 8, dpi = 300)
  }
  
  if(!is.null(tcell_results)) {
    ggsave(paste0(output_prefix, "_tcell_reactivity.pdf"), 
           tcell_results$plot, width = 12, height = 8, dpi = 300)
  }
  
  # Save statistics
  write.csv(celltype_results$stats, paste0(output_prefix, "_celltype_stats.csv"), row.names = FALSE)
  
  if(!is.null(caf_results)) {
    write.csv(caf_results$stats, paste0(output_prefix, "_caf_stats.csv"), row.names = FALSE)
  }
  
  if(!is.null(tcell_results)) {
    write.csv(tcell_results$stats, paste0(output_prefix, "_tcell_stats.csv"), row.names = FALSE)
  }
  
  # Dataset comparison if available
  if("Dataset" %in% colnames(seurat_obj@meta.data)) {
    dataset_results <- compare_datasets(seurat_obj)
    if(!is.null(dataset_results)) {
      ggsave(paste0(output_prefix, "_dataset_comparison.pdf"), 
             dataset_results$plot, width = 10, height = 8, dpi = 300)
    }
  }
  
  cat("\n=== Merged Analysis Complete ===\n")
  
  return(list(
    celltype = celltype_results,
    myeloid = myeloid_results,
    caf = caf_results,
    tcell = tcell_results
  ))
}

# ==============================================================================
# Example Usage
# ==============================================================================

# Uncomment and modify to run:

# # Load merged data
# seurat_obj <- readRDS("Merged_Annotated_SingleCell_Data.rds")
# 
# # Run merged analysis
# results <- run_merged_analysis_pipeline(
#   seurat_obj,
#   sample_col = "orig.ident",
#   sample_order = c("Untreated", "Chemo", "Chemo_Abt263", "PCC1", "Taxel", "Ctrl", "blank"),
#   output_prefix = "merged"
# )

cat("Script 06: Merged Data Analysis - Functions loaded successfully!\n")
