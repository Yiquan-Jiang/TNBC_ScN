# ==============================================================================
# Script 01: Data Processing and CAF Subtyping
# ScN Mouse Single-Cell RNA-seq Analysis Pipeline
# ==============================================================================
# 
# Description:
#   This script performs initial data loading, sample filtering, treatment group
#   assignment, UMAP recalculation, and CAF subtype classification using marker
#   gene module scores.
#
# Input:
#   - Seurat object file (e.g., "your_seurat_object.rds")
#
# Output:
#   - Processed Seurat object with CAF subtypes
#   - UMAP visualization plots
#
# Author: Biomedical Research Assistant
# Date: 2025
# ==============================================================================

# Load required packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(patchwork)
library(wesanderson)
library(ggsci)

# Set high-quality plot parameters
options(repr.plot.width = 12, repr.plot.height = 8)

# ==============================================================================
# Color Scheme Configuration
# ==============================================================================

elegant_colors <- list(
  treatment = c("Untreated" = "#2E86AB", "Chemo" = "#A23B72", 
                "Chemo_Abt263" = "#F18F01", "Chemo_Abt263_IL1B" = "#457B9D"),
  
  caf_subtype = c(
    "iCAF" = "#E63946",      # Inflammatory CAF - Red
    "myCAF" = "#457B9D",     # Myofibroblastic CAF - Blue  
    "apCAF" = "#F77F00",     # Antigen-presenting CAF - Orange
    "matCAF" = "#2A9D8F"     # Matrix CAF - Green
  ),
  
  cell_types = c(
    "Fibroblast" = "#E31A1C",
    "NK_T_cell" = "#1F78B4",
    "Tumor_cell" = "#33A02C",
    "Myeloid" = "#FF7F00",
    "B_Plasma" = "#6A3D9A",
    "Neutrophil" = "#FB9A99",
    "Endothelial" = "#A6CEE3"
  )
)

# Journal-quality theme
journal_theme <- theme_classic() +
  theme(
    text = element_text(family = "sans", size = 15),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 15, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.line = element_line(colour = "black", linewidth = 0.5),
    plot.margin = margin(20, 20, 20, 20)
  )

# ==============================================================================
# Data Loading
# ==============================================================================

cat("Loading data...\n")

# Load your Seurat object - modify path as needed
# seurat_obj <- readRDS("your_seurat_object.rds")

# Example: Check data structure
# cat("Dataset information:\n")
# cat("Number of cells:", ncol(seurat_obj), "\n")
# cat("Number of genes:", nrow(seurat_obj), "\n")
# cat("Treatment groups:", unique(seurat_obj$orig.ident), "\n")
# cat("Cell types:", unique(seurat_obj$Celltype_major_man), "\n")

# ==============================================================================
# Sample Filtering and Treatment Group Assignment
# ==============================================================================

# Define treatment groups to keep
# keep_groups <- c("Untreated", "Chemo", "Chemo_Abt263", "Chemo_Abt263_IL1B")
# seurat_filtered <- subset(seurat_obj, subset = orig.ident %in% keep_groups)

# Create Treatment column
# seurat_filtered$Treatment <- seurat_filtered$orig.ident

# ==============================================================================
# UMAP Recalculation
# ==============================================================================

recalculate_umap <- function(seurat_obj, dims = 30, verbose = FALSE) {
  cat("Recalculating UMAP for all cells...\n")
  
  # Normalize data
  seurat_obj <- NormalizeData(seurat_obj, verbose = verbose)
  
  # Find variable features
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", 
                                     nfeatures = 2000, verbose = verbose)
  
  # Scale data
  seurat_obj <- ScaleData(seurat_obj, verbose = verbose)
  
  # PCA
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj), 
                       verbose = verbose)
  
  # Find neighbors
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:dims, verbose = verbose)
  
  # UMAP
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:dims, verbose = verbose)
  
  cat("UMAP calculation completed\n")
  return(seurat_obj)
}

# ==============================================================================
# CAF Subtype Classification
# ==============================================================================

# CAF marker gene sets
caf_markers_refined <- list(
  # Inflammatory CAF - Pro-inflammatory cytokines and chemokines
  iCAF = c("Il6", "Il1b", "Cxcl1", "Cxcl2", "Ccl2", "Ccl7", "Tnf", "Il1a", 
           "Ptgs2", "Nfkb1", "Stat3", "Cxcl12", "Cxcl14"),
  
  # Myofibroblastic CAF - Contractile and ECM remodeling
  myCAF = c("Acta2", "Tagln", "Myh11", "Cnn1", "Postn", "Lox", "Tgfb1", 
            "Pdgfra", "Pdgfrb", "Vim", "Fn1"),
  
  # Antigen-presenting CAF - MHC class II and co-stimulatory molecules
  apCAF = c("H2-Aa", "H2-Ab1", "H2-Eb1", "Cd74", "Ciita", "H2-DMa", 
            "H2-DMb1", "Cd80", "Cd86", "Ccl19", "Ccl21a"),
  
  # Matrix CAF - Collagens and ECM components
  matCAF = c("Col1a1", "Col1a2", "Col3a1", "Col5a1", "Fn1", "Sparc", 
             "Tnc", "Loxl2", "Mmp2", "Mmp14", "Fbn1", "Dpt")
)

# Function to calculate CAF module scores
calculate_caf_scores <- function(fibroblast_cells, caf_markers) {
  for(caf_type in names(caf_markers)) {
    available_genes <- caf_markers[[caf_type]][caf_markers[[caf_type]] %in% rownames(fibroblast_cells)]
    cat("Available", caf_type, "genes:", length(available_genes), "out of", 
        length(caf_markers[[caf_type]]), "\n")
    
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
        cat("Error in AddModuleScore for", caf_type, ":", e$message, "\n")
        expr_data <- GetAssayData(fibroblast_cells, assay = "RNA", layer = "data")
        caf_expr <- expr_data[available_genes, , drop = FALSE]
        fibroblast_cells@meta.data[[paste0(caf_type, "_Score")]] <<- colMeans(caf_expr, na.rm = TRUE)
      })
    } else {
      fibroblast_cells@meta.data[[paste0(caf_type, "_Score")]] <- rnorm(ncol(fibroblast_cells), mean = 0, sd = 0.1)
      cat("Warning: Too few", caf_type, "genes, created dummy scores\n")
    }
  }
  return(fibroblast_cells)
}

# CAF classification function
classify_caf_refined <- function(scores) {
  # Normalize scores
  scores_norm <- (scores - min(scores, na.rm = TRUE)) / 
                 (max(scores, na.rm = TRUE) - min(scores, na.rm = TRUE))
  
  # Find maximum score
  max_score_idx <- which.max(scores_norm)
  max_score <- scores_norm[max_score_idx]
  score_names <- names(scores_norm)
  
  # If max score is clearly higher (difference > 0.1), classify as that subtype
  other_scores <- scores_norm[-max_score_idx]
  if(max_score - max(other_scores, na.rm = TRUE) > 0.1) {
    return(gsub("_Score", "", score_names[max_score_idx]))
  } else {
    # If scores are close, choose based on biological priority
    # Priority: myCAF > iCAF > apCAF > matCAF
    priority_order <- c("myCAF_Score", "iCAF_Score", "apCAF_Score", "matCAF_Score")
    for(priority_type in priority_order) {
      if(priority_type %in% score_names && scores_norm[priority_type] > 0.3) {
        return(gsub("_Score", "", priority_type))
      }
    }
    return(gsub("_Score", "", score_names[max_score_idx]))
  }
}

# ==============================================================================
# Main Processing Pipeline (Example)
# ==============================================================================

run_processing_pipeline <- function(seurat_obj, treatment_groups = NULL) {
  
  # Filter samples if specified
  if(!is.null(treatment_groups)) {
    seurat_obj <- subset(seurat_obj, subset = orig.ident %in% treatment_groups)
    seurat_obj$Treatment <- seurat_obj$orig.ident
  }
  
  # Recalculate UMAP
  seurat_obj <- recalculate_umap(seurat_obj, dims = 30)
  
  # Extract fibroblasts
  fibroblast_cells <- subset(seurat_obj, subset = Celltype_major_man == "Fibroblast")
  
  if(ncol(fibroblast_cells) > 0) {
    cat("Number of Fibroblast cells:", ncol(fibroblast_cells), "\n")
    
    # Calculate CAF scores
    fibroblast_cells <- calculate_caf_scores(fibroblast_cells, caf_markers_refined)
    
    # Classify CAF subtypes
    score_data <- fibroblast_cells@meta.data[, c("iCAF_Score", "myCAF_Score", "apCAF_Score", "matCAF_Score")]
    fibroblast_cells$CAF_Subtype <- apply(score_data, 1, classify_caf_refined)
    
    cat("CAF subtype distribution:\n")
    print(table(fibroblast_cells$CAF_Subtype))
    
    # Recalculate UMAP for fibroblasts
    fibroblast_cells <- recalculate_umap(fibroblast_cells, dims = 20)
    
    # Add CAF subtype to main object
    seurat_obj$CAF_Subtype <- NA
    seurat_obj$CAF_Subtype[colnames(seurat_obj) %in% colnames(fibroblast_cells)] <- 
      fibroblast_cells$CAF_Subtype[match(colnames(seurat_obj)[colnames(seurat_obj) %in% colnames(fibroblast_cells)], 
                                         colnames(fibroblast_cells))]
    seurat_obj$CAF_Subtype[is.na(seurat_obj$CAF_Subtype)] <- "Non-Fibroblast"
  }
  
  return(list(seurat_obj = seurat_obj, fibroblast_cells = fibroblast_cells))
}

# ==============================================================================
# UMAP Plotting Functions
# ==============================================================================

plot_umap_celltypes <- function(seurat_obj, colors = elegant_colors$cell_types) {
  available_cell_types <- unique(seurat_obj$Celltype_major_man)
  available_cell_types <- available_cell_types[!is.na(available_cell_types)]
  
  celltype_colors_final <- colors[available_cell_types]
  missing_types <- available_cell_types[!available_cell_types %in% names(colors)]
  if(length(missing_types) > 0) {
    additional_colors <- rainbow(length(missing_types), alpha = 0.8)
    names(additional_colors) <- missing_types
    celltype_colors_final <- c(celltype_colors_final, additional_colors)
  }
  
  p <- DimPlot(seurat_obj, 
               reduction = "umap", 
               group.by = "Celltype_major_man",
               cols = celltype_colors_final,
               pt.size = 0.5,
               raster = TRUE) +
    ggtitle("UMAP: All Cell Types") +
    xlab("UMAP 1") + ylab("UMAP 2") +
    journal_theme +
    guides(color = guide_legend(title = "Cell Type", 
                                override.aes = list(size = 4), ncol = 1))
  
  return(p)
}

plot_umap_treatment <- function(seurat_obj, colors = elegant_colors$treatment) {
  p <- DimPlot(seurat_obj, 
               reduction = "umap", 
               group.by = "Treatment",
               cols = colors,
               pt.size = 0.5,
               raster = TRUE) +
    ggtitle("UMAP: Treatment Groups") +
    xlab("UMAP 1") + ylab("UMAP 2") +
    journal_theme +
    guides(color = guide_legend(title = "Treatment", 
                                override.aes = list(size = 4), ncol = 1))
  
  return(p)
}

plot_umap_caf_subtypes <- function(fibroblast_cells, colors = elegant_colors$caf_subtype) {
  available_subtypes <- unique(fibroblast_cells$CAF_Subtype)
  available_subtypes <- available_subtypes[!is.na(available_subtypes)]
  
  caf_colors_final <- colors[available_subtypes]
  missing_subtypes <- available_subtypes[!available_subtypes %in% names(colors)]
  if(length(missing_subtypes) > 0) {
    additional_colors <- viridis::viridis(length(missing_subtypes), alpha = 0.8, option = "plasma")
    names(additional_colors) <- missing_subtypes
    caf_colors_final <- c(caf_colors_final, additional_colors)
  }
  
  p <- DimPlot(fibroblast_cells, 
               reduction = "umap", 
               group.by = "CAF_Subtype",
               cols = caf_colors_final,
               pt.size = 1,
               raster = TRUE) +
    ggtitle("UMAP: CAF Subtypes") +
    xlab("UMAP 1") + ylab("UMAP 2") +
    journal_theme +
    guides(color = guide_legend(title = "CAF Subtype", 
                                override.aes = list(size = 4), ncol = 1))
  
  return(p)
}

# ==============================================================================
# Example Usage
# ==============================================================================

# Uncomment and modify the following lines to run the pipeline:

# # Load data
# seurat_obj <- readRDS("your_seurat_object.rds")
# 
# # Run processing pipeline
# results <- run_processing_pipeline(seurat_obj, 
#                                    treatment_groups = c("Untreated", "Chemo", "Chemo_Abt263"))
# 
# # Extract processed data
# processed_obj <- results$seurat_obj
# fibroblast_cells <- results$fibroblast_cells
# 
# # Generate plots
# p1 <- plot_umap_celltypes(processed_obj)
# p2 <- plot_umap_treatment(processed_obj)
# p3 <- plot_umap_caf_subtypes(fibroblast_cells)
# 
# # Save plots
# ggsave("umap_all_celltypes.pdf", p1, width = 10, height = 8, dpi = 300)
# ggsave("umap_treatment_groups.pdf", p2, width = 10, height = 8, dpi = 300)
# ggsave("umap_caf_subtypes.pdf", p3, width = 10, height = 8, dpi = 300)
# 
# # Save processed data
# saveRDS(processed_obj, "processed_seurat_object.rds")

cat("Script 01: Data Processing - Functions loaded successfully!\n")
