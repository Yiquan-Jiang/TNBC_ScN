# ==============================================================================
# Script 02: UMAP Visualization using SCP Package
# ScN Mouse Single-Cell RNA-seq Analysis Pipeline
# ==============================================================================
# 
# Description:
#   This script generates high-quality UMAP visualizations using the SCP 
#   (Single Cell Pipeline) package, including cell type distributions,
#   treatment groups, and CAF subtype visualizations.
#
# Input:
#   - Processed Seurat object from Script 01
#
# Output:
#   - Publication-quality UMAP plots (PDF and PNG)
#   - Statistical summary plots
#
# Author: Biomedical Research Assistant
# Date: 2025
# ==============================================================================

# Set renv environment
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}
renv::activate()

# Load required packages
suppressMessages({
  library(SCP)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# Set high-quality plot parameters
options(repr.plot.width = 12, repr.plot.height = 8)

# Set PDF font type for editable format
pdf.options(useDingbats = FALSE)
options(device = function(...) {
  pdf(..., useDingbats = FALSE)
})

cat("=== SCP Package UMAP Visualization ===\n")

# ==============================================================================
# Load Data
# ==============================================================================

# Modify the path to your processed data file
# seurat_obj <- readRDS("processed_seurat_object.rds")

# cat("Dataset information:\n")
# cat("Number of cells:", ncol(seurat_obj), "\n")
# cat("Number of genes:", nrow(seurat_obj), "\n")
# cat("Treatment groups:", paste(unique(seurat_obj$Treatment), collapse = ", "), "\n")
# cat("Cell types:", paste(unique(seurat_obj$Celltype_major_man), collapse = ", "), "\n")

# ==============================================================================
# SCP UMAP Plotting Functions
# ==============================================================================

generate_scp_celldimplot <- function(seurat_obj, group_by, theme_use = "theme_blank", pt.size = 1.2) {
  p <- CellDimPlot(
    srt = seurat_obj,
    group.by = group_by,
    reduction = "UMAP",
    theme_use = theme_use,
    pt.size = pt.size
  )
  return(p)
}

generate_scp_cellstatplot <- function(seurat_obj, stat_by, group_by, theme_use = "theme_blank") {
  p <- CellStatPlot(
    srt = seurat_obj,
    stat.by = stat_by,
    group.by = group_by,
    theme_use = theme_use
  )
  return(p)
}

# ==============================================================================
# Main Visualization Pipeline
# ==============================================================================

run_scp_visualization <- function(seurat_obj, output_prefix = "scp_umap") {
  
  cat("Generating SCP UMAP visualizations...\n")
  
  # Cell type UMAP
  cat("Creating cell type UMAP...\n")
  p1_scp <- generate_scp_celldimplot(seurat_obj, "Celltype_major_man")
  
  # Treatment group UMAP
  cat("Creating treatment group UMAP...\n")
  p2_scp <- generate_scp_celldimplot(seurat_obj, "Treatment")
  
  # Extract Fibroblast cells for CAF analysis
  fibroblast_cells <- subset(seurat_obj, subset = Celltype_major_man == "Fibroblast")
  
  if(ncol(fibroblast_cells) > 0) {
    cat("Processing CAF cells and recalculating UMAP...\n")
    cat("Fibroblast cell count:", ncol(fibroblast_cells), "\n")
    
    # Recalculate UMAP for CAF cells
    fibroblast_cells <- NormalizeData(fibroblast_cells, verbose = FALSE)
    fibroblast_cells <- FindVariableFeatures(fibroblast_cells, 
                                             selection.method = "vst", 
                                             nfeatures = 2000, 
                                             verbose = FALSE)
    fibroblast_cells <- ScaleData(fibroblast_cells, verbose = FALSE)
    fibroblast_cells <- RunPCA(fibroblast_cells, 
                               features = VariableFeatures(object = fibroblast_cells), 
                               verbose = FALSE)
    fibroblast_cells <- FindNeighbors(fibroblast_cells, dims = 1:20, verbose = FALSE)
    fibroblast_cells <- RunUMAP(fibroblast_cells, dims = 1:20, verbose = FALSE)
    
    cat("CAF UMAP recalculation completed\n")
    
    # CAF subtype visualization
    if("CAF_Subtype" %in% colnames(fibroblast_cells@meta.data)) {
      cat("CAF subtype distribution:\n")
      print(table(fibroblast_cells$CAF_Subtype))
      
      p3_scp <- generate_scp_celldimplot(fibroblast_cells, "CAF_Subtype")
      p4_scp <- generate_scp_celldimplot(fibroblast_cells, "Treatment")
      
    } else {
      cat("Warning: CAF_Subtype not found\n")
      p3_scp <- ggplot() + ggtitle("CAF subtype not found") + theme_void()
      p4_scp <- ggplot() + ggtitle("CAF subtype not found") + theme_void()
    }
  } else {
    cat("Warning: No Fibroblast cells found\n")
    p3_scp <- ggplot() + ggtitle("No Fibroblast cells") + theme_void()
    p4_scp <- ggplot() + ggtitle("No Fibroblast cells") + theme_void()
    fibroblast_cells <- NULL
  }
  
  # Generate statistical plots
  cat("Generating statistical plots...\n")
  p5_scp <- generate_scp_cellstatplot(seurat_obj, "Celltype_major_man", "Treatment")
  
  if(!is.null(fibroblast_cells) && "CAF_Subtype" %in% colnames(fibroblast_cells@meta.data)) {
    p6_scp <- generate_scp_cellstatplot(fibroblast_cells, "CAF_Subtype", "Treatment")
  } else {
    p6_scp <- ggplot() + ggtitle("CAF statistics not available") + theme_void()
  }
  
  # Combine plots
  cat("Saving plots...\n")
  
  # Comprehensive PDF
  pdf(paste0(output_prefix, "_comprehensive_analysis.pdf"), 
      width = 16, height = 20, useDingbats = FALSE)
  
  if(exists("p3_scp") && !inherits(p3_scp, "try-error")) {
    combined_plot_scp <- (p1_scp | p2_scp) / (p3_scp | p4_scp) / (p5_scp | p6_scp)
  } else {
    combined_plot_scp <- (p1_scp | p2_scp) / p5_scp
  }
  
  print(combined_plot_scp)
  dev.off()
  
  # Individual plots
  pdf(paste0(output_prefix, "_all_celltypes.pdf"), width = 10, height = 8, useDingbats = FALSE)
  print(p1_scp)
  dev.off()
  
  pdf(paste0(output_prefix, "_treatment_groups.pdf"), width = 10, height = 8, useDingbats = FALSE)
  print(p2_scp)
  dev.off()
  
  if(exists("p3_scp") && !inherits(p3_scp, "try-error")) {
    pdf(paste0(output_prefix, "_caf_subtypes.pdf"), width = 10, height = 8, useDingbats = FALSE)
    print(p3_scp)
    dev.off()
    
    pdf(paste0(output_prefix, "_caf_treatments.pdf"), width = 10, height = 8, useDingbats = FALSE)
    print(p4_scp)
    dev.off()
  }
  
  pdf(paste0(output_prefix, "_celltype_statistics.pdf"), width = 12, height = 8, useDingbats = FALSE)
  print(p5_scp)
  dev.off()
  
  if(exists("p6_scp") && !inherits(p6_scp, "try-error")) {
    pdf(paste0(output_prefix, "_caf_statistics.pdf"), width = 12, height = 8, useDingbats = FALSE)
    print(p6_scp)
    dev.off()
  }
  
  # High-resolution PNG
  png(paste0(output_prefix, "_comprehensive_analysis.png"), 
      width = 4800, height = 6000, res = 300)
  print(combined_plot_scp)
  dev.off()
  
  # Return results
  return(list(
    p1 = p1_scp,
    p2 = p2_scp,
    p3 = p3_scp,
    p4 = p4_scp,
    p5 = p5_scp,
    p6 = p6_scp,
    fibroblast_cells = fibroblast_cells
  ))
}

# ==============================================================================
# Statistical Report
# ==============================================================================

generate_statistical_report <- function(seurat_obj, fibroblast_cells = NULL) {
  cat("\n=== Statistical Report ===\n")
  cat("Total cells:", ncol(seurat_obj), "\n")
  cat("Total genes:", nrow(seurat_obj), "\n")
  
  cat("\nTreatment group distribution:\n")
  treatment_table <- table(seurat_obj$Treatment)
  print(treatment_table)
  cat("Treatment group proportions (%):\n")
  print(round(prop.table(treatment_table) * 100, 1))
  
  cat("\nCell type distribution:\n")
  celltype_table <- table(seurat_obj$Celltype_major_man)
  print(celltype_table)
  cat("Cell type proportions (%):\n")
  print(round(prop.table(celltype_table) * 100, 1))
  
  if(!is.null(fibroblast_cells)) {
    cat("\nFibroblast statistics:\n")
    cat("Total Fibroblasts:", ncol(fibroblast_cells), "\n")
    
    if("CAF_Subtype" %in% colnames(fibroblast_cells@meta.data)) {
      cat("CAF subtype distribution:\n")
      caf_table <- table(fibroblast_cells$CAF_Subtype)
      print(caf_table)
      cat("CAF subtype proportions (%):\n")
      print(round(prop.table(caf_table) * 100, 1))
      
      cat("\nCAF subtype distribution by treatment:\n")
      caf_treatment_table <- table(fibroblast_cells$Treatment, fibroblast_cells$CAF_Subtype)
      print(caf_treatment_table)
      
      cat("\nCAF subtype proportions by treatment (%):\n")
      caf_treatment_prop <- prop.table(caf_treatment_table, 1) * 100
      print(round(caf_treatment_prop, 1))
    }
  }
}

# ==============================================================================
# 3D UMAP (Optional)
# ==============================================================================

generate_3d_umap <- function(seurat_obj, output_file = "scp_umap_3d_celltypes.html") {
  tryCatch({
    if("umap_3d" %in% names(seurat_obj@reductions)) {
      p7_scp <- CellDimPlot3D(
        srt = seurat_obj,
        group.by = "Celltype_major_man",
        reduction = "umap_3d",
        theme_use = "theme_blank"
      )
      
      htmlwidgets::saveWidget(p7_scp, output_file)
      cat("3D UMAP saved as:", output_file, "\n")
      return(p7_scp)
    } else {
      cat("3D UMAP coordinates not found, skipping 3D visualization\n")
      return(NULL)
    }
  }, error = function(e) {
    cat("3D UMAP generation failed:", e$message, "\n")
    return(NULL)
  })
}

# ==============================================================================
# Example Usage
# ==============================================================================

# Uncomment and modify the following lines to run:

# # Load data
# seurat_obj <- readRDS("processed_seurat_object.rds")
# 
# # Run SCP visualization
# results <- run_scp_visualization(seurat_obj, output_prefix = "scp_umap")
# 
# # Generate statistical report
# generate_statistical_report(seurat_obj, results$fibroblast_cells)
# 
# # Save recalculated CAF data
# if(!is.null(results$fibroblast_cells)) {
#   saveRDS(results$fibroblast_cells, "CAF_cells_recalculated_UMAP.rds")
# }

cat("Script 02: UMAP Visualization - Functions loaded successfully!\n")
cat("Output files will include:\n")
cat("1. Comprehensive analysis: scp_umap_comprehensive_analysis.pdf\n")
cat("2. All cell types: scp_umap_all_celltypes.pdf\n")
cat("3. Treatment groups: scp_umap_treatment_groups.pdf\n")
cat("4. CAF subtypes: scp_umap_caf_subtypes.pdf\n")
cat("5. CAF treatments: scp_umap_caf_treatments.pdf\n")
cat("6. Cell type statistics: scp_umap_celltype_statistics.pdf\n")
cat("7. CAF statistics: scp_umap_caf_statistics.pdf\n")
