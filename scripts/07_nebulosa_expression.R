# ==============================================================================
# Script 07: Nebulosa Expression Density Visualization
# ScN Mouse Single-Cell RNA-seq Analysis Pipeline
# ==============================================================================
# 
# Description:
#   This script generates expression density visualizations using the Nebulosa
#   package for:
#   - IL1R1 and COL1A1 expression patterns
#   - Co-expression analysis in CAF subtypes
#   - Comparison with traditional feature plots
#
# Input:
#   - CAF Seurat object with UMAP coordinates
#
# Output:
#   - Expression density plots (PDF)
#   - Co-expression visualizations
#   - Statistical summaries (CSV)
#
# Author: Biomedical Research Assistant
# Date: 2025
# ==============================================================================

cat("=== Nebulosa Expression Density Analysis ===\n")

# Set renv environment
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}
renv::activate()

# Load required packages
suppressMessages({
  library(Nebulosa)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(viridis)
})

# Set PDF parameters
pdf.options(useDingbats = FALSE)

cat("Using renv environment, R version:", R.version.string, "\n")

# ==============================================================================
# Color Configuration
# ==============================================================================

elegant_caf_colors <- c(
  "iCAF" = "#E63946",
  "myCAF" = "#457B9D",
  "apCAF" = "#F77F00",
  "matCAF" = "#2A9D8F"
)

# ==============================================================================
# Gene Availability Check
# ==============================================================================

check_gene_availability <- function(seurat_obj, gene_list) {
  available_genes <- list()
  
  for(gene in gene_list) {
    # Check different case variations
    variations <- c(gene, toupper(gene), tolower(gene), 
                   paste0(toupper(substr(gene, 1, 1)), tolower(substr(gene, 2, nchar(gene)))))
    
    for(var in variations) {
      if(var %in% rownames(seurat_obj)) {
        available_genes[[gene]] <- var
        break
      }
    }
  }
  
  return(available_genes)
}

# ==============================================================================
# Single Gene Density Plot
# ==============================================================================

plot_gene_density <- function(seurat_obj, gene, size = 2, 
                              reduction = "umap", title = NULL) {
  
  if(!gene %in% rownames(seurat_obj)) {
    cat("Gene", gene, "not found\n")
    return(NULL)
  }
  
  if(is.null(title)) {
    title <- paste(gene, "Expression Density")
  }
  
  p <- plot_density(seurat_obj, c(gene), size = size) + 
    ggtitle(title) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 18)
    )
  
  return(p)
}

# ==============================================================================
# Traditional Feature Plot (for comparison)
# ==============================================================================

plot_gene_traditional <- function(seurat_obj, gene, pt.size = 2, 
                                  reduction = "umap", title = NULL) {
  
  if(!gene %in% rownames(seurat_obj)) {
    cat("Gene", gene, "not found\n")
    return(NULL)
  }
  
  if(is.null(title)) {
    title <- paste(gene, "Expression (Traditional)")
  }
  
  p <- FeaturePlot(
    seurat_obj, 
    features = gene,
    reduction = reduction,
    pt.size = pt.size
  ) + 
    ggtitle(title) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 18)
    )
  
  return(p)
}

# ==============================================================================
# Co-expression Density Plot
# ==============================================================================

plot_coexpression_density <- function(seurat_obj, gene1, gene2, 
                                      size = 2, joint = TRUE, combine = FALSE) {
  
  if(!all(c(gene1, gene2) %in% rownames(seurat_obj))) {
    missing <- c(gene1, gene2)[!c(gene1, gene2) %in% rownames(seurat_obj)]
    cat("Missing genes:", paste(missing, collapse = ", "), "\n")
    return(NULL)
  }
  
  p <- plot_density(seurat_obj, c(gene1, gene2), 
                   joint = joint, combine = combine, size = size)
  
  return(p)
}

# ==============================================================================
# CAF Subtype-specific Density
# ==============================================================================

plot_subtype_density <- function(seurat_obj, gene, subtype, 
                                 subtype_col = "CAF_Subtype", size = 2) {
  
  if(!gene %in% rownames(seurat_obj)) {
    cat("Gene", gene, "not found\n")
    return(NULL)
  }
  
  if(!subtype %in% seurat_obj@meta.data[[subtype_col]]) {
    cat("Subtype", subtype, "not found\n")
    return(NULL)
  }
  
  # Subset to specific subtype
  subtype_cells <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[subtype_col]] == subtype]
  
  if(length(subtype_cells) < 20) {
    cat("Insufficient cells for", subtype, "(", length(subtype_cells), ")\n")
    return(NULL)
  }
  
  subtype_data <- subset(seurat_obj, cells = subtype_cells)
  
  p <- plot_density(subtype_data, c(gene), size = size) + 
    ggtitle(paste(subtype, gene, "Expression Density")) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 18)
    )
  
  return(p)
}

# ==============================================================================
# Expression Statistics
# ==============================================================================

calculate_expression_stats <- function(seurat_obj, gene, group_col = "CAF_Subtype") {
  
  if(!gene %in% rownames(seurat_obj)) {
    return(NULL)
  }
  
  # Get expression
  tryCatch({
    expr <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")[gene, ]
  }, error = function(e) {
    expr <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")[gene, ]
  })
  
  expr_df <- data.frame(
    Cell = colnames(seurat_obj),
    Expression = as.numeric(expr),
    Group = seurat_obj@meta.data[[group_col]]
  )
  
  stats <- expr_df %>%
    group_by(Group) %>%
    summarise(
      Cell_Count = n(),
      Expressing_Cells = sum(Expression > 0),
      Expression_Percentage = round(sum(Expression > 0) / n() * 100, 2),
      Mean_Expression = round(mean(Expression), 4),
      Median_Expression = round(median(Expression), 4),
      .groups = 'drop'
    )
  
  return(stats)
}

# ==============================================================================
# Main Analysis Pipeline
# ==============================================================================

run_nebulosa_analysis <- function(caf_obj, target_genes = c("Il1r1", "Col1a1"),
                                  subtype_col = "CAF_Subtype",
                                  output_prefix = "nebulosa") {
  
  cat("=== Starting Nebulosa Analysis ===\n")
  cat("CAF cells:", ncol(caf_obj), "\n")
  cat("Genes:", nrow(caf_obj), "\n")
  
  # Check reduction
  available_reductions <- names(caf_obj@reductions)
  cat("Available reductions:", paste(available_reductions, collapse = ", "), "\n")
  
  reduction_name <- "umap"
  if(!"umap" %in% available_reductions) {
    if("UMAP" %in% available_reductions) {
      reduction_name <- "UMAP"
    } else {
      cat("Error: No UMAP found\n")
      return(NULL)
    }
  }
  
  # Check genes
  available_genes <- check_gene_availability(caf_obj, target_genes)
  cat("Available genes:", paste(names(available_genes), "->", unlist(available_genes), collapse = ", "), "\n")
  
  results <- list()
  
  # Process each gene
  for(gene_name in names(available_genes)) {
    gene <- available_genes[[gene_name]]
    cat("\nProcessing", gene, "...\n")
    
    # Expression statistics
    stats <- calculate_expression_stats(caf_obj, gene, subtype_col)
    if(!is.null(stats)) {
      cat("Expression statistics:\n")
      print(stats)
      write.csv(stats, paste0(output_prefix, "_", gene_name, "_stats.csv"), row.names = FALSE)
    }
    
    # Density plot
    p_density <- plot_gene_density(caf_obj, gene, size = 2)
    if(!is.null(p_density)) {
      grDevices::pdf(paste0(output_prefix, "_", gene_name, "_density.pdf"), 
                    width = 5, height = 4, useDingbats = FALSE)
      print(p_density)
      grDevices::dev.off()
      results[[paste0(gene_name, "_density")]] <- p_density
    }
    
    # Traditional plot
    p_traditional <- plot_gene_traditional(caf_obj, gene, pt.size = 2)
    if(!is.null(p_traditional)) {
      grDevices::pdf(paste0(output_prefix, "_", gene_name, "_traditional.pdf"), 
                    width = 5, height = 4, useDingbats = FALSE)
      print(p_traditional)
      grDevices::dev.off()
      results[[paste0(gene_name, "_traditional")]] <- p_traditional
    }
    
    # Subtype-specific plots
    if(subtype_col %in% colnames(caf_obj@meta.data)) {
      main_subtypes <- c("matCAF", "myCAF")
      
      for(subtype in main_subtypes) {
        p_subtype <- plot_subtype_density(caf_obj, gene, subtype, subtype_col, size = 2)
        if(!is.null(p_subtype)) {
          grDevices::pdf(paste0(output_prefix, "_", gene_name, "_", subtype, ".pdf"), 
                        width = 10, height = 8, useDingbats = FALSE)
          print(p_subtype)
          grDevices::dev.off()
          results[[paste0(gene_name, "_", subtype)]] <- p_subtype
        }
      }
    }
  }
  
  # Co-expression analysis
  if(length(available_genes) >= 2) {
    gene1 <- available_genes[[1]]
    gene2 <- available_genes[[2]]
    
    cat("\nCo-expression analysis:", gene1, "and", gene2, "\n")
    
    p_coexpr <- plot_coexpression_density(caf_obj, gene1, gene2, 
                                          size = 2, joint = TRUE, combine = FALSE)
    
    if(!is.null(p_coexpr) && is.list(p_coexpr) && length(p_coexpr) > 0) {
      grDevices::pdf(paste0(output_prefix, "_coexpression.pdf"), 
                    width = 10, height = 8, useDingbats = FALSE)
      if(length(p_coexpr) >= 2) {
        combined <- p_coexpr[[1]] | p_coexpr[[2]]
        print(combined)
      } else {
        print(p_coexpr[[1]])
      }
      grDevices::dev.off()
      results$coexpression <- p_coexpr
    }
  }
  
  # CAF subtype distribution plot
  if(subtype_col %in% colnames(caf_obj@meta.data)) {
    p_subtype_dist <- DimPlot(
      caf_obj, 
      reduction = reduction_name, 
      group.by = subtype_col,
      pt.size = 2,
      cols = elegant_caf_colors
    ) + 
      ggtitle("CAF Subtype Distribution") +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16)
      )
    
    grDevices::pdf(paste0(output_prefix, "_subtype_distribution.pdf"), 
                  width = 5, height = 4, useDingbats = FALSE)
    print(p_subtype_dist)
    grDevices::dev.off()
    results$subtype_distribution <- p_subtype_dist
  }
  
  # Comprehensive comparison plot
  if(length(results) > 0) {
    cat("\nCreating comprehensive plot...\n")
    
    plots_to_combine <- list()
    for(gene_name in names(available_genes)) {
      if(paste0(gene_name, "_density") %in% names(results)) {
        plots_to_combine <- c(plots_to_combine, list(results[[paste0(gene_name, "_density")]]))
      }
      if(paste0(gene_name, "_traditional") %in% names(results)) {
        plots_to_combine <- c(plots_to_combine, list(results[[paste0(gene_name, "_traditional")]]))
      }
    }
    
    if(length(plots_to_combine) >= 2) {
      grDevices::pdf(paste0(output_prefix, "_comprehensive.pdf"), 
                    width = 20, height = 12, useDingbats = FALSE)
      comprehensive <- wrap_plots(plots_to_combine, ncol = 2)
      print(comprehensive)
      grDevices::dev.off()
    }
  }
  
  cat("\n=== Nebulosa Analysis Complete ===\n")
  cat("Output files:\n")
  cat("- Density plots:", paste0(output_prefix, "_*_density.pdf\n"))
  cat("- Traditional plots:", paste0(output_prefix, "_*_traditional.pdf\n"))
  cat("- Subtype-specific plots:", paste0(output_prefix, "_*_subtype.pdf\n"))
  cat("- Co-expression plot:", paste0(output_prefix, "_coexpression.pdf\n"))
  cat("- Statistics:", paste0(output_prefix, "_*_stats.csv\n"))
  cat("- Comprehensive plot:", paste0(output_prefix, "_comprehensive.pdf\n"))
  
  return(results)
}

# ==============================================================================
# Example Usage
# ==============================================================================

# Uncomment and modify to run:

# # Load CAF data
# caf_obj <- readRDS("CAF_cells_recalculated_UMAP.rds")
# 
# # Run Nebulosa analysis
# results <- run_nebulosa_analysis(
#   caf_obj,
#   target_genes = c("Il1r1", "Col1a1"),
#   subtype_col = "CAF_Subtype",
#   output_prefix = "caf_nebulosa"
# )

cat("Script 07: Nebulosa Expression Density - Functions loaded successfully!\n")
cat("\nKey parameters:\n")
cat("- plot_density size = 2\n")
cat("- Co-expression: joint = TRUE, combine = FALSE\n")
cat("- All plots are 300 DPI high-quality PDFs\n")
