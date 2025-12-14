# ==============================================================================
# Script 04: Differential Expression and Enrichment Analysis
# ScN Mouse Single-Cell RNA-seq Analysis Pipeline
# ==============================================================================
# 
# Description:
#   This script performs differential expression analysis and pathway 
#   enrichment analysis for CAF subtypes using the SCP package, including:
#   - Differential expression between CAF subtypes
#   - GO and KEGG pathway enrichment
#   - Feature heatmaps with functional annotations
#   - Volcano plots and network visualizations
#
# Input:
#   - CAF Seurat object with subtype annotations
#
# Output:
#   - Differential gene lists (CSV)
#   - Enrichment plots (PDF)
#   - Feature heatmaps (PDF)
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
  library(ComplexHeatmap)
})

# Set high-quality plot parameters
options(repr.plot.width = 12, repr.plot.height = 8)

# Set PDF font type for editable format
pdf.options(useDingbats = FALSE)
plt.rcParams <- list(
  'pdf.fonttype' = 42,
  'ps.fonttype' = 42
)

cat("=== CAF Differential Expression and Enrichment Analysis ===\n")

# ==============================================================================
# Differential Expression Analysis
# ==============================================================================

run_differential_expression <- function(caf_obj, group_by = "CAF_Subtype", 
                                        fc_threshold = 1, min_cells = 10) {
  
  cat("\n=== Running Differential Expression Analysis ===\n")
  
  # Check CAF subtype distribution
  caf_subtype_table <- table(caf_obj[[group_by]])
  cat("CAF subtype distribution:\n")
  print(caf_subtype_table)
  
  # Check for sufficient cells
  sufficient_groups <- names(caf_subtype_table)[caf_subtype_table >= min_cells]
  cat("Groups with sufficient cells (>=", min_cells, "):", 
      paste(sufficient_groups, collapse = ", "), "\n")
  
  if(length(sufficient_groups) < 2) {
    cat("Warning: Not enough CAF subtypes for differential analysis\n")
    return(NULL)
  }
  
  # Run differential expression
  cat("Running differential expression analysis...\n")
  caf_obj <- RunDEtest(
    srt = caf_obj, 
    group_by = group_by, 
    fc.threshold = fc_threshold, 
    only.pos = FALSE,
    verbose = TRUE
  )
  
  # Check results
  de_tool_name <- paste0("DEtest_", group_by)
  if(!de_tool_name %in% names(caf_obj@tools)) {
    cat("Error: Differential expression analysis failed\n")
    return(NULL)
  }
  
  cat("Differential expression analysis completed\n")
  return(caf_obj)
}

# ==============================================================================
# Extract and Filter DEGs
# ==============================================================================

extract_degs <- function(caf_obj, group_by = "CAF_Subtype", 
                        log2fc_threshold = 1, padj_threshold = 0.05) {
  
  de_tool_name <- paste0("DEtest_", group_by)
  
  # Extract all markers
  DEGs <- caf_obj@tools[[de_tool_name]]$AllMarkers_wilcox
  cat("Total differential genes:", nrow(DEGs), "\n")
  
  # Filter significant DEGs
  DEGs_significant <- DEGs[with(DEGs, avg_log2FC > log2fc_threshold & 
                                       p_val_adj < padj_threshold), ]
  cat("Significant upregulated genes (FC>", 2^log2fc_threshold, 
      ", padj<", padj_threshold, "):", nrow(DEGs_significant), "\n")
  
  # If too few, relax threshold
  if(nrow(DEGs_significant) < 10) {
    cat("Relaxing threshold...\n")
    DEGs_significant <- DEGs[with(DEGs, avg_log2FC > 0.5 & p_val_adj < 0.05), ]
    cat("With relaxed threshold:", nrow(DEGs_significant), "\n")
  }
  
  # Show top genes per subtype
  cat("\nTop genes per CAF subtype:\n")
  for(subtype in unique(DEGs_significant$group1)) {
    subtype_genes <- DEGs_significant[DEGs_significant$group1 == subtype, ]
    if(nrow(subtype_genes) > 0) {
      top_genes <- head(subtype_genes[order(subtype_genes$avg_log2FC, 
                                             decreasing = TRUE), ], 5)
      cat(subtype, "top 5:", paste(top_genes$gene, collapse = ", "), "\n")
    }
  }
  
  return(list(all = DEGs, significant = DEGs_significant))
}

# ==============================================================================
# Feature Annotation
# ==============================================================================

annotate_features <- function(caf_obj, species = "Mus_musculus", 
                             db = c("TF", "CSPA")) {
  
  cat("\n=== Annotating Features ===\n")
  cat("Species:", species, "\n")
  cat("Databases:", paste(db, collapse = ", "), "\n")
  
  tryCatch({
    caf_obj <- AnnotateFeatures(caf_obj, species = species, db = db)
    cat("Feature annotation completed\n")
  }, error = function(e) {
    cat("Feature annotation failed:", e$message, "\n")
  })
  
  return(caf_obj)
}

# ==============================================================================
# Volcano Plot
# ==============================================================================

generate_volcano_plot <- function(caf_obj, group_by = "CAF_Subtype", 
                                  output_file = "volcano_plot.pdf") {
  
  cat("\n=== Generating Volcano Plot ===\n")
  
  volcano_plot <- VolcanoPlot(
    srt = caf_obj, 
    group_by = group_by
  )
  
  pdf(output_file, width = 12, height = 10, useDingbats = FALSE)
  print(volcano_plot)
  dev.off()
  
  cat("Volcano plot saved:", output_file, "\n")
  return(volcano_plot)
}

# ==============================================================================
# Feature Heatmap
# ==============================================================================

generate_feature_heatmap <- function(caf_obj, DEGs_significant, 
                                     group_by = "CAF_Subtype",
                                     max_genes_per_group = 25,
                                     species = "Mus_musculus",
                                     output_file = "feature_heatmap.pdf") {
  
  cat("\n=== Generating Feature Heatmap ===\n")
  
  # Select genes for heatmap
  if(nrow(DEGs_significant) > max_genes_per_group * length(unique(DEGs_significant$group1))) {
    top_genes_per_group <- DEGs_significant %>%
      group_by(group1) %>%
      slice_max(order_by = avg_log2FC, n = max_genes_per_group) %>%
      pull(gene) %>%
      unique()
  } else {
    top_genes_per_group <- DEGs_significant$gene
  }
  
  cat("Genes for heatmap:", length(top_genes_per_group), "\n")
  
  if(length(top_genes_per_group) == 0) {
    cat("Warning: No genes available for heatmap\n")
    return(NULL)
  }
  
  tryCatch({
    ht <- FeatureHeatmap(
      srt = caf_obj, 
      group.by = group_by, 
      features = top_genes_per_group, 
      feature_split = DEGs_significant$group1[match(top_genes_per_group, 
                                                     DEGs_significant$gene)],
      species = species, 
      db = c("GO_BP", "KEGG", "WikiPathway"), 
      anno_terms = TRUE,
      feature_annotation = c("TF"),
      feature_annotation_palcolor = list(c("gold", "steelblue")),
      height = 8, 
      width = 6
    )
    
    pdf(output_file, width = 12, height = 16, useDingbats = FALSE)
    print(ht$plot)
    dev.off()
    
    cat("Feature heatmap saved:", output_file, "\n")
    return(ht)
    
  }, error = function(e) {
    cat("Feature heatmap failed:", e$message, "\n")
    return(NULL)
  })
}

# ==============================================================================
# Enrichment Analysis
# ==============================================================================

run_enrichment_analysis <- function(caf_obj, group_by = "CAF_Subtype",
                                    db = "GO_BP", species = "Mus_musculus",
                                    de_threshold = "avg_log2FC > log2(1.5) & p_val_adj < 0.05") {
  
  cat("\n=== Running Enrichment Analysis ===\n")
  cat("Database:", db, "\n")
  
  tryCatch({
    caf_obj <- RunEnrichment(
      srt = caf_obj, 
      group_by = group_by, 
      db = db, 
      species = species,
      DE_threshold = de_threshold
    )
    
    cat("Enrichment analysis completed\n")
    return(caf_obj)
    
  }, error = function(e) {
    cat("Enrichment analysis failed:", e$message, "\n")
    return(caf_obj)
  })
}

# ==============================================================================
# Enrichment Visualization
# ==============================================================================

generate_enrichment_plots <- function(caf_obj, group_by = "CAF_Subtype",
                                      output_prefix = "enrichment") {
  
  cat("\n=== Generating Enrichment Plots ===\n")
  
  enrichment_tool_name <- paste0("Enrichment_", group_by)
  
  if(!enrichment_tool_name %in% names(caf_obj@tools)) {
    cat("Warning: No enrichment results found\n")
    return(NULL)
  }
  
  enrichment_results <- caf_obj@tools[[enrichment_tool_name]]
  available_groups <- names(enrichment_results)
  cat("Available groups:", paste(available_groups, collapse = ", "), "\n")
  
  if(length(available_groups) < 2) {
    cat("Warning: Need at least 2 groups for comparison\n")
    return(NULL)
  }
  
  compare_groups <- available_groups[1:min(2, length(available_groups))]
  cat("Comparing:", paste(compare_groups, collapse = " vs "), "\n")
  
  plots <- list()
  
  # Bar plot
  tryCatch({
    cat("Creating bar plot...\n")
    bar_plot <- EnrichmentPlot(
      srt = caf_obj, 
      group_by = group_by, 
      group_use = compare_groups,
      plot_type = "bar"
    )
    
    pdf(paste0(output_prefix, "_bar_plot.pdf"), width = 14, height = 10, useDingbats = FALSE)
    print(bar_plot)
    dev.off()
    
    plots$bar <- bar_plot
  }, error = function(e) {
    cat("Bar plot failed:", e$message, "\n")
  })
  
  # Word cloud
  tryCatch({
    cat("Creating word cloud...\n")
    wordcloud_plot <- EnrichmentPlot(
      srt = caf_obj, 
      group_by = group_by, 
      group_use = compare_groups,
      plot_type = "wordcloud"
    )
    
    pdf(paste0(output_prefix, "_wordcloud.pdf"), width = 12, height = 8, useDingbats = FALSE)
    print(wordcloud_plot)
    dev.off()
    
    plots$wordcloud <- wordcloud_plot
  }, error = function(e) {
    cat("Word cloud failed:", e$message, "\n")
  })
  
  # Network plot
  if(length(available_groups) >= 1) {
    tryCatch({
      cat("Creating network plot...\n")
      network_plot <- EnrichmentPlot(
        srt = caf_obj, 
        group_by = group_by, 
        group_use = compare_groups[1],
        plot_type = "network"
      )
      
      pdf(paste0(output_prefix, "_network.pdf"), width = 12, height = 10, useDingbats = FALSE)
      print(network_plot)
      dev.off()
      
      plots$network <- network_plot
    }, error = function(e) {
      cat("Network plot failed:", e$message, "\n")
    })
  }
  
  # Enrichment map
  tryCatch({
    cat("Creating enrichment map...\n")
    enrichmap_plot <- EnrichmentPlot(
      srt = caf_obj, 
      group_by = group_by, 
      group_use = compare_groups[1],
      plot_type = "enrichmap"
    )
    
    pdf(paste0(output_prefix, "_map.pdf"), width = 12, height = 10, useDingbats = FALSE)
    print(enrichmap_plot)
    dev.off()
    
    plots$enrichmap <- enrichmap_plot
  }, error = function(e) {
    cat("Enrichment map failed:", e$message, "\n")
  })
  
  # Comparison plot
  if(length(available_groups) > 1) {
    tryCatch({
      cat("Creating comparison plot...\n")
      comparison_plot <- EnrichmentPlot(
        srt = caf_obj, 
        group_by = group_by, 
        plot_type = "comparison"
      )
      
      pdf(paste0(output_prefix, "_comparison.pdf"), width = 14, height = 10, useDingbats = FALSE)
      print(comparison_plot)
      dev.off()
      
      plots$comparison <- comparison_plot
    }, error = function(e) {
      cat("Comparison plot failed:", e$message, "\n")
    })
  }
  
  return(plots)
}

# ==============================================================================
# Main Analysis Pipeline
# ==============================================================================

run_de_enrichment_pipeline <- function(caf_obj, group_by = "CAF_Subtype",
                                       output_prefix = "caf_de",
                                       species = "Mus_musculus") {
  
  cat("=== Starting DE and Enrichment Pipeline ===\n")
  
  # Step 1: Differential expression
  caf_obj <- run_differential_expression(caf_obj, group_by)
  if(is.null(caf_obj)) return(NULL)
  
  # Step 2: Extract DEGs
  degs <- extract_degs(caf_obj, group_by)
  
  # Save DEG tables
  write.csv(degs$all, paste0(output_prefix, "_all_genes.csv"), row.names = FALSE)
  write.csv(degs$significant, paste0(output_prefix, "_significant_genes.csv"), row.names = FALSE)
  
  # Step 3: Volcano plot
  generate_volcano_plot(caf_obj, group_by, paste0(output_prefix, "_volcano.pdf"))
  
  # Step 4: Feature annotation
  caf_obj <- annotate_features(caf_obj, species)
  
  # Step 5: Feature heatmap
  generate_feature_heatmap(caf_obj, degs$significant, group_by, 
                          output_file = paste0(output_prefix, "_heatmap.pdf"))
  
  # Step 6: GO enrichment
  caf_obj <- run_enrichment_analysis(caf_obj, group_by, db = "GO_BP", species = species)
  generate_enrichment_plots(caf_obj, group_by, paste0(output_prefix, "_GO"))
  
  # Step 7: KEGG enrichment
  caf_obj_kegg <- run_enrichment_analysis(caf_obj, group_by, db = "KEGG", species = species)
  generate_enrichment_plots(caf_obj_kegg, group_by, paste0(output_prefix, "_KEGG"))
  
  # Save updated object
  saveRDS(caf_obj, paste0(output_prefix, "_analysis.rds"))
  
  cat("\n=== Pipeline Complete ===\n")
  cat("Output files:\n")
  cat("- Volcano plot:", paste0(output_prefix, "_volcano.pdf\n"))
  cat("- Feature heatmap:", paste0(output_prefix, "_heatmap.pdf\n"))
  cat("- GO enrichment plots:", paste0(output_prefix, "_GO_*.pdf\n"))
  cat("- KEGG enrichment plots:", paste0(output_prefix, "_KEGG_*.pdf\n"))
  cat("- DEG tables:", paste0(output_prefix, "_*.csv\n"))
  cat("- Analysis object:", paste0(output_prefix, "_analysis.rds\n"))
  
  return(caf_obj)
}

# ==============================================================================
# Example Usage
# ==============================================================================

# Uncomment and modify to run:

# # Load CAF data
# caf_obj <- readRDS("CAF_cells_recalculated_UMAP.rds")
# 
# # Run pipeline
# caf_obj <- run_de_enrichment_pipeline(
#   caf_obj,
#   group_by = "CAF_Subtype",
#   output_prefix = "caf_de",
#   species = "Mus_musculus"
# )

cat("Script 04: Differential Expression & Enrichment - Functions loaded successfully!\n")
