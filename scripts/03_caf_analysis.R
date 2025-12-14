# ==============================================================================
# Script 03: CAF Subtype Analysis and Characterization
# ScN Mouse Single-Cell RNA-seq Analysis Pipeline
# ==============================================================================
# 
# Description:
#   This script performs comprehensive CAF (Cancer-Associated Fibroblast) 
#   analysis including:
#   - CAF subtype classification and validation
#   - IL1R1 expression analysis across CAF subtypes
#   - ECM remodeling and inflammatory capacity scoring
#   - Matrix CAF and CXCL13+ CAF proportion analysis
#
# Input:
#   - Processed Seurat object with CAF subtype annotations
#
# Output:
#   - CAF analysis visualizations (PDF)
#   - Statistical summary tables (CSV)
#
# Author: Biomedical Research Assistant
# Date: 2025
# ==============================================================================

# Load required packages
suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(RColorBrewer)
    library(viridis)
    library(tibble)
    library(tidyr)
    library(ggpubr)
    library(rstatix)
    library(ggsignif)
    library(ggridges)
    library(ggbeeswarm)
})

# Set random seed
set.seed(42)

# PDF font settings
options(pdf.fonttype = 42)

# ==============================================================================
# Theme and Color Configuration
# ==============================================================================

# Elegant modern theme
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

# Color schemes
treatment_colors <- c(
    "Untreated" = "#2E86AB",
    "Chemo" = "#A23B72",
    "Chemo_Abt263" = "#F18F01",
    "Chemo_Abt263_IL1B" = "#6A994E"
)

caf_subtype_colors <- c(
    "iCAF" = "#E63946",
    "myCAF" = "#457B9D",
    "apCAF" = "#F77F00",
    "matCAF" = "#2A9D8F"
)

# ==============================================================================
# Statistical Test Functions
# ==============================================================================

perform_stat_test <- function(data, group_col, value_col, method = "wilcox.test") {
    tryCatch({
        valid_data <- data[!is.na(data[[value_col]]) & !is.na(data[[group_col]]), ]
        
        if(nrow(valid_data) < 4) {
            cat("Insufficient data for statistical test\n")
            return(NULL)
        }
        
        groups <- unique(valid_data[[group_col]])
        if(length(groups) < 2) {
            cat("Need at least 2 groups for comparison\n")
            return(NULL)
        }
        
        if(method == "wilcox.test") {
            stat_test <- valid_data %>%
                wilcox_test(as.formula(paste(value_col, "~", group_col))) %>%
                add_significance() %>%
                add_xy_position(x = group_col)
        } else if(method == "t.test") {
            stat_test <- valid_data %>%
                t_test(as.formula(paste(value_col, "~", group_col))) %>%
                add_significance() %>%
                add_xy_position(x = group_col)
        }
        return(stat_test)
    }, error = function(e) {
        cat("Error in statistical test:", e$message, "\n")
        return(NULL)
    })
}

# ==============================================================================
# CAF Score Calculation
# ==============================================================================

calculate_caf_marker_scores <- function(fibroblast_cells) {
    
    # CAF marker genes
    icaf_genes <- c("Il1b", "Il6", "Cxcl1", "Cxcl2", "Ccl2")
    mycaf_genes <- c("Acta2", "Tagln", "Myh11", "Cnn1")
    apcaf_genes <- c("H2-Aa", "H2-Ab1", "Cd74", "Ciita")
    matcaf_genes <- c("Col1a1", "Col1a2", "Col3a1", "Fn1", "Lum", "Dcn")
    
    # Check available genes
    available_icaf <- intersect(icaf_genes, rownames(fibroblast_cells))
    available_mycaf <- intersect(mycaf_genes, rownames(fibroblast_cells))
    available_apcaf <- intersect(apcaf_genes, rownames(fibroblast_cells))
    available_matcaf <- intersect(matcaf_genes, rownames(fibroblast_cells))
    
    cat("Available CAF genes - iCAF:", length(available_icaf), 
        "myCAF:", length(available_mycaf), 
        "apCAF:", length(available_apcaf), 
        "matCAF:", length(available_matcaf), "\n")
    
    # Calculate scores
    suppressWarnings({
        if(length(available_icaf) > 0) {
            icaf_data <- GetAssayData(fibroblast_cells, slot = "data")[available_icaf, , drop = FALSE]
            fibroblast_cells@meta.data$iCAF_Score <- colMeans(icaf_data)
        } else {
            fibroblast_cells@meta.data$iCAF_Score <- 0
        }
        
        if(length(available_mycaf) > 0) {
            mycaf_data <- GetAssayData(fibroblast_cells, slot = "data")[available_mycaf, , drop = FALSE]
            fibroblast_cells@meta.data$myCAF_Score <- colMeans(mycaf_data)
        } else {
            fibroblast_cells@meta.data$myCAF_Score <- 0
        }
        
        if(length(available_apcaf) > 0) {
            apcaf_data <- GetAssayData(fibroblast_cells, slot = "data")[available_apcaf, , drop = FALSE]
            fibroblast_cells@meta.data$apCAF_Score <- colMeans(apcaf_data)
        } else {
            fibroblast_cells@meta.data$apCAF_Score <- 0
        }
        
        if(length(available_matcaf) > 0) {
            matcaf_data <- GetAssayData(fibroblast_cells, slot = "data")[available_matcaf, , drop = FALSE]
            fibroblast_cells@meta.data$matCAF_Score <- colMeans(matcaf_data)
        } else {
            fibroblast_cells@meta.data$matCAF_Score <- 0
        }
    })
    
    # Classify CAF subtypes
    caf_scores <- fibroblast_cells@meta.data[, c("iCAF_Score", "myCAF_Score", "apCAF_Score", "matCAF_Score")]
    fibroblast_cells@meta.data$CAF_Subtype <- apply(caf_scores, 1, function(x) {
        max_score <- which.max(x)
        c("iCAF", "myCAF", "apCAF", "matCAF")[max_score]
    })
    
    return(fibroblast_cells)
}

# ==============================================================================
# Analysis 1: Matrix CAF Proportion
# ==============================================================================

analyze_matcaf_proportion <- function(fibroblast_cells, treatment_col = "Treatment") {
    
    cat("\n=== Analysis: Matrix CAF Proportion ===\n")
    
    # Calculate Matrix CAF proportion
    matcaf_summary <- fibroblast_cells@meta.data %>%
        group_by(.data[[treatment_col]]) %>%
        summarise(
            total_caf = n(),
            matcaf_count = sum(CAF_Subtype == "matCAF"),
            matcaf_proportion = matcaf_count / total_caf * 100,
            .groups = "drop"
        )
    
    cat("Matrix CAF proportions:\n")
    print(matcaf_summary)
    
    # Create plot
    p <- ggplot(matcaf_summary, aes_string(x = treatment_col, y = "matcaf_proportion", 
                                            fill = treatment_col)) +
        geom_col(alpha = 0.85, width = 0.65, color = "white", linewidth = 1.5) +
        scale_fill_manual(values = treatment_colors) +
        labs(title = "Matrix CAF Proportion Analysis",
             subtitle = "ECM-producing fibroblasts across treatment groups",
             x = "Treatment Group", y = "Matrix CAF Proportion (%)") +
        guides(fill = "none") +
        geom_text(aes(label = paste0(round(matcaf_proportion, 1), "%\n(", 
                                     matcaf_count, "/", total_caf, ")")), 
                  vjust = -0.5, size = 5, fontface = "bold", color = "#2C3E50") +
        coord_cartesian(ylim = c(0, max(matcaf_summary$matcaf_proportion) * 1.15)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"))
    
    return(list(plot = p, summary = matcaf_summary))
}

# ==============================================================================
# Analysis 2: IL1R1 Expression in CAF Subtypes
# ==============================================================================

analyze_il1r1_expression <- function(fibroblast_cells) {
    
    cat("\n=== Analysis: IL1R1 Expression in CAF Subtypes ===\n")
    
    # Check IL1R1 expression
    if(!"Il1r1" %in% rownames(fibroblast_cells)) {
        cat("IL1R1 gene not found\n")
        return(NULL)
    }
    
    # Get IL1R1 expression
    suppressWarnings({
        il1r1_expr <- GetAssayData(fibroblast_cells, slot = "data")["Il1r1", ]
        fibroblast_cells@meta.data$Il1r1_Expression <- il1r1_expr
        fibroblast_cells@meta.data$Il1r1_positive <- il1r1_expr > 0
    })
    
    # Calculate statistics
    il1r1_caf_summary <- fibroblast_cells@meta.data %>%
        group_by(CAF_Subtype) %>%
        summarise(
            n_cells = n(),
            mean_il1r1 = mean(Il1r1_Expression, na.rm = TRUE),
            median_il1r1 = median(Il1r1_Expression, na.rm = TRUE),
            positive_cells = sum(Il1r1_positive, na.rm = TRUE),
            positive_rate = positive_cells / n_cells * 100,
            se_il1r1 = sd(Il1r1_Expression, na.rm = TRUE) / sqrt(n_cells),
            .groups = "drop"
        ) %>%
        arrange(desc(mean_il1r1))
    
    cat("IL1R1 expression in CAF subtypes:\n")
    print(il1r1_caf_summary)
    
    # Violin plot
    p_violin <- ggplot(fibroblast_cells@meta.data, 
                       aes(x = CAF_Subtype, y = Il1r1_Expression, fill = CAF_Subtype)) +
        geom_violin(alpha = 0.75, width = 0.8, color = "white", linewidth = 1, trim = FALSE) +
        geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.9, 
                     color = "#2C3E50", linewidth = 0.8) +
        stat_summary(fun = median, geom = "point", size = 2.5, color = "#2C3E50") +
        scale_fill_manual(values = caf_subtype_colors) +
        labs(title = "IL1R1 Expression Levels in CAF Subtypes",
             subtitle = "Expression distribution across CAF functional phenotypes",
             x = "CAF Subtype", y = "IL1R1 Expression Level") +
        guides(fill = "none") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"))
    
    # Proportion plot
    p_rate <- ggplot(il1r1_caf_summary, 
                     aes(x = reorder(CAF_Subtype, positive_rate), y = positive_rate, 
                         fill = CAF_Subtype)) +
        geom_col(alpha = 0.85, width = 0.65, color = "white", linewidth = 1.5) +
        scale_fill_manual(values = caf_subtype_colors) +
        labs(title = "IL1R1 Positive Rate in CAF Subtypes",
             subtitle = "Percentage of IL1R1+ cells across functional CAF phenotypes",
             x = "CAF Subtype", y = "IL1R1 Positive Rate (%)") +
        guides(fill = "none") +
        geom_text(aes(label = paste0(round(positive_rate, 1), "%\n(", 
                                     positive_cells, "/", n_cells, ")")), 
                  vjust = -0.3, size = 4.5, fontface = "bold", color = "#2C3E50") +
        coord_cartesian(ylim = c(0, max(il1r1_caf_summary$positive_rate) * 1.15)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"))
    
    return(list(
        violin_plot = p_violin,
        rate_plot = p_rate,
        summary = il1r1_caf_summary,
        data = fibroblast_cells
    ))
}

# ==============================================================================
# Analysis 3: CAF Functional Scores (ECM & Inflammatory)
# ==============================================================================

analyze_caf_functional_scores <- function(fibroblast_cells, treatment_col = "Treatment") {
    
    cat("\n=== Analysis: CAF Functional Scores ===\n")
    
    # Define functional gene sets
    caf_function_genes <- list(
        "ECM_Remodeling" = c("Col1a1", "Col1a2", "Col3a1", "Col5a1", "Fn1", "Lox", 
                            "Loxl2", "Mmp2", "Mmp9", "Mmp14", "Sparc", "Tnc", 
                            "Postn", "Dpt", "Dcn", "Bgn"),
        "Inflammatory" = c("Il6", "Il1b", "Il1a", "Tnf", "Cxcl1", "Cxcl2", "Ccl2", 
                          "Ccl7", "Ptgs2", "Nfkb1", "Stat3", "Cxcl12", "Cxcl14", "Ccl11")
    )
    
    # Calculate functional scores
    for(function_name in names(caf_function_genes)) {
        available_genes <- caf_function_genes[[function_name]][
            caf_function_genes[[function_name]] %in% rownames(fibroblast_cells)
        ]
        
        cat(function_name, "genes available:", length(available_genes), "out of", 
            length(caf_function_genes[[function_name]]), "\n")
        
        if(length(available_genes) > 0) {
            suppressWarnings({
                function_data <- GetAssayData(fibroblast_cells, slot = "data")[
                    available_genes, , drop = FALSE
                ]
                fibroblast_cells@meta.data[[paste0(function_name, "_Score")]] <- 
                    colMeans(function_data, na.rm = TRUE)
            })
        } else {
            fibroblast_cells@meta.data[[paste0(function_name, "_Score")]] <- 0
        }
    }
    
    # Calculate summary statistics
    caf_function_summary <- fibroblast_cells@meta.data %>%
        group_by(.data[[treatment_col]]) %>%
        summarise(
            n_cells = n(),
            ECM_Remodeling_mean = mean(ECM_Remodeling_Score, na.rm = TRUE),
            ECM_Remodeling_se = sd(ECM_Remodeling_Score, na.rm = TRUE) / sqrt(n()),
            Inflammatory_mean = mean(Inflammatory_Score, na.rm = TRUE),
            Inflammatory_se = sd(Inflammatory_Score, na.rm = TRUE) / sqrt(n()),
            .groups = "drop"
        )
    
    cat("CAF functional scores by treatment:\n")
    print(caf_function_summary)
    
    # ECM remodeling violin plot
    p_ecm <- ggplot(fibroblast_cells@meta.data, 
                    aes_string(x = treatment_col, y = "ECM_Remodeling_Score", 
                               fill = treatment_col)) +
        geom_violin(alpha = 0.75, width = 0.8, color = "white", linewidth = 1, trim = FALSE) +
        geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.9, 
                     color = "#2C3E50", linewidth = 0.8) +
        stat_summary(fun = median, geom = "point", size = 2.5, color = "#2C3E50") +
        scale_fill_manual(values = treatment_colors) +
        labs(title = "CAF ECM Remodeling Capacity",
             subtitle = "Extracellular matrix remodeling functional score",
             x = "Treatment Group", y = "ECM Remodeling Score") +
        guides(fill = "none") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold"))
    
    # Inflammatory violin plot
    p_inflammatory <- ggplot(fibroblast_cells@meta.data, 
                             aes_string(x = treatment_col, y = "Inflammatory_Score", 
                                        fill = treatment_col)) +
        geom_violin(alpha = 0.75, width = 0.8, color = "white", linewidth = 1, trim = FALSE) +
        geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.9, 
                     color = "#2C3E50", linewidth = 0.8) +
        stat_summary(fun = median, geom = "point", size = 2.5, color = "#2C3E50") +
        scale_fill_manual(values = treatment_colors) +
        labs(title = "CAF Inflammatory Capacity",
             subtitle = "Inflammatory cytokine production functional score",
             x = "Treatment Group", y = "Inflammatory Score") +
        guides(fill = "none") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold"))
    
    # Scatter plot showing correlation
    p_scatter <- ggplot(fibroblast_cells@meta.data, 
                        aes_string(x = "ECM_Remodeling_Score", y = "Inflammatory_Score", 
                                   color = treatment_col)) +
        geom_point(alpha = 0.6, size = 1.5) +
        scale_color_manual(values = treatment_colors) +
        stat_ellipse(aes_string(fill = treatment_col), alpha = 0.2, geom = "polygon") +
        scale_fill_manual(values = treatment_colors) +
        labs(title = "CAF Functional Score Correlation",
             subtitle = "ECM remodeling vs inflammatory capacity",
             x = "ECM Remodeling Score", y = "Inflammatory Score",
             color = "Treatment", fill = "Treatment") +
        theme(legend.position = "right")
    
    return(list(
        ecm_plot = p_ecm,
        inflammatory_plot = p_inflammatory,
        scatter_plot = p_scatter,
        summary = caf_function_summary,
        data = fibroblast_cells
    ))
}

# ==============================================================================
# Analysis 4: CXCL13+ CAF Proportion
# ==============================================================================

analyze_cxcl13_caf <- function(fibroblast_cells, treatment_col = "Treatment") {
    
    cat("\n=== Analysis: CXCL13+ CAF Proportion ===\n")
    
    if(!"Cxcl13" %in% rownames(fibroblast_cells)) {
        cat("CXCL13 gene not found\n")
        return(NULL)
    }
    
    # Get CXCL13 expression
    suppressWarnings({
        cxcl13_expr <- GetAssayData(fibroblast_cells, slot = "data")["Cxcl13", ]
        fibroblast_cells@meta.data$Cxcl13_Expression <- cxcl13_expr
        fibroblast_cells@meta.data$Cxcl13_positive <- cxcl13_expr > 0
    })
    
    # Calculate proportions
    cxcl13_summary <- fibroblast_cells@meta.data %>%
        group_by(.data[[treatment_col]]) %>%
        summarise(
            total_caf = n(),
            cxcl13_count = sum(Cxcl13_positive, na.rm = TRUE),
            cxcl13_proportion = cxcl13_count / total_caf * 100,
            .groups = "drop"
        )
    
    cat("CXCL13+ CAF proportions:\n")
    print(cxcl13_summary)
    
    # Create plot
    p <- ggplot(cxcl13_summary, aes_string(x = treatment_col, y = "cxcl13_proportion", 
                                            fill = treatment_col)) +
        geom_col(alpha = 0.8, width = 0.7, color = "white", linewidth = 1.2) +
        scale_fill_manual(values = treatment_colors) +
        labs(title = "CXCL13+ CAF Proportion Analysis",
             subtitle = "Chemokine-expressing fibroblasts across treatment groups",
             x = "Treatment Group", y = "CXCL13+ CAF Proportion (%)") +
        guides(fill = "none") +
        geom_text(aes(label = paste0(round(cxcl13_proportion, 1), "%\n(", 
                                     cxcl13_count, "/", total_caf, ")")), 
                  vjust = -0.3, size = 4.5, fontface = "bold", color = "#2C3E50") +
        coord_cartesian(ylim = c(0, max(cxcl13_summary$cxcl13_proportion) * 1.2)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold"))
    
    return(list(plot = p, summary = cxcl13_summary))
}

# ==============================================================================
# Main Analysis Pipeline
# ==============================================================================

run_caf_analysis_pipeline <- function(seurat_obj, treatment_col = "Treatment", 
                                      output_prefix = "caf") {
    
    cat("=== Starting CAF Analysis Pipeline ===\n")
    
    # Extract fibroblasts
    fibroblast_cells <- subset(seurat_obj, subset = Celltype_major_man == "Fibroblast")
    
    if(ncol(fibroblast_cells) == 0) {
        cat("No fibroblast cells found!\n")
        return(NULL)
    }
    
    cat("Number of CAF cells:", ncol(fibroblast_cells), "\n")
    
    # Calculate CAF scores if not present
    if(!"CAF_Subtype" %in% colnames(fibroblast_cells@meta.data)) {
        fibroblast_cells <- calculate_caf_marker_scores(fibroblast_cells)
    }
    
    # Run analyses
    matcaf_results <- analyze_matcaf_proportion(fibroblast_cells, treatment_col)
    il1r1_results <- analyze_il1r1_expression(fibroblast_cells)
    functional_results <- analyze_caf_functional_scores(fibroblast_cells, treatment_col)
    cxcl13_results <- analyze_cxcl13_caf(fibroblast_cells, treatment_col)
    
    # Save individual plots
    if(!is.null(matcaf_results)) {
        ggsave(paste0(output_prefix, "_matcaf_proportion.pdf"), 
               matcaf_results$plot, width = 4, height = 7, dpi = 300)
        write.csv(matcaf_results$summary, 
                  paste0(output_prefix, "_matcaf_stats.csv"), row.names = FALSE)
    }
    
    if(!is.null(il1r1_results)) {
        ggsave(paste0(output_prefix, "_il1r1_violin.pdf"), 
               il1r1_results$violin_plot, width = 4, height = 7, dpi = 300)
        ggsave(paste0(output_prefix, "_il1r1_proportion.pdf"), 
               il1r1_results$rate_plot, width = 4, height = 7, dpi = 300)
        write.csv(il1r1_results$summary, 
                  paste0(output_prefix, "_il1r1_stats.csv"), row.names = FALSE)
    }
    
    if(!is.null(functional_results)) {
        ggsave(paste0(output_prefix, "_ecm_score.pdf"), 
               functional_results$ecm_plot, width = 4, height = 7, dpi = 300)
        ggsave(paste0(output_prefix, "_inflammatory_score.pdf"), 
               functional_results$inflammatory_plot, width = 4, height = 7, dpi = 300)
        ggsave(paste0(output_prefix, "_functional_correlation.pdf"), 
               functional_results$scatter_plot, width = 7, height = 7, dpi = 300)
        write.csv(functional_results$summary, 
                  paste0(output_prefix, "_functional_stats.csv"), row.names = FALSE)
    }
    
    if(!is.null(cxcl13_results)) {
        ggsave(paste0(output_prefix, "_cxcl13_proportion.pdf"), 
               cxcl13_results$plot, width = 4, height = 7, dpi = 300)
        write.csv(cxcl13_results$summary, 
                  paste0(output_prefix, "_cxcl13_stats.csv"), row.names = FALSE)
    }
    
    # Create comprehensive plot
    if(!is.null(matcaf_results) && !is.null(il1r1_results) && 
       !is.null(functional_results) && !is.null(cxcl13_results)) {
        
        comprehensive_plot <- (matcaf_results$plot + cxcl13_results$plot) / 
                              (il1r1_results$rate_plot + functional_results$ecm_plot) / 
                              functional_results$inflammatory_plot +
            plot_layout(heights = c(1, 1, 0.8)) +
            plot_annotation(
                title = "Comprehensive CAF Analysis",
                subtitle = "CAF heterogeneity and functional characterization",
                theme = theme(
                    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
                    plot.subtitle = element_text(size = 16, hjust = 0.5)
                )
            )
        
        ggsave(paste0(output_prefix, "_comprehensive_analysis.pdf"), 
               comprehensive_plot, width = 8, height = 21, dpi = 300)
    }
    
    cat("\n=== CAF Analysis Complete ===\n")
    
    return(list(
        matcaf = matcaf_results,
        il1r1 = il1r1_results,
        functional = functional_results,
        cxcl13 = cxcl13_results
    ))
}

# ==============================================================================
# Example Usage
# ==============================================================================

# Uncomment and modify to run:

# # Load data
# seurat_obj <- readRDS("processed_seurat_object.rds")
# 
# # Run CAF analysis
# results <- run_caf_analysis_pipeline(seurat_obj, 
#                                      treatment_col = "Treatment",
#                                      output_prefix = "caf")

cat("Script 03: CAF Analysis - Functions loaded successfully!\n")
