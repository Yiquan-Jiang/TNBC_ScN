#!/usr/bin/env Rscript

# ====================================================================
# Breast Cancer Spatial Transcriptomics - CellChat Analysis
# ====================================================================
# Cell-cell communication analysis using CellChat for:
# - IL1R1+ Fibroblasts (CAF) and CXCL13+ T cells
# - Ligand-receptor interaction identification
# - Signaling pathway analysis
#
# Author: Research Team
# Date: 2024
# Version: 1.0
# ====================================================================

cat("====================================================================\n")
cat("CellChat Analysis: IL1R1+ CAF and CXCL13+ T Cell Communication\n")
cat("====================================================================\n\n")

# =========================================================================
# Load Required Libraries
# =========================================================================

cat("Loading required libraries...\n")

suppressMessages({
    library(CellChat)
    library(Seurat)
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(ComplexHeatmap)
    library(circlize)
    library(future)
})

# Configure parallel computing and memory
options(future.globals.maxSize = 12 * 1024^3)  # 12 GB
plan("multisession", workers = 2)

# =========================================================================
# Setup Output Directory
# =========================================================================

output_dir <- "cellchat_results"
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# =========================================================================
# Step 1: Load and Preprocess Data
# =========================================================================

cat("\n=== Step 1: Loading Seurat Data ===\n")

rds_file <- "BRCA_GSE176078_TS.rds"
if (!file.exists(rds_file)) {
    stop("Error: RDS file not found")
}

cat("Loading RDS file...\n")
seurat_obj <- readRDS(rds_file)
cat(sprintf("✓ Loaded data: %d cells, %d genes\n", 
           ncol(seurat_obj), nrow(seurat_obj)))

# Data overview
cat("Data overview:\n")
print(seurat_obj)

# =========================================================================
# Step 2: Cell Type Annotation
# =========================================================================

cat("\n=== Step 2: Cell Type Annotation ===\n")

# Check for existing cell type annotations
if ("cell_type" %in% colnames(seurat_obj@meta.data)) {
    cat("Using existing cell type annotations\n")
    cell_types <- seurat_obj@meta.data$cell_type
} else if ("seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
    cat("Using Seurat clusters as cell types\n")
    cell_types <- paste0("Cluster_", seurat_obj@meta.data$seurat_clusters)
} else {
    cat("Performing basic clustering...\n")
    plan("sequential")
    
    seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
    seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
    seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
    seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
    seurat_obj <- FindNeighbors(seurat_obj, verbose = FALSE)
    seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, verbose = FALSE)
    cell_types <- paste0("Cluster_", seurat_obj@meta.data$seurat_clusters)
    
    plan("multisession", workers = 2)
}

seurat_obj@meta.data$cell_type <- cell_types

cat("Cell type distribution:\n")
print(table(cell_types))

# =========================================================================
# Step 3: Identify Target Cell Populations
# =========================================================================

cat("\n=== Step 3: Identifying Target Cell Populations ===\n")

# Check gene availability
genes_to_check <- c("IL1R1", "CXCL13")
available_genes <- genes_to_check[genes_to_check %in% rownames(seurat_obj)]
cat(sprintf("Available genes: %s\n", paste(available_genes, collapse = ", ")))

if (length(available_genes) < 2) {
    stop("Error: IL1R1 or CXCL13 not found in data")
}

# Get expression matrix
expr_matrix <- GetAssayData(seurat_obj, slot = "data")

# Identify fibroblasts
fibroblast_patterns <- c("Fibroblast", "CAF", "Stromal", "Mesenchymal")
fibroblast_cells <- which(grepl(paste(fibroblast_patterns, collapse = "|"), 
                               cell_types, ignore.case = TRUE))

if (length(fibroblast_cells) == 0) {
    cat("Warning: No explicit fibroblast annotation found\n")
    fibroblast_cells <- 1:ncol(seurat_obj)
}

# Identify T cells
t_cell_patterns <- c("T cell", "T-cell", "Tcell", "CD4", "CD8", "Helper", "Regulatory")
t_cells <- which(grepl(paste(t_cell_patterns, collapse = "|"), 
                      cell_types, ignore.case = TRUE))

if (length(t_cells) == 0) {
    cat("Warning: No explicit T cell annotation found\n")
    t_cells <- 1:ncol(seurat_obj)
}

# Identify IL1R1+ fibroblasts
il1r1_expr <- expr_matrix["IL1R1", ]
il1r1_threshold <- quantile(il1r1_expr[fibroblast_cells], 0.75, na.rm = TRUE)
il1r1_pos_fibroblasts <- intersect(fibroblast_cells, 
                                  which(il1r1_expr > il1r1_threshold))

# Identify CXCL13+ T cells
cxcl13_expr <- expr_matrix["CXCL13", ]
cxcl13_threshold <- quantile(cxcl13_expr[t_cells], 0.75, na.rm = TRUE)
cxcl13_pos_t_cells <- intersect(t_cells, 
                               which(cxcl13_expr > cxcl13_threshold))

cat("Identified cell populations:\n")
cat(sprintf("- Total fibroblasts: %d\n", length(fibroblast_cells)))
cat(sprintf("- IL1R1+ fibroblasts: %d\n", length(il1r1_pos_fibroblasts)))
cat(sprintf("- Total T cells: %d\n", length(t_cells)))
cat(sprintf("- CXCL13+ T cells: %d\n", length(cxcl13_pos_t_cells)))

# Create CellChat labels
cellchat_labels <- as.character(cell_types)
cellchat_labels[il1r1_pos_fibroblasts] <- "IL1R1_pos_Fibroblasts"
cellchat_labels[cxcl13_pos_t_cells] <- "CXCL13_pos_T_cells"

# Filter to target cells
target_cells <- unique(c(il1r1_pos_fibroblasts, cxcl13_pos_t_cells))
if (length(target_cells) < 50) {
    stop("Error: Too few target cells for CellChat analysis")
}

cat(sprintf("Final analysis cells: %d\n", length(target_cells)))

# =========================================================================
# Step 4: Create CellChat Object
# =========================================================================

cat("\n=== Step 4: Creating CellChat Object ===\n")

seurat_subset <- seurat_obj[, target_cells]
expr_data <- GetAssayData(seurat_subset, slot = "data")

target_labels <- cellchat_labels[target_cells]
meta_data <- data.frame(
    labels = target_labels,
    row.names = colnames(seurat_subset)
)

cat(sprintf("Expression matrix: %d cells\n", ncol(expr_data)))
cat(sprintf("Metadata rows: %d\n", nrow(meta_data)))
cat("Label distribution:\n")
print(table(target_labels))

# Create CellChat object
cellchat <- createCellChat(object = expr_data, meta = meta_data, group.by = "labels")

# Load CellChatDB database
cat("Loading CellChatDB database...\n")
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)

CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use

# =========================================================================
# Step 5: Data Preprocessing
# =========================================================================

cat("\n=== Step 5: Data Preprocessing ===\n")

cellchat <- subsetData(cellchat)
cat(sprintf("After filtering: %d cells, %d genes\n", 
           length(cellchat@idents), nrow(cellchat@data)))

# Identify overexpressed genes and interactions
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# =========================================================================
# Step 6: Cell Communication Analysis
# =========================================================================

cat("\n=== Step 6: Cell Communication Analysis ===\n")

cat("Computing communication probability...\n")
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)

# Filter low probability communications
cellchat <- filterCommunication(cellchat, min.cells = 10)

cat("Inferring cell communication network...\n")
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# =========================================================================
# Step 7: Visualization
# =========================================================================

cat("\n=== Step 7: Generating Visualizations ===\n")

options(repr.plot.width=12, repr.plot.height=8)

# 1. Network overview
cat("Generating network overview...\n")
pdf(file.path(output_dir, "cellchat_network_overview.pdf"), width = 12, height = 8)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                weight.scale = TRUE, label.edge = FALSE, 
                title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                weight.scale = TRUE, label.edge = FALSE, 
                title.name = "Interaction weights/strength")

dev.off()

# 2. Signaling pathway heatmap
cat("Generating pathway heatmap...\n")
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

pdf(file.path(output_dir, "cellchat_pathway_heatmap.pdf"), width = 12, height = 10)
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
dev.off()

# 3. Top signaling pathways
pathways.show <- cellchat@netP$pathways[1:min(5, length(cellchat@netP$pathways))]
if (length(pathways.show) > 0) {
    cat(sprintf("Analyzing top %d signaling pathways...\n", length(pathways.show)))
    
    for (i in 1:length(pathways.show)) {
        pathway <- pathways.show[i]
        cat(sprintf("  Pathway: %s\n", pathway))
        
        # Pathway network
        pdf(file.path(output_dir, sprintf("cellchat_pathway_%s.pdf", pathway)), 
            width = 10, height = 8)
        netVisual_aggregate(cellchat, signaling = pathway, layout = "circle")
        dev.off()
        
        # Chord diagram
        pdf(file.path(output_dir, sprintf("cellchat_chord_%s.pdf", pathway)), 
            width = 8, height = 8)
        netVisual_chord_gene(cellchat, signaling = pathway)
        dev.off()
    }
}

# 4. Ligand-receptor pairs
cat("Generating ligand-receptor pair analysis...\n")
if (length(pathways.show) > 0) {
    pdf(file.path(output_dir, "cellchat_lr_pairs.pdf"), width = 12, height = 10)
    netAnalysis_contribution(cellchat, signaling = pathways.show[1])
    dev.off()
}

# =========================================================================
# Step 8: Save Results
# =========================================================================

cat("\n=== Step 8: Saving Results ===\n")

# Save CellChat object
saveRDS(cellchat, file = file.path(output_dir, "cellchat_object.rds"))

# Communication summary
comm_summary <- data.frame(
    source = cellchat@net$count %>% rownames(),
    target = cellchat@net$count %>% colnames(),
    count = as.vector(cellchat@net$count),
    weight = as.vector(cellchat@net$weight)
)
write.csv(comm_summary, file.path(output_dir, "communication_summary.csv"), row.names = FALSE)

# Signaling pathways
if (length(cellchat@netP$pathways) > 0) {
    pathway_df <- data.frame(
        pathway = cellchat@netP$pathways,
        prob = sapply(cellchat@netP$pathways, function(x) {
            if (x %in% names(cellchat@netP$prob)) {
                sum(cellchat@netP$prob[[x]], na.rm = TRUE)
            } else {
                0
            }
        })
    )
    write.csv(pathway_df, file.path(output_dir, "signaling_pathways.csv"), row.names = FALSE)
}

# Ligand-receptor pairs
if (length(pathways.show) > 0) {
    tryCatch({
        lr_pairs <- extractEnrichedLR(cellchat, signaling = pathways.show[1:min(3, length(pathways.show))])
        if (nrow(lr_pairs) > 0) {
            write.csv(lr_pairs, file.path(output_dir, "ligand_receptor_pairs.csv"), row.names = FALSE)
        }
    }, error = function(e) {
        cat("Warning extracting L-R pairs:", e$message, "\n")
    })
}

# Generate analysis report
report_file <- file.path(output_dir, "cellchat_analysis_report.txt")
cat("CellChat Analysis Report\n", file = report_file)
cat("========================\n\n", file = report_file, append = TRUE)
cat(sprintf("Analysis time: %s\n", Sys.time()), file = report_file, append = TRUE)
cat(sprintf("Data file: %s\n", rds_file), file = report_file, append = TRUE)
cat(sprintf("Cells analyzed: %d\n", length(cellchat@idents)), file = report_file, append = TRUE)
cat(sprintf("Cell types: %s\n", paste(unique(cellchat@idents), collapse = ", ")), 
    file = report_file, append = TRUE)
cat(sprintf("Pathways detected: %d\n", length(cellchat@netP$pathways)), 
    file = report_file, append = TRUE)
cat(sprintf("Top pathways: %s\n", paste(pathways.show, collapse = ", ")), 
    file = report_file, append = TRUE)

# =========================================================================
# Completion
# =========================================================================

cat("\n====================================================================\n")
cat("✓ CellChat Analysis Complete!\n")
cat("====================================================================\n\n")
cat(sprintf("Results saved in: %s\n", output_dir))
cat("Output files:\n")
cat("  - cellchat_object.rds: Complete CellChat object\n")
cat("  - communication_summary.csv: Cell communication summary\n")
cat("  - signaling_pathways.csv: Signaling pathway analysis\n")
cat("  - ligand_receptor_pairs.csv: Ligand-receptor pairs\n")
cat("  - cellchat_*.pdf: Visualization figures\n")
cat("  - cellchat_analysis_report.txt: Analysis report\n")

