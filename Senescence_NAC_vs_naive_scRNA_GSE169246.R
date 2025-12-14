####----env----####
library(limma)
library(Seurat)
library(dplyr)
library(tidyr)
library(data.table)
library(magrittr)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer)
library(SummarizedExperiment)
library(scuttle)
#library(metap)
#library(lessR)
library(ggrepel)
library(patchwork)
#library(scCancer)
library(extrafont)
library(Cairo)
library(showtext)
library(ggrepel)
library(ggpubr)
library(cowplot)
library(ggExtra)
library(scales)
library(cowplot)
library(ggpubr)
library(ggthemes)
#library(monocle3)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
#library(scCancer)
library(extrafont)
library(Cairo)
library(showtext)
library(ggrepel)
library(ggpubr)
library(cowplot)
library(ggExtra)
library(scales)
library(cowplot)
library("ggplot2")
library("ggsci")
library(scDblFinder)
library(broom)

library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(Azimuth)
library(ggplot2)
library(patchwork)

options(future.globals.maxSize = 30000000000)

allcolour=c("grey","grey","#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#9400D3","#7B68EE","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
allcolour1=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
             "#808000","#FF00FF","#FA8072","#9400D3","#7B68EE","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
             "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
             "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")



####----data input and integration----####
#set9 breast cancer PD-1 or chemotherapy
setwd("C:/tgscRNA/ICI scRNA/pd1_breast_cancer")
setwd("~/ICB_scRNA/pd1_breast_cancer")
#cohort2
`1867-counts_cells_cohort2` <- readRDS("1867-counts_cells_cohort2.rds")
set9Seurat <- CreateSeuratObject(counts = `1867-counts_cells_cohort2`, project = "set9")
set9Seurat[["percent.mt"]] <- PercentageFeatureSet(object = set9Seurat, pattern = "^MT-")


library(readr)
set9metadata <- read_csv("1871-BIOKEY_metaData_cohort2_web.csv")
cell_cycle_metadata_pd1_breast_cancer_scrna <- read_csv("cell_cycle_metadata_pd1_breast_cancer_scrna.csv")
cell_cycle_metadata_pd1_breast_cancer_scrna<-cell_cycle_metadata_pd1_breast_cancer_scrna[,-1]
colnames(cell_cycle_metadata_pd1_breast_cancer_scrna)[1]<-"cell.id"
colnames(set9metadata)[1]<-"cell.id"
set9metadata<-dplyr::left_join(set9metadata,cell_cycle_metadata_pd1_breast_cancer_scrna,by="cell.id")
cell.id<-colnames(set9Seurat)
cell.id<-as.data.frame(cell.id)
cell.id$set<-"set9"
metadata9<-dplyr::left_join(cell.id,set9metadata,by="cell.id")
set9Seurat[["patient"]]<-metadata9$patient_id
set9Seurat[["origin_celltype"]]<-metadata9$cellType
set9Seurat[["Treatment.Status"]]<-metadata9$timepoint
set9Seurat[["treatment"]]<-metadata9$cohort
set9Seurat[["set"]]<-metadata9$set
set9Seurat[["cell_phase"]]<-metadata9$assigned.phases
set9Seurat[["G1_score"]]<-metadata9$G1
set9Seurat[["S_score"]]<-metadata9$S
set9Seurat[["G2M_score"]]<-metadata9$G2M

set9Seurat<-subset(set9Seurat,subset=Treatment.Status=="Pre")
set9Seurat <- NormalizeData(set9Seurat, normalization.method = "LogNormalize", scale.factor = 10000)

rm(`1867-counts_cells_cohort2`)
gc()

#cohort1
`1863-counts_cells_cohort1` <- readRDS("1863-counts_cells_cohort1.rds")
set9Seurat2 <- CreateSeuratObject(counts = `1863-counts_cells_cohort1`, project = "set9")
set9Seurat2[["percent.mt"]] <- PercentageFeatureSet(object = set9Seurat2, pattern = "^MT-")

set9metadata2 <- read_csv("1872-BIOKEY_metaData_cohort1_web.csv")

cell_cycle_metadata_pd1_breast_cancer_scrna_cohort1 <- read_csv("cell_cycle_metadata_pd1_breast_cancer_scrna_cohort1.csv")
cell_cycle_metadata_pd1_breast_cancer_scrna_cohort1<-cell_cycle_metadata_pd1_breast_cancer_scrna_cohort1[,-1]
colnames(cell_cycle_metadata_pd1_breast_cancer_scrna_cohort1)[1]<-"cell.id"
colnames(set9metadata2)[1]<-"cell.id"
set9metadata2<-dplyr::left_join(set9metadata2,cell_cycle_metadata_pd1_breast_cancer_scrna_cohort1,by="cell.id")

cell.id<-colnames(set9Seurat2)
cell.id<-as.data.frame(cell.id)
cell.id$set<-"set9"
metadata9_1<-dplyr::left_join(cell.id,set9metadata2,by="cell.id")

set9Seurat2[["patient"]]<-metadata9_1$patient_id
set9Seurat2[["origin_celltype"]]<-metadata9_1$cellType
set9Seurat2[["Treatment.Status"]]<-metadata9_1$timepoint
set9Seurat2[["treatment"]]<-metadata9_1$cohort
set9Seurat2[["set"]]<-metadata9_1$set

set9Seurat2[["cell_phase"]]<-metadata9_1$assigned.phases
set9Seurat2[["G1_score"]]<-metadata9_1$G1
set9Seurat2[["S_score"]]<-metadata9_1$S
set9Seurat2[["G2M_score"]]<-metadata9_1$G2M

rm(`1863-counts_cells_cohort1`)
gc()
set9Seurat2<-subset(set9Seurat2,subset=Treatment.Status=="Pre")

####----cell cycle estimation----####
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("scran")
library(scran)
library(Seurat)
library(SingleCellExperiment)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
##
#‘scuttle’, ‘BiocParallel’, ‘edgeR’, ‘limma’, ‘BiocSingular’, ‘bluster’
#cohort2
`1867-counts_cells_cohort2` <- readRDS("C:/tgscRNA/ICI scRNA/pd1_breast_cancer/1867-counts_cells_cohort2.rds")
sce <- SingleCellExperiment(list(counts=`1867-counts_cells_cohort2`)) 

# hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
str(mm.pairs)
head(mm.pairs$G1)

ensembl <- mapIds(org.Hs.eg.db, keys=rownames(sce), keytype="SYMBOL", column="ENSEMBL")
mm.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
assigned <- cyclone(sce, pairs=mm.pairs, gene.names=ensembl)
head(assigned$scores)
table(assigned$phases)
plot(assigned$score$G1, assigned$score$G2M, 
     xlab="G1 score", ylab="G2/M score", pch=16)

metadata_pd1_breast_cancer_scrna<-data.frame(colnames(`1867-counts_cells_cohort2`),assigned$phases,assigned$normalized.scores)
write.csv(metadata_pd1_breast_cancer_scrna,file = "cell_cycle_metadata_pd1_breast_cancer_scrna.csv")

####----integration----####

set9Seurat <- SCTransform(set9Seurat,  vars.to.regress = c("percent.mt","nCount_RNA"), verbose = FALSE,vst.flavor = "v2")
set9Seurat2 <- SCTransform(set9Seurat2,  vars.to.regress = c("percent.mt","nCount_RNA"), verbose = FALSE,vst.flavor = "v2")


##integration using RPCA
pancreas.list <- list(set9Seurat2,set9Seurat)
features <- SelectIntegrationFeatures(object.list = pancreas.list, nfeatures = 3000)
pancreas.list <- PrepSCTIntegration(object.list = pancreas.list, anchor.features = features)
pancreas.list <- lapply(X = pancreas.list, FUN = RunPCA, features = features)

immune.anchors <- FindIntegrationAnchors(object.list = pancreas.list, normalization.method = "SCT",
                                         anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)
pancreas.integrated_RPCA <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT", dims = 1:30)

pancreas.integrated_RPCA <- RunPCA(pancreas.integrated_RPCA, verbose = FALSE)
pancreas.integrated_RPCA <- RunUMAP(pancreas.integrated_RPCA, reduction = "pca", dims = 1:30)
pancreas.integrated_RPCA <- RunTSNE(pancreas.integrated_RPCA, dims.use = 1:50, verbose=TRUE)


pancreas.integrated_RPCA <- FindNeighbors(pancreas.integrated_RPCA, dims = 1:30, verbose = FALSE)
pancreas.integrated_RPCA <- FindClusters(pancreas.integrated_RPCA, resolution=0.5,verbose = TRUE)

rm(pancreas.integrated)
gc()
#plot for a look
p1 <- DimPlot(pancreas.integrated_RPCA, reduction = "umap", group.by = "treatment")
p2 <- DimPlot(pancreas.integrated_RPCA, reduction = "umap", group.by = "origin_celltype", label = TRUE,
              repel = TRUE)
p1 + p2

#correct b cells idents
select.cells <- CellSelector(plot = p2)
Idents(pancreas.integrated_RPCA)<-"origin_celltype"
Idents(pancreas.integrated_RPCA, cells = select.cells) <- "Unknown"
pancreas.integrated_RPCA[["origin_celltype"]] <- Idents(object = pancreas.integrated_RPCA)

#rename ident of treatment
Idents(pancreas.integrated_RPCA)<-"treatment"
pancreas.integrated_RPCA <- RenameIdents(object = pancreas.integrated_RPCA, `neoadjuvant_chemo` = "NAC")
pancreas.integrated_RPCA <- RenameIdents(object = pancreas.integrated_RPCA, `treatment_naive` = "Ctrl")
pancreas.integrated_RPCA[["treatment"]] <- Idents(object = pancreas.integrated_RPCA)

#rename ident of ScN
Idents(pancreas.integrated)<-"ScN"
pancreas.integrated <- RenameIdents(object = pancreas.integrated, `Non_ScN` = "NonScN")
pancreas.integrated[["ScN"]] <- Idents(object = pancreas.integrated)




####----annotation using scATOMIC----####
library(scATOMIC)
library(plyr)
library(dplyr)
library(data.table)
library(randomForest)
library(caret)
library(parallel)
library(reticulate)
library(Rmagic)
library(Matrix)
library(Seurat)
library(agrmt)
library(cutoff.scATOMIC)
library(copykat)
library(ggplot2)

#where input_matrix is any non sparse count matrix
sparse_matrix <- as(as.matrix(input_matrix), "sparseMatrix")


annotations_scATOMIC<-NULL
#where seurat_object is any scRNA-seq Seurat object

pancreas.integrated_RPCA_BIOKEY_9 <- subset(pancreas.integrated_RPCA,subset = patient %in% c("BIOKEY_9"))
pancreas.integrated_RPCA_BIOKEY_36 <- subset(pancreas.integrated_RPCA,subset = patient %in% c("BIOKEY_36"))

sparse_matrix <- pancreas.integrated_RPCA_BIOKEY_36@assays$RNA@counts
cell_predictions <- run_scATOMIC(sparse_matrix, breast_mode = T, mc.cores =1)
results_lung_CNV <- create_summary_matrix(prediction_list = cell_predictions, use_CNVs = T, modify_results = T, mc.cores = 1, raw_counts = sparse_matrix, min_prop = 0.5)

# Load required libraries
library(foreach)
library(doParallel)

# Register the number of cores to use for parallel computing
no_cores <- parallel::detectCores() - 10  # Leave one core free
registerDoParallel(cores = no_cores)

# Use foreach for parallel computation
foreach(i = unique(pancreas.integrated_RPCA$patient)) %dopar% {
  patient1 <- subset(pancreas.integrated_RPCA, subset = patient == i)
  sparse_matrix <- patient1@assays$RNA@counts
  cell_predictions <- run_scATOMIC(sparse_matrix, breast_mode = F, mc.cores =1,confidence_cutoff = FALSE)
  results_lung_CNV <- create_summary_matrix(prediction_list = cell_predictions, use_CNVs = T, modify_results = T, mc.cores = 1, raw_counts = sparse_matrix, min_prop = 0.5,confidence_cutoff = FALSE,breast_mode = F)
  # Save the result as an RDS file
  saveRDS(results_lung_CNV, file = paste0(i, "_ATOMIC_results.rds"))
}


patient1<-subset(pancreas.integrated_RPCA,subset=patient=="BIOKEY_13")
sparse_matrix <- patient1@assays$RNA@counts

cell_predictions <- run_scATOMIC(sparse_matrix,breast_mode = T,mc.cores = 10)

results_lung_CNV <- create_summary_matrix(prediction_list = cell_predictions, use_CNVs = T, modify_results = T, mc.cores = 10, raw_counts = sparse_matrix, min_prop = 0.5 )


# Get a list of all the RDS files in your working directory
setwd("~/pd1_breast_cancer/ATOMIC_1")
rds_files <- list.files(pattern = "*.rds")

# Read them all into a list
list_of_dfs <- lapply(rds_files, readRDS)

# Find the common columns
common_cols <- Reduce(intersect, lapply(list_of_dfs, colnames))

# Subset each data frame to only include the common columns
list_of_dfs <- lapply(list_of_dfs, function(df) df[, common_cols, drop = FALSE])

# Bind all the data frames together
combined_df <- do.call(rbind, list_of_dfs)

# add metadata
colnames(combined_df)[4]<-"cell.id"

cell.id<-colnames(pancreas.integrated_RPCA)
cell.id<-as.data.frame(cell.id)

metadata<-dplyr::left_join(cell.id,combined_df,by="cell.id")
pancreas.integrated_RPCA[["celltype_ATOMIC"]]<-metadata$scATOMIC_pred
pancreas.integrated_RPCA[["celltype_ATOMIC_layer2"]]<-metadata$layer_2
pancreas.integrated_RPCA[["celltype_ATOMIC_layer3"]]<-metadata$layer_3
pancreas.integrated_RPCA[["celltype_ATOMIC_s_score"]]<-metadata$S.Score
pancreas.integrated_RPCA[["celltype_ATOMIC_g2m_score"]]<-metadata$G2M.Score

DimPlot(pancreas.integrated_RPCA,group.by = "celltype_ATOMIC",cols = allcolour1,label = T,pt.size = 1)

Idents(pancreas.integrated_RPCA)<-"celltype_ATOMIC"
pancreas.integrated_RPCA<-PrepSCTFindMarkers(pancreas.integrated_RPCA)

ATOMIC_annot_celltype_markers<-FindAllMarkers(pancreas.integrated_RPCA,logfc.threshold = 0.5,only.pos = T,verbose = T)




####----DEG analyses across conditions(chemo vs. naive)----####
#Identify differential expressed genes across conditions(chemo vs. naive), with integration using RPCA
Idents(pancreas.integrated_RPCA)<-"origin_celltype"
pancreas.integrated_RPCA$celltype.stim <- paste(Idents(pancreas.integrated_RPCA), pancreas.integrated_RPCA$treatment, sep = "_")
pancreas.integrated_RPCA$celltype <- Idents(pancreas.integrated_RPCA)
Idents(pancreas.integrated_RPCA) <- "celltype.stim"

DefaultAssay(pancreas.integrated_RPCA)<-"SCT"

#using SCT assay
pancreas.integrated_RPCA <- PrepSCTFindMarkers(pancreas.integrated_RPCA)
DEGChemovsNaivetotal_rpca <- NULL
for (i in unique(pancreas.integrated_RPCA@meta.data[["origin_celltype"]])){
  DEGChemovsNaive<-FindMarkers(pancreas.integrated_RPCA,ident.2 = paste0(i,"_Ctrl"),ident.1 = paste0(i,"_NAC"), verbose = TRUE,min.pct = 0.05,logfc.threshold = 0.05,only.pos = F)
  colnames(DEGChemovsNaive)[1]<-"p.value"
  colnames(DEGChemovsNaive)[5]<-"q.value"
  colnames(DEGChemovsNaive)[2]<-"logFC"
  DEGChemovsNaive$cell_type<-i
  DEGChemovsNaive$gene<-row.names(DEGChemovsNaive)
  DEGChemovsNaivetotal_rpca<-rbind(DEGChemovsNaivetotal_rpca,DEGChemovsNaive)
}
write.csv(DEGChemovsNaivetotal_rpca,file = "DEGChemovsNaivetotal_rpca.csv")

#using RNA assay with normalize
pancreas.integrated_RPCA<-NormalizeData(pancreas.integrated_RPCA)
DEGChemovsNaivetotal_rpca_with_normalizeRNA <- NULL
for (i in unique(pancreas.integrated_RPCA@meta.data[["celltype"]])){
  DEGChemovsNaive<-FindMarkers(pancreas.integrated_RPCA,ident.1 = paste0(i,"_treatment_naive"),ident.2 = paste0(i,"_neoadjuvant_chemo"), verbose = TRUE,min.pct = 0.05,logfc.threshold = log(1.5))
  colnames(DEGChemovsNaive)[1]<-"p.value"
  colnames(DEGChemovsNaive)[5]<-"q.value"
  colnames(DEGChemovsNaive)[2]<-"logFC"
  DEGChemovsNaive$cell_type<-i
  DEGChemovsNaive$gene<-row.names(DEGChemovsNaive)
  DEGChemovsNaivetotal_rpca_with_normalizeRNA<-rbind(DEGChemovsNaivetotal_rpca_with_normalizeRNA,DEGChemovsNaive)
}


cell_age_gene<-read.csv(file = "cell_aging_genes.csv")
DEGChemovsNaivetotal_RPCA_senescence<-DEGChemovsNaivetotal[DEGChemovsNaivetotal$gene %in% cell_age_gene$gene_name,]

##Identify differential expressed genes across conditions(chemo vs. naive), with integration using CCA
Idents(pancreas.integrated)<-"origin_celltype"
pancreas.integrated$celltype.stim <- paste(Idents(pancreas.integrated), pancreas.integrated$treatment, sep = "_")
pancreas.integrated$celltype <- Idents(pancreas.integrated)
Idents(pancreas.integrated) <- "celltype.stim"

DefaultAssay(pancreas.integrated)<-"RNA"

#using RNA assay without normalize
DEGChemovsNaivetotal_cca <- NULL
for (i in unique(pancreas.integrated@meta.data[["celltype"]])){
  DEGChemovsNaive<-FindMarkers(pancreas.integrated,ident.2 = paste0(i,"_Ctrl"),ident.1 = paste0(i,"_NAC"), verbose = TRUE,only.pos = F,min.pct = 0.05,logfc.threshold = 0)
  colnames(DEGChemovsNaive)[1]<-"p.value"
  colnames(DEGChemovsNaive)[5]<-"q.value"
  colnames(DEGChemovsNaive)[2]<-"logFC"
  DEGChemovsNaive$cell_type<-i
  DEGChemovsNaive$gene<-row.names(DEGChemovsNaive)
  DEGChemovsNaive$logFC[DEGChemovsNaive$q.value>0.05]<-0.0001
  DEGChemovsNaivetotal_cca<-rbind(DEGChemovsNaivetotal_cca,DEGChemovsNaive)
}

#using RNA assay with normalize
DEGChemovsNaivetotal_cca_RNAnormalized <- NULL
for (i in unique(pancreas.integrated@meta.data[["celltype"]])){
  DEGChemovsNaive<-FindMarkers(pancreas.integrated,ident.2 = paste0(i,"_treatment_naive"),ident.1 = paste0(i,"_neoadjuvant_chemo"), verbose = TRUE,only.pos = F,min.pct = 0.05,logfc.threshold = 0,assay = "RNA")
  colnames(DEGChemovsNaive)[1]<-"p.value"
  colnames(DEGChemovsNaive)[5]<-"q.value"
  colnames(DEGChemovsNaive)[2]<-"logFC"
  DEGChemovsNaive$cell_type<-i
  DEGChemovsNaive$gene<-row.names(DEGChemovsNaive)
  DEGChemovsNaivetotal_cca_RNAnormalized<-rbind(DEGChemovsNaivetotal_cca_RNAnormalized,DEGChemovsNaive)
}

#using SCT assay 
DEGChemovsNaivetotal_cca_SCT <- NULL
for (i in unique(pancreas.integrated@meta.data[["celltype"]])){
  DEGChemovsNaive<-FindMarkers(pancreas.integrated,ident.2 = paste0(i,"_treatment_naive"),ident.1 = paste0(i,"_neoadjuvant_chemo"), verbose = TRUE,only.pos = F,min.pct = 0.05,logfc.threshold = 0,assay = "SCT")
  colnames(DEGChemovsNaive)[1]<-"p.value"
  colnames(DEGChemovsNaive)[5]<-"q.value"
  colnames(DEGChemovsNaive)[2]<-"logFC"
  DEGChemovsNaive$cell_type<-i
  DEGChemovsNaive$gene<-row.names(DEGChemovsNaive)
  DEGChemovsNaivetotal_cca_SCT<-rbind(DEGChemovsNaivetotal_cca_SCT,DEGChemovsNaive)
}

cell_age_gene<-read.csv(file = "cell_aging_genes.csv")
DEGChemovsNaivetotal_cca_senescence<-DEGChemovsNaivetotal_cca[DEGChemovsNaivetotal_cca$gene %in% cell_age_gene$gene_name,]
write.csv(DEGChemovsNaivetotal_cca_senescence,file = "scRNA_DEG_ChemovsNaive_senescence.csv")


#different integration methods visualize
p1 <- DimPlot(pancreas.integrated, reduction = "umap", label = TRUE)
p2 <- DimPlot(pancreas.integrated_RPCA, reduction = "umap", label = TRUE,repel = TRUE)
p1 + p2

#Venn plot for different cell type senescence genes
vennlist<-list('Fibroblast'=as.character(DEGChemovsNaivetotal_cca_senescence$gene[DEGChemovsNaivetotal_cca_senescence$cell_type=="Fibroblast"]),
               'Endothelial_cell'=as.character(DEGChemovsNaivetotal_cca_senescence$gene[DEGChemovsNaivetotal_cca_senescence$cell_type=="Endothelial_cell"]),
               'pDC'=as.character(DEGChemovsNaivetotal_cca_senescence$gene[DEGChemovsNaivetotal_cca_senescence$cell_type=="pDC"]),
               'T_cell'=as.character(DEGChemovsNaivetotal_cca_senescence$gene[DEGChemovsNaivetotal_cca_senescence$cell_type=="T_cell"]),
               'Myeloid_cell'=as.character(DEGChemovsNaivetotal_cca_senescence$gene[DEGChemovsNaivetotal_cca_senescence$cell_type=="Myeloid_cell"]),
               'B_cell'=as.character(DEGChemovsNaivetotal_cca_senescence$gene[DEGChemovsNaivetotal_cca_senescence$cell_type=="B_cell"]),
               'Cancer_cell'=as.character(DEGChemovsNaivetotal_cca_senescence$gene[DEGChemovsNaivetotal_cca_senescence$cell_type=="Cancer_cell"]))
install.packages("ggvenn") # install via CRAN
library(ggvenn)
ggvenn(vennlist,c("Fibroblast","Endothelial_cell","Cancer_cell","T_cell"), show_elements = TRUE,text_size = 1)
ggvenn(vennlist,c("pDC","Myeloid_cell","Cancer_cell","B_cell"), show_elements = TRUE,text_size = 1)


####----Identify DEGs among cell types----####
Idents(pancreas.integrated)<-"celltype"
cell_type_markers<-FindAllMarkers(pancreas.integrated,verbose = T,only.pos = T)
DotPlot(pancreas.integrated,features = c("FCGR3A","S100A8","S100A9","CD14",
                                         "C1QA","C1QB","C1QC","ITGAX","CD163",
                                         "CD206","TGFB1","VEGFA","VEGFB","ARG1",
                                         "CD19","CD20","CD79A","SDC1","XBP1",
                                         "FAP","ACTA2","CXCL12",
                                         "CD3D","CD8A","CD4",
                                         "FOXP3",
                                         "NKG7","KLRB1",
                                         "HAVCR2","PDCD1"),group.by = "RNA_snn_res.0.5")
#rename ident of celltype
Idents(pancreas.integrated)<-"RNA_snn_res.0.5"
pancreas.integrated <- RenameIdents(object = pancreas.integrated, `6` = "Myeloid_Pro_Tumor")
pancreas.integrated <- RenameIdents(object = pancreas.integrated, `7` = "Myeloid")
pancreas.integrated <- RenameIdents(object = pancreas.integrated, `14` = "Neutrophil")
pancreas.integrated <- RenameIdents(object = pancreas.integrated, `1` = "Cancer_cell")
pancreas.integrated <- RenameIdents(object = pancreas.integrated, `11` = "Cancer_cell")
pancreas.integrated <- RenameIdents(object = pancreas.integrated, `19` = "Cancer_cell")
pancreas.integrated <- RenameIdents(object = pancreas.integrated, `3` = "Cancer_cell")
pancreas.integrated <- RenameIdents(object = pancreas.integrated, `18` = "Cancer_cell")
pancreas.integrated <- RenameIdents(object = pancreas.integrated, `8` = "Cancer_cell")
pancreas.integrated <- RenameIdents(object = pancreas.integrated, `17` = "Cancer_cell")
pancreas.integrated <- RenameIdents(object = pancreas.integrated, `2` = "Endothelial")
pancreas.integrated <- RenameIdents(object = pancreas.integrated, `4` = "iCAF")
pancreas.integrated <- RenameIdents(object = pancreas.integrated, `10` = "mCAF")
pancreas.integrated <- RenameIdents(object = pancreas.integrated, `20` = "pDC")
pancreas.integrated <- RenameIdents(object = pancreas.integrated, `15` = "Treg")
pancreas.integrated <- RenameIdents(object = pancreas.integrated, `16` = "CD8_Tex")
pancreas.integrated <- RenameIdents(object = pancreas.integrated, `5` = "CD8_Tem")
pancreas.integrated <- RenameIdents(object = pancreas.integrated, `0` = "T_naive")
pancreas.integrated <- RenameIdents(object = pancreas.integrated, `13` = "NK_cell")
pancreas.integrated <- RenameIdents(object = pancreas.integrated, `12` = "B_cell")
pancreas.integrated <- RenameIdents(object = pancreas.integrated, `9` = "undefined")


pancreas.integrated[["celltype_for_cibersort"]] <- Idents(object = pancreas.integrated)

pancreas.integrated_cibersort<-subset(pancreas.integrated,subset=celltype_for_cibersort!="undefined")
Idents(pancreas.integrated_cibersort)<-"celltype_for_cibersort"
pancreas.integrated_cibersort <- FindVariableFeatures(pancreas.integrated_cibersort, selection.method = "vst", nfeatures = 500)

table(pancreas.integrated_cibersort@meta.data[["celltype_for_cibersort"]])
Idents(pancreas.integrated_cibersort) <- "celltype_for_cibersort"
DefaultAssay(pancreas.integrated_cibersort)<-"RNA"
cluster.averages <- AverageExpression(pancreas.integrated_cibersort)
cluster.averages_sct<-cluster.averages[["RNA"]]
cluster.averages_sct<-as.data.frame(cluster.averages_sct)
cluster.averages_sct$Gene<-row.names(cluster.averages_sct)
cluster.averages_sct<-select(cluster.averages_sct,15,1:14,everything())
cluster.averages_sct<-cluster.averages_sct[cluster.averages_sct$Gene %in% VariableFeatures(pancreas.integrated_cibersort),]
write.table(cluster.averages_sct,file = "scRNA_breastcancer_celltype_expression.txt",sep = "\t",row.names = F)

#myeloid
Myeloid_subset<-subset(pancreas.integrated,subset=celltype=="Myeloid_cell")
Myeloid_subset <- FindNeighbors(Myeloid_subset, dims = 1:30, verbose = FALSE)
Myeloid_subset <- FindClusters(Myeloid_subset, resolution=0.05,verbose = TRUE)
Myeloid_subset_markers<-FindAllMarkers(Myeloid_subset,only.pos = T,verbose = T)
DotPlot(Myeloid_subset,features = c("FCGR3A","S100A8","S100A9","CD14","C1QA","C1QB","C1QC","ITGAX"))


####----Senescence annotation and DEG analyses----####
##define senescence cells in each cell type
ScN_gene_status<-FetchData(pancreas.integrated_RPCA,vars = c("celltype","CDKN2A","CDKN1A","cell_phase"))


ScN_gene_status$celltype_ScN[ScN_gene_status$CDKN2A>0 & ScN_gene_status$cell_phase=="G1"]<-"ScN"
ScN_gene_status$celltype_ScN[ScN_gene_status$CDKN1A>0 & ScN_gene_status$cell_phase=="G1"]<-"ScN"
ScN_gene_status$celltype_ScN[is.na(ScN_gene_status$celltype_ScN)]<-"Non_ScN"

pancreas.integrated_RPCA[["ScN"]]<-ScN_gene_status$celltype_ScN

DimPlot(pancreas.integrated_RPCA,group.by = "ScN",label = T)

##
ScN_gene_status$celltype_ScN[ScN_gene_status$CDKN2A>quantile(ScN_gene_status$CDKN2A,0.67)&
                               ScN_gene_status$CDKN1A>quantile(ScN_gene_status$CDKN1A,0.67)]<-"ScN"
ScN_gene_status$celltype_ScN[is.na(ScN_gene_status$celltype_ScN)]<-"Non_ScN"

pancreas.integrated[["ScN"]]<-ScN_gene_status$celltype_ScN

#senescence genes expression output
ScN_gene_status<-FetchData(set9Seurat,vars = c("celltype","CDKN2A","CDKN1A","cell_phase","G1_score","S_score","G2M_score"))



##Identify differential expressed genes across conditions(ScN vs. Non_ScN), with integration using CCA
Idents(pancreas.integrated_RPCA)<-"celltype"
pancreas.integrated_RPCA$celltype_ScN <- paste(Idents(pancreas.integrated_RPCA), pancreas.integrated_RPCA$ScN, sep = "_")
Idents(pancreas.integrated_RPCA) <- "celltype_ScN"

DefaultAssay(pancreas.integrated_RPCA)<-"SCT"

DEG_ScN_vs_NonScN_total <- NULL
for (i in unique(pancreas.integrated_RPCA@meta.data[["celltype"]])){
  DEG_ScN_vs_NonScN<-FindMarkers(pancreas.integrated_RPCA,ident.1 = paste0(i,"_ScN"),ident.2 = paste0(i,"_Non_ScN"), verbose = TRUE,only.pos = F,min.pct = 0.05,logfc.threshold = 0.02)
  colnames(DEG_ScN_vs_NonScN)[1]<-"p.value"
  colnames(DEG_ScN_vs_NonScN)[5]<-"q.value"
  colnames(DEG_ScN_vs_NonScN)[2]<-"logFC"
  DEG_ScN_vs_NonScN$cell_type<-i
  DEG_ScN_vs_NonScN$gene<-row.names(DEG_ScN_vs_NonScN)
  DEG_ScN_vs_NonScN_total<-rbind(DEG_ScN_vs_NonScN_total,DEG_ScN_vs_NonScN)
}
write.csv(DEG_ScN_vs_NonScN_total,file = "DEG_ScN_vs_NonScN_total.csv")

#ScN vs. Non_ScN Myeloid cells DEG
Myeloid_ScN_vs_NonScN_DEG<-FindMarkers(pancreas.integrated,ident.1 = "Myeloid_cell_ScN",ident.2 ="Myeloid_cell_Non_ScN",verbose = TRUE,only.pos = F,min.pct = 0.05,logfc.threshold = 0 )
Myeloid_ScN_vs_NonScN_DEG$gene<-row.names(Myeloid_ScN_vs_NonScN_DEG)

#subset myeloid and visualize
VlnPlot(pancreas.integrated,idents = c("Myeloid_cell_ScN","Myeloid_cell_Non_ScN"),
        features = c("CDKN1A","CDKN2A","IL6","CXCL1","CXCL2"),pt.size = 0,stack = T,flip = T,assay = "RNA")
RidgePlot(pancreas.integrated,idents = c("Myeloid_cell_ScN","Myeloid_cell_Non_ScN"),
          features = c("CDKN1A","CDKN2A","IL6","CXCL1","CXCL2"),log = T)
DotPlot(pancreas.integrated,idents = c("Myeloid_cell_ScN","Myeloid_cell_Non_ScN"),
        features = c("CDKN1A","CDKN2A","IL6","CXCL1","CXCL2"),assay ="RNA",dot.scale = 20)

myeloid_subset<-subset(pancreas.integrated,subset=celltype=="Myeloid_cell")

FeatureScatter(myeloid_subset, feature1 = "CDKN2A", feature2 = "IL6")

Myeloid_ScN_vs_NonScN_SASP<-FetchData(myeloid_subset,vars = c("celltype_ScN","CDKN1A","CDKN2A","IL6","CXCL1","CXCL2"))

library(ggplot2)

ggplot(Myeloid_ScN_vs_NonScN_SASP) +
 aes(x = celltype_ScN, y = IL6, fill = celltype_ScN) +
 geom_violin(adjust = 4.6, 
 scale = "area") +
 scale_fill_brewer(palette = "Pastel1", direction = 1) +
 scale_y_continuous(trans = "log1p") +
 theme_classic() +
 theme(legend.position = "top", plot.caption = element_text(size = 12L), axis.title.y = element_text(size = 16L, 
 face = "bold"), axis.title.x = element_text(size = 17L, face = "bold"))

##gsea analyses for ScN vs. Non_ScN myeloid
library(msigdbr) #?ṩMSigdb???ݿ???????
library(fgsea)
library(data.table)
library(tibble)

mdb_h<-msigdbr(species = "Homo sapiens", category = "H")
mdb_c2 <- msigdbr(species = "Homo sapiens", category = "C2")
mdb_c5 <- msigdbr(species = "Homo sapiens", category = "C5")
mdb_kegg = mdb_c2 [grep("^KEGG",mdb_c2 $gs_name),]
mdb_GO =mdb_c5 [grep("^GO",mdb_c2 $gs_name),]

fgsea_setsKEGG<- mdb_kegg %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_setsGO<- mdb_GO %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_setsH<- mdb_h %>% split(x = .$gene_symbol, f = .$gs_name)

#??????logFC????
Myeloid_ScN_vs_NonScN_DEG_adj<-Myeloid_ScN_vs_NonScN_DEG[Myeloid_ScN_vs_NonScN_DEG$p_val_adj<0.01,]
Myeloid_ScN_vs_NonScN_DEG_adj$genes = rownames(Myeloid_ScN_vs_NonScN_DEG_adj)

cluster0.genes<- Myeloid_ScN_vs_NonScN_DEG_adj %>% arrange(desc(avg_log2FC)) %>% dplyr::select(genes,avg_log2FC)
ranks<- deframe(cluster0.genes)

#fgsea
fgseaResKEGG<- fgsea(fgsea_setsKEGG, stats = ranks, nperm = 1000)
fgseaResGO<- fgsea(fgsea_setsGO, stats = ranks, nperm = 1000)
fgseaResH<- fgsea(fgsea_setsH, stats = ranks, nperm = 1000)
fgseaRes<-rbind(fgseaResKEGG,fgseaResGO)

fgseaRes<-as.data.frame(fgseaRes)
fgseaRes$annotation<-as.character(fgseaRes$leadingEdge)
write.csv(fgseaRes[,c(1:7,9)],file = "Myeloid_ScN_vs_NonScN_DEG_GSVA.csv")

library(Hmisc)
fgseaRes$pathway<-tolower(fgseaRes$pathway)
fgseaRes$pathway<-capitalize(fgseaRes$pathway)


ggplot(fgseaRes[fgseaRes$padj<0.05,],aes(x=NES,y= reorder(pathway,NES)))+
  geom_point(aes(size=-1*log10(padj),color=NES))+
  scale_color_gradient(low="gray",high = "red")+ #
  labs(size="-1*log10(adj.P.Val)",  ##expression???????庯????ʽ []?????±꣬^?????ϱ?
       x="Normalized Enrichment Score",      ##?Զ???????
       y="Pathway name",
       title="Pathway enrichment")+ ##?Զ?????????
  theme_bw() #ȥ?����?




####----MDSC annotation----####
MDSC_markers<-FetchData(pancreas.integrated_RPCA,vars = c("celltype","CD84","CD33","ITGAM","HLA-DRB1","HLA-DRA","ARG1","treatment"))
quantile(MDSC_markers$`CD84`,0.9)
quantile(MDSC_markers$`HLA-DRB1`,0.5)

MDSC_markers$celltype<-as.character(MDSC_markers$celltype)
MDSC_markers$celltype[MDSC_markers$CD84>1 &
                            MDSC_markers$celltype=="Myeloid_cell"]<-"MDSC"
mdsc_detail<-MDSC_markers[MDSC_markers$celltype=="MDSC",]


pancreas.integrated_RPCA[["celltype_with_mdsc"]]<-MDSC_markers$celltype

##calculate the MDSC proportion using bulk-seq data
table(pancreas.integrated_RPCA@meta.data[["celltype_with_mdsc"]])
Idents(pancreas.integrated_RPCA) <- "celltype_with_mdsc"
DefaultAssay(pancreas.integrated_RPCA)<-"RNA"
cluster.averages <- AverageExpression(pancreas.integrated_RPCA)
cluster.averages_sct<-cluster.averages[["RNA"]]
cluster.averages_sct<-as.data.frame(cluster.averages_sct)
cluster.averages_sct$Gene<-row.names(cluster.averages_sct)
cluster.averages_sct<-select(cluster.averages_sct,11,1:10,everything())
DefaultAssay(pancreas.integrated_RPCA)<-"RNA"
pancreas.integrated_RPCA <- FindVariableFeatures(pancreas.integrated_RPCA, selection.method = "vst", nfeatures = 2000)
cluster.averages_sct<-cluster.averages_sct[cluster.averages_sct$Gene %in% VariableFeatures(pancreas.integrated_RPCA),]
write.table(cluster.averages_sct,file = "scRNA_breastcancer_MDSC_for_cibersort.txt",sep = "\t",row.names = F)
write.csv(cluster.averages_sct,file = "scRNA_breastcancer_ScN_celltype_expression.csv",row.names = T)


#mdsc proportion analysis
mdsc_metadata<-FetchData(pancreas.integrated,vars = c("celltype_with_mdsc","treatment","patient"))
library(reshape2)
mdsc_metadata_cast<-dcast(mdsc_metadata,treatment+patient~celltype_with_mdsc)
mdsc_metadata_cast$mdsc_myeloid_ratio<-mdsc_metadata_cast$MDSC/mdsc_metadata_cast$Myeloid_cell
mdsc_metadata_cast$mdsc_CD45_ratio<-mdsc_metadata_cast$MDSC/(mdsc_metadata_cast$Myeloid_cell+mdsc_metadata_cast$T_cell+mdsc_metadata_cast$B_cell+mdsc_metadata_cast$pDC)
mdsc_metadata_cast$myeloid_CD45_ratio<-mdsc_metadata_cast$Myeloid_cell/(mdsc_metadata_cast$Myeloid_cell+mdsc_metadata_cast$T_cell+mdsc_metadata_cast$B_cell+mdsc_metadata_cast$pDC)


####----plot for fund----####
pancreas.integrated<-AddMetaData(pancreas.integrated,pancreas.integrated@reductions[["umap"]]@cell.embeddings,col.name = colnames(pancreas.integrated@reductions[["umap"]]@cell.embeddings))
allcolour=c("#DC143C","#ffffb3","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#8dd3c7","#bc80bd","#bebada","#fb8072","#80b1d3","#fdb462",
            "#b3de69","#fccde5","#d9d9d9","#bc80bd","#FF6347","#6A5ACD","#9932CC","#8B008B","#7B68EE","#8B4513","#9400D3",
            "#DEB887","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C",
            "#FFFFE0","#EE82EE")
class_avg <- pancreas.integrated@meta.data %>%
  group_by(origin_celltype) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )
CairoPDF(file="UMAP.pdf",width=12.5,height=9) 
ggplot(pancreas.integrated@meta.data ,aes(x=UMAP_1,y=UMAP_2))+
  geom_point(size=1,aes(color=origin_celltype))+
  scale_color_manual(values = allcolour)+
  geom_text(aes(label = origin_celltype), data = class_avg)+
  theme(text=element_text(family="A",size=18)) +
  theme(panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(), axis.title = element_text(color='black',
                                                              family="A",size=18),axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(1,"lines"),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=18),
        axis.title.y=element_text(colour='black', size=18),
        axis.text=element_text(colour='black',size=18),
        legend.title=element_blank(),
        legend.text=element_text(family="A", size=18),
        legend.key=element_blank())+
  theme(plot.title = element_text(size=22,colour = "black",family = "A",face = "bold"))  + 
  guides(colour = guide_legend(override.aes = list(size=6)))

dev.off()

DefaultAssay(pancreas.integrated)<-"RNA"
Idents(pancreas.integrated_RPCA)<-"origin_celltype"

#To prove chemotherapy would induce senescence
p1<-DimPlot(pancreas.integrated,group.by = "celltype_modified_by_selector",shape.by = "treatment",cols = allcolour,label = T,pt.size = 1)
p2<-FeaturePlot(pancreas.integrated,split.by = "treatment",features = c("CDKN1A","CDKN2A"),raster = T,min.cutoff = "q10", max.cutoff = "q60",label = T)
p3<-DotPlot(pancreas.integrated_RPCA,features = c("CDKN1A","CDKN2A"),split.by = "treatment",dot.scale = 8,cols = c("gray","red"))

#To prove myeloid SASP
mdsc_chemokines<-FetchData(pancreas.integrated,vars =c("patient","celltype","CXCL1","IL6",
                                                       "CXCL2", "CCL2","CXCL12") )
mdsc_chemokines<-melt(mdsc_chemokines,id=c("patient","celltype"))

mdsc_chemokines<-cast(mdsc_chemokines,celltype~variable,mean)

write.csv(mdsc_chemokines,file = "mdsc_recruit_chemokines_bulk.csv")
Idents(pancreas.integrated)<-"origin_celltype"
p4<-VlnPlot(pancreas.integrated,features = c("CXCL1","IL6",
                                             "CXCL2", "CCL2","CXCL12"),
            pt.size = 0,stack = T,flip = T)

p5<-DotPlot(pancreas.integrated,features = c("CXCL1","IL6",
                                         "CXCL2",  "CCL2","CSF1R","CXCL10","CXCL12", "CCL8", "CCL4", "TNF",
                                         "CCL4L2", "IL1B", "CCL3", "CXCL9",
                                         "CCL20", "CCL5", "CXCL11", "CCL3L1", "IFNGR2",
                                          "OSM", "CCR1", 
                                         "TNFSF13B", "TNFSF10", "IL2RG", "IL22RA2", 
                                         "TNFRSF9", "CXCL16", "CD40", "TNFSF13", 
                                         "LTBR", "IL15RA", "IL7R", "TGFB1", "IL1A",
                                         "IL21R", "IFNAR2", "TNFRSF1A", "RELT", 
                                         "IL10RA", "CSF2RB", "IL17RA", "CXCL3"),
        split.by = "ScN",dot.scale = 8,cols = c("grey","red"))



select.cells <- CellSelector(plot = p1)
pancreas.integrated[["celltype_modified_by_selector"]]<-Idents(pancreas.integrated)
Idents(pancreas.integrated)<-"celltype_modified_by_selector"
Idents(pancreas.integrated, cells = select.cells) <- "Cancer_cell"


#To prove co-expression of senescence genes and SASP chemokines
FeatureScatter(myeloid_subset, feature1 = "CDKN2A", feature2 = "CXCL2")


####----calculate the ScN cell proportion using bulk-seq data----####
table(pancreas.integrated@meta.data[["celltype_ScN"]])
Idents(pancreas.integrated) <- "celltype_ScN"
DefaultAssay(pancreas.integrated)<-"RNA"
cluster.averages <- AverageExpression(pancreas.integrated)
cluster.averages_sct<-cluster.averages[["RNA"]]
cluster.averages_sct<-as.data.frame(cluster.averages_sct)
cluster.averages_sct$Gene<-row.names(cluster.averages_sct)
cluster.averages_sct<-select(cluster.averages_sct,17,1:16,everything())
cluster.averages_sct<-cluster.averages_sct[cluster.averages_sct$Gene %in% VariableFeatures(pancreas.integrated),]
write.table(cluster.averages_sct,file = "scRNA_breastcancer_ScN_celltype_expression.txt",sep = "\t",row.names = F)
write.csv(cluster.averages_sct,file = "scRNA_breastcancer_ScN_celltype_expression.csv",row.names = T)



####----ScN vs. Non_ScN immune_modulators genes----####
immune_modulators<-c("ADORA2A","ARG1","BTLA","BTN3A1","BTN3A2","CCL5","CD27","CD274","CD276","CD28",
                     "CD40","CD40LG","CD70","CD80","CTLA4","CX3CL1","CXCL10","CXCL9","EDNRB","ENTPD1",
                     "GZMA","HAVCR2","HLA-A","HLA-B","HLA-C","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQA2",
                     "HLA-DQB1","HLA-DQB2","HLA-DRA","HLA-DRB1","HLA-DRB3","HLA-DRB4","HLA-DRB5","HMGB1",
                     "ICAM1","ICOS","ICOSLG","IDO1","IFNA1","IFNA2","IFNG","IL10","IL12A","IL13","IL1A",
                     "IL1B","IL2","IL2RA","IL4","ITGB2","KIR2DL1","KIR2DL2","KIR2DL3","LAG3","MICA","MICB",
                     "PDCD1","PDCD1LG2","PRF1","SELP","SLAMF7","TGFB1","TIGIT","TLR4","TNF","TNFRSF14",
                     "TNFRSF18","TNFRSF4","TNFRSF9","TNFSF4","TNFSF9","VEGFA","VEGFB","C10orf54","VTCN1")
DEG_ScN_vs_NonScN_total_immune_modulators1<-DEG_ScN_vs_NonScN_total[DEG_ScN_vs_NonScN_total$gene %in% immune_modulators,]

write.csv(DEG_ScN_vs_NonScN_total_immune_modulators1,file = "pd1bc_DEG_ScN_vs_NonScN_immune_modulators1.csv")

VlnPlot(setTNBC,features = c("CTLA4","HAVCR2","CD276","TGFB1","VEGFA","VEGFB","TIGIT","IL10"),group.by = "Major_celltype",
        split.by = "ScN",pt.size = 0)

####----CXCL13 T cells in ctrl vs. NAC
CXCL13_T<-FetchData(pancreas.integrated_RPCA,vars = c("celltype","CXCL13","patient","treatment"))
write.csv(CXCL13_T,file = "pd1bc_CXCL13_T_ctrl_NAC.csv")
CXCL13_T %>%
  filter(celltype %in% "T_cell") %>%
  ggplot() +
  aes(
    x = patient,
    y = CXCL13,
    fill = treatment
  ) +
  geom_violin(adjust = 2.6, scale = "width") +
  scale_fill_brewer(palette = "Set1", direction = 1) +
  scale_color_distiller(palette = "Set1", direction = 1) +
  theme_minimal() +
  facet_wrap(vars(treatment))


####----ScN proportion per patient----####
ScN_proportion<-FetchData(pancreas.integrated_RPCA,vars = c("celltype","celltype_ScN","patient","ScN","treatment"))
ScN_proportion_dcast<-reshape2::dcast(ScN_proportion[ScN_proportion$celltype=="Myeloid_cell",],treatment+patient~ScN)
ScN_proportion_dcast$ScN_ratio<-ScN_proportion_dcast$ScN/(ScN_proportion_dcast$ScN+ScN_proportion_dcast$Non_ScN)
wilcox.test(ScN_ratio~treatment,data=ScN_proportion_dcast)




####----pySCENIC----####
Idents(pancreas.integrated_RPCA)<-"celltype"
myeloid_subset<-subset(pancreas.integrated_RPCA,subset = celltype=="Myeloid_cell")
Idents(myeloid_subset)<-"treatment"
myeloid_subset_pd1_TNBC<-subset(myeloid_subset,downsample=2000)
myeloid_subset_pd1_TNBC@assays$integrated<-NULL
myeloid_subset_pd1_TNBC@assays$SCT<-NULL
DefaultAssay(myeloid_subset_pd1_TNBC)<-"RNA"
myeloid_subset_pd1_TNBC<-FindVariableFeatures(myeloid_subset_pd1_TNBC,nfeatures = 4000,selection.method = "vst")
myeloid_subset_pd1_TNBC@assays$RNA<-myeloid_subset_pd1_TNBC@assays$RNA[myeloid_subset_pd1_TNBC@assays$RNA@var.features,]
DEG_myeloid_subset_pd1_TNBC<-FindAllMarkers(myeloid_subset_pd1_TNBC,verbose = T,min.pct = 0.01)
saveRDS(myeloid_subset,file = "~/pyscenic/myeloid_subset_pd1_TNBC.rds")


library(igraph)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)

library(tidyverse)
library(igraph)
library(ggraph)
library(colormap)
library(wesanderson)
library(reshape2)


#regulonAUC plot
load("myeloid_subset_pd1_tnbc_regulon_RSS.Rdata")
AUCmatrix <- regulonAUC@assays@data@listData$AUC
AUCmatrix <- data.frame(t(AUCmatrix), check.names=F)
RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(','_',RegulonName_AUC)
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC)
colnames(AUCmatrix) <- RegulonName_AUC
scRNAauc <- AddMetaData(PRO, AUCmatrix)
scRNAauc@assays$integrated <- NULL
saveRDS(scRNAauc,'scRNAauc.rds')

#plot RSS
rss<-readRDS("myeloid_subset_pd1_tnbc_rss.rds")
rssPlot$plot

#regulonAUC data
dat <- tibble::enframe(regulons) %>% data.frame() %>% rename("regulon"=1,"targets"=2)  ### 将list转换成dataframe 导出
dat$targets <- vapply(dat$targets, paste, collapse = ", ", character(1L))

#network igraph
regul_network<-fread("myeloid_subset_pd1_tnbc_grn.tsv")
regul_network<-regul_network[regul_network$importance>50,]
regul_network_graph<-graph_from_data_frame(regul_network)
regul_network_graph

plot(regul_network_graph,#vertex.label="",
     rescale = T,
     seed=1,
     layout=layout_nicely, 
     edge.curved=0
)


dataUU <-fread("AdjacencyUndirectedUnweighted.txt",header=TRUE,check.names = FALSE)

connect<-fread("myeloid_subset_pd1_tnbc_grn.tsv")
tf_cdkn1a<-connect[connect$target=="CDKN1A",]
tf_cdkn1a<-tf_cdkn1a[tf_cdkn1a$importance>5,]
tf_cdkn1a<-as.character(tf_cdkn1a$TF)
connect<-connect[connect$TF %in% tf_cdkn1a,]
connect<-connect[connect$importance>15,]
colnames(connect)<-c("from","to","value")
# Transform the adjacency matrix in a long format
#connect <- dataUU %>% gather(key="to", value="value", -1) %>% mutate(to = gsub("\\.", " ",to)) %>% na.omit()
# Number of connection per person
c( as.character(connect$from), as.character(connect$to)) %>%
  as.tibble() %>%
  group_by(value) %>%
  summarize(n=n()) -> vertices
colnames(vertices) <- c("name", "n")

# Create a graph object with igraph
mygraph <- graph_from_data_frame( connect, vertices = vertices, directed = FALSE )

# Find community
com <- walktrap.community(mygraph)

#Reorder dataset and make the graph
vertices <- vertices %>% 
  mutate( group = com$membership) %>%
  mutate(group=as.numeric(factor(group,
                                 levels=sort(summary (as.factor(group)),index.return=TRUE,decreasing = T)$ix,
                                 order=TRUE)))%>%
  filter( group<10) %>%
  arrange(group,desc(n)) %>%
  mutate(name=factor(name, name))

# keep only this people in edges
connect <- connect %>%
  filter(from %in% vertices$name) %>%
  filter(to %in% vertices$name)%>%
  left_join(vertices,by=c('from'='name'))

# Create a graph object with igraph
mygraph <- graph_from_data_frame( connect, vertices = vertices, directed = FALSE )

mycolor <- wes_palette("Darjeeling1", max(vertices$group), type = "continuous")
mycolor <- sample(mycolor, length(mycolor))


#-----------------------------------MDS(多尺度分析布局)-------------------------------------------
ggraph(mygraph,layout='mds') + 
  geom_edge_link(edge_colour="black", edge_alpha=0.2, edge_width=0.3) +
  geom_node_point(aes(size=n, fill=as.factor(group)), shape=21,color='black',alpha=0.9) +
  scale_size_continuous(range=c(0.5,10)) +
  scale_fill_manual(values=mycolor) +
  geom_node_text(aes(label=ifelse(n>5, as.character(name), "")), size=3, color="black") +
  expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2))+
  theme_minimal() +
  theme(
    legend.position="none",
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.ticks =element_blank(),
    axis.text =element_blank(),
    axis.title = element_blank()
    #plot.margin=unit(c(0,0,0,0), "null"),
    #panel.spacing=unit(c(0,0,0,0), "null")
  )


#----------------------------------fr(Fruchterman-Reingold 算法的力导向布局)-----------------------
ggraph(mygraph,layout='fr') + 
  geom_edge_link(edge_colour="black", edge_alpha=0.2, edge_width=0.3) +
  geom_node_point(aes(size=n, fill=as.factor(group)), shape=21,color='black',alpha=0.9) +
  scale_size_continuous(range=c(0.5,10)) +
  scale_fill_manual(values=mycolor) +
  geom_node_text(aes(label=ifelse(n>5, as.character(name), "")), size=3, color="black") +
  geom_node_text(aes(label=ifelse(name=="Arg2", as.character(name), "")), size=3, color="black") +
  expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2))+
  theme_minimal() +
  theme(
    legend.position="none",
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.ticks =element_blank(),
    axis.text =element_blank(),
    axis.title = element_blank()
    #plot.margin=unit(c(0,0,0,0), "null"),
    #panel.spacing=unit(c(0,0,0,0), "null")
  )



##plot TFs regulating CDKN1A in NAC vs. Ctrl
devtools::install_github("Simon-Leonard/FlexDotPlot")
library(FlexDotPlot)
myeloid_NAC_subset<-subset(myeloid_subset,subset = treatment=="NAC")
dp<-DotPlot(myeloid_NAC_subset,features =c("ATF3","MDM2","BAX","SPIB","CD83","RELB","KDM2B"),group.by = "ScN" )
dot_plot(dp$data[,c(4,3,1,2,5)],shape_var = "pct.exp",col_var = "avg.exp.scaled",shape_legend = "Percent Expressed",
         col_legend = "Average Expression",x.lab.pos = "bottom",
         hclust_method = "ward.D2")

DEG_myeloid_subset_pd1_TNBC_NAC_vs_ctrl<-FindMarkers(myeloid_subset,ident.1 = "NAC",ident.2 = "Ctrl",only.pos = T,verbose = T,assay = "RNA",logfc.threshold = 0.2)
DEG_myeloid_subset_pd1_TNBC_NAC_vs_ctrl$gene<-row.names(DEG_myeloid_subset_pd1_TNBC_NAC_vs_ctrl)
write.csv(DEG_myeloid_subset_pd1_TNBC_NAC_vs_ctrl,file = "~/ICB_scRNA/pd1_breast_cancer/DEG_myeloid_subset_pd1_TNBC_NAC_vs_ctrl.csv")


####----myeloid subtypes ScN status----####
pancreas.integrated_RPCA_myeloid<-subset(pancreas.integrated_RPCA,subset=celltype_ATOMIC_layer3=="Macrophage or Monocyte")
pancreas.integrated_RPCA_myeloid<-subset(pancreas.integrated_RPCA_myeloid,subset=celltype_ATOMIC_layer3!="Macrophage")
pancreas.integrated_RPCA_myeloid<-subset(pancreas.integrated_RPCA_myeloid,subset=celltype_ATOMIC_layer3!="Normal Tissue Cell")
saveRDS(pancreas.integrated_RPCA_myeloid,file = "pancreas.integrated_RPCA_myeloid.RDS")

Idents(pancreas.integrated_RPCA_myeloid)<-"IL1B_annot"
pancreas.integrated_RPCA_myeloid$IL1B_annot_ScN <- paste(Idents(pancreas.integrated_RPCA_myeloid), pancreas.integrated_RPCA_myeloid$ScN, sep = "_")
Idents(pancreas.integrated_RPCA_myeloid) <- "IL1B_annot_ScN"

Idents(pancreas.integrated_RPCA_myeloid)<-"IL1B_annot"
pancreas.integrated_RPCA_myeloid$IL1B_annot_treat <- paste(Idents(pancreas.integrated_RPCA_myeloid), pancreas.integrated_RPCA_myeloid$treatment, sep = "_")
Idents(pancreas.integrated_RPCA_myeloid) <- "IL1B_annot_treat"

DefaultAssay(pancreas.integrated_RPCA_myeloid)<-"integrated"
pancreas.integrated_RPCA_myeloid <- RunPCA(pancreas.integrated_RPCA_myeloid, verbose = FALSE)
pancreas.integrated_RPCA_myeloid <- RunUMAP(pancreas.integrated_RPCA_myeloid, reduction = "pca", dims = 1:30)
pancreas.integrated_RPCA_myeloid <- RunTSNE(pancreas.integrated_RPCA_myeloid, dims.use = 1:50, verbose=TRUE)

set.seed(1234)
pancreas.integrated_RPCA_myeloid <- FindNeighbors(pancreas.integrated_RPCA_myeloid, dims = 1:30, verbose = FALSE)
pancreas.integrated_RPCA_myeloid <- FindClusters(pancreas.integrated_RPCA_myeloid, resolution=0.5,verbose = TRUE)

DimPlot(pancreas.integrated_RPCA_myeloid,group.by = c("seurat_clusters","celltype_ATOMIC"),label = T,cols = allcolour1)
DefaultAssay(pancreas.integrated_RPCA_myeloid)<-"SCT"
FeaturePlot(pancreas.integrated_RPCA_myeloid,features = c("CDKN1A","IL1B","RGS2"))

pancreas.integrated_RPCA_myeloid <- RenameIdents(
  object = pancreas.integrated_RPCA_myeloid,
  '5' = 'Pro_angio_infla_Mac',
  '6' = 'Pro_angio_infla_Mac',
  '4' = 'Pro_angio_infla_Mac',
  '9' = 'Pro_angio_infla_Mac',
  '2' = 'Pro_angio_infla_Mac',
  '10' = 'Pro_angio_infla_Mac',
  '7' = 'Phagocytosis_Mac',
  '3' = 'Phagocytosis_Mac',
  '1' = 'Phagocytosis_Mac',
  '6' = 'Phagocytosis_Mac',
  '11' = 'Phagocytosis_Mac',
  '0' = 'Pro_angio_infla_Mac',
  '8' = 'Phagocytosis_Mac'
)


pancreas.integrated_RPCA_myeloid[["IL1B_annot"]]<-Idents(pancreas.integrated_RPCA_myeloid)

pancreas.integrated_RPCA_myeloid<-PrepSCTFindMarkers(pancreas.integrated_RPCA_myeloid)
IL1B_annot_markers<-FindAllMarkers(pancreas.integrated_RPCA_myeloid,only.pos = T,logfc.threshold = 0.5,assay = "SCT")

##Find IL1B co-expression genes
# Extract the count matrix from the Seurat object
data_matrix <- pancreas.integrated_RPCA_myeloid@assays$RNA@data

# Get the expression values of IL1B
il1b_exp <- data_matrix["IL1B", ]

# Register a parallel backend for foreach
num_cores <- detectCores()
registerDoParallel(cores = num_cores-10)

# Calculate correlations in parallel
corr_values <- foreach(gene = rownames(data_matrix), .combine = 'c', .packages = 'stats') %dopar% {
  gene_exp <- data_matrix[gene, ]
  cor(il1b_exp, gene_exp, method = "pearson")
}

# Associate each correlation value with its respective gene
names(corr_values) <- rownames(data_matrix)

# Sort the genes by the absolute value of their correlation with IL1B
sorted_genes <- sort(corr_values, decreasing = TRUE, index.return = TRUE)

# Display the top 10 co-expressed genes with IL1B
print(names(sorted_genes$ix[1:10]))

rm(data_matrix)




DotPlot(pancreas.integrated_RPCA_myeloid,features = c("CDKN1A","IL1B","RGS2","ATF3","STAT1"),dot.scale = 8,cols = c("gray","red"),group.by = "IL1B_annot_ScN")
DotPlot(pancreas.integrated_RPCA_myeloid,features = c("CDKN1A","PDGFB","IL1A","IL1B","IL10","IL6"),dot.scale = 8,cols = c("gray","red"),group.by = "IL1B_annot_ScN")+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

DotPlot(pancreas.integrated_RPCA_myeloid,features = c("CDKN1A","PDGFB","IL1A","IL1B","IL10","IL6"),dot.scale = 8,cols = c("gray","red"),group.by = "IL1B_annot_treat")+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

VlnPlot(pancreas.integrated_RPCA_myeloid,features = c("s_score","g2m_score"),group.by = "celltype_ATOMIC",pt.size = 0)

#IL1B CDKN1A RGS2 coexpression
install.packages("scCustomize")
library(scCustomize) # 需要Seurat版本4.3.0
library(viridis)
library(RColorBrewer)
library(gridExtra)
#单基因
p000 <- Plot_Density_Custom(seurat_object = pancreas.integrated_RPCA_myeloid, features = "IL1B")
#双基因
p111 <- Plot_Density_Joint_Only(seurat_object = pancreas.integrated_RPCA_myeloid, 
                                features = c("IL1B", "CDKN1A"))
#多基因
p222 <- Plot_Density_Joint_Only(seurat_object = pancreas.integrated_RPCA_myeloid, 
                                features = c("IL1B", "CDKN1A","ATF3"),
                                custom_palette = BlueAndRed())
p000 + p111 + p222

#plot pro-inflam mac heatmap
library(tidyverse)
library(Scillus)
library(Seurat)
library(magrittr)
library(purrr)

plot_stat(pancreas.integrated_RPCA_myeloid, plot_type = "prop_multi", pal_setup = "Set3",group_by = "celltype_ATOMIC")

#plot ATF3 violin plot
library(dittoSeq)
library(dplyr)
pancreas.integrated_RPCA_myeloid_IL1B<-subset(pancreas.integrated_RPCA_myeloid,subset=IL1B_annot=="Pro_angio_infla_Mac")
dittoPlot(pancreas.integrated_RPCA_myeloid_IL1B,
          "ATF3", 
          group.by = "treatment",
          plots=c("vlnplot","boxplot"),
          boxplot.fill=F,
          boxplot.color='white',
          color.panel = dittoColors(),
          colors = c(1,3),
          theme = theme(axis.text = element_text(size = 12, color = 'black'),
                        axis.line = element_line(size = 1),
                        axis.title.y = element_text(size = 15, color = 'black'),
                        plot.title = element_text(size=15,hjust=0.5, color = 'black')),
          ylab = 'Expression',
          y.breaks = seq(0,8,1),
          xlab = '',
          x.labels = c("CTRL","NAC"),
          x.labels.rotate =F,
          max=3,
          min=0,
          main = "ATF3",
          legend.show = F)



####----pyscenic----####
DefaultAssay(pancreas.integrated_RPCA_myeloid)<-"RNA"
# Normalize the data
pancreas.integrated_RPCA_myeloid <- NormalizeData(pancreas.integrated_RPCA_myeloid)

# Identify the top 3000 highly variable genes
pancreas.integrated_RPCA_myeloid <- FindVariableFeatures(pancreas.integrated_RPCA_myeloid, nfeatures = 3000)

# Retain only the top 3000 highly variable genes
pancreas.integrated_RPCA_myeloid_pyscenic <- subset(pancreas.integrated_RPCA_myeloid, features = VariableFeatures(pancreas.integrated_RPCA_myeloid))

saveRDS(pancreas.integrated_RPCA_myeloid_pyscenic,file = "/home/shpc_101793/pyscenic/bc_pd1_scData_myeloid.RDS")

#run in linux
Rscript get_count_from_seurat.R -i bc_pd1_scData_myeloid.RDS -d treatment -s 2000 -l bc_pd1_scData_myeloid
python create_loom_file_from_scanpy.py -i bc_pd1_scData_myeloid.csv

sudo docker run -it --rm \
-v /home/shpc_101793/pyscenic:/pyscenic aertslab/pyscenic:0.12.1 \
pyscenic grn \
--num_workers 10 \
--method grnboost2 \
--output /pyscenic/bc_pd1_scData_myeloid_grn.csv \
--sparse \
/pyscenic/bc_pd1_scData_myeloid.loom \
/pyscenic/allTFs_hg38.txt

sudo docker run -it --rm \
-v /home/shpc_101793/pyscenic:/pyscenic aertslab/pyscenic:0.12.1 \
pyscenic ctx \
/pyscenic/grn.csv \
/pyscenic/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
--annotations_fname /pyscenic/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname /pyscenic/bc_pd1_scData_myeloid.loom \
--mode "custom_multiprocessing" \
--output /pyscenic/bc_pd1_scData_myeloid_regulons.csv \
--num_workers 10


sudo docker run -it --rm \
-v /home/shpc_101793/pyscenic:/pyscenic aertslab/pyscenic:0.12.1 \
pyscenic aucell \
/pyscenic/bc_pd1_scData_myeloid \
/pyscenic/bc_pd1_scData_myeloid_regulons.csv \
-o /pyscenic/bc_pd1_scData_myeloid_sample_SCENIC.loom \
--num_workers 10


library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)

