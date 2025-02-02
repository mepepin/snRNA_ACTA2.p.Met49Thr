---
title: "Single Nuclear Resolution of ACTA2 Mutation in Thoracic Ascending Aortic Tissue"
author: "Mark E. Pepin, MD, PhD, MS"
date: "02/01/2025"
output:
  html_document: 
    code_folding: hide
    keep_md: true
    toc: true
    toc_float: true
    fig_caption: true
  pdf_document:
    toc: yes
header-includes:
- \usepackage{booktabs}
- \usepackage{longtable}
- \usepackage{array}
- \usepackage{multirow}
- \usepackage[table]{xcolor}
- \usepackage{wrapfig}
- \usepackage{float}
- \usepackage{colortbl}
- \usepackage{pdflscape}
- \usepackage{tabu}
- \usepackage{threeparttable}
mainfont: Times
fontsize: 10pt
always_allow_html: yes
editor_options:
  markdown:
    wrap: 72
  chunk_output_type: inline
---



**Code Author**: Mark E. Pepin, MD, PhD, MS **Contact**:
[pepinme\@gmail.com](mailto:pepin@broadinstitute.org){.email}\
**Institution**: Brigham and Women's Hospital \| Harvard Medical School \| Stanford Cardiovascular Institute
\| Broad Institute of Harvard and MIT\
**Location**: Stanford, CA, USA

# Environment Setup


``` r
# Directories
ifelse(!dir.exists(file.path(paste0("../2_Output/"))), dir.create(file.path(paste0("../2_Output/"))), FALSE)
ifelse(!dir.exists(file.path(paste0("../2_Output/1_ACTA/"))), dir.create(file.path(paste0("../2_Output/1_ACTA/"))), FALSE)
ifelse(!dir.exists(file.path(paste0("../2_Output/1_ACTA/1_QC/"))), dir.create(file.path(paste0("../2_Output/1_ACTA/1_QC/"))), FALSE)
ifelse(!dir.exists(file.path(paste0("../2_Output/1_ACTA/2_Clustering/"))), dir.create(file.path(paste0("../2_Output/1_ACTA/2_Clustering/"))), FALSE)
ifelse(!dir.exists(file.path(paste0("../2_Output/1_ACTA/3_Differential.Expression/"))), dir.create(file.path(paste0("../2_Output/1_ACTA/3_Differential.Expression/"))), FALSE)
ifelse(!dir.exists(file.path(paste0("../2_Output/1_ACTA/4_Trajectory/"))), dir.create(file.path(paste0("../2_Output/1_ACTA/4_Trajectory/"))), FALSE)
ifelse(!dir.exists(file.path(paste0("../2_Output/2_TAA/"))), dir.create(file.path(paste0("../2_Output/2_TAA/"))), FALSE)
ifelse(!dir.exists(file.path(paste0("../2_Output/2_TAA/1_QC/"))), dir.create(file.path(paste0("../2_Output/2_TAA/1_QC/"))), FALSE)
ifelse(!dir.exists(file.path(paste0("../2_Output/2_TAA/2_Clustering/"))), dir.create(file.path(paste0("../2_Output/2_TAA/2_Clustering/"))), FALSE)
ifelse(!dir.exists(file.path(paste0("../2_Output/2_TAA/3_Differential.Expression/"))), dir.create(file.path(paste0("../2_Output/2_TAA/3_Differential.Expression/"))), FALSE)
ifelse(!dir.exists(file.path(paste0("../2_Output/2_TAA/4_Trajectory/"))), dir.create(file.path(paste0("../2_Output/2_TAA/4_Trajectory/"))), FALSE)
ifelse(!dir.exists(file.path(paste0("../2_Output/3_Merged/"))), dir.create(file.path(paste0("../2_Output/3_Merged/"))), FALSE)
ifelse(!dir.exists(file.path(paste0("../2_Output/3_Merged/1_QC/"))), dir.create(file.path(paste0("../2_Output/3_Merged/1_QC/"))), FALSE)
ifelse(!dir.exists(file.path(paste0("../2_Output/3_Merged/2_Clustering/"))), dir.create(file.path(paste0("../2_Output/3_Merged/2_Clustering/"))), FALSE)
ifelse(!dir.exists(file.path(paste0("../2_Output/3_Merged/3_Differential.Expression/"))), dir.create(file.path(paste0("../2_Output/3_Merged/3_Differential.Expression/"))), FALSE)
# Parallelization
library(future)
plan("multicore", workers = 8)
# set this option when analyzing large datasets
options(future.globals.maxSize = 3e+09)
options(Seurat.object.assay.version = "v5")
```

#ACTA2 Mutant Sample \## Uploading 10X Data

Data were analyzed using the Harvard single-cell sequencing core.


``` r
library(dplyr)
library(Seurat)
library(patchwork)
library(openxlsx)
# Load the dataset from the 10X sequenced data
taa.dat <- Read10X(data.dir = "../1_Input/Sequencing/Data_raw")
# Initialize the Seurat object with the raw (non-normalized data).
taa <- CreateSeuratObject(counts = taa.dat, project = "taa_acta2", min.cells = 3, min.features = 200)
taa
```

## Quality Control and Filtering

Standard filtering and QC tehcniques were used to remove poor-quality
reads and nuclei.


``` r
# Create a feature that quantifies the number of mitochondrial counts (impurity)
taa[["percent.mt"]] <- PercentageFeatureSet(taa, pattern = "^MT-")
head(taa@meta.data, 5)
# Visualize QC metrics as a violin plot
VlnPlot(taa, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# We can visualize feature-feature comparisons using the following 'FeatureScatter'
plot1 <- FeatureScatter(taa, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(taa, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
# The data were filtered to remove nuclei containing >5% mitochondrial counts, as well as unique feature counts greater than 2,500 or less than 200 to minimize the effects of multiplets and empty droplets, respectively.
taa_filtered <- subset(taa, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
taa_filtered
VlnPlot(taa, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
```

## Data Normalization and Scaling

A standard global-scaling normalization method "LogNormalize" was used
that normalizes the feature expression measurements for each cell by the
total expression, multiplies this by a scale factor (10,000 by default),
and log-transforms the result.


``` r
taa_normalized <- NormalizeData(taa_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
taa_normalized <- FindVariableFeatures(taa_normalized, selection.method = "vst", nfeatures = 10000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(taa_normalized), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(taa_normalized)
plot2 <- LabelPoints(plot = plot1, points = c(top10,"ACTA2"), repel = TRUE)
plot2
# The next task is to perform linear transformation needed for dimensionality reduction algorithms (PCA, UMAP, etc).
all.genes <- rownames(taa_normalized)
acta2_scaled <- ScaleData(taa_normalized, features = all.genes) # note: the PCA only requires the top X number of most variable genes, as defined by the 'VariableFeatures' command. To make it run faster, we can run the scaling on these alone (which it does by default, i.e. without the 'features' parameter)
```

## Dimensionality Reduction

We can first visualize the primary features/genes driving dimensionality
in our dataset; specifically, the principal components with largest
inter-cellular variability can be seen. Once this is accomplished, we
identify the "dimensionality" of the dataset, or the number of principal
components. From these heuristics, the first 10 dimensions reflect most
of the variabilitiy present in this dataset.


``` r
# to control for any confounders, we can add the 'vars.to.regress' parameter, including any variable in the metadata <- ScaleData(pbmc, vars.to.regress = "percent.mt")
acta2_scaled <- RunPCA(acta2_scaled, features = VariableFeatures(object = acta2_scaled))
acta2_scaled <- RunTSNE(acta2_scaled, assay = "RNA")
# Visualize the features/genes
VizDimLoadings(acta2_scaled, dims = 1:5, reduction = "pca")
DimPlot(acta2_scaled, reduction = "pca")
DimPlot(acta2_scaled, reduction = "tsne")
pdf(file = "../2_Output/1_ACTA/2_Clustering/PC_Heatmaps.pdf")
DimHeatmap(acta2_scaled, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()
```

## Unsupervised Clustering

Standard methods were used to cluster nuclei based on Uniform Manifold
Approximation and Projection (UMAP), a non-linear dimensionality
reduction technique that projects high-dimensional data into lower
dimensions while preserving its local structure and relationships. It
works by constructing a weighted graph representing the data’s local
neighborhood and then optimizing its layout in a low-dimensional space
to capture both local and global data patterns.


``` r
# Create a K-nearest neighbor (KNN) graph based on the euclidean distance in PCA space
acta2_scaled <- FindNeighbors(acta2_scaled, dims = 1:10) # uses the number of dimensions identified in the previous step ("dimensionality")
# Now  modularity optimization techniques such as the Louvain algorithm
acta2_scaled <- FindClusters(acta2_scaled, resolution = c(0.3, 0.4, 0.5, 0.6))
# Non-linear dimensional reductino (UMAP/tSNE)
acta2_scaled <- RunUMAP(acta2_scaled, dims = 1:30)
# Create the UMAP:
UMAP_Clustered<-DimPlot(acta2_scaled, reduction = "umap", group.by = "RNA_snn_res.0.3", label = T) + NoLegend()
tSNA_Clustered <- DimPlot(acta2_scaled, reduction = "tsne", group.by = "RNA_snn_res.0.3", label = T) + NoLegend()
UMAP_Clustered+tSNA_Clustered
pdf(file = "../2_Output/1_ACTA/2_Clustering/UMAP_clustered.pdf")
UMAP_Clustered+tSNA_Clustered
dev.off()
Idents(acta2_scaled) <- "RNA_snn_res.0.3"
# Save this instance to avoid needing to re-run:
saveRDS(acta2_scaled, file = "ACTA2_clustered.rds")
```

## Annotated UMAP, Density Plots

Annotation of clusters was accomplished by comparing the expression
patters of canonical cell type marker genes, shown below.


``` r
library(ggplot2)
library(dplyr)
library(Seurat)
acta2_scaled<-readRDS(file = "ACTA2_clustered.rds")
# Find all gene markers
DEGs_Clusters<-FindAllMarkers(acta2_scaled)
write.csv(DEGs_Clusters, "../2_Output/1_ACTA/DEGs_Clusters.csv")
#Identify Clusters corresponding with known gene markers:
VSMC1_genes<-c("PKD1", "COL6A2", "PDLIM7", "FLNA", "SMTN")
VSMC2_genes<-c("MYH11", "ACTA2", "ITGA8", "PRKG1", "CRISPLD1")
Fibroblast1_genes<-c("ABCA10", "C3", "ADGRD1", "FBLN1", "DCN")
Fibroblast2_genes<-c("UACA", "NFASC", "PRSS23", "SAMD5") # "TEX41",
Fibromyocyte_genes<-c("DGKG", "ADAMTS1", "RGS6", "TNC", "GRIP2") #, "ANGPT2"
EC1_genes<-c("DIPK2B", "ARHGEF15", "STC1", "FLT1") # "NOTCH4",
EC2_genes<-c("VWF", "BMPER", "BMX", "NOS1")
NKT_genes<-c("SKAP1", "RIPOR2", "FYN", "ITGAL", "CD96")
Macrophage_genes<-c("MRC1", "LGMN", "RBPJ", "F13A1", "RBM47")
Dendritic_genes<-c("ITGAX", "S100A9", "CSF3R", "CXCL8")
chondromyocytes <- c("COL2A1", "ACAN", "SOX9")
# Plot density function
library(Nebulosa)
VSMC1_density<-plot_density(acta2_scaled, VSMC1_genes, joint = TRUE, combine = FALSE, pal = "magma")
VSMC2_density<-plot_density(acta2_scaled, VSMC2_genes, joint = TRUE, combine = FALSE, pal = "magma")
Fibroblast1_density<-plot_density(acta2_scaled, Fibroblast1_genes, joint = TRUE, combine = FALSE, pal = "magma")
Fibroblast2_density<-plot_density(acta2_scaled, Fibroblast2_genes, joint = TRUE, combine = FALSE, pal = "magma")
Fibromyocyte_density<-plot_density(acta2_scaled, Fibromyocyte_genes, joint = TRUE, combine = FALSE, pal = "magma")
EC1_density<-plot_density(acta2_scaled, EC1_genes, joint = TRUE, combine = FALSE, pal = "magma")
EC2_density<-plot_density(acta2_scaled, EC2_genes, joint = TRUE, combine = FALSE, pal = "magma")
Macrophage_density<-plot_density(acta2_scaled, Macrophage_genes, joint = TRUE, combine = FALSE, pal = "magma")
NKT_density<-plot_density(acta2_scaled, NKT_genes, joint = TRUE, combine = FALSE, pal = "magma")
Dendritic_density<-plot_density(acta2_scaled, Dendritic_genes, joint = TRUE, combine = FALSE, pal = "magma")
Chondro_density<-plot_density(acta2_scaled, chondromyocytes, joint = TRUE, combine = FALSE, pal = "magma")
Chondro_density
# Print
pdf(file = "../2_Output/1_ACTA/2_Clustering/Density_plots.pdf", height = 3, width = 3)
VSMC1_density
VSMC2_density
Fibroblast1_density
Fibroblast2_density
Fibromyocyte_density
EC1_density
EC2_density
Macrophage_density
NKT_density
Dendritic_density
dev.off()
# Overlay these gene markers onto the UMAP to identify clusters
FeaturePlot(acta2_scaled, reduction = "umap", label = T, features = Fibroblast1_genes)
FeaturePlot(acta2_scaled, reduction = "umap", label = T, features = Fibroblast2_genes)
FeaturePlot(acta2_scaled, reduction = "umap", label = T, features = EC1_genes)
FeaturePlot(acta2_scaled, reduction = "umap", label = T, features = EC2_genes)
FeaturePlot(acta2_scaled, reduction = "umap", label = T, features = NKT_genes)
FeaturePlot(acta2_scaled, reduction = "umap", label = T, features = Macrophage_genes)
FeaturePlot(acta2_scaled, reduction = "umap", label = T, features = VSMC1_genes)
FeaturePlot(acta2_scaled, reduction = "umap", label = T, features = VSMC2_genes)
FeaturePlot(acta2_scaled, reduction = "umap", label = T, features = Dendritic_genes)
#Create a figure of cell-type specific markers overlying the UMAP
pdf(file = "../2_Output/1_ACTA/2_Clustering/Cell.Type_UMAP.pdf", height = 10, width = 15)
FeaturePlot(acta2_scaled, reduction = "umap", label = T, features = c("ACTA2", "FBLN1", "F13A1", "VWF", "STC1", "NFASC", "ITGAL", "ITGAX"))
dev.off()
# Create cluster labels based on the feature plot delineations of cellular phenotypes
acta2_scaled$umap_cluster <- Idents(acta2_scaled) # save old cluster names (for trajectory analysis)
acta2_scaled <- RenameIdents(object = acta2_scaled, # Rename cluster names
                          `0` = "VSMC_1", 
                          `1` = "VSMC_2", 
                          `2` = "VSMC_3",
                          `3` = "Fibroblast",
                          `4` = "Macrophage",
                          `5` = "EC2",
                          `6` = "EC1",
                          `7` = "NKT",
                          `8` = "Dendritic")
acta2_scaled$cell_type <- Idents(acta2_scaled)
acta2_scaled$cell_type <- factor(acta2_scaled$cell_type, 
                     levels = c("VSMC_1",
                                "VSMC_2",
                                "VSMC_3",
                                "Fibroblast", 
                                "EC1",
                                "EC2",
                                "Macrophage",
                                "NKT",
                                "Dendritic"))
acta2_scaled <- SetIdent(acta2_scaled, value = "cell_type")
UMAP_CellTypes<-DimPlot(acta2_scaled, reduction = "umap", label = T, repel = T,label.size = 4,
                        cols = c("coral2", 
                                 "coral4",
                                 "coral3",
                                 "wheat", 
                                 "steelblue4",
                                 "deepskyblue3",
                                 "azure4",
                                 "gray",
                                 "tan2")) + NoLegend()
tSNE_CellTypes<-DimPlot(acta2_scaled, reduction = "tsne", label = T, repel = T,label.size = 4,
                        cols = c("coral2", 
                                 "coral4",
                                 "coral3",
                                 "wheat", 
                                 "steelblue4",
                                 "deepskyblue3",
                                 "azure4",
                                 "gray",
                                 "tan2")) + NoLegend()
UMAP_CellTypes+tSNE_CellTypes
pdf(file = "../2_Output/1_ACTA/2_Clustering/UMAP_Annotated.pdf", height = 5, width = 5)
UMAP_CellTypes
dev.off()
## Use ridgeplots to identify bimodal gene marker distributions (enriched clusters)
pdf(file = "../2_Output/1_ACTA/2_Clustering/Celltype_RidgePlots.pdf", height = 10, width = 15)
RidgePlot(acta2_scaled, features = Fibroblast1_genes, ncol = 2)
RidgePlot(acta2_scaled, features = Fibroblast2_genes, ncol = 2)
RidgePlot(acta2_scaled, features = EC1_genes, ncol = 2)
RidgePlot(acta2_scaled, features = EC2_genes, ncol = 2)
RidgePlot(acta2_scaled, features = NKT_genes, ncol = 2)
RidgePlot(acta2_scaled, features = Macrophage_genes, ncol = 2)
RidgePlot(acta2_scaled, features = VSMC1_genes, ncol = 2)
RidgePlot(acta2_scaled, features = VSMC2_genes, ncol = 2)
RidgePlot(acta2_scaled, features = Dendritic_genes, ncol = 2)
dev.off()
# Differential Expression
# Export DEGs using cell-type clusters
DEGs_CellTypes<-FindAllMarkers(acta2_scaled)
write.csv(DEGs_CellTypes, "../2_Output/1_ACTA/DEGs_Clusters.csv")
#For each cluster
EC1.markers <- FindMarkers(acta2_scaled, ident.1 = "EC1", min.pct = 0.25)
EC2.markers <- FindMarkers(acta2_scaled, ident.1 = "EC2", min.pct = 0.25)
Dendritic.markers <- FindMarkers(acta2_scaled, ident.1 = "Dendritic", min.pct = 0.25)
VSMC1.markers <- FindMarkers(acta2_scaled, ident.1 = "VSMC_1", min.pct = 0.25)
VSMC2.markers <- FindMarkers(acta2_scaled, ident.1 = "VSMC_2", min.pct = 0.25)
Fibroblast.markers <- FindMarkers(acta2_scaled, ident.1 = "Fibroblast", min.pct = 0.25)
Macrophage.markers <- FindMarkers(acta2_scaled, ident.1 = "Macrophage", min.pct = 0.25)
NKT.markers <- FindMarkers(acta2_scaled, ident.1 = "NKT", min.pct = 0.25)
# Create a dot-bplot of the top 5 markers for each 
library(scCustomize)
paletteLength <- 100
myColor <- colorRampPalette(c("dodgerblue4", "white", "coral2"))(paletteLength)
top5_markers <- Extract_Top_Markers(marker_dataframe = DEGs_CellTypes, num_genes = 5, named_vector = FALSE,
    make_unique = TRUE)
pdf(file = "../2_Output/1_ACTA/2_Clustering/Dot.plot_top5markers.pdf", height = 10, width = 7)
Clustered_DotPlot(seurat_object = acta2_scaled, features = top5_markers, k = 10, colors_use_exp = myColor)
dev.off()
# Create excel sheet
library(openxlsx)
wb_DESeq<-createWorkbook()
  addWorksheet(wb_DESeq, "DEGs_all")
  writeData(wb_DESeq, "DEGs_all", DEGs_CellTypes, startCol = 1, rowNames = T)
  addWorksheet(wb_DESeq, "EC1")
  writeData(wb_DESeq, "EC1", EC1.markers, startCol = 1, rowNames = T)
    addWorksheet(wb_DESeq, "EC2")
  writeData(wb_DESeq, "EC2", EC2.markers, startCol = 1, rowNames = T)
    addWorksheet(wb_DESeq, "VSMC1")
  writeData(wb_DESeq, "VSMC1", VSMC1.markers, startCo = 1, rowNames = T)
    addWorksheet(wb_DESeq, "VSMC2")
  writeData(wb_DESeq, "VSMC2", VSMC2.markers, startCo = 1, rowNames = T)
      addWorksheet(wb_DESeq, "Fibroblast")
  writeData(wb_DESeq, "Fibroblast", Fibroblast.markers, startCol = 1, rowNames = T)
    addWorksheet(wb_DESeq, "Macrophage")
  writeData(wb_DESeq, "Macrophage", Macrophage.markers, startCol = 1, rowNames = T)
    addWorksheet(wb_DESeq, "Dendritic")
  writeData(wb_DESeq, "Dendritic", Dendritic.markers, startCol = 1, rowNames = T)
      addWorksheet(wb_DESeq, "NKT")
  writeData(wb_DESeq, "NKT", NKT.markers, startCol = 1, rowNames = T)
saveWorkbook(wb_DESeq, file = "../2_Output/1_ACTA/3_Differential.Expression/DEGs_CellType.Specific.xlsx", overwrite = T)
# Heatmap of Clusters
DEGs_CellTypes %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
pdf(file = "../2_Output/1_ACTA/2_Clustering/Celltype_Heatmap.pdf", height = 10, width = 15)
DoHeatmap(acta2_scaled, features = top10$gene, size = 3, disp.min = -2, disp.max = 2) + scale_fill_gradientn(colors = c("dodgerblue4", "white", "coral2")) + labs(title = "Heatmap of Top-10 most variable genes within each cluster")
dev.off()
#
# Contractility score
contractility_genes <- c("CNN1", "TAGLN", "ACTA2", "MYH11")
gene_list <- contractility_genes[contractility_genes %in% rownames(acta2_scaled)] 
gene_expr <- GetAssayData(object = acta2_scaled, assay = "RNA", slot = "data")[gene_list, ]
sum_expr <- colSums(gene_expr)
acta2_scaled$contractility <- scale(sum_expr)
#osteogenic score
gene_list <- c("CBFA1", "MSX2", "RUNX2", "SOX9", "SPP1", "BGLAP", "ALPL", "COL2A1")
gene_list <- gene_list[gene_list %in% rownames(acta2_scaled)]
gene_expr <- GetAssayData(object = acta2_scaled, assay = "RNA", slot = "data")[gene_list, ]
sum_expr <- colSums(gene_expr)
acta2_scaled$osteogenic <- scale(sum_expr)
# Synthetic score
gene_list <-c("CCND1","PCNA", "MKI67")
gene_list <- gene_list[gene_list %in% rownames(acta2_scaled)]
gene_expr <- GetAssayData(object = acta2_scaled, assay = "RNA", slot = "data")[gene_list, ]
sum_expr <- colSums(gene_expr)
acta2_scaled$Synthetic <- scale(sum_expr)
# EZH2 Target Score
# Load the package
library(GSEABase)
# Specify the file path to your GMT file
gmt_file <- "../1_Input/EZH2/EZH2.gmt"
# Read the GMT file
gene_sets <- getGmt(gmt_file)
# Explore the data
gene_sets_list <- geneIds(gene_sets) # Extract genes for each set
names(gene_sets_list) # View gene set names
EZH2_Targets <- gene_sets_list[[3]] # View genes in the first set
gene_list <- intersect(EZH2_Targets, rownames(acta2_scaled))
gene_expr <- GetAssayData(object = acta2_scaled, assay = "RNA", slot = "data")[gene_list, ]
sum_expr <- colSums(gene_expr)
acta2_scaled$EZH2_Score <- scale(sum_expr)
# Save file for downstream analysis
saveRDS(acta2_scaled, file = "../1_Input/ACTA2_annotated.rds")
```

## VSMC - Proliferation scores

To explore whether differences in VSMC subpopulations reflects
differential regulation of proliferative pathways, as previously
reported by Kwartler et al. 2023, we used the Gene Ontology term
GO0008283.


``` r
# Load the necessary libraries
acta2_scaled<-readRDS(file = "../1_Input/ACTA2_annotated.rds")
acta2_scaled <- subset(acta2_scaled, idents = c("VSMC_1", "VSMC_2", "VSMC_3")) # select only the VSMCs
acta2_scaled <- RunPCA(acta2_scaled) # re-cluster
acta2_scaled <- RunUMAP(acta2_scaled, assay = "RNA", dims = 1:40) # re-cluster
Idents(acta2_scaled) <- acta2_scaled$RNA_snn_res.0.3
library(clusterProfiler)
library(org.Mm.eg.db)  # Mouse gene database
# Define the mouse GO term of interest, e.g., GO:0008150 (biological process)
go_term <- "GO:0008283" # used in the original 2023 paper
# Extract gene symbols associated with the specified GO term
genes_vector <- AnnotationDbi::select(
  org.Mm.eg.db,                 # Mouse gene database
  keys = go_term,               # GO term of interest
  columns = "SYMBOL",           # Retrieve the gene symbols
  keytype = "GOALL"             # Define the key type as GO terms
)$SYMBOL
# Proliferation Score
genes_vector <- toupper(genes_vector)
gene_list <- genes_vector[genes_vector %in% rownames(acta2_scaled)] 
gene_expr <- GetAssayData(object = acta2_scaled, slot = "data")[gene_list, ]
sum_expr <- colSums(gene_expr)
acta2_scaled$proliferation <- scale(sum_expr)
pdf(file = "../2_Output/1_ACTA/4_Trajectory/Proliferation.Score_UMAP.pdf", height = 4, width = 4)
FeaturePlot(acta2_scaled, features = "proliferation", reduction = "umap") + ggtitle("Proliferation Score") +
scale_color_gradient2(low = "darkcyan", mid = "white", high = "darkred", midpoint = median(acta2_scaled$proliferation, na.rm = FALSE))
dev.off()
pdf(file = "../2_Output/1_ACTA/4_Trajectory/EZH2.Score_UMAP.pdf", height = 4, width = 4)
FeaturePlot(acta2_scaled, features = "EZH2_Score", reduction = "umap") + ggtitle("EZH2 Target Score") +
scale_color_gradient2(low = "darkcyan", mid = "white", high = "darkred", midpoint = median(acta2_scaled$EZH2_Score, na.rm = FALSE))
dev.off()
# Contractility
go_term <- "GO:0006940" # used in the original 2023 paper
genes_vector <- AnnotationDbi::select(
  org.Mm.eg.db,                 # Mouse gene database
  keys = go_term,               # GO term of interest
  columns = "SYMBOL",           # Retrieve the gene symbols
  keytype = "GOALL"             # Define the key type as GO terms
)$SYMBOL
genes_vector <- toupper(genes_vector)
gene_list <- genes_vector[genes_vector %in% rownames(acta2_scaled)] 
gene_expr <- GetAssayData(object = acta2_scaled, slot = "data")[gene_list, ]
sum_expr <- colSums(gene_expr)
acta2_scaled$contractility <- scale(sum_expr)
pdf(file = "../2_Output/1_ACTA/4_Trajectory/Contractility.Score_UMAP.pdf", height = 4, width = 4)
FeaturePlot(acta2_scaled, features = "contractility", reduction = "umap") + ggtitle("Contractility Score") +
scale_color_gradient2(low = "darkcyan", mid = "white", high = "darkred", midpoint = median(acta2_scaled$contractility, na.rm = FALSE))
dev.off()
```

## VSMC - Trajectory Analysis (Monocle 3)


``` r
library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(ggtrace)
ann_colors = list(CellType = c(VSMC_1="coral2",
              VSMC_2="coral4",
              VSMC_3="coral3",
                Fibroblast="wheat", 
                EC1="steelblue4",
                EC2="deepskyblue3",
                Macrophage="azure4",
                NKT="gray",
                Dendritic="tan2"))
acta2_scaled<-readRDS(file = "../1_Input/ACTA2_annotated.rds")
# Subsetting and re-normalization
acta2_scaled <- subset(acta2_scaled, idents = c("VSMC_1", "VSMC_2", "VSMC_3")) # select only the VSMCs
acta2_scaled <- RunPCA(acta2_scaled) # re-cluster
acta2_scaled <- RunUMAP(acta2_scaled, assay = "RNA", dims = 1:40) # re-cluster
Idents(acta2_scaled) <- acta2_scaled$RNA_snn_res.0.3
# Clustering
pdf(file = "../2_Output/1_ACTA/4_Trajectory/UMAP_VSMC.pdf", height = 4, width = 4)
DimPlot(acta2_scaled, group.by = "RNA_snn_res.0.3", reduction = "umap", label = T, repel = T,label.size = 4, cols = c("coral2", "coral4", "coral3"))+NoLegend()
dev.off()
# 
DefaultAssay(acta2_scaled) = "RNA" # This changed in seurat V5
# acta2_scaled <- ScaleData(object = acta2_scaled)
cds <- as.cell_data_set(acta2_scaled)
cds <- preprocess_cds(cds, num_dim = 100)
cds <- estimate_size_factors(cds)
# Include gene names (not done by default by the seurat conversion)
rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name
#  Run Monocle
cds <- cluster_cells(cds) # 0.000000002 This step creates "partitions" that are used in the trajectory inference
plot_cells(cds, show_trajectory_graph = FALSE, color_cells_by = "partition") # this shows the partitions overlain on the UMAP
cds <- learn_graph(cds, use_partition = TRUE) # creating trajectory inference within each partition
# cds <- order_cells(cds) # Used when setting the nodes; if already known, use the next line
root_cell <- "CAGATTCGTAACGTGA-1" # names(which(cds@principal_graph_aux$UMAP$pseudotime==0))
cds<-order_cells(cds, root_cells = root_cell)
# Plot the pseudotime on UMAP
plot_cells(cds, 
           color_cells_by = "pseudotime", 
           label_branch_points = FALSE, 
           label_leaves = FALSE)
pdf(file = "../2_Output/1_ACTA/4_Trajectory/ACTA2_UMAP_Trajectory_Partition.pdf", height = 3, width = 3)
plot_cells(cds, 
           color_cells_by = "partition", 
           label_branch_points = FALSE, 
           label_leaves = F,
           show_trajectory_graph = F,
           label_roots = F)
dev.off()
pdf(file = "../2_Output/1_ACTA/4_Trajectory/ACTA2_UMAP_Trajectory_Pseudotime.pdf", height = 3, width = 4)
plot_cells(cds, 
           color_cells_by = "pseudotime", 
           label_branch_points = FALSE, 
           label_leaves = F,
           show_trajectory_graph = F,
           label_roots = F)
dev.off()
# Examine specific genes
plot_cells(cds, 
           genes = c("ACTA2"),
           label_cell_groups = FALSE,
           show_trajectory_graph = FALSE,
           min_expr = 3)
# Identify pseudotime
modulated_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=8) # Identify differentially-expressed genes with pseudotime
modulated_genes <- na.omit(modulated_genes) # remove NA's
modulated_genes <- modulated_genes %>% filter(modulated_genes$q_value < 0.05 & modulated_genes$status =="OK") # filter cds results down
modulated_genes <- modulated_genes[order(-modulated_genes$morans_test_statistic), ] # order by moran's test
modulated_genes <- top_n(modulated_genes, 1000, -q_value)
#### Create a heatmap of genes with similar pseudotime kinetics
genes <- row.names(subset(modulated_genes, q_value < 0.05))
openxlsx::write.xlsx(modulated_genes, "../2_Output/1_ACTA2/Trajectory/Pseudotime_DEGs.xlsx")
library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(circlize)
library(monocle3)
pt.matrix <- as.data.frame(exprs(cds)[match(genes,rownames(rowData(cds))),order(pseudotime(cds))])
#
cell_names <- colnames(pt.matrix)
Index<-as.data.frame(cds@colData) %>% dplyr::select(cell_type)
Index<-subset(Index, row.names(Index) %in% cell_names)
Index$CellType <- factor(Index$cell_type,
                         levels = c("VSMC_1",
                                    "VSMC_2",
                                    "VSMC_3",
                                    "Fibroblast", 
                                    "EC1",
                                    "EC2",
                                    "Macrophage",
                                    "NKT",
                                    "Dendritic"))
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=6)$y})) # Create a spline that smooths the pseudotime-based expression along 6 degrees of freedom.
filtered_matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes
colnames(pt.matrix) <- cell_names
###########
paletteLength <- 20
myColor <- colorRampPalette(c("dodgerblue4", "white", "goldenrod1"))(paletteLength)
heatmap_DMC<-pheatmap::pheatmap(pt.matrix, scale="row", 
      cluster_cols = F, 
      cluster_rows = TRUE,
      cutree_rows = 4,
      fontsize_col = 8,
      color = myColor,
      annotation_col = Index,
      annotation_colors = ann_colors,
      show_colnames = F,
      show_rownames = F,
      # right_annotation = ha,
      border_color = NA)
################################
hc <-heatmap_DMC$tree_row
lbl <- cutree(hc, 4)
cluster1<-which(lbl==1)
cluster2<-which(lbl==2)
cluster3<-which(lbl==3)
cluster4<-which(lbl==4)
Cluster1_data<-pt.matrix[cluster1,]
Cluster2_data<-pt.matrix[cluster2,]
write.csv(Cluster2_data, "Cluster2.csv")
Cluster3_data<-pt.matrix[cluster3,]
Cluster4_data<-pt.matrix[cluster4,]
Cluster1_GENES <- rownames(Cluster1_data)
Cluster2_GENES <- rownames(Cluster2_data)
Cluster3_GENES <- rownames(Cluster3_data)
Cluster4_GENES <- rownames(Cluster4_data)
##Enrichr
library(enrichR)
library(dplyr)
dbs <- c("WikiPathway_2023_Human")
enriched_path__1 <- enrichr(Cluster1_GENES, dbs)
enrich_path_1<-enriched_path__1[[dbs]] %>% filter(Adjusted.P.value < 0.05)
head(enrich_path_1)
enriched_path__2 <- enrichr(Cluster2_GENES, dbs)
enrich_path_2<-enriched_path__2[[dbs]] %>% filter(Adjusted.P.value < 0.05)
head(enrich_path_2)
enriched_path__3 <- enrichr(Cluster3_GENES, dbs)
enrich_path_3<-enriched_path__3[[dbs]] %>% filter(Adjusted.P.value < 0.05)
head(enrich_path_3)
enriched_path__4 <- enrichr(Cluster4_GENES, dbs)
enrich_path_4<-enriched_path__4[[dbs]] %>% filter(Adjusted.P.value < 0.05)
head(enrich_path_4)
library(openxlsx)
wb_DESeq<-createWorkbook()
  addWorksheet(wb_DESeq, "Cluster 1")
  writeData(wb_DESeq, "Cluster 1", enrich_path_1, startCol = 1)
  addWorksheet(wb_DESeq, "Cluster 2")
  writeData(wb_DESeq, "Cluster 2", enrich_path_2, startCol = 1)
  addWorksheet(wb_DESeq, "Cluster 3")
  writeData(wb_DESeq, "Cluster 3", enrich_path_3, startCol = 1)
  addWorksheet(wb_DESeq, "Cluster 4")
  writeData(wb_DESeq, "Cluster 4", enrich_path_4, startCol = 1)
saveWorkbook(wb_DESeq, file = paste0("../2_Output/1_ACTA/4_Trajectory/Trajectory_Pathway.xlsx"), overwrite = TRUE)
dbs <- c("GWAS_Catalog_2023")
enriched_1 <- enrichr(Cluster1_GENES, dbs)
enrich_1<-enriched_1[[dbs]] %>% filter(Adjusted.P.value < 0.05)
head(enrich_1)
enriched_2 <- enrichr(Cluster2_GENES, dbs)
enrich_2<-enriched_2[[dbs]] %>% filter(Adjusted.P.value < 0.05)
head(enrich_2)
enriched_3 <- enrichr(Cluster3_GENES, dbs)
enrich_3<-enriched_3[[dbs]] %>% filter(Adjusted.P.value < 0.05)
head(enrich_3)
enriched_4 <- enrichr(Cluster4_GENES, dbs)
enrich_4<-enriched_4[[dbs]] %>% filter(Adjusted.P.value < 0.05)
head(enrich_4)
library(openxlsx)
wb_DESeq<-createWorkbook()
  addWorksheet(wb_DESeq, "Cluster 1")
  writeData(wb_DESeq, "Cluster 1", enrich_1, startCol = 1)
  addWorksheet(wb_DESeq, "Cluster 2")
  writeData(wb_DESeq, "Cluster 2", enrich_2, startCol = 1)
  addWorksheet(wb_DESeq, "Cluster 3")
  writeData(wb_DESeq, "Cluster 3", enrich_3, startCol = 1)
  addWorksheet(wb_DESeq, "Cluster 4")
  writeData(wb_DESeq, "Cluster 4", enrich_4, startCol = 1)
saveWorkbook(wb_DESeq, file = paste0("../2_Output/1_ACTA/4_Trajectory/Trajectory_GWAS.Traits.xlsx"), overwrite = TRUE)
library(stringr)
Aortic_genes <- enrich_2 %>% filter(str_detect(Term,"Aort"))
Aortic_genes <- Aortic_genes$Genes
Aortic_genes <- unique(unlist(strsplit(Aortic_genes, ";", fixed = T)))
# proliferative score
gene_list <-Aortic_genes
gene_list <- gene_list[gene_list %in% rownames(acta2_scaled)]
gene_expr <- GetAssayData(object = acta2_scaled, assay = "RNA", slot = "data")[gene_list, ]
sum_expr <- colSums(gene_expr)
acta2_scaled$Aortopathy <- scale(sum_expr)
Cluster1_Names<-unlist(strsplit(enrich_1$Genes[1], ";", fixed = T))
Cluster2_Names<-unlist(strsplit(enrich_2$Genes[1], ";", fixed = T))
Cluster3_Names<-unlist(strsplit(enrich_3$Genes[1], ";", fixed = T))
Cluster4_Names<-unlist(strsplit(enrich_4$Genes[1], ";", fixed = T))
pdf(file = "../2_Output/1_ACTA/4_Trajectory/Aortic_genes.pdf")
plot_cells(cds, 
           genes = Aortic_genes,
           label_cell_groups = F,
           label_roots = F,
           label_branch_points = F,
           label_leaves = F,
           show_trajectory_graph = F,
           min_expr = 1,
           alpha = 0.8,
           trajectory_graph_color = "grey28",
           ) +
  theme(
  axis.text = element_text(size = 6),    # Adjust the size as needed
  axis.title = element_text(size = 8),  # Adjust the size as needed
  legend.text = element_text(size = 4),  # Adjust the size as needed
  legend.title = element_text(size = 6),
  legend.key.size = unit(2, 'mm')) +
  scale_colour_viridis_c(option = "inferno")
dev.off()

TAA.combined <- readRDS(file = "../1_Input/TAA_snRNA.rds")
library(patchwork)
  plots_up <- VlnPlot(TAA.combined, features = Aortic_genes, split.by = "Disease", group.by = "cell_type",
      pt.size = 0, combine = FALSE)
  plots_up <- lapply(seq_along(plots_up), function(i) {
    plot <- plots_up[[i]] + theme(axis.title.x = element_blank())  # Remove x-axis label
    if (i == length(plots_up)) {
      return(plot)  # Keep the legend in the last plot
    } else {
      return(plot + NoLegend())  # Remove legend in other plots
    }
  })
  pdf(paste0("../2_Output/3_Merged/3_Differential.Expression/Conserved.pdf"), height = 6, width = 10)
  print(wrap_plots(plots = plots_up, ncol = 6))
  dev.off()

############################################################################
GENES_HM<-c(Aortic_genes) #
pt.df<-as.data.frame(pt.matrix)
pt.df$gene_name<-rownames(pt.matrix)
ha = rowAnnotation(link = anno_mark(at = which(pt.df$gene_name %in% GENES_HM),
                   labels = as.character(pt.df[which(pt.df$gene_name %in% GENES_HM), "gene_name"]),
                   labels_gp = gpar(fontsize = 8),
                   padding = unit(1, "mm"))
                   )
heatmap_combination<-ComplexHeatmap::pheatmap(pt.matrix, scale="row",
                    cluster_cols = F,
                    cluster_rows = T,
                    cutree_rows = 4,
                    # cutree_cols = 3,
                     fontsize_col = 5,
                     color = myColor,
                    annotation_names_col = FALSE,
                    show_colnames = F,
                     show_rownames = F,
                     border_color = NA,
                    annotation_colors = NA,
                    right_annotation = ha,
                    annotation_col = Index,
                    border = TRUE)
pdf(file = "../2_Output/1_ACTA/4_Trajectory/Trajectory_Heatmap.pdf", height = 6, width = 6)
heatmap_combination
dev.off()
pdf(file = "../2_Output/1_ACTA/4_Trajectory/Aortopathy.Score_UMAP.pdf", height = 4, width = 4)
FeaturePlot(acta2_scaled, features = "Aortopathy", reduction = "umap") + ggtitle("Aortopathy Score") +
scale_color_gradient2(low = "dodgerblue3", mid = "white", high = "goldenrod2", midpoint = median(acta2_scaled$Aortopathy, na.rm = FALSE))
dev.off()
```

# TAA Control FFPE Sample

Data were analyzed using the Harvard single-cell sequencing core.


``` r
library(dplyr)
library(Seurat)
library(patchwork)
library(openxlsx)
# Load the dataset from the 10X sequenced data
taa.dat <- Read10X(data.dir = "../1_Input/Sequencing/bri2327_TAA/seurat")
# Initialize the Seurat object with the raw (non-normalized data).
taa <- CreateSeuratObject(counts = taa.dat, project = "taa_acta2", min.cells = 3, min.features = 200)
taa
```

## Quality Control and Filtering

Standard filtering and QC tehcniques were used to remove poor-quality
reads and nuclei.


``` r
# Create a feature that quantifies the number of mitochondrial counts (impurity)
taa[["percent.mt"]] <- PercentageFeatureSet(taa, pattern = "^MT-")
head(taa@meta.data, 5)
# Visualize QC metrics as a violin plot
VlnPlot(taa, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# We can visualize feature-feature comparisons using the following 'FeatureScatter'
plot1 <- FeatureScatter(taa, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(taa, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
# The data were filtered to remove nuclei containing >5% mitochondrial counts, as well as unique feature counts greater than 2,500 or less than 200 to minimize the effects of multiplets and empty droplets, respectively.
taa_filtered <- subset(taa, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
taa_filtered
VlnPlot(taa, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
```

## Data Normalization

A standard global-scaling normalization method "LogNormalize" was used
that normalizes the feature expression measurements for each cell by the
total expression, multiplies this by a scale factor (10,000 by default),
and log-transforms the result.


``` r
taa_normalized <- NormalizeData(taa_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
```

## Cellular Differential Expression

We sought to identify features/genes exhibiting the highest cell-cell
variability.


``` r
taa_normalized <- FindVariableFeatures(taa_normalized, selection.method = "vst", nfeatures = 10000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(taa_normalized), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(taa_normalized)
plot2 <- LabelPoints(plot = plot1, points = c(top10,"ACTA2"), repel = TRUE)
plot2
# The next task is to perform linear transformation needed for dimensionality reduction algorithms (PCA, UMAP, etc).
all.genes <- rownames(taa_normalized)
con_scaled <- ScaleData(taa_normalized, features = all.genes) # note: the PCA only requires the top X number of most variable genes, as defined by the 'VariableFeatures' command. To make it run faster, we can run the scaling on these alone (which it does by default, i.e. without the 'features' parameter)
```

## Dimensionality

We can first visualize the primary features/genes driving dimensionality
in our dataset; specifically, the principal components with largest
inter-cellular variability can be seen. Once this is accomplished, we
identify the "dimensionality" of the dataset, or the number of principal
components. From these heuristics, the first 10 dimensions reflect most
of the variabilitiy present in this dataset.


``` r
# to TAAng for any confounders, we can add the 'vars.to.regress' parameter, including any variable in the metadata <- ScaleData(pbmc, vars.to.regress = "percent.mt")
TAAng_scaled <- RunPCA(con_scaled, features = VariableFeatures(object = con_scaled))
TAAng_scaled <- RunTSNE(TAAng_scaled, assay = "RNA")
# Visualize the features/genes
VizDimLoadings(TAAng_scaled, dims = 1:5, reduction = "pca")
DimPlot(TAAng_scaled, reduction = "pca")
DimPlot(TAAng_scaled, reduction = "tsne")
pdf(file = "../2_Output/2_TAA/2_Clustering/PC_Heatmaps.pdf")
DimHeatmap(TAAng_scaled, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()
```

## Clustering


``` r
# Create a K-nearest neighbor (KNN) graph based on the euclidean distance in PCA space
TAAng_scaled <- FindNeighbors(TAAng_scaled, dims = 1:10) # uses the number of dimensions identified in the previous step ("dimensionality")
# Now  modularity optimization techniques such as the Louvain algorithm
TAAng_scaled <- FindClusters(TAAng_scaled, resolution = c(0.3, 0.4, 0.5, 0.6))
# Non-linear dimensional reductino (UMAP/tSNE)
TAAng_scaled <- RunUMAP(TAAng_scaled, dims = 1:30)
# Create the UMAP:
UMAP_Clustered<-DimPlot(TAAng_scaled, reduction = "umap", group.by = "RNA_snn_res.0.3", label = T) + NoLegend()
tSNA_Clustered <- DimPlot(TAAng_scaled, reduction = "tsne", group.by = "RNA_snn_res.0.3", label = T) + NoLegend()
UMAP_Clustered+tSNA_Clustered
pdf(file = "../2_Output/2_TAA/2_Clustering/UMAP_clustered.pdf")
UMAP_Clustered+tSNA_Clustered
dev.off()
Idents(TAAng_scaled) <- "RNA_snn_res.0.3"
# Save this instance to avoid needing to re-run:
saveRDS(TAAng_scaled, file = "../1_Input/TAA_nongenomic_clustered.rds")
```

## Gene Markers


``` r
library(ggplot2)
library(dplyr)
library(Seurat)
TAAng_scaled<-readRDS(file = "../1_Input/TAA_nongenomic_clustered.rds")
# Find all gene markers
DEGs_Clusters<-FindAllMarkers(TAAng_scaled)
write.csv(DEGs_Clusters, "../2_Output/2_TAA/DEGs_Clusters.csv")
#Identify Clusters corresponding with known gene markers:
VSMC1_genes<-c("PKD1", "COL6A2", "PDLIM7", "FLNA", "SMTN")
VSMC2_genes<-c("MYH11", "ACTA2", "ITGA8", "PRKG1", "CRISPLD1")
Fibroblast1_genes<-c("ABCA10", "C3", "ADGRD1", "FBLN1", "DCN")
Fibroblast2_genes<-c("UACA", "NFASC", "PRSS23", "SAMD5") # "TEX41",
Fibromyocyte_genes<-c("DGKG", "ADAMTS1", "RGS6", "TNC", "GRIP2") #, "ANGPT2"
EC1_genes<-c("DIPK2B", "ARHGEF15", "STC1", "FLT1") # "NOTCH4",
EC2_genes<-c("VWF", "BMPER", "BMX", "NOS1")
NKT_genes<-c("SKAP1", "RIPOR2", "FYN", "ITGAL", "CD96")
Macrophage_genes<-c("MRC1", "LGMN", "RBPJ", "F13A1", "RBM47")
Dendritic_genes<-c("ITGAX", "S100A9", "CSF3R", "CXCL8")
# Plot density function
library(Nebulosa)
VSMC1_density<-plot_density(TAAng_scaled, VSMC1_genes, joint = TRUE, combine = FALSE, pal = "magma")
VSMC2_density<-plot_density(TAAng_scaled, VSMC2_genes, joint = TRUE, combine = FALSE, pal = "magma")
Fibroblast1_density<-plot_density(TAAng_scaled, Fibroblast1_genes, joint = TRUE, combine = FALSE, pal = "magma")
Fibroblast2_density<-plot_density(TAAng_scaled, Fibroblast2_genes, joint = TRUE, combine = FALSE, pal = "magma")
Fibromyocyte_density<-plot_density(TAAng_scaled, Fibromyocyte_genes, joint = TRUE, combine = FALSE, pal = "magma")
EC1_density<-plot_density(TAAng_scaled, EC1_genes, joint = TRUE, combine = FALSE, pal = "magma")
EC2_density<-plot_density(TAAng_scaled, EC2_genes, joint = TRUE, combine = FALSE, pal = "magma")
Macrophage_density<-plot_density(TAAng_scaled, Macrophage_genes, joint = TRUE, combine = FALSE, pal = "magma")
NKT_density<-plot_density(TAAng_scaled, NKT_genes, joint = TRUE, combine = FALSE, pal = "magma")
Dendritic_density<-plot_density(TAAng_scaled, Dendritic_genes, joint = TRUE, combine = FALSE, pal = "magma")
# Print
pdf(file = "../2_Output/2_TAA/2_Clustering/Density_plots.pdf")
VSMC1_density
VSMC2_density
Fibroblast1_density
Fibroblast2_density
Fibromyocyte_density
EC1_density
EC2_density
Macrophage_density
NKT_density
Dendritic_density
dev.off()
# Overlay these gene markers onto the UMAP to identify clusters
FeaturePlot(TAAng_scaled, reduction = "umap", label = T, features = Fibroblast1_genes)
FeaturePlot(TAAng_scaled, reduction = "umap", label = T, features = Fibroblast2_genes)
FeaturePlot(TAAng_scaled, reduction = "umap", label = T, features = EC1_genes)
FeaturePlot(TAAng_scaled, reduction = "umap", label = T, features = EC2_genes)
FeaturePlot(TAAng_scaled, reduction = "umap", label = T, features = NKT_genes)
FeaturePlot(TAAng_scaled, reduction = "umap", label = T, features = Macrophage_genes)
FeaturePlot(TAAng_scaled, reduction = "umap", label = T, features = VSMC1_genes)
FeaturePlot(TAAng_scaled, reduction = "umap", label = T, features = VSMC2_genes)
FeaturePlot(TAAng_scaled, reduction = "umap", label = T, features = Dendritic_genes)
#Create a figure of cell-type specific markers overlying the UMAP
pdf(file = "../2_Output/2_TAA/2_Clustering/Cell.Type_UMAP.pdf", height = 10, width = 15)
FeaturePlot(TAAng_scaled, reduction = "umap", label = T, features = c("TAAng", "FBLN1", "F13A1", "VWF", "STC1", "NFASC", "ITGAL", "ITGAX"))
dev.off()
# Create cluster labels based on the feature plot delineations of cellular phenotypes
TAAng_scaled$umap_cluster <- Idents(TAAng_scaled) # save old cluster names (for trajectory analysis)
TAAng_scaled <- RenameIdents(object = TAAng_scaled, # Rename cluster names
                          `0` = "VSMC", 
                          `1` = "VSMC", 
                          `2` = "Macrophage",
                          `3` = "Macrophage",
                          `4` = "Fibromyocyte",
                          `5` = "Fibroblast",
                          `6` = "NKT",
                          `7` = "EC",
                          `8` = "Dendritic",
                          `9` = "Lymphoid",
                          `10` = "10")
TAAng_scaled$cell_type <- Idents(TAAng_scaled)
TAAng_scaled$cell_type <- factor(TAAng_scaled$cell_type, levels = c("VSMC", 
                                                                    "Fibroblast", 
                                                                    "EC",
                                                                    "Macrophage",
                                                                    "NKT",
                                                                    "Dendritic",
                                                                    "Fibromyocyte",
                                                                    "Lymphoid",
                                                                    "10"))
TAAng_scaled <- SetIdent(TAAng_scaled, value = "cell_type")
UMAP_CellTypes<-DimPlot(TAAng_scaled, reduction = "umap", label = T, repel = T,label.size = 4,
                        cols = c("coral2", 
                                 "wheat", 
                                 "steelblue4",
                                 "azure4",
                                 "gray",
                                 "tan2",
                                 "darkcyan",
                                 "darkgray",
                                 "black")) + NoLegend()
tSNE_CellTypes<-DimPlot(TAAng_scaled, reduction = "tsne", label = T, repel = T,label.size = 4,
                        cols = c("coral2", 
                                 "wheat", 
                                 "steelblue4",
                                 "azure4",
                                 "gray",
                                 "tan2",
                                 "darkcyan",
                                 "darkgray",
                                 "black")) + NoLegend()

UMAP_CellTypes+tSNE_CellTypes
pdf(file = "../2_Output/2_TAA/2_Clustering/UMAP_Annotated.pdf")
UMAP_CellTypes
dev.off()
## Use ridgeplots to identify bimodal gene marker distributions (enriched clusters)
pdf(file = "../2_Output/2_TAA/2_Clustering/Celltype_RidgePlots.pdf", height = 10, width = 15)
RidgePlot(TAAng_scaled, features = Fibroblast1_genes, ncol = 2)
RidgePlot(TAAng_scaled, features = Fibroblast2_genes, ncol = 2)
RidgePlot(TAAng_scaled, features = EC1_genes, ncol = 2)
RidgePlot(TAAng_scaled, features = EC2_genes, ncol = 2)
RidgePlot(TAAng_scaled, features = NKT_genes, ncol = 2)
RidgePlot(TAAng_scaled, features = Macrophage_genes, ncol = 2)
RidgePlot(TAAng_scaled, features = VSMC1_genes, ncol = 2)
RidgePlot(TAAng_scaled, features = VSMC2_genes, ncol = 2)
RidgePlot(TAAng_scaled, features = Dendritic_genes, ncol = 2)
dev.off()
# Differential Expression
# Export DEGs using cell-type clusters
DEGs_CellTypes<-FindAllMarkers(TAAng_scaled)
write.csv(DEGs_CellTypes, "../2_Output/2_TAA/DEGs_Clusters.csv")
#For each cluster
VSMC.markers <- FindMarkers(TAAng_scaled, ident.1 = "VSMC", min.pct = 0.25)
EC.markers <- FindMarkers(TAAng_scaled, ident.1 = "EC", min.pct = 0.25)
Dendritic.markers <- FindMarkers(TAAng_scaled, ident.1 = "Dendritic", min.pct = 0.25)
VSMC.markers <- FindMarkers(TAAng_scaled, ident.1 = "VSMC", min.pct = 0.25)
Fibroblast.markers <- FindMarkers(TAAng_scaled, ident.1 = "Fibroblast", min.pct = 0.25)
Macrophage.markers <- FindMarkers(TAAng_scaled, ident.1 = "Macrophage", min.pct = 0.25)
NKT.markers <- FindMarkers(TAAng_scaled, ident.1 = "NKT", min.pct = 0.25)
# Create a dot-bplot of the top 5 markers for each 
library(scCustomize)
paletteLength <- 100
myColor <- colorRampPalette(c("dodgerblue4", "white", "coral2"))(paletteLength)
top5_markers <- Extract_Top_Markers(marker_dataframe = DEGs_CellTypes, num_genes = 5, named_vector = FALSE,
    make_unique = TRUE)
pdf(file = "../2_Output/2_TAA/2_Clustering/Dot.plot_top5markers.pdf", height = 10, width = 7)
Clustered_DotPlot(seurat_object = TAAng_scaled, features = top5_markers, k = 10, colors_use_exp = myColor)
dev.off()
# Create excel sheet
library(openxlsx)
wb_DESeq<-createWorkbook()
  addWorksheet(wb_DESeq, "DEGs_all")
  writeData(wb_DESeq, "DEGs_all", DEGs_CellTypes, startCol = 1, rowNames = T)
  addWorksheet(wb_DESeq, "EC1")
  writeData(wb_DESeq, "EC1", EC1.markers, startCol = 1, rowNames = T)
    addWorksheet(wb_DESeq, "EC2")
  writeData(wb_DESeq, "EC2", EC2.markers, startCol = 1, rowNames = T)
    addWorksheet(wb_DESeq, "VSMC")
  writeData(wb_DESeq, "VSMC", VSMC.markers, startCo = 1, rowNames = T)
    addWorksheet(wb_DESeq, "Fibroblast")
  writeData(wb_DESeq, "Fibroblast", Fibroblast.markers, startCol = 1, rowNames = T)
    addWorksheet(wb_DESeq, "Macrophage")
  writeData(wb_DESeq, "Macrophage", Macrophage.markers, startCol = 1, rowNames = T)
    addWorksheet(wb_DESeq, "Dendritic")
  writeData(wb_DESeq, "Dendritic", Dendritic.markers, startCol = 1, rowNames = T)
      addWorksheet(wb_DESeq, "NKT")
  writeData(wb_DESeq, "NKT", NKT.markers, startCol = 1, rowNames = T)
saveWorkbook(wb_DESeq, file = "../2_Output/2_TAA/3_Differential.Expression/DEGs_CellType.Specific.xlsx", overwrite = T)
# Heatmap of Clusters
DEGs_CellTypes %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
pdf(file = "../2_Output/2_TAA/2_Clustering/Celltype_Heatmap.pdf", height = 10, width = 15)
DoHeatmap(TAAng_scaled, features = top10$gene, size = 3, disp.min = -2, disp.max = 2) + scale_fill_gradientn(colors = c("dodgerblue4", "white", "coral2")) + labs(title = "Heatmap of Top-10 most variable genes within each cluster")
dev.off()
#
# Contractility score
contractility_genes <- c("CNN1", "TAGLN", "TAAng", "MYH11")
gene_list <- contractility_genes[contractility_genes %in% rownames(TAAng_scaled)] 
gene_expr <- GetAssayData(object = TAAng_scaled, assay = "RNA", slot = "data")[gene_list, ]
sum_expr <- colSums(gene_expr)
TAAng_scaled$contractility <- scale(sum_expr)
#osteogenic score
gene_list <- c("CBFA1", "MSX2", "RUNX2", "SOX9", "SPP1", "BGLAP", "ALPL", "COL2A1")
gene_list <- gene_list[gene_list %in% rownames(TAAng_scaled)]
gene_expr <- GetAssayData(object = TAAng_scaled, assay = "RNA", slot = "data")[gene_list, ]
sum_expr <- colSums(gene_expr)
TAAng_scaled$osteogenic <- scale(sum_expr)
# Synthetic score
gene_list <-c("CCND1","PCNA", "MKI67")
gene_list <- gene_list[gene_list %in% rownames(TAAng_scaled)]
gene_expr <- GetAssayData(object = TAAng_scaled, assay = "RNA", slot = "data")[gene_list, ]
sum_expr <- colSums(gene_expr)
TAAng_scaled$Synthetic <- scale(sum_expr)
# Save file for downstream analysis
saveRDS(TAAng_scaled, file = "../1_Input/TAAng_annotated.rds")
```

## VSMC-specific Trajectory Analysis (Monocle 3)


``` r
library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(ggtrace)
ann_colors = list(CellType = c(VSMC="coral2", 
                             Fibroblast="wheat", 
                             EC="steelblue4",
                             Macrophage="azure4",
                             NKT="gray",
                             Dendritic="tan2",
                             Fibromyocyte="darkcyan",
                             Lymphoid="darkgray",
                             `10`="black"))
TAAng_scaled<-readRDS(file = "../1_Input/TAAng_annotated.rds")
# Subsetting and re-normalization
TAAng_scaled <- subset(TAAng_scaled, idents = "VSMC") # select only the VSMCs
TAAng_scaled <- RunPCA(TAAng_scaled) # re-cluster
TAAng_scaled <- RunUMAP(TAAng_scaled, assay = "RNA", dims = 1:40) # re-cluster
Idents(TAAng_scaled) <- TAAng_scaled$RNA_snn_res.0.3
# Clustering
pdf(file = "../2_Output/2_TAA/4_Trajectory/UMAP_VSMC.pdf", height = 4, width = 4)
DimPlot(TAAng_scaled, group.by = "RNA_snn_res.0.3", reduction = "umap", label = T, repel = T,label.size = 4, cols = c("coral2", "coral4", "coral3"))+NoLegend()
dev.off()
#Contractility
library(ggplot2)
pdf(file = "../2_Output/2_TAA/4_Trajectory/Contractility.Score_UMAP.pdf", height = 4, width = 4)
FeaturePlot(TAAng_scaled, features = "contractility", reduction = "umap") + ggtitle("Contractility Score") +
scale_color_gradient2(low = "darkcyan", mid = "white", high = "darkred", midpoint = median(TAAng_scaled$contractility, na.rm = FALSE))
dev.off()
#Osteogenic
pdf(file = "../2_Output/2_TAA/4_Trajectory/Osteogenicity.Score_UMAP.pdf", height = 4, width = 4)
FeaturePlot(TAAng_scaled, features = "osteogenic", reduction = "umap") + ggtitle("Osteogenicity Score") +
scale_colour_gradient2(low = "darkcyan", mid = "white", high = "darkred", midpoint = median(TAAng_scaled$osteogenic, na.rm = FALSE))
dev.off()
#Synthetic
pdf(file = "../2_Output/2_TAA/4_Trajectory/Synthetic.Score_UMAP.pdf", height = 4, width = 4)
FeaturePlot(TAAng_scaled, features = "Synthetic", reduction = "umap", alpha = 0.6) + ggtitle("Synthetic Score") +
scale_colour_gradient2(low = "darkcyan", mid = "white", high = "darkred", midpoint = mean(TAAng_scaled$Synthetic, na.rm = FALSE))
dev.off()
# 
DefaultAssay(TAAng_scaled) = "RNA" # This changed in seurat V5
# TAAng_scaled <- ScaleData(object = TAAng_scaled)
cds <- as.cell_data_set(TAAng_scaled)
cds <- preprocess_cds(cds, num_dim = 100)
cds <- estimate_size_factors(cds)
# Include gene names (not done by default by the seurat conversion)
rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name
#  Run Monocle
cds <- cluster_cells(cds) # 0.000000002 This step creates "partitions" that are used in the trajectory inference
plot_cells(cds, show_trajectory_graph = FALSE, color_cells_by = "partition") # this shows the partitions overlain on the UMAP
cds <- learn_graph(cds, use_partition = TRUE) # creating trajectory inference within each partition
cds <- order_cells(cds) # Used when setting the nodes; if already known, use the next line
root_cell <- "TCATCAGTCAAGCTCC-1" # names(which(cds@principal_graph_aux$UMAP$pseudotime==0))
cds<-order_cells(cds, root_cells = root_cell)
# Plot the pseudotime on UMAP
plot_cells(cds, 
           color_cells_by = "pseudotime", 
           label_branch_points = FALSE, 
           label_leaves = FALSE)
pdf(file = "../2_Output/2_TAA/4_Trajectory/TAAng_UMAP_Trajectory_Partition.pdf", height = 3, width = 3)
plot_cells(cds, 
           color_cells_by = "partition", 
           label_branch_points = FALSE, 
           label_leaves = F,
           show_trajectory_graph = F,
           label_roots = F)
dev.off()
pdf(file = "../2_Output/2_TAA/4_Trajectory/TAAng_UMAP_Trajectory_Pseudotime.pdf", height = 3, width = 4)
plot_cells(cds, 
           color_cells_by = "pseudotime", 
           label_branch_points = FALSE, 
           label_leaves = F,
           show_trajectory_graph = F,
           label_roots = F)
dev.off()
# Examine specific genes
plot_cells(cds, 
           genes = c("ACTA2"),
           label_cell_groups = FALSE,
           show_trajectory_graph = FALSE,
           min_expr = 3)
# Identify pseudotime
modulated_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=8) # Identify differentially-expressed genes with pseudotime
modulated_genes <- na.omit(modulated_genes) # remove NA's
modulated_genes <- modulated_genes %>% filter(modulated_genes$q_value < 0.05 & modulated_genes$status =="OK") # filter cds results down
modulated_genes <- modulated_genes[order(-modulated_genes$morans_test_statistic), ] # order by moran's test
modulated_genes <- top_n(modulated_genes, 1000, -q_value)
#### Create a heatmap of genes with similar pseudotime kinetics
genes <- row.names(subset(modulated_genes, q_value < 0.05))
openxlsx::write.xlsx(modulated_genes, "../2_Output/1_TAAng/Trajectory/Pseudotime_DEGs.xlsx")
library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(circlize)
library(monocle3)
pt.matrix <- as.data.frame(exprs(cds)[match(genes,rownames(rowData(cds))),order(pseudotime(cds))])
#
cell_names <- colnames(pt.matrix)
Index<-as.data.frame(cds@colData) %>% dplyr::select(cell_type)
Index<-subset(Index, row.names(Index) %in% cell_names)
Index$CellType <- factor(Index$cell_type,
                         levels = c("VSMC", 
                                    "Fibroblast", 
                                    "EC1",
                                    "EC2",
                                    "Macrophage",
                                    "NKT",
                                    "Dendritic"))
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=6)$y})) # Create a spline that smooths the pseudotime-based expression along 6 degrees of freedom.
filtered_matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes
colnames(pt.matrix) <- cell_names
###########
paletteLength <- 20
myColor <- colorRampPalette(c("dodgerblue4", "white", "goldenrod1"))(paletteLength)
heatmap_DMC<-pheatmap::pheatmap(pt.matrix, scale="row", 
      cluster_cols = F, 
      cluster_rows = TRUE,
      cutree_rows = 4,
      fontsize_col = 8,
      color = myColor,
      annotation_col = Index,
      annotation_colors = ann_colors,
      show_colnames = F,
      show_rownames = F,
      # right_annotation = ha,
      border_color = NA)
################################
hc <-heatmap_DMC$tree_row
lbl <- cutree(hc, 4)
cluster1<-which(lbl==1)
cluster2<-which(lbl==2)
cluster3<-which(lbl==3)
cluster4<-which(lbl==4)
Cluster1_data<-pt.matrix[cluster1,]
Cluster2_data<-pt.matrix[cluster2,]
write.csv(Cluster2_data, "Cluster2.csv")
Cluster3_data<-pt.matrix[cluster3,]
Cluster4_data<-pt.matrix[cluster4,]
Cluster1_GENES <- rownames(Cluster1_data)
Cluster2_GENES <- rownames(Cluster2_data)
Cluster3_GENES <- rownames(Cluster3_data)
Cluster4_GENES <- rownames(Cluster4_data)
##Enrichr
library(enrichR)
library(dplyr)
dbs <- c("WikiPathway_2023_Human")
enriched_path__1 <- enrichr(Cluster1_GENES, dbs)
enrich_path_1<-enriched_path__1[[dbs]] %>% filter(Adjusted.P.value < 0.05)
head(enrich_path_1)
enriched_path__2 <- enrichr(Cluster2_GENES, dbs)
enrich_path_2<-enriched_path__2[[dbs]] %>% filter(Adjusted.P.value < 0.05)
head(enrich_path_2)
enriched_path__3 <- enrichr(Cluster3_GENES, dbs)
enrich_path_3<-enriched_path__3[[dbs]] %>% filter(Adjusted.P.value < 0.05)
head(enrich_path_3)
enriched_path__4 <- enrichr(Cluster4_GENES, dbs)
enrich_path_4<-enriched_path__4[[dbs]] %>% filter(Adjusted.P.value < 0.05)
head(enrich_path_4)
library(openxlsx)
wb_DESeq<-createWorkbook()
  addWorksheet(wb_DESeq, "Cluster 1")
  writeData(wb_DESeq, "Cluster 1", enrich_path_1, startCol = 1)
  addWorksheet(wb_DESeq, "Cluster 2")
  writeData(wb_DESeq, "Cluster 2", enrich_path_2, startCol = 1)
  addWorksheet(wb_DESeq, "Cluster 3")
  writeData(wb_DESeq, "Cluster 3", enrich_path_3, startCol = 1)
  addWorksheet(wb_DESeq, "Cluster 4")
  writeData(wb_DESeq, "Cluster 4", enrich_path_4, startCol = 1)
saveWorkbook(wb_DESeq, file = paste0("../2_Output/2_TAA/4_Trajectory/Trajectory_Pathway.xlsx"), overwrite = TRUE)
dbs <- c("GWAS_Catalog_2023")
enriched_1 <- enrichr(Cluster1_GENES, dbs)
enrich_1<-enriched_1[[dbs]] %>% filter(Adjusted.P.value < 0.05)
head(enrich_1)
enriched_2 <- enrichr(Cluster2_GENES, dbs)
enrich_2<-enriched_2[[dbs]] %>% filter(Adjusted.P.value < 0.05)
head(enrich_2)
enriched_3 <- enrichr(Cluster3_GENES, dbs)
enrich_3<-enriched_3[[dbs]] %>% filter(Adjusted.P.value < 0.05)
head(enrich_3)
enriched_4 <- enrichr(Cluster4_GENES, dbs)
enrich_4<-enriched_4[[dbs]] %>% filter(Adjusted.P.value < 0.05)
head(enrich_4)
library(openxlsx)
wb_DESeq<-createWorkbook()
  addWorksheet(wb_DESeq, "Cluster 1")
  writeData(wb_DESeq, "Cluster 1", enrich_1, startCol = 1)
  addWorksheet(wb_DESeq, "Cluster 2")
  writeData(wb_DESeq, "Cluster 2", enrich_2, startCol = 1)
  addWorksheet(wb_DESeq, "Cluster 3")
  writeData(wb_DESeq, "Cluster 3", enrich_3, startCol = 1)
  addWorksheet(wb_DESeq, "Cluster 4")
  writeData(wb_DESeq, "Cluster 4", enrich_4, startCol = 1)
saveWorkbook(wb_DESeq, file = paste0("../2_Output/2_TAA/4_Trajectory/Trajectory_GWAS.Traits.xlsx"), overwrite = TRUE)
library(stringr)
Aortic_genes <- enrich_1 %>% filter(str_detect(Term,"Aort"))
Aortic_genes <- Aortic_genes$Genes
Aortic_genes <- unique(unlist(strsplit(Aortic_genes, ";", fixed = T)))
# proliferative score
gene_list <-Aortic_genes
gene_list <- gene_list[gene_list %in% rownames(TAAng_scaled)]
gene_expr <- GetAssayData(object = TAAng_scaled, assay = "RNA", slot = "data")[gene_list, ]
sum_expr <- colSums(gene_expr)
TAAng_scaled$Aortic <- scale(sum_expr)
Cluster1_Names<-unlist(strsplit(enrich_1$Genes[1], ";", fixed = T))
Cluster2_Names<-unlist(strsplit(enrich_2$Genes[1], ";", fixed = T))
Cluster3_Names<-unlist(strsplit(enrich_3$Genes[1], ";", fixed = T))
Cluster4_Names<-unlist(strsplit(enrich_4$Genes[1], ";", fixed = T))
pdf(file = "../2_Output/2_TAA/4_Trajectory/Aortic_genes.pdf")
plot_cells(cds, 
           genes = Aortic_genes,
           label_cell_groups = F,
           label_roots = F,
           label_branch_points = F,
           label_leaves = F,
           show_trajectory_graph = F,
           min_expr = 1,
           alpha = 0.8,
           trajectory_graph_color = "grey28",
           ) +
  theme(
  axis.text = element_text(size = 6),    # Adjust the size as needed
  axis.title = element_text(size = 8),  # Adjust the size as needed
  legend.text = element_text(size = 4),  # Adjust the size as needed
  legend.title = element_text(size = 6),
  legend.key.size = unit(2, 'mm')) +
  scale_colour_viridis_c(option = "inferno")
dev.off()
###############################
GENES_HM<-c(Aortic_genes) #
pt.df<-as.data.frame(pt.matrix)
pt.df$gene_name<-rownames(pt.matrix)
ha = rowAnnotation(link = anno_mark(at = which(pt.df$gene_name %in% GENES_HM),
                   labels = as.character(pt.df[which(pt.df$gene_name %in% GENES_HM), "gene_name"]),
                   labels_gp = gpar(fontsize = 8),
                   padding = unit(1, "mm"))
                   )
heatmap_combination<-ComplexHeatmap::pheatmap(pt.matrix, scale="row",
                    cluster_cols = F,
                    cluster_rows = T,
                    cutree_rows = 4,
                    # cutree_cols = 3,
                     fontsize_col = 5,
                     color = myColor,
                    annotation_names_col = FALSE,
                    show_colnames = F,
                     show_rownames = F,
                     border_color = NA,
                    annotation_colors = NA,
                    right_annotation = ha,
                    annotation_col = Index,
                    border = TRUE)
pdf(file = "../2_Output/2_TAA/4_Trajectory/Trajectory_Heatmap.pdf", height = 8, width = 5)
heatmap_combination
dev.off()
```

# Integrated and harmonized Analysis (FFPE - TAA and ACTA2)

Once the individual snRNA-seq datasets were interrogated for quality and
normalized, they were fully integrated using standard protocol within
Seurat. To accomplish this, shared cell populations were matched across
datasets ("anchors") were used to correct for technical differences
between datasets (i.e. batch effects), as well as to perform comparative
sn-RNA-seq analysis across experimental conditions (ACTA2-mut vs. TAA
vs. CON).

Data from Chou et al. 2022 were incorporated to provide a more
generalizeable understanding of cellular heterogeneity at single-nuclear
resolution. Once the individual snRNA-seq datasets were interrogated for
quality and normalized, they were fully integrated using standard
protocol within Seurat. To accomplish this, shared cell populations were
matched across datasets ("anchors") were used to correct for technical
differences between datasets (i.e. batch effects), as well as to perform
comparative sn-RNA-seq analysis across experimental conditions.


``` r
options(future.globals.maxSize= 891289600000000000)
library(Seurat)
ACTA2_scaled <- readRDS("../1_Input/ACTA2_annotated.rds")
ACTA2_scaled$Disease <- "TAA_ACTA2"
control_scaled <- readRDS("../1_Input/control_annotated.rds")
control_scaled$Disease <- "TAA_NG"
ifnb <- merge(ACTA2_scaled, y = control_scaled, project = "ACTA2", assay = "RNA")
ifnb <- NormalizeData(ifnb)
ifnb <- FindVariableFeatures(ifnb)
ifnb <- ScaleData(ifnb)
ifnb <- RunPCA(ifnb)
ifnb <- FindNeighbors(ifnb, dims = 1:30, reduction = "pca")
ifnb <- FindClusters(ifnb, resolution = 2, cluster.name = "unintegrated_clusters")
ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
pdf("../2_Output/3_Merged/2_Clustering/Unharmonized_UMAP.pdf", height = 5, width = 10)
DimPlot(ifnb, reduction = "umap.unintegrated", group.by = c("Disease", "cell_type"))
dev.off()
# Harmonize
ifnb <- IntegrateLayers(object = ifnb, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
    verbose = FALSE)
ifnb[["RNA"]] <- JoinLayers(ifnb[["RNA"]]) # re-join layers after integration
ifnb <- FindNeighbors(ifnb, reduction = "integrated.cca", dims = 1:30)
ifnb <- FindClusters(ifnb, resolution = c(0.1,0.2,0.3,0.4,0.5))
ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "integrated.cca")
pdf("../2_Output/3_Merged/2_Clustering/Harmonized_UMAP.pdf", height = 4, width = 12)
DimPlot(ifnb, reduction = "umap", group.by = c("Disease", "cell_type", "RNA_snn_res.0.4"))
dev.off()
saveRDS(ifnb, "../1_Input/FFPE_integrated.rds")
```

### Cell Annotation


``` r
library(Seurat)
library(ggrepel)
library(RColorBrewer)
TAA_integrated<-readRDS(file = "../1_Input/FFPE_integrated.rds")
TAA_integrated$Disease <- factor(TAA_integrated$Disease, levels = c("TAA_NG", "TAA_ACTA2"))
colors = c("coral2",
           "coral3",
           "azure4", 
           "steelblue4",
           "wheat",
           "cadetblue3",
           "gray",
           "azure2",
           "cadetblue",
           "dodgerblue2",
           "black")
TAA_integrated[["RNA"]] <- JoinLayers(TAA_integrated[["RNA"]])
set2_colors <- brewer.pal(n = 8, name = "Set2")
UMAP_0.1<-DimPlot(TAA_integrated, reduction = "umap", label = T, repel = T,label.size = 4, group.by = "RNA_snn_res.0.1", split.by = "orig.ident", cols = ) + NoLegend()
UMAP_0.2<-DimPlot(TAA_integrated, reduction = "umap", label = T, repel = T,label.size = 4, group.by = "RNA_snn_res.0.2", split.by = "orig.ident", cols = colors) + NoLegend()
UMAP_0.3<-DimPlot(TAA_integrated, reduction = "umap", label = T, repel = T,label.size = 4, group.by = "RNA_snn_res.0.3", split.by = "orig.ident", cols = colors) + NoLegend()
UMAP_0.4<-DimPlot(TAA_integrated, reduction = "umap", label = T, repel = T,label.size = 4, group.by = "RNA_snn_res.0.4", split.by = "orig.ident", cols = colors) + NoLegend()
UMAP_0.5<-DimPlot(TAA_integrated, reduction = "umap", label = T, repel = T,label.size = 4, group.by = "RNA_snn_res.0.5", split.by = "orig.ident", cols = colors) + NoLegend()
pdf(file = "../2_Output/3_Merged/2_Clustering/UMAP.pdf", height = 5, width = 4)
UMAP_0.1
UMAP_0.2
UMAP_0.3
UMAP_0.4
UMAP_0.5
dev.off()
pdf(file = "../2_Output/3_Merged/2_Clustering/UMAP_pub.pdf")
DimPlot(TAA_integrated, reduction = "umap", label = T, repel = T,label.size = 4, group.by = "RNA_snn_res.0.4", split.by = "orig.ident", cols = colors) + NoLegend()
dev.off()
#Identify Clusters corresponding with known gene markers:
VSMC1_genes<-c("PKD1", "COL6A2", "PDLIM7", "FLNA", "SMTN")
VSMC1_available <- intersect(VSMC1_genes, rownames(TAA_integrated))
VSMC2_genes<-c("MYH11", "ACTA2", "ITGA8", "PRKG1", "CRISPLD1")
VSMC2_available <- intersect(VSMC2_genes, rownames(TAA_integrated))
Fibroblast1_genes<-c("ABCA10", "C3", "ADGRD1", "FBLN1", "DCN")
Fibroblast1_available <- intersect(Fibroblast1_genes, rownames(TAA_integrated))
Fibroblast2_genes<-c("NFASC",  "SAMD5", "PRSS23","UACA","TEX41") #
Fibroblast2_available <- intersect(Fibroblast2_genes, rownames(TAA_integrated))
Fibromyocyte_genes<-c("ADAMTS1", "RGS6", "TNC") # , "ANGPT2", "DGKG", "GRIP2"
Fibromyocyte_available <- intersect(Fibromyocyte_genes, rownames(TAA_integrated))
EC1_genes<-c("DIPK2B", "ARHGEF15", "STC1", "FLT1") # , "NOTCH4"
EC1_available <- intersect(EC1_genes, rownames(TAA_integrated))
EC2_genes<-c("VWF", "BMPER", "BMX", "NOS1")
EC2_available <- intersect(EC2_genes, rownames(TAA_integrated))
NKT_genes<-c("SKAP1", "RIPOR2", "ITGAL", "CD96") #  "RBPJ", "FYN",
NKT_available <- intersect(NKT_genes, rownames(TAA_integrated))
Macrophage_genes<-c("MRC1", "LGMN", "F13A1", "RBM47") #
Macrophage_available <- intersect(Macrophage_genes, rownames(TAA_integrated))
Dendritic_genes<-c("ITGAX", "S100A9", "CSF3R", "CXCL8") #
Dendritic_available <- intersect(Dendritic_genes, rownames(TAA_integrated))
B_GENES<-c("IGHM", "IGHD") #
B_GENES <- intersect(B_GENES, rownames(TAA_integrated))
T_GENES<-c("CD3","CD200","CCR7") #
T_GENES <- intersect(T_GENES, rownames(TAA_integrated))

# Plot density function
library(Nebulosa)
# VSMC1_density<-plot_density(TAA_integrated, VSMC1_available, joint = TRUE, combine = FALSE, pal = "magma")
VSMC1_density<-plot_density(TAA_integrated, VSMC1_available, joint = TRUE, combine = FALSE, pal = "magma")
VSMC2_density<-plot_density(TAA_integrated, VSMC2_available, joint = TRUE, combine = FALSE, pal = "magma")
Fibroblast1_density<-plot_density(TAA_integrated, Fibroblast1_available, joint = TRUE, combine = FALSE, pal = "magma")
Fibroblast2_density<-plot_density(TAA_integrated, Fibroblast2_available, joint = TRUE, combine = FALSE, pal = "magma")
Fibromyocyte_density<-plot_density(TAA_integrated, Fibromyocyte_available, joint = TRUE, combine = FALSE, pal = "magma")
EC1_density<-plot_density(TAA_integrated, EC1_available, joint = TRUE, combine = FALSE, pal = "magma")
EC2_density<-plot_density(TAA_integrated, EC2_available, joint = TRUE, combine = FALSE, pal = "magma")
Macrophage_density<-plot_density(TAA_integrated, Macrophage_available, joint = TRUE, combine = FALSE, pal = "magma")
NKT_density<-plot_density(TAA_integrated, NKT_available, joint = TRUE, combine = FALSE, pal = "magma")
Dendritic_density<-plot_density(TAA_integrated, Dendritic_available, joint = TRUE, combine = FALSE, pal = "magma")
B_density<-plot_density(TAA_integrated, B_GENES, joint = TRUE, combine = FALSE, pal = "magma")
# T_density<-plot_density(TAA_integrated, T_GENES, joint = TRUE, combine = FALSE, pal = "magma")
pdf(file = "../2_Output/3_Merged/2_Clustering/FFPE_Density_plots.pdf", height = 3, width = 3)
# VSMC1_density
VSMC1_density
VSMC2_density
Fibroblast1_density
Fibroblast2_density
EC1_density
EC2_density
Macrophage_density
NKT_density
Dendritic_density
Fibromyocyte_density
B_density
T_density
dev.off()
# Overlay these gene markers onto the UMAP to identify clusters
pdf(file = "../2_Output/3_Merged/2_Clustering/FFPE_FeaturePlots.pdf")
FeaturePlot(TAA_integrated, reduction = "umap", label = T, features = Fibroblast1_genes)
FeaturePlot(TAA_integrated, reduction = "umap", label = T, features = Fibroblast2_genes)
FeaturePlot(TAA_integrated, reduction = "umap", label = T, features = EC1_genes)
FeaturePlot(TAA_integrated, reduction = "umap", label = T, features = EC2_genes)
FeaturePlot(TAA_integrated, reduction = "umap", label = T, features = NKT_genes)
FeaturePlot(TAA_integrated, reduction = "umap", label = T, features = Macrophage_genes)
# FeaturePlot(TAA_integrated, reduction = "umap", label = T, features = VSMC1_genes)
FeaturePlot(TAA_integrated, reduction = "umap", label = T, features = VSMC2_genes)
FeaturePlot(TAA_integrated, reduction = "umap", label = T, features = Dendritic_genes)
dev.off()
#Create a figure of cell-type specific markers overlying the UMAP
pdf(file = "../2_Output/3_Merged/2_Clustering/CellType_Differentiation.pdf")
plot_density(TAA_integrated, c("ACTA2", "FBLN1", "F13A1", "VWF", "STC1", "NFASC", "ITGAL", "ITGAX"), joint = TRUE, combine = FALSE, pal = "magma")
dev.off()
# Create cluster labels based on the feature plot delineations of cellular phenotypes
Idents(TAA_integrated) <- TAA_integrated$RNA_snn_res.0.4
TAA_integrated <- RenameIdents(object = TAA_integrated, 
                          `0` = "VSMC_2", 
                          `1` = "VSMC_1", 
                          `2` = "Macrophage",
                          `3` = "Fibroblast",
                          `4` = "EC",
                          `5` = "Fibromyocyte",
                          `6` = "NKT",
                          `7` = "Dendritic",
                          `8` = "Fibromyocyte",
                          `9` = "B Lymph",
                          `10` = "T Lymph")
TAA_integrated$cell_type <- Idents(TAA_integrated) # redefine the cell-types based on the clustering of the integration, above
TAA_integrated$cell_type <- factor(TAA_integrated$cell_type, 
                                   levels = c(
                                    "VSMC_1", 
                                    "VSMC_2",
                                    "Fibromyocyte",
                                    "Fibroblast", 
                                    "EC",
                                    "Macrophage",
                                    "NKT",
                                    "Dendritic",
                                    "B Lymph",
                                    "T Lymph")
                                   )
Idents(TAA_integrated) <- TAA_integrated$cell_type
colors_cell = c("coral2",
           "coral3",
          "cadetblue3",
           "steelblue4",
           "wheat",
           "azure4", 
           "gray",
           "azure2",
           "dodgerblue2",
           "black")
# Create the UMAP of cell-type specific clusters
UMAP_Disease<-DimPlot(TAA_integrated,  reduction = "umap",label = F, group.by = "Disease",label.size = 4,
                        cols = c("steelblue4",
                                 "goldenrod2"))
UMAP_CellTypes<-DimPlot(TAA_integrated,  reduction = "umap",label = T, label.size = 4,
                        cols = colors_cell) + NoLegend()
pdf("../2_Output/3_Merged/2_Clustering/UMAP_Disease.pdf", height = 5, width = 5)
UMAP_Disease
dev.off()
pdf("../2_Output/3_Merged/2_Clustering/UMAP_CellTypes.pdf", height = 4, width = 4)
UMAP_CellTypes
dev.off()
# Differential Expression
DEGs_CellTypes<-FindAllMarkers(TAA_integrated)
write.csv(DEGs_CellTypes, "../2_Output/3_Merged/3_Differential.Expression/DEGs_Clusters.csv")
# Create a dot-bplot of the top 5 markers for each 
library(scCustomize)
paletteLength <- 100
myColor <- colorRampPalette(c("dodgerblue4", "white", "coral2"))(paletteLength)
top5_markers <- Extract_Top_Markers(marker_dataframe = DEGs_CellTypes, num_genes = 10, named_vector = FALSE,
    make_unique = TRUE)
pdf(file = "../2_Output/3_Merged/2_Clustering/Dot.plot_top5markers.pdf", height = 10, width = 7)
Clustered_DotPlot(seurat_object = TAA_integrated, features = top5_markers, k = 10, colors_use_exp = myColor, assay = "RNA")
dev.off()
# Heatmap of Clusters
library(dplyr)
DEGs_CellTypes %>%
    group_by(cluster) %>%
    top_n(n = 50, wt = avg_log2FC) -> top10
write.csv(top10, "../2_Output/3_Merged/3_Differential.Expression/top10.csv")
pdf(file = "../2_Output/3_Merged/2_Clustering/Celltype_Heatmap.pdf", height = 10, width = 15)
DoHeatmap(TAA_integrated, features = top10$gene, size = 3, disp.min = -2, disp.max = 2) + scale_fill_gradientn(colors = c("dodgerblue4", "white", "coral2")) + labs(title = "Heatmap of Top-10 most variable genes within each cluster")
dev.off()
# Probability distribution of specific genes across all clusters/cell types
VlnPlot(TAA_integrated, features = c("VWF", "HDAC9", "ACTA2"))
DotPlot(TAA_integrated, features = VSMC1_genes) #+ RotatedAxis()
saveRDS(TAA_integrated, file = "../1_Input/TAA_snRNA.rds")
```

## Proportional Comparison

## UMAP and Cellular Proportions


``` r
# Cell-type Specific Differential Expression
library(ggplot2)
library(dplyr)
library(Seurat)
library(dittoSeq)
TAA.combined <- readRDS(file = "../1_Input/TAA_snRNA.rds")
#Figure 2C -  Proportional Graph
pdf(file = "../2_Output/3_Merged/2_Clustering/Proportional.Bar_Treatment.pdf", height = 4, width = 4)
dittoBarPlot(
  object = TAA_integrated,
  var = "cell_type",
  group.by = "Disease"
) + 
  ggtitle(NULL) +
  scale_fill_manual(values = 
    c(VSMC_1="coral2",
     VSMC_2="coral3",
     Fibroblast="steelblue4", 
     EC="wheat",
     Macrophage="azure4",
     NKT="gray",
     Dendritic="tan2",
     Fibromyocyte="cadetblue3",
     `B Lymph`="darkgray",
     `T Lymph`="black")
    )
dev.off()

BarPlot <- dittoBarPlot(
    object = TAA.combined,
    var = "cell_type",
    group.by = "Disease", 
    data.out = T)
# Extract Data for export
BarPlot_sampledata <- BarPlot$data %>% 
  select(-count, -label.count.total.per.facet) %>% 
  tidyr::pivot_wider(., names_from = "grouping", values_from = "percent")
openxlsx::write.xlsx(BarPlot_sampledata, "../2_Output/BarPlot_CellTypes.xlsx")
# Visualize the features/genes
VizDimLoadings(TAA.combined, dims = 1:5, reduction = "pca")
pdf(file = "../2_Output/PC_Heatmaps.pdf")
DimHeatmap(TAA.combined, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()
cell_type_disease_counts <- TAA.combined@meta.data %>%
  group_by(Disease, cell_type = Idents(TAA.combined)) %>%
  summarise(cell_count = n())
write.csv(cell_type_disease_counts, "cell_types_Disease.csv")
```

# Differential Analysis


``` r
library(Seurat)
library(dplyr)
library(ggtrace)
library(ggplot2)
library(ggrepel)
TAA.combined<-readRDS(file = "../1_Input/TAA_snRNA.rds")
Idents(TAA.combined) <- "Disease"

#Volcano Plot
# Load packages
library(dplyr)
library(ggplot2)
library(ggrepel)
library(openxlsx)
library(Nebulosa)
# Read data from the web
options(ggrepel.max.overlaps = Inf)
cell_types <- unique(TAA.combined$cell_type)
# Loop through each cell type
for (CELL in cell_types) {
  # Subset the Seurat object based on the cell type
  subset_seurat <- subset(TAA.combined, subset = cell_type == CELL)
  # Generate a UMAP plot for the subsetted cell type
  DimPlot(subset_seurat, reduction = "umap", label = TRUE, repel = TRUE, label.size = 4) + 
    NoLegend()
  # Differential Expression Analysis for the current cell type
  de_markers <- FindMarkers(subset_seurat, ident.1 = "TAA_ACTA2", ident.2 = "TAA_NG", min.pct = 0.25)
  # Save the results to an Excel file
  output_file <- paste0("../2_Output/3_Merged/3_Differential.Expression/", CELL, "_DEGs.xlsx")
  openxlsx::write.xlsx(de_markers, output_file, rowNames = TRUE)
    ##############
  # Violin Plots (top 5)
  ###############
  top_degs.up <- de_markers %>% filter(pct.1>0.5 | pct.2>0.5) %>% top_n(., n = 5, avg_log2FC)
  top_degs.down <- de_markers %>% filter(pct.1>0.5 | pct.2>0.5) %>% top_n(., n = 5, -avg_log2FC)
  plots_up <- VlnPlot(subset_seurat, features = rownames(top_degs.up), split.by = "Disease",
      pt.size = 0, combine = FALSE)
  plots_up <- lapply(seq_along(plots_up), function(i) {
    plot <- plots_up[[i]] + theme(axis.title.x = element_blank())  # Remove x-axis label
    if (i == length(plots_up)) {
      return(plot)  # Keep the legend in the last plot
    } else {
      return(plot + NoLegend())  # Remove legend in other plots
    }
  })
  plots_down <- VlnPlot(subset_seurat, features = rownames(top_degs.down), split.by = "Disease",
      pt.size = 0, combine = FALSE)
  plots_down <- lapply(seq_along(plots_down), function(i) {
    plot <- plots_down[[i]] + theme(axis.title.x = element_blank())  # Remove x-axis label
    if (i == length(plots_down)) {
      return(plot)  # Keep the legend in the last plot
    } else {
      return(plot + NoLegend())  # Remove legend in other plots
    }
  })
  pdf(paste0("../2_Output/3_Merged/3_Differential.Expression/", CELL, "_TAA_ACTA2.v.NG_UP.pdf"), height = 3, width = 10)
  print(wrap_plots(plots = plots_up, ncol = 5))
  dev.off()
  pdf(paste0("../2_Output/3_Merged/3_Differential.Expression/", CELL, "_TAA_ACTA2.v.NG_DOWN.pdf"), height = 3, width = 10)
  print(wrap_plots(plots = plots_down, ncol = 5))
  dev.off()
  # Prepare data for Volcano Plot
    max_pval <- min(de_markers$p_val[de_markers$p_val != 0])  # Calculate the highest non-zero p_val value
    results <- de_markers %>%
      mutate(p_val = ifelse(p_val == 0, max_pval, p_val)) %>%  # Replace p_val == 0 with max_pval
      mutate(minuslogpvalue = -log(p_val), log2FC = avg_log2FC) %>%
      filter(p_val != 0) %>%
      mutate(sig = ifelse(p_val < 10^-8 & log2FC > 2, "P < 10^-8 and Fold-Change > 2", ifelse(p_val < 10^-8 & log2FC < -2, "P < 10^-8 and Fold-Change < -2", "Not Sig"))) %>%
      mutate(gene_name = rownames(.))  # Add gene_name as a new column
  #############
  #### Volcano Plot (TAA_ACTA2 vs. TAA_NG)
  ############
p <- ggplot(results, aes(log2FC, minuslogpvalue)) + 
  theme_classic() +
  geom_point(aes(fill = sig, size = minuslogpvalue), colour = "black", shape = 21, stroke = 0, alpha = 7/10) +
  geom_segment(aes(x = 2, xend = 2, y = 0, yend = max(results$minuslogpvalue, na.rm = TRUE)), size = 0.5, linetype = "dashed", lineend = "round") +  # Vertical line at x = 2
  geom_segment(aes(x = -2, xend = -2, y = 0, yend = max(results$minuslogpvalue, na.rm = TRUE)), size = 0.5, linetype = "dashed", lineend = "round") +  # Vertical line at x = -2
  geom_segment(aes(x = -5, xend = 5, y = 0 - log(10^-8), yend = 0 - log(10^-8)), size = 0.5, linetype = "dashed", lineend = "round") +
  geom_segment(aes(x = -5, xend = 5, y = 0, yend = 0), size = 1, linetype = "solid", lineend = "round") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(results$minuslogpvalue, na.rm = TRUE)), size = 1, linetype = "solid", lineend = "round") +  # Vertical line at x = 0
  labs(x = expression(Log[2](Fold-Change)), y = expression(-Log[10](P-value))) + 
  scale_fill_manual(values = c("darkgray", "dodgerblue4", "goldenrod2")) +
  scale_size_continuous(range = c(0.1, 3)) +
  theme(axis.title.x = element_text(color = "black", size = 12),  # X-axis label black and larger
        axis.title.y = element_text(color = "black", size = 12),  # Y-axis label black and larger
        axis.text = element_text(color = "black", size = 10),  # Hide axis text
        axis.ticks = element_blank(),  # Hide axis ticks
        legend.position = "none")  # Hide legend
  # Annotate significant genes
  p <- p + 
    geom_text_repel(data = top_n(filter(results, log2FC < -2), 10, minuslogpvalue), aes(label = gene_name)) +
    geom_text_repel(data = top_n(filter(results, log2FC > 2), 10, minuslogpvalue), aes(label = gene_name))
  # Save the Volcano Plot
  pdf_file <- paste0("../2_Output/3_Merged/3_Differential.Expression/", CELL, "_VolcanoPlot.pdf")
  pdf(pdf_file, height = 5, width = 8)
  print(p)
  dev.off()
  ##Enrichr
  #### Pathways
  library(enrichR)
  library(dplyr)
  DEGs_up <- de_markers %>% filter(avg_log2FC > 0) %>% top_n(500, abs(avg_log2FC))
  DEGs_down <- de_markers %>% filter(avg_log2FC < 0) %>% top_n(500, abs(avg_log2FC))  
  #############
  #### Wikipathways Enrichment
  ############
  dbs <- c("WikiPathway_2023_Human")
  enriched_up <- enrichr(rownames(DEGs_up), dbs)
  enrich_up<-enriched_up[[dbs]] %>% filter(Adjusted.P.value < 0.05)
  head(enrich_up)
  enriched_down <- enrichr(rownames(DEGs_down), dbs)
  enrich_down<-enriched_down[[dbs]] %>% filter(Adjusted.P.value < 0.05)
  head(enrich_down)
  library(openxlsx)
  file_path <- "../2_Output/3_Merged/3_Differential.Expression/Wikipathways_All.Celltypes_TAAng.vs.ACTA2.xlsx"
      if (file.exists(file_path)) { # Check if the file exists and load or create the workbook accordingly
        wb <- loadWorkbook(file_path)  # Load existing workbook
        message("Workbook loaded successfully.")
        existing_sheets <- getSheetNames(file_path)  # Get existing sheet names using getSheetNames() only if the file exists
      } else {
        wb <- createWorkbook()  # Create new workbook
        message("New workbook created.")
        existing_sheets <- character(0)  # Empty vector indicating no sheets exist
      }
      if (!(paste0(CELL, "_UP") %in% existing_sheets)) { # Adding sheets conditionally based on their existence
        addWorksheet(wb, sheetName = paste0(CELL, "_UP"))
        writeData(wb, sheet = paste0(CELL, "_UP"), enrich_up, startCol = 1)
      } else {
        message(paste("Sheet", paste0(CELL, "_UP"), "already exists. Skipping."))
      }
      if (!(paste0(CELL, "_DOWN") %in% existing_sheets)) {
        addWorksheet(wb, sheetName = paste0(CELL, "_DOWN"))
        writeData(wb, sheet = paste0(CELL, "_DOWN"), enrich_down, startCol = 1)
      } else {
        message(paste("Sheet", paste0(CELL, "_DOWN"), "already exists. Skipping."))
      }
      saveWorkbook(wb, file_path, overwrite = TRUE) # Save the workbook
      message("Workbook saved successfully.")
  #############
  #### GWAS Enrichment
  ############
  dbs <- c("GWAS_Catalog_2023")
  enriched_up <- enrichr(rownames(DEGs_up), dbs)
  enrich_up<-enriched_up[[dbs]] %>% filter(Adjusted.P.value < 0.05)
  head(enrich_down)
  enriched_down <- enrichr(rownames(DEGs_down), dbs)
  enrich_down<-enriched_down[[dbs]] %>% filter(Adjusted.P.value < 0.05)
  head(enrich_up)
  file_path <- "../2_Output/3_Merged/3_Differential.Expression/GWAS_All.Celltypes_TAAng.vs.ACTA2.xlsx"
      if (file.exists(file_path)) { # Check if the file exists and load or create the workbook accordingly
        wb <- loadWorkbook(file_path)  # Load existing workbook
        message("Workbook loaded successfully.")
        existing_sheets <- getSheetNames(file_path)  # Get existing sheet names using getSheetNames() only if the file exists
      } else {
        wb <- createWorkbook()  # Create new workbook
        message("New workbook created.")
        existing_sheets <- character(0)  # Empty vector indicating no sheets exist
      }
      if (!(paste0(CELL, "_UP") %in% existing_sheets)) { # Adding sheets conditionally based on their existence
        addWorksheet(wb, sheetName = paste0(CELL, "_UP"))
        writeData(wb, sheet = paste0(CELL, "_UP"), enrich_up, startCol = 1)
      } else {
        message(paste("Sheet", paste0(CELL, "_UP"), "already exists. Skipping."))
      }
      if (!(paste0(CELL, "_DOWN") %in% existing_sheets)) {
        addWorksheet(wb, sheetName = paste0(CELL, "_DOWN"))
        writeData(wb, sheet = paste0(CELL, "_DOWN"), enrich_down, startCol = 1)
      } else {
        message(paste("Sheet", paste0(CELL, "_DOWN"), "already exists. Skipping."))
      }
      saveWorkbook(wb, file_path, overwrite = TRUE) # Save the workbook
      message("Workbook saved successfully.")
}
```

# Comparing DEGs across Cell Types


``` r
library(Seurat)
library(dplyr)
TAA.combined<-readRDS(file = "../1_Input/TAA_snRNA.rds")
cell_states <- unique(TAA.combined$cell_type)  # Get unique disease states
deg_list <- list()
# Loop through each disease state to find DEGs
for (CELL in cell_states) {
  markers <- openxlsx::read.xlsx(paste0("../2_Output/3_Merged/3_Differential.Expression/", CELL, "_DEGs.xlsx"), rowNames = T) %>% filter(p_val < 0.05)
  deg_list[[CELL]] <- rownames(markers[markers$p_val_adj < 0.05, ])  # Extract DEGs with p-value < 0.05
}
### Cell Type Comparison
# Create a data frame from the list for UpSetR
deg_matrix <- data.frame(
  gene = unique(unlist(deg_list)),  # Unique genes
  stringsAsFactors = FALSE
)
# Add binary presence/absence indicators for each cell type
for (cell_type in names(deg_list)) {
  deg_matrix[[cell_type]] <- ifelse(deg_matrix$gene %in% deg_list[[cell_type]], 1, 0)
}
desired_order <- rev(c("VSMC_1", "VSMC_2", "Fibromyocyte", "Fibroblast", "Macrophage", "Dendritic", "NKT", "EC"))
t <- upset(deg_matrix, 
      sets = desired_order,  # Names of the cell types
      # text.scale = c(2, 1.5, 2, 2, 2, 2),
      mb.ratio = c(0.5, 0.5),
      keep.order = TRUE,  # Keep the order of sets as provided
      main.bar.color = "#56B4E9",  # Color for main bars
      sets.bar.color = "#D55E00",  # Color for sets bars
      order.by = "freq")  # Order by frequency of intersections
t
pdf("../2_Output/3_Merged/3_Differential.Expression/Venn_UpSetR.pdf", height = 4, width = 7)
print(t)
dev.off()
# Load libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(ReactomePA)
for (CELL in cell_states) {
  markers <- openxlsx::read.xlsx(paste0("../2_Output/3_Merged/3_Differential.Expression/", CELL, "_DEGs.xlsx"), rowNames = T) %>% top_n(n = 100, wt = abs(avg_log2FC))
  deg_list[[CELL]] <- rownames(markers[markers$p_val_adj < 0.05, ])  # Extract DEGs with p-value < 0.05
}
# Prepare DEG list (Example DEG list structure)
deg_list <- deg_list[!names(deg_list) %in% c("B Lymph","T Lymph")]
common_genes <- Reduce(intersect, deg_list)
# Convert gene symbols to Entrez IDs
deg_list_entrez <- lapply(deg_list, function(genes) {
  bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
})
# Perform Reactome pathway enrichment analysis for each cell type
reactome_results <- lapply(deg_list_entrez, function(entrez_genes) {
  enrichPathway(gene = entrez_genes, organism = "human", qvalueCutoff = 0.05)
})
# Assign cell type names
names(reactome_results) <- names(deg_list)
# Create a compareCluster object for Reactome pathway enrichment
compare_cluster_results <- compareCluster(
  geneCluster = deg_list_entrez,
  fun = "enrichPathway",
  organism = "human",
  pvalueCutoff = 0.05
)
# compare_cluster_results <- compareCluster(
#   geneCluster = deg_list_entrez,  # List of genes for each cell type in Entrez ID format
#   fun = "enrichKEGG",             # Use "enrichKEGG" for KEGG pathway analysis
#   organism = "hsa",               # Specify organism as "hsa" for human in KEGG
#   pvalueCutoff = 0.05             # p-value cutoff for pathway enrichment
# )
# compare_cluster_results <- compareCluster(
#   geneCluster = deg_list_entrez,  # List of genes for each cell type in Entrez ID format
#   fun = "enrichDO",             # Use "enrichKEGG" for KEGG pathway analysis
#   organism = "hsa",               # Specify organism as "hsa" for human in KEGG
#   pvalueCutoff = 0.05             # p-value cutoff for pathway enrichment
# )
# compare_cluster_results <- compareCluster(
#   geneCluster = deg_list_entrez,
#   fun = "enrichWP",
#   organism = "Homo sapiens",
#   pvalueCutoff = 0.05
# )
# compare_cluster_results <- compareCluster(
#   geneCluster = deg_list_entrez,
#   fun = "enrichDGN",
#   pvalueCutoff = 0.01
# )

###############

compare_cluster_results@compareClusterResult$cell_type <- compare_cluster_results@compareClusterResult$Cluster
pdf("../2_Output/3_Merged/3_Differential.Expression/DotPlot_SPLIT.pdf", height = 4, width = 11)
dotplot(compare_cluster_results, showCategory = 3, x = "GeneRatio", split = "cell_type", label_format = 50) +
  facet_grid(. ~ cell_type) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10))
dev.off()

# Visualize using dotplot
pdf("../2_Output/3_Merged/3_Differential.Expression/DotPlot.pdf", height = 4, width = 9)
dotplot(compare_cluster_results, color = "p.adjust", showCategory = 2) + 
  scale_x_discrete(labels = function(x) gsub("\\s*\\(\\d+\\)$", "", x)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        xis.text.y = element_text(size = 6))
dev.off()

pathway_genes <- compare_cluster_results@compareClusterResult %>%
  filter(Description == "Extracellular matrix organization") %>%
  pull(geneID) %>%
  strsplit("/") %>%
  unlist()
# Step 2: Convert Entrez IDs to gene symbols
converted_genes <- bitr(
  pathway_genes,
  fromType = "ENTREZID",
  toType = "SYMBOL",
  OrgDb = org.Hs.eg.db
)
gene_symbols <- converted_genes$SYMBOL
genes_curated <- c("TIMP1", "COL8A1", "EFEMP1", "COL14A1", "SERPINE1", "ADAMTS1", "VCAN", "ITGA7", "A2M", "ITGA8", "MMP2", "HSPG2")
 # Violin Plots (top 5)
  ###############
library(Seurat)
library(patchwork)
  plots_up <- VlnPlot(TAA.combined, features = common_genes, split.by = "Disease", group.by = "cell_type",
      pt.size = 0, combine = FALSE)
  plots_up <- lapply(seq_along(plots_up), function(i) {
    plot <- plots_up[[i]] + theme(axis.title.x = element_blank())  # Remove x-axis label
    if (i == length(plots_up)) {
      return(plot)  # Keep the legend in the last plot
    } else {
      return(plot + NoLegend())  # Remove legend in other plots
    }
  })
  pdf(paste0("../2_Output/3_Merged/3_Differential.Expression/Conserved.pdf"), height = 6, width = 10)
  print(wrap_plots(plots = plots_up, ncol = 3))
  dev.off()
```

# Cell-cell Interactions


``` r
# Load necessary libraries
library(Seurat)
library(tidyverse)
library(magrittr)
library(liana)
library(openxlsx)
library(circlize)
library(ComplexHeatmap)
# Create folder structure for output if it does not already exist
ifelse(!dir.exists(file.path("../2_Output/Regulation/")), dir.create(file.path("../2_Output/Regulation/")), FALSE)

# Import the Seurat object containing the ACTA2-mutant sample
acta2_scaled <- readRDS(file = "../1_Input/ACTA2_annotated.rds")

# Subset data to include only relevant cell types if needed, e.g., VSMCs, ECs, Myeloid cells, etc.
# Uncomment and modify the line below if specific cell types need to be isolated
# acta2_scaled <- subset(acta2_scaled, idents = c("VSMC", "EC", "Myeloid", "Fibroblast"))

# Run LIANA to identify potential cell-cell interactions
liana_results <- liana_wrap(acta2_scaled) # Use the appropriate resource if working with human data
liana_results <- liana_aggregate(liana_results)

# Filter LIANA results for significant interactions (e.g., p-value <= 0.05)
liana_filtered <- liana_results %>% filter(cellphonedb.pvalue <= 0.05)

# Save filtered results to an Excel file
openxlsx::write.xlsx(liana_filtered, "../2_Output/Regulation/LIANA_ACTA2_Mutant_Significant.xlsx")

# Check unique targets
unique_targets <- unique(liana_results$target)
print(unique_targets)  # Print to ensure it's populated correctly

# Create a heatmap of interaction frequencies
pdf(file = "../2_Output/Regulation/LIANA_Heatmap_ACTA2_Mutant.pdf", height = 4, width = 5)
heat_freq(liana_filtered)
dev.off()

# Visualize cell-cell interactions using a dot plot
pdf(file = "../2_Output/Regulation/LIANA_DotPlot_ACTA2_Mutant.pdf", height = 7, width = 7)
liana_results %>%
  liana_dotplot(source_groups = c("VSMC_2"),  # Define the source cell type(s)
                target_groups = unique_targets[!is.na(unique_targets)],  # Remove NA if any
                ntop = 20) + 
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1), 
        plot.title = element_text(size = 0), 
        axis.title.x = element_text(size = 0))
dev.off()

# Define color scheme for each cell type, including all targets in unique_targets
grid.col <-c(VSMC_1="#b95e4b",
             VSMC_2="coral4",
             VSMC_3="coral3",
             Fibroblast="wheat", 
             EC1="steelblue4",
             EC2="deepskyblue3",
             Macrophage="azure4",
             NKT="gray",
             Dendritic="tan2")

# Aggregate interactions to create a weight for each source-target pair
liana_data_aggregated <- liana_filtered %>%
  mutate(weight = -log10(aggregate_rank)) %>%  # Convert rank to weight for visualization
  group_by(source, target) %>%
  summarise(weight = sum(weight, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(!is.infinite(weight) & weight > 0)  # Exclude infinite or zero weights

# Convert to data frame for chord diagram compatibility
liana_data_aggregated <- as.data.frame(liana_data_aggregated)

# Save the chord diagram to a PDF
pdf(file = "../2_Output/Regulation/LIANA_ChordDiagram_ACTA2_Mutant.pdf", height = 6, width = 9)
circlize::chordDiagram(
  x = liana_data_aggregated,
  grid.col = grid.col,
  transparency = 0.5,
  directional = 1,
  direction.type = c("diffHeight", "arrows"),
  link.arr.type = "big.arrow",
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.1)
)

# Add cell type labels with custom spacing
circlize::circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1] + 0.3, CELL_META$sector.index, facing = "bending", niceFacing = TRUE, adj = c(0.5, 0.5))
}, bg.border = NA)

# Create and draw a legend for cell types
legend <- ComplexHeatmap::Legend(
  at = names(grid.col), 
  title = "Cell Type", 
  legend_gp = gpar(fill = grid.col)
)
ComplexHeatmap::draw(legend, x = unit(1, "npc") - unit(5, "mm"), just = "right")

# Close the PDF device
dev.off()
library(Nebulosa)
plot_density(acta2_scaled, "MYL9", joint = TRUE, combine = FALSE, pal = "magma")
FeaturePlot(acta2_scaled, features = "LGALS1", reduction = "umap") +
    scale_color_viridis_c(option = "magma") +  # Apply the "mma" color scale
    ggtitle("UMAP Plot Showing Expression of LGALS1") +
    theme_minimal()  # Optional: Apply a minimal theme for a cleaner appearance
```

# ShinyR App Development


``` r
library(Seurat)
library(ShinyCell)
TAA.combined<-readRDS(file = "../1_Input/TAA_snRNA.rds")
library(clusterProfiler)
library(org.Mm.eg.db)  # Mouse gene database
# Define the mouse GO term of interest, e.g., GO:0008150 (biological process)
go_term <- "GO:0008283" # used in the original 2023 paper
# Extract gene symbols associated with the specified GO term
genes_vector <- AnnotationDbi::select(
  org.Mm.eg.db,                 # Mouse gene database
  keys = go_term,               # GO term of interest
  columns = "SYMBOL",           # Retrieve the gene symbols
  keytype = "GOALL"             # Define the key type as GO terms
)$SYMBOL
# Proliferation Score
genes_vector <- toupper(genes_vector)
gene_list <- genes_vector[genes_vector %in% rownames(TAA.combined)] 
gene_expr <- GetAssayData(object = TAA.combined, slot = "data")[gene_list, ]
sum_expr <- colSums(gene_expr)
TAA.combined$Proliferation <- scale(sum_expr)
# Load required libraries
library(Seurat)
library(ShinyCell)

# Make a config file
scConf <- createConfig(
  obj = TAA.combined,                # Your Seurat object
  meta.to.include = c("cell_type", "Disease", "Proliferation", "RNA_snn_res.0.4"),  # Metadata columns to include
  legendCols = 4,                    # Number of columns in categorical metadata legends
  maxLevels = 50                     # Max number of levels for categorical metadata
)
makeShinyApp(
  obj = TAA.combined,               # Seurat object or file path to the dataset
  scConf = scConf,                  # Configuration data.table generated in Step 1
  default.dimred = c("umap_1", "umap_2"),
  shiny.footnotes = "Pepin et al. 2024",
  gex.assay = "RNA",                # Assay in Seurat object to use (e.g., "RNA", "integrated")
  gex.slot = "data",                # Slot in the Seurat assay to use ("data" is default for normalized expression)
  shiny.title = "ACTA2-p.Met49Thr: snRNA-Sequencing Analysis",  # Title for the Shiny app
  shiny.dir = "ShinyCellApp",       # Directory to create the app files
  enableSubset = TRUE,              # Enable subsetting cells functionality in the app
  defPtSiz = 1.25,                  # Default point size for single cells in plots
  default.gene1 = "MYH11",          # Primary default gene to show in feature plots
  default.gene2 = "KLF4"            # Secondary default gene for comparison
)
```

# Supplemental Table: R Session Information

All packages and setting are acquired using the following command.


``` r
# # unloading system variables after session
# homer_tE<-Sys.time()
# homer_t<-homer_tE - homer_tS
# homer_t
# 
# end_t<-Sys.time()
# Total_time<-end_t - start_t
# Total_time
# Sys.unsetenv("motifLoc_UP")
# Sys.unsetenv("motifLoc_DOWN")
# Sys.unsetenv("motifLoc_WD")
# Write 
options(kableExtra.latex.load_packages = FALSE)
library(kableExtra)
sinfo<-devtools::session_info()
sinfo$platform
```

```
##  setting  value
##  version  R version 4.4.2 (2024-10-31)
##  os       macOS Sequoia 15.2
##  system   aarch64, darwin20
##  ui       X11
##  language (EN)
##  collate  en_US.UTF-8
##  ctype    en_US.UTF-8
##  tz       America/Los_Angeles
##  date     2025-02-01
##  pandoc   3.4 @ /opt/anaconda3/bin/ (via rmarkdown)
```

``` r
sinfo$packages %>% kable(
                         align="c",
                         longtable=T,
                         booktabs=T,
                         caption="Packages and Required Dependencies") %>%
    kable_styling(latex_options=c("striped", "repeat_header", "condensed"))
```

<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
<caption>Packages and Required Dependencies</caption>
 <thead>
  <tr>
   <th style="text-align:left;">  </th>
   <th style="text-align:center;"> package </th>
   <th style="text-align:center;"> ondiskversion </th>
   <th style="text-align:center;"> loadedversion </th>
   <th style="text-align:center;"> path </th>
   <th style="text-align:center;"> loadedpath </th>
   <th style="text-align:center;"> attached </th>
   <th style="text-align:center;"> is_base </th>
   <th style="text-align:center;"> date </th>
   <th style="text-align:center;"> source </th>
   <th style="text-align:center;"> md5ok </th>
   <th style="text-align:center;"> library </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> abind </td>
   <td style="text-align:center;"> abind </td>
   <td style="text-align:center;"> 1.4.8 </td>
   <td style="text-align:center;"> 1.4-8 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/abind </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/abind </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-12 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AnnotationDbi </td>
   <td style="text-align:center;"> AnnotationDbi </td>
   <td style="text-align:center;"> 1.68.0 </td>
   <td style="text-align:center;"> 1.68.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/AnnotationDbi </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/AnnotationDbi </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-11-08 </td>
   <td style="text-align:center;"> Bioconductor 3.20 (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ape </td>
   <td style="text-align:center;"> ape </td>
   <td style="text-align:center;"> 5.8 </td>
   <td style="text-align:center;"> 5.8 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ape </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ape </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-04-11 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> aplot </td>
   <td style="text-align:center;"> aplot </td>
   <td style="text-align:center;"> 0.2.3 </td>
   <td style="text-align:center;"> 0.2.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/aplot </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/aplot </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-06-17 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> beeswarm </td>
   <td style="text-align:center;"> beeswarm </td>
   <td style="text-align:center;"> 0.4.0 </td>
   <td style="text-align:center;"> 0.4.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/beeswarm </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/beeswarm </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-06-01 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Biobase </td>
   <td style="text-align:center;"> Biobase </td>
   <td style="text-align:center;"> 2.65.1 </td>
   <td style="text-align:center;"> 2.65.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/Biobase </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/Biobase </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-08-28 </td>
   <td style="text-align:center;"> Bioconductor 3.20 (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> BiocGenerics </td>
   <td style="text-align:center;"> BiocGenerics </td>
   <td style="text-align:center;"> 0.51.3 </td>
   <td style="text-align:center;"> 0.51.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/BiocGenerics </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/BiocGenerics </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-10-02 </td>
   <td style="text-align:center;"> Bioconductor 3.20 (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> BiocParallel </td>
   <td style="text-align:center;"> BiocParallel </td>
   <td style="text-align:center;"> 1.36.0 </td>
   <td style="text-align:center;"> 1.36.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/BiocParallel </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/BiocParallel </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-10-24 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Biostrings </td>
   <td style="text-align:center;"> Biostrings </td>
   <td style="text-align:center;"> 2.74.0 </td>
   <td style="text-align:center;"> 2.74.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/Biostrings </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/Biostrings </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-11-08 </td>
   <td style="text-align:center;"> Bioconductor 3.20 (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> bit </td>
   <td style="text-align:center;"> bit </td>
   <td style="text-align:center;"> 4.5.0 </td>
   <td style="text-align:center;"> 4.5.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/bit </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/bit </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-20 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> bit64 </td>
   <td style="text-align:center;"> bit64 </td>
   <td style="text-align:center;"> 4.5.2 </td>
   <td style="text-align:center;"> 4.5.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/bit64 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/bit64 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-22 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> blob </td>
   <td style="text-align:center;"> blob </td>
   <td style="text-align:center;"> 1.2.4 </td>
   <td style="text-align:center;"> 1.2.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/blob </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/blob </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-03-17 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> bslib </td>
   <td style="text-align:center;"> bslib </td>
   <td style="text-align:center;"> 0.8.0 </td>
   <td style="text-align:center;"> 0.8.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/bslib </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/bslib </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-29 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cachem </td>
   <td style="text-align:center;"> cachem </td>
   <td style="text-align:center;"> 1.1.0 </td>
   <td style="text-align:center;"> 1.1.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/cachem </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/cachem </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-05-16 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> circlize </td>
   <td style="text-align:center;"> circlize </td>
   <td style="text-align:center;"> 0.4.16 </td>
   <td style="text-align:center;"> 0.4.16 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/circlize </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/circlize </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-02-20 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cli </td>
   <td style="text-align:center;"> cli </td>
   <td style="text-align:center;"> 3.6.3 </td>
   <td style="text-align:center;"> 3.6.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/cli </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/cli </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-06-21 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cluster </td>
   <td style="text-align:center;"> cluster </td>
   <td style="text-align:center;"> 2.1.6 </td>
   <td style="text-align:center;"> 2.1.6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/cluster </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/cluster </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-12-01 </td>
   <td style="text-align:center;"> CRAN (R 4.4.2) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> clusterProfiler </td>
   <td style="text-align:center;"> clusterProfiler </td>
   <td style="text-align:center;"> 4.14.3 </td>
   <td style="text-align:center;"> 4.14.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/clusterProfiler </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/clusterProfiler </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-11-12 </td>
   <td style="text-align:center;"> Bioconductor 3.20 (R 4.4.2) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> codetools </td>
   <td style="text-align:center;"> codetools </td>
   <td style="text-align:center;"> 0.2.20 </td>
   <td style="text-align:center;"> 0.2-20 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/codetools </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/codetools </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-03-31 </td>
   <td style="text-align:center;"> CRAN (R 4.4.2) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> colorspace </td>
   <td style="text-align:center;"> colorspace </td>
   <td style="text-align:center;"> 2.1.1 </td>
   <td style="text-align:center;"> 2.1-1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/colorspace </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/colorspace </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-26 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cowplot </td>
   <td style="text-align:center;"> cowplot </td>
   <td style="text-align:center;"> 1.1.3 </td>
   <td style="text-align:center;"> 1.1.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/cowplot </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/cowplot </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-01-22 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> crayon </td>
   <td style="text-align:center;"> crayon </td>
   <td style="text-align:center;"> 1.5.3 </td>
   <td style="text-align:center;"> 1.5.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/crayon </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/crayon </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-06-20 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:center;"> data.table </td>
   <td style="text-align:center;"> 1.16.2 </td>
   <td style="text-align:center;"> 1.16.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/data.table </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/data.table </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-10-10 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> DBI </td>
   <td style="text-align:center;"> DBI </td>
   <td style="text-align:center;"> 1.2.3 </td>
   <td style="text-align:center;"> 1.2.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/DBI </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/DBI </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-06-02 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> DelayedArray </td>
   <td style="text-align:center;"> DelayedArray </td>
   <td style="text-align:center;"> 0.32.0 </td>
   <td style="text-align:center;"> 0.32.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/DelayedArray </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/DelayedArray </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-11-08 </td>
   <td style="text-align:center;"> Bioconductor 3.20 (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> deldir </td>
   <td style="text-align:center;"> deldir </td>
   <td style="text-align:center;"> 2.0.4 </td>
   <td style="text-align:center;"> 2.0-4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/deldir </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/deldir </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-02-28 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> devtools </td>
   <td style="text-align:center;"> devtools </td>
   <td style="text-align:center;"> 2.4.5 </td>
   <td style="text-align:center;"> 2.4.5 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/devtools </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/devtools </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-10-11 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> digest </td>
   <td style="text-align:center;"> digest </td>
   <td style="text-align:center;"> 0.6.37 </td>
   <td style="text-align:center;"> 0.6.37 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/digest </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/digest </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-08-19 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> DOSE </td>
   <td style="text-align:center;"> DOSE </td>
   <td style="text-align:center;"> 4.0.0 </td>
   <td style="text-align:center;"> 4.0.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/DOSE </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/DOSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-11-08 </td>
   <td style="text-align:center;"> Bioconductor 3.20 (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> dotCall64 </td>
   <td style="text-align:center;"> dotCall64 </td>
   <td style="text-align:center;"> 1.2 </td>
   <td style="text-align:center;"> 1.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/dotCall64 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/dotCall64 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-10-04 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> dplyr </td>
   <td style="text-align:center;"> dplyr </td>
   <td style="text-align:center;"> 1.1.4 </td>
   <td style="text-align:center;"> 1.1.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/dplyr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/dplyr </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-11-17 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ellipsis </td>
   <td style="text-align:center;"> ellipsis </td>
   <td style="text-align:center;"> 0.3.2 </td>
   <td style="text-align:center;"> 0.3.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ellipsis </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ellipsis </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-04-29 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> enrichplot </td>
   <td style="text-align:center;"> enrichplot </td>
   <td style="text-align:center;"> 1.26.2 </td>
   <td style="text-align:center;"> 1.26.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/enrichplot </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/enrichplot </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-11-08 </td>
   <td style="text-align:center;"> Bioconductor 3.20 (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> evaluate </td>
   <td style="text-align:center;"> evaluate </td>
   <td style="text-align:center;"> 1.0.1 </td>
   <td style="text-align:center;"> 1.0.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/evaluate </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/evaluate </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-10-10 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fansi </td>
   <td style="text-align:center;"> fansi </td>
   <td style="text-align:center;"> 1.0.6 </td>
   <td style="text-align:center;"> 1.0.6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fansi </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fansi </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-12-08 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> farver </td>
   <td style="text-align:center;"> farver </td>
   <td style="text-align:center;"> 2.1.2 </td>
   <td style="text-align:center;"> 2.1.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/farver </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/farver </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-05-13 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fastDummies </td>
   <td style="text-align:center;"> fastDummies </td>
   <td style="text-align:center;"> 1.7.4 </td>
   <td style="text-align:center;"> 1.7.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fastDummies </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fastDummies </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-08-16 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fastmap </td>
   <td style="text-align:center;"> fastmap </td>
   <td style="text-align:center;"> 1.2.0 </td>
   <td style="text-align:center;"> 1.2.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fastmap </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fastmap </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-05-15 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fastmatch </td>
   <td style="text-align:center;"> fastmatch </td>
   <td style="text-align:center;"> 1.1.4 </td>
   <td style="text-align:center;"> 1.1-4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fastmatch </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fastmatch </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-08-18 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fgsea </td>
   <td style="text-align:center;"> fgsea </td>
   <td style="text-align:center;"> 1.32.0 </td>
   <td style="text-align:center;"> 1.32.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fgsea </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fgsea </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-11-08 </td>
   <td style="text-align:center;"> Bioconductor 3.20 (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fitdistrplus </td>
   <td style="text-align:center;"> fitdistrplus </td>
   <td style="text-align:center;"> 1.2.1 </td>
   <td style="text-align:center;"> 1.2-1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fitdistrplus </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fitdistrplus </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-12 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> forcats </td>
   <td style="text-align:center;"> forcats </td>
   <td style="text-align:center;"> 1.0.0 </td>
   <td style="text-align:center;"> 1.0.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/forcats </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/forcats </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-01-29 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fs </td>
   <td style="text-align:center;"> fs </td>
   <td style="text-align:center;"> 1.6.5 </td>
   <td style="text-align:center;"> 1.6.5 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fs </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fs </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-10-30 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> future </td>
   <td style="text-align:center;"> future </td>
   <td style="text-align:center;"> 1.34.0 </td>
   <td style="text-align:center;"> 1.34.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/future </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/future </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-29 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> future.apply </td>
   <td style="text-align:center;"> future.apply </td>
   <td style="text-align:center;"> 1.11.3 </td>
   <td style="text-align:center;"> 1.11.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/future.apply </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/future.apply </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-10-27 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> generics </td>
   <td style="text-align:center;"> generics </td>
   <td style="text-align:center;"> 0.1.3 </td>
   <td style="text-align:center;"> 0.1.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/generics </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/generics </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-07-05 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GenomeInfoDb </td>
   <td style="text-align:center;"> GenomeInfoDb </td>
   <td style="text-align:center;"> 1.42.0 </td>
   <td style="text-align:center;"> 1.42.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/GenomeInfoDb </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/GenomeInfoDb </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-11-08 </td>
   <td style="text-align:center;"> Bioconductor 3.20 (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GenomeInfoDbData </td>
   <td style="text-align:center;"> GenomeInfoDbData </td>
   <td style="text-align:center;"> 1.2.13 </td>
   <td style="text-align:center;"> 1.2.13 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/GenomeInfoDbData </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/GenomeInfoDbData </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-11-10 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GenomicRanges </td>
   <td style="text-align:center;"> GenomicRanges </td>
   <td style="text-align:center;"> 1.58.0 </td>
   <td style="text-align:center;"> 1.58.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/GenomicRanges </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/GenomicRanges </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-11-08 </td>
   <td style="text-align:center;"> Bioconductor 3.20 (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ggbeeswarm </td>
   <td style="text-align:center;"> ggbeeswarm </td>
   <td style="text-align:center;"> 0.7.2 </td>
   <td style="text-align:center;"> 0.7.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggbeeswarm </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggbeeswarm </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-04-29 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ggfun </td>
   <td style="text-align:center;"> ggfun </td>
   <td style="text-align:center;"> 0.1.7 </td>
   <td style="text-align:center;"> 0.1.7 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggfun </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggfun </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-10-24 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ggplot2 </td>
   <td style="text-align:center;"> ggplot2 </td>
   <td style="text-align:center;"> 3.5.1 </td>
   <td style="text-align:center;"> 3.5.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggplot2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggplot2 </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-04-23 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ggplotify </td>
   <td style="text-align:center;"> ggplotify </td>
   <td style="text-align:center;"> 0.1.2 </td>
   <td style="text-align:center;"> 0.1.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggplotify </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggplotify </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-08-09 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ggprism </td>
   <td style="text-align:center;"> ggprism </td>
   <td style="text-align:center;"> 1.0.5 </td>
   <td style="text-align:center;"> 1.0.5 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggprism </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggprism </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-03-21 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ggrastr </td>
   <td style="text-align:center;"> ggrastr </td>
   <td style="text-align:center;"> 1.0.2 </td>
   <td style="text-align:center;"> 1.0.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggrastr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggrastr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-06-01 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ggrepel </td>
   <td style="text-align:center;"> ggrepel </td>
   <td style="text-align:center;"> 0.9.6 </td>
   <td style="text-align:center;"> 0.9.6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggrepel </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggrepel </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-07 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ggridges </td>
   <td style="text-align:center;"> ggridges </td>
   <td style="text-align:center;"> 0.5.6 </td>
   <td style="text-align:center;"> 0.5.6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggridges </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggridges </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-01-23 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ggtangle </td>
   <td style="text-align:center;"> ggtangle </td>
   <td style="text-align:center;"> 0.0.4 </td>
   <td style="text-align:center;"> 0.0.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggtangle </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggtangle </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-10-30 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ggtree </td>
   <td style="text-align:center;"> ggtree </td>
   <td style="text-align:center;"> 3.14.0 </td>
   <td style="text-align:center;"> 3.14.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggtree </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggtree </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-11-08 </td>
   <td style="text-align:center;"> Bioconductor 3.20 (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GlobalOptions </td>
   <td style="text-align:center;"> GlobalOptions </td>
   <td style="text-align:center;"> 0.1.2 </td>
   <td style="text-align:center;"> 0.1.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/GlobalOptions </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/GlobalOptions </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-06-10 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> globals </td>
   <td style="text-align:center;"> globals </td>
   <td style="text-align:center;"> 0.16.3 </td>
   <td style="text-align:center;"> 0.16.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/globals </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/globals </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-03-08 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> glue </td>
   <td style="text-align:center;"> glue </td>
   <td style="text-align:center;"> 1.8.0 </td>
   <td style="text-align:center;"> 1.8.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/glue </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/glue </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-30 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO.db </td>
   <td style="text-align:center;"> GO.db </td>
   <td style="text-align:center;"> 3.20.0 </td>
   <td style="text-align:center;"> 3.20.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/GO.db </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/GO.db </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-11-17 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> goftest </td>
   <td style="text-align:center;"> goftest </td>
   <td style="text-align:center;"> 1.2.3 </td>
   <td style="text-align:center;"> 1.2-3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/goftest </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/goftest </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-10-07 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GOSemSim </td>
   <td style="text-align:center;"> GOSemSim </td>
   <td style="text-align:center;"> 2.32.0 </td>
   <td style="text-align:center;"> 2.32.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/GOSemSim </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/GOSemSim </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-11-08 </td>
   <td style="text-align:center;"> Bioconductor 3.20 (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> gridExtra </td>
   <td style="text-align:center;"> gridExtra </td>
   <td style="text-align:center;"> 2.3 </td>
   <td style="text-align:center;"> 2.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/gridExtra </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/gridExtra </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2017-09-09 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> gridGraphics </td>
   <td style="text-align:center;"> gridGraphics </td>
   <td style="text-align:center;"> 0.5.1 </td>
   <td style="text-align:center;"> 0.5-1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/gridGraphics </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/gridGraphics </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-12-13 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> gson </td>
   <td style="text-align:center;"> gson </td>
   <td style="text-align:center;"> 0.1.0 </td>
   <td style="text-align:center;"> 0.1.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/gson </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/gson </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-03-07 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> gtable </td>
   <td style="text-align:center;"> gtable </td>
   <td style="text-align:center;"> 0.3.6 </td>
   <td style="text-align:center;"> 0.3.6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/gtable </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/gtable </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-10-25 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> htmltools </td>
   <td style="text-align:center;"> htmltools </td>
   <td style="text-align:center;"> 0.5.8.1 </td>
   <td style="text-align:center;"> 0.5.8.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/htmltools </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/htmltools </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-04-04 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> htmlwidgets </td>
   <td style="text-align:center;"> htmlwidgets </td>
   <td style="text-align:center;"> 1.6.4 </td>
   <td style="text-align:center;"> 1.6.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/htmlwidgets </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/htmlwidgets </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-12-06 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> httpuv </td>
   <td style="text-align:center;"> httpuv </td>
   <td style="text-align:center;"> 1.6.15 </td>
   <td style="text-align:center;"> 1.6.15 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/httpuv </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/httpuv </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-03-26 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> httr </td>
   <td style="text-align:center;"> httr </td>
   <td style="text-align:center;"> 1.4.7 </td>
   <td style="text-align:center;"> 1.4.7 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/httr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/httr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-08-15 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ica </td>
   <td style="text-align:center;"> ica </td>
   <td style="text-align:center;"> 1.0.3 </td>
   <td style="text-align:center;"> 1.0-3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ica </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ica </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-07-08 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> igraph </td>
   <td style="text-align:center;"> igraph </td>
   <td style="text-align:center;"> 2.1.1 </td>
   <td style="text-align:center;"> 2.1.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/igraph </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/igraph </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-10-19 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IRanges </td>
   <td style="text-align:center;"> IRanges </td>
   <td style="text-align:center;"> 2.39.2 </td>
   <td style="text-align:center;"> 2.39.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/IRanges </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/IRanges </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-17 </td>
   <td style="text-align:center;"> Bioconductor 3.20 (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> irlba </td>
   <td style="text-align:center;"> irlba </td>
   <td style="text-align:center;"> 2.3.5.1 </td>
   <td style="text-align:center;"> 2.3.5.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/irlba </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/irlba </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-10-03 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> janitor </td>
   <td style="text-align:center;"> janitor </td>
   <td style="text-align:center;"> 2.2.0 </td>
   <td style="text-align:center;"> 2.2.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/janitor </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/janitor </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-02-02 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> jquerylib </td>
   <td style="text-align:center;"> jquerylib </td>
   <td style="text-align:center;"> 0.1.4 </td>
   <td style="text-align:center;"> 0.1.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/jquerylib </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/jquerylib </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-04-26 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> jsonlite </td>
   <td style="text-align:center;"> jsonlite </td>
   <td style="text-align:center;"> 1.8.9 </td>
   <td style="text-align:center;"> 1.8.9 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/jsonlite </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/jsonlite </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-20 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> kableExtra </td>
   <td style="text-align:center;"> kableExtra </td>
   <td style="text-align:center;"> 1.4.0 </td>
   <td style="text-align:center;"> 1.4.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/kableExtra </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/kableExtra </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-01-24 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KEGGREST </td>
   <td style="text-align:center;"> KEGGREST </td>
   <td style="text-align:center;"> 1.46.0 </td>
   <td style="text-align:center;"> 1.46.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/KEGGREST </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/KEGGREST </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-11-08 </td>
   <td style="text-align:center;"> Bioconductor 3.20 (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KernSmooth </td>
   <td style="text-align:center;"> KernSmooth </td>
   <td style="text-align:center;"> 2.23.24 </td>
   <td style="text-align:center;"> 2.23-24 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/KernSmooth </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/KernSmooth </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-05-17 </td>
   <td style="text-align:center;"> CRAN (R 4.4.2) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> knitr </td>
   <td style="text-align:center;"> knitr </td>
   <td style="text-align:center;"> 1.49 </td>
   <td style="text-align:center;"> 1.49 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/knitr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/knitr </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-11-08 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ks </td>
   <td style="text-align:center;"> ks </td>
   <td style="text-align:center;"> 1.14.3 </td>
   <td style="text-align:center;"> 1.14.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ks </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ks </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-20 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> later </td>
   <td style="text-align:center;"> later </td>
   <td style="text-align:center;"> 1.3.2 </td>
   <td style="text-align:center;"> 1.3.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/later </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/later </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-12-06 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lattice </td>
   <td style="text-align:center;"> lattice </td>
   <td style="text-align:center;"> 0.22.6 </td>
   <td style="text-align:center;"> 0.22-6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/lattice </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/lattice </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-03-20 </td>
   <td style="text-align:center;"> CRAN (R 4.4.2) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lazyeval </td>
   <td style="text-align:center;"> lazyeval </td>
   <td style="text-align:center;"> 0.2.2 </td>
   <td style="text-align:center;"> 0.2.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/lazyeval </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/lazyeval </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-03-15 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> leiden </td>
   <td style="text-align:center;"> leiden </td>
   <td style="text-align:center;"> 0.4.3.1 </td>
   <td style="text-align:center;"> 0.4.3.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/leiden </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/leiden </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-11-17 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lifecycle </td>
   <td style="text-align:center;"> lifecycle </td>
   <td style="text-align:center;"> 1.0.4 </td>
   <td style="text-align:center;"> 1.0.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/lifecycle </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/lifecycle </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-11-07 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> listenv </td>
   <td style="text-align:center;"> listenv </td>
   <td style="text-align:center;"> 0.9.1 </td>
   <td style="text-align:center;"> 0.9.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/listenv </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/listenv </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-01-29 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lmtest </td>
   <td style="text-align:center;"> lmtest </td>
   <td style="text-align:center;"> 0.9.40 </td>
   <td style="text-align:center;"> 0.9-40 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/lmtest </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/lmtest </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-03-21 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lubridate </td>
   <td style="text-align:center;"> lubridate </td>
   <td style="text-align:center;"> 1.9.3 </td>
   <td style="text-align:center;"> 1.9.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/lubridate </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/lubridate </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-09-27 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> magrittr </td>
   <td style="text-align:center;"> magrittr </td>
   <td style="text-align:center;"> 2.0.3 </td>
   <td style="text-align:center;"> 2.0.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/magrittr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/magrittr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-03-30 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MASS </td>
   <td style="text-align:center;"> MASS </td>
   <td style="text-align:center;"> 7.3.61 </td>
   <td style="text-align:center;"> 7.3-61 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/MASS </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/MASS </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-06-13 </td>
   <td style="text-align:center;"> CRAN (R 4.4.2) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Matrix </td>
   <td style="text-align:center;"> Matrix </td>
   <td style="text-align:center;"> 1.7.1 </td>
   <td style="text-align:center;"> 1.7-1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/Matrix </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/Matrix </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-10-18 </td>
   <td style="text-align:center;"> CRAN (R 4.4.2) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MatrixGenerics </td>
   <td style="text-align:center;"> MatrixGenerics </td>
   <td style="text-align:center;"> 1.14.0 </td>
   <td style="text-align:center;"> 1.14.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/MatrixGenerics </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/MatrixGenerics </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-10-24 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> matrixStats </td>
   <td style="text-align:center;"> matrixStats </td>
   <td style="text-align:center;"> 1.4.1 </td>
   <td style="text-align:center;"> 1.4.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/matrixStats </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/matrixStats </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-08 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mclust </td>
   <td style="text-align:center;"> mclust </td>
   <td style="text-align:center;"> 6.1.1 </td>
   <td style="text-align:center;"> 6.1.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/mclust </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/mclust </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-04-29 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> memoise </td>
   <td style="text-align:center;"> memoise </td>
   <td style="text-align:center;"> 2.0.1 </td>
   <td style="text-align:center;"> 2.0.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/memoise </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/memoise </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-11-26 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mime </td>
   <td style="text-align:center;"> mime </td>
   <td style="text-align:center;"> 0.12 </td>
   <td style="text-align:center;"> 0.12 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/mime </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/mime </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-09-28 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miniUI </td>
   <td style="text-align:center;"> miniUI </td>
   <td style="text-align:center;"> 0.1.1.1 </td>
   <td style="text-align:center;"> 0.1.1.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/miniUI </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/miniUI </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2018-05-18 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> munsell </td>
   <td style="text-align:center;"> munsell </td>
   <td style="text-align:center;"> 0.5.1 </td>
   <td style="text-align:center;"> 0.5.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/munsell </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/munsell </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-04-01 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mvtnorm </td>
   <td style="text-align:center;"> mvtnorm </td>
   <td style="text-align:center;"> 1.3.2 </td>
   <td style="text-align:center;"> 1.3-2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/mvtnorm </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/mvtnorm </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-11-04 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Nebulosa </td>
   <td style="text-align:center;"> Nebulosa </td>
   <td style="text-align:center;"> 1.16.0 </td>
   <td style="text-align:center;"> 1.16.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/Nebulosa </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/Nebulosa </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-11-08 </td>
   <td style="text-align:center;"> Bioconductor 3.20 (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> nlme </td>
   <td style="text-align:center;"> nlme </td>
   <td style="text-align:center;"> 3.1.166 </td>
   <td style="text-align:center;"> 3.1-166 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/nlme </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/nlme </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-08-14 </td>
   <td style="text-align:center;"> CRAN (R 4.4.2) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> openxlsx </td>
   <td style="text-align:center;"> openxlsx </td>
   <td style="text-align:center;"> 4.2.7.1 </td>
   <td style="text-align:center;"> 4.2.7.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/openxlsx </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/openxlsx </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-20 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> org.Mm.eg.db </td>
   <td style="text-align:center;"> org.Mm.eg.db </td>
   <td style="text-align:center;"> 3.20.0 </td>
   <td style="text-align:center;"> 3.20.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/org.Mm.eg.db </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/org.Mm.eg.db </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-11-17 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> paletteer </td>
   <td style="text-align:center;"> paletteer </td>
   <td style="text-align:center;"> 1.6.0 </td>
   <td style="text-align:center;"> 1.6.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/paletteer </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/paletteer </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-01-21 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> parallelly </td>
   <td style="text-align:center;"> parallelly </td>
   <td style="text-align:center;"> 1.39.0 </td>
   <td style="text-align:center;"> 1.39.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/parallelly </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/parallelly </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-11-07 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> patchwork </td>
   <td style="text-align:center;"> patchwork </td>
   <td style="text-align:center;"> 1.3.0 </td>
   <td style="text-align:center;"> 1.3.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/patchwork </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/patchwork </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-16 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pbapply </td>
   <td style="text-align:center;"> pbapply </td>
   <td style="text-align:center;"> 1.7.2 </td>
   <td style="text-align:center;"> 1.7-2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pbapply </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pbapply </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-06-27 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pillar </td>
   <td style="text-align:center;"> pillar </td>
   <td style="text-align:center;"> 1.9.0 </td>
   <td style="text-align:center;"> 1.9.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pillar </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pillar </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-03-22 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pkgbuild </td>
   <td style="text-align:center;"> pkgbuild </td>
   <td style="text-align:center;"> 1.4.5 </td>
   <td style="text-align:center;"> 1.4.5 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pkgbuild </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pkgbuild </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-10-28 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pkgconfig </td>
   <td style="text-align:center;"> pkgconfig </td>
   <td style="text-align:center;"> 2.0.3 </td>
   <td style="text-align:center;"> 2.0.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pkgconfig </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pkgconfig </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-09-22 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pkgload </td>
   <td style="text-align:center;"> pkgload </td>
   <td style="text-align:center;"> 1.4.0 </td>
   <td style="text-align:center;"> 1.4.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pkgload </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pkgload </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-06-28 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> plotly </td>
   <td style="text-align:center;"> plotly </td>
   <td style="text-align:center;"> 4.10.4 </td>
   <td style="text-align:center;"> 4.10.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/plotly </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/plotly </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-01-13 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> plyr </td>
   <td style="text-align:center;"> plyr </td>
   <td style="text-align:center;"> 1.8.9 </td>
   <td style="text-align:center;"> 1.8.9 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/plyr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/plyr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-10-02 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> png </td>
   <td style="text-align:center;"> png </td>
   <td style="text-align:center;"> 0.1.8 </td>
   <td style="text-align:center;"> 0.1-8 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/png </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/png </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-11-29 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polyclip </td>
   <td style="text-align:center;"> polyclip </td>
   <td style="text-align:center;"> 1.10.7 </td>
   <td style="text-align:center;"> 1.10-7 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/polyclip </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/polyclip </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-23 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pracma </td>
   <td style="text-align:center;"> pracma </td>
   <td style="text-align:center;"> 2.4.4 </td>
   <td style="text-align:center;"> 2.4.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pracma </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pracma </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-11-10 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> profvis </td>
   <td style="text-align:center;"> profvis </td>
   <td style="text-align:center;"> 0.4.0 </td>
   <td style="text-align:center;"> 0.4.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/profvis </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/profvis </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-20 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> progressr </td>
   <td style="text-align:center;"> progressr </td>
   <td style="text-align:center;"> 0.15.0 </td>
   <td style="text-align:center;"> 0.15.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/progressr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/progressr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-10-29 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> promises </td>
   <td style="text-align:center;"> promises </td>
   <td style="text-align:center;"> 1.3.0 </td>
   <td style="text-align:center;"> 1.3.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/promises </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/promises </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-04-05 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> purrr </td>
   <td style="text-align:center;"> purrr </td>
   <td style="text-align:center;"> 1.0.2 </td>
   <td style="text-align:center;"> 1.0.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/purrr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/purrr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-08-10 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> qvalue </td>
   <td style="text-align:center;"> qvalue </td>
   <td style="text-align:center;"> 2.38.0 </td>
   <td style="text-align:center;"> 2.38.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/qvalue </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/qvalue </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-11-08 </td>
   <td style="text-align:center;"> Bioconductor 3.20 (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> R.methodsS3 </td>
   <td style="text-align:center;"> R.methodsS3 </td>
   <td style="text-align:center;"> 1.8.2 </td>
   <td style="text-align:center;"> 1.8.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/R.methodsS3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/R.methodsS3 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-06-13 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> R.oo </td>
   <td style="text-align:center;"> R.oo </td>
   <td style="text-align:center;"> 1.27.0 </td>
   <td style="text-align:center;"> 1.27.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/R.oo </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/R.oo </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-11-01 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> R.utils </td>
   <td style="text-align:center;"> R.utils </td>
   <td style="text-align:center;"> 2.12.3 </td>
   <td style="text-align:center;"> 2.12.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/R.utils </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/R.utils </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-11-18 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> R6 </td>
   <td style="text-align:center;"> R6 </td>
   <td style="text-align:center;"> 2.5.1 </td>
   <td style="text-align:center;"> 2.5.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/R6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/R6 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-08-19 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RANN </td>
   <td style="text-align:center;"> RANN </td>
   <td style="text-align:center;"> 2.6.2 </td>
   <td style="text-align:center;"> 2.6.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/RANN </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/RANN </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-08-25 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RColorBrewer </td>
   <td style="text-align:center;"> RColorBrewer </td>
   <td style="text-align:center;"> 1.1.3 </td>
   <td style="text-align:center;"> 1.1-3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/RColorBrewer </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/RColorBrewer </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-04-03 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Rcpp </td>
   <td style="text-align:center;"> Rcpp </td>
   <td style="text-align:center;"> 1.0.13.1 </td>
   <td style="text-align:center;"> 1.0.13-1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/Rcpp </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/Rcpp </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-11-02 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RcppAnnoy </td>
   <td style="text-align:center;"> RcppAnnoy </td>
   <td style="text-align:center;"> 0.0.22 </td>
   <td style="text-align:center;"> 0.0.22 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/RcppAnnoy </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/RcppAnnoy </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-01-23 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RcppHNSW </td>
   <td style="text-align:center;"> RcppHNSW </td>
   <td style="text-align:center;"> 0.6.0 </td>
   <td style="text-align:center;"> 0.6.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/RcppHNSW </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/RcppHNSW </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-02-04 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rematch2 </td>
   <td style="text-align:center;"> rematch2 </td>
   <td style="text-align:center;"> 2.1.2 </td>
   <td style="text-align:center;"> 2.1.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/rematch2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/rematch2 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-05-01 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> remotes </td>
   <td style="text-align:center;"> remotes </td>
   <td style="text-align:center;"> 2.5.0 </td>
   <td style="text-align:center;"> 2.5.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/remotes </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/remotes </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-03-17 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> reshape2 </td>
   <td style="text-align:center;"> reshape2 </td>
   <td style="text-align:center;"> 1.4.4 </td>
   <td style="text-align:center;"> 1.4.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/reshape2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/reshape2 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-04-09 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> reticulate </td>
   <td style="text-align:center;"> reticulate </td>
   <td style="text-align:center;"> 1.39.0 </td>
   <td style="text-align:center;"> 1.39.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/reticulate </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/reticulate </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-05 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rlang </td>
   <td style="text-align:center;"> rlang </td>
   <td style="text-align:center;"> 1.1.4 </td>
   <td style="text-align:center;"> 1.1.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/rlang </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/rlang </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-06-04 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rmarkdown </td>
   <td style="text-align:center;"> rmarkdown </td>
   <td style="text-align:center;"> 2.29 </td>
   <td style="text-align:center;"> 2.29 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/rmarkdown </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/rmarkdown </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-11-04 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ROCR </td>
   <td style="text-align:center;"> ROCR </td>
   <td style="text-align:center;"> 1.0.11 </td>
   <td style="text-align:center;"> 1.0-11 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ROCR </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ROCR </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-05-02 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RSpectra </td>
   <td style="text-align:center;"> RSpectra </td>
   <td style="text-align:center;"> 0.16.2 </td>
   <td style="text-align:center;"> 0.16-2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/RSpectra </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/RSpectra </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-18 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RSQLite </td>
   <td style="text-align:center;"> RSQLite </td>
   <td style="text-align:center;"> 2.3.7 </td>
   <td style="text-align:center;"> 2.3.7 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/RSQLite </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/RSQLite </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-05-27 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rstudioapi </td>
   <td style="text-align:center;"> rstudioapi </td>
   <td style="text-align:center;"> 0.17.1 </td>
   <td style="text-align:center;"> 0.17.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/rstudioapi </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/rstudioapi </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-10-22 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Rtsne </td>
   <td style="text-align:center;"> Rtsne </td>
   <td style="text-align:center;"> 0.17 </td>
   <td style="text-align:center;"> 0.17 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/Rtsne </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/Rtsne </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-12-07 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S4Arrays </td>
   <td style="text-align:center;"> S4Arrays </td>
   <td style="text-align:center;"> 1.6.0 </td>
   <td style="text-align:center;"> 1.6.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/S4Arrays </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/S4Arrays </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-11-08 </td>
   <td style="text-align:center;"> Bioconductor 3.20 (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S4Vectors </td>
   <td style="text-align:center;"> S4Vectors </td>
   <td style="text-align:center;"> 0.43.2 </td>
   <td style="text-align:center;"> 0.43.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/S4Vectors </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/S4Vectors </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-17 </td>
   <td style="text-align:center;"> Bioconductor 3.20 (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sass </td>
   <td style="text-align:center;"> sass </td>
   <td style="text-align:center;"> 0.4.9 </td>
   <td style="text-align:center;"> 0.4.9 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/sass </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/sass </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-03-15 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> scales </td>
   <td style="text-align:center;"> scales </td>
   <td style="text-align:center;"> 1.3.0 </td>
   <td style="text-align:center;"> 1.3.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/scales </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/scales </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-11-28 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> scattermore </td>
   <td style="text-align:center;"> scattermore </td>
   <td style="text-align:center;"> 1.2 </td>
   <td style="text-align:center;"> 1.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/scattermore </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/scattermore </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-06-12 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> scCustomize </td>
   <td style="text-align:center;"> scCustomize </td>
   <td style="text-align:center;"> 2.1.2 </td>
   <td style="text-align:center;"> 2.1.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/scCustomize </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/scCustomize </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-02-28 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sctransform </td>
   <td style="text-align:center;"> sctransform </td>
   <td style="text-align:center;"> 0.4.1 </td>
   <td style="text-align:center;"> 0.4.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/sctransform </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/sctransform </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-10-19 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sessioninfo </td>
   <td style="text-align:center;"> sessioninfo </td>
   <td style="text-align:center;"> 1.2.2 </td>
   <td style="text-align:center;"> 1.2.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/sessioninfo </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/sessioninfo </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-12-06 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Seurat </td>
   <td style="text-align:center;"> Seurat </td>
   <td style="text-align:center;"> 5.1.0 </td>
   <td style="text-align:center;"> 5.1.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/Seurat </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/Seurat </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-05-10 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SeuratObject </td>
   <td style="text-align:center;"> SeuratObject </td>
   <td style="text-align:center;"> 5.0.2 </td>
   <td style="text-align:center;"> 5.0.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/SeuratObject </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/SeuratObject </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-05-08 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> shape </td>
   <td style="text-align:center;"> shape </td>
   <td style="text-align:center;"> 1.4.6.1 </td>
   <td style="text-align:center;"> 1.4.6.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/shape </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/shape </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-02-23 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> shiny </td>
   <td style="text-align:center;"> shiny </td>
   <td style="text-align:center;"> 1.9.1 </td>
   <td style="text-align:center;"> 1.9.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/shiny </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/shiny </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-08-01 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SingleCellExperiment </td>
   <td style="text-align:center;"> SingleCellExperiment </td>
   <td style="text-align:center;"> 1.24.0 </td>
   <td style="text-align:center;"> 1.24.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/SingleCellExperiment </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/SingleCellExperiment </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-10-24 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> snakecase </td>
   <td style="text-align:center;"> snakecase </td>
   <td style="text-align:center;"> 0.11.1 </td>
   <td style="text-align:center;"> 0.11.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/snakecase </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/snakecase </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-08-27 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sp </td>
   <td style="text-align:center;"> sp </td>
   <td style="text-align:center;"> 2.1.4 </td>
   <td style="text-align:center;"> 2.1-4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/sp </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/sp </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-04-30 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spam </td>
   <td style="text-align:center;"> spam </td>
   <td style="text-align:center;"> 2.11.0 </td>
   <td style="text-align:center;"> 2.11-0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spam </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spam </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-10-03 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SparseArray </td>
   <td style="text-align:center;"> SparseArray </td>
   <td style="text-align:center;"> 1.6.0 </td>
   <td style="text-align:center;"> 1.6.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/SparseArray </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/SparseArray </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-11-08 </td>
   <td style="text-align:center;"> Bioconductor 3.20 (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spatstat.data </td>
   <td style="text-align:center;"> spatstat.data </td>
   <td style="text-align:center;"> 3.1.2 </td>
   <td style="text-align:center;"> 3.1-2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.data </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.data </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-06-21 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spatstat.explore </td>
   <td style="text-align:center;"> spatstat.explore </td>
   <td style="text-align:center;"> 3.3.3 </td>
   <td style="text-align:center;"> 3.3-3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.explore </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.explore </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-10-22 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spatstat.geom </td>
   <td style="text-align:center;"> spatstat.geom </td>
   <td style="text-align:center;"> 3.3.3 </td>
   <td style="text-align:center;"> 3.3-3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.geom </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.geom </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-18 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spatstat.random </td>
   <td style="text-align:center;"> spatstat.random </td>
   <td style="text-align:center;"> 3.3.2 </td>
   <td style="text-align:center;"> 3.3-2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.random </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.random </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-18 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spatstat.sparse </td>
   <td style="text-align:center;"> spatstat.sparse </td>
   <td style="text-align:center;"> 3.1.0 </td>
   <td style="text-align:center;"> 3.1-0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.sparse </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.sparse </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-06-21 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spatstat.univar </td>
   <td style="text-align:center;"> spatstat.univar </td>
   <td style="text-align:center;"> 3.1.1 </td>
   <td style="text-align:center;"> 3.1-1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.univar </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.univar </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-11-05 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spatstat.utils </td>
   <td style="text-align:center;"> spatstat.utils </td>
   <td style="text-align:center;"> 3.1.1 </td>
   <td style="text-align:center;"> 3.1-1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.utils </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.utils </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-11-03 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> stringi </td>
   <td style="text-align:center;"> stringi </td>
   <td style="text-align:center;"> 1.8.4 </td>
   <td style="text-align:center;"> 1.8.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/stringi </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/stringi </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-05-06 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> stringr </td>
   <td style="text-align:center;"> stringr </td>
   <td style="text-align:center;"> 1.5.1 </td>
   <td style="text-align:center;"> 1.5.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/stringr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/stringr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-11-14 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SummarizedExperiment </td>
   <td style="text-align:center;"> SummarizedExperiment </td>
   <td style="text-align:center;"> 1.32.0 </td>
   <td style="text-align:center;"> 1.32.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/SummarizedExperiment </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/SummarizedExperiment </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-10-24 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> survival </td>
   <td style="text-align:center;"> survival </td>
   <td style="text-align:center;"> 3.7.0 </td>
   <td style="text-align:center;"> 3.7-0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/survival </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/survival </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-06-05 </td>
   <td style="text-align:center;"> CRAN (R 4.4.2) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> svglite </td>
   <td style="text-align:center;"> svglite </td>
   <td style="text-align:center;"> 2.1.3 </td>
   <td style="text-align:center;"> 2.1.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/svglite </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/svglite </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-12-08 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> systemfonts </td>
   <td style="text-align:center;"> systemfonts </td>
   <td style="text-align:center;"> 1.1.0 </td>
   <td style="text-align:center;"> 1.1.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/systemfonts </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/systemfonts </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-05-15 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> tensor </td>
   <td style="text-align:center;"> tensor </td>
   <td style="text-align:center;"> 1.5 </td>
   <td style="text-align:center;"> 1.5 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/tensor </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/tensor </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2012-05-05 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> tibble </td>
   <td style="text-align:center;"> tibble </td>
   <td style="text-align:center;"> 3.2.1 </td>
   <td style="text-align:center;"> 3.2.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/tibble </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/tibble </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-03-20 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> tidyr </td>
   <td style="text-align:center;"> tidyr </td>
   <td style="text-align:center;"> 1.3.1 </td>
   <td style="text-align:center;"> 1.3.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/tidyr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/tidyr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-01-24 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> tidyselect </td>
   <td style="text-align:center;"> tidyselect </td>
   <td style="text-align:center;"> 1.2.1 </td>
   <td style="text-align:center;"> 1.2.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/tidyselect </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/tidyselect </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-03-11 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> tidytree </td>
   <td style="text-align:center;"> tidytree </td>
   <td style="text-align:center;"> 0.4.6 </td>
   <td style="text-align:center;"> 0.4.6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/tidytree </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/tidytree </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-12-12 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> timechange </td>
   <td style="text-align:center;"> timechange </td>
   <td style="text-align:center;"> 0.3.0 </td>
   <td style="text-align:center;"> 0.3.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/timechange </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/timechange </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-01-18 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> treeio </td>
   <td style="text-align:center;"> treeio </td>
   <td style="text-align:center;"> 1.30.0 </td>
   <td style="text-align:center;"> 1.30.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/treeio </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/treeio </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-11-08 </td>
   <td style="text-align:center;"> Bioconductor 3.20 (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UCSC.utils </td>
   <td style="text-align:center;"> UCSC.utils </td>
   <td style="text-align:center;"> 1.2.0 </td>
   <td style="text-align:center;"> 1.2.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/UCSC.utils </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/UCSC.utils </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-11-08 </td>
   <td style="text-align:center;"> Bioconductor 3.20 (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> urlchecker </td>
   <td style="text-align:center;"> urlchecker </td>
   <td style="text-align:center;"> 1.0.1 </td>
   <td style="text-align:center;"> 1.0.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/urlchecker </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/urlchecker </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-11-30 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> usethis </td>
   <td style="text-align:center;"> usethis </td>
   <td style="text-align:center;"> 3.0.0 </td>
   <td style="text-align:center;"> 3.0.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/usethis </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/usethis </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-29 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> utf8 </td>
   <td style="text-align:center;"> utf8 </td>
   <td style="text-align:center;"> 1.2.4 </td>
   <td style="text-align:center;"> 1.2.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/utf8 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/utf8 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-10-22 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> uwot </td>
   <td style="text-align:center;"> uwot </td>
   <td style="text-align:center;"> 0.2.2 </td>
   <td style="text-align:center;"> 0.2.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/uwot </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/uwot </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-04-21 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> vctrs </td>
   <td style="text-align:center;"> vctrs </td>
   <td style="text-align:center;"> 0.6.5 </td>
   <td style="text-align:center;"> 0.6.5 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/vctrs </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/vctrs </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-12-01 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> vipor </td>
   <td style="text-align:center;"> vipor </td>
   <td style="text-align:center;"> 0.4.7 </td>
   <td style="text-align:center;"> 0.4.7 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/vipor </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/vipor </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-12-18 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> viridisLite </td>
   <td style="text-align:center;"> viridisLite </td>
   <td style="text-align:center;"> 0.4.2 </td>
   <td style="text-align:center;"> 0.4.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/viridisLite </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/viridisLite </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-05-02 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> withr </td>
   <td style="text-align:center;"> withr </td>
   <td style="text-align:center;"> 3.0.2 </td>
   <td style="text-align:center;"> 3.0.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/withr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/withr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-10-28 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> xfun </td>
   <td style="text-align:center;"> xfun </td>
   <td style="text-align:center;"> 0.49 </td>
   <td style="text-align:center;"> 0.49 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/xfun </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/xfun </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-10-31 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> xml2 </td>
   <td style="text-align:center;"> xml2 </td>
   <td style="text-align:center;"> 1.3.6 </td>
   <td style="text-align:center;"> 1.3.6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/xml2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/xml2 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-12-04 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> xtable </td>
   <td style="text-align:center;"> xtable </td>
   <td style="text-align:center;"> 1.8.4 </td>
   <td style="text-align:center;"> 1.8-4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/xtable </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/xtable </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-04-21 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> XVector </td>
   <td style="text-align:center;"> XVector </td>
   <td style="text-align:center;"> 0.46.0 </td>
   <td style="text-align:center;"> 0.46.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/XVector </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/XVector </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-11-08 </td>
   <td style="text-align:center;"> Bioconductor 3.20 (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> yaml </td>
   <td style="text-align:center;"> yaml </td>
   <td style="text-align:center;"> 2.3.10 </td>
   <td style="text-align:center;"> 2.3.10 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/yaml </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/yaml </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-26 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> yulab.utils </td>
   <td style="text-align:center;"> yulab.utils </td>
   <td style="text-align:center;"> 0.1.8 </td>
   <td style="text-align:center;"> 0.1.8 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/yulab.utils </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/yulab.utils </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-11-07 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> zip </td>
   <td style="text-align:center;"> zip </td>
   <td style="text-align:center;"> 2.3.1 </td>
   <td style="text-align:center;"> 2.3.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/zip </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/zip </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-01-27 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> zlibbioc </td>
   <td style="text-align:center;"> zlibbioc </td>
   <td style="text-align:center;"> 1.52.0 </td>
   <td style="text-align:center;"> 1.52.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/zlibbioc </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/zlibbioc </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-11-08 </td>
   <td style="text-align:center;"> Bioconductor 3.20 (R 4.4.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> zoo </td>
   <td style="text-align:center;"> zoo </td>
   <td style="text-align:center;"> 1.8.12 </td>
   <td style="text-align:center;"> 1.8-12 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/zoo </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/zoo </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-04-13 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
</tbody>
</table>
