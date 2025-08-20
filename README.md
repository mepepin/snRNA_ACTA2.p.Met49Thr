---
title: "Single Nuclear Resolution of ACTA2 Mutation in Thoracic Ascending Aortic Tissue"
author: "Mark E. Pepin, MD, PhD, MS"
date: "07/16/2023"
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



**Code Author**: Mark E. Pepin, MD, PhD, MS 
**Contact**: [pepinme\@gmail.com](mailto:mepepin@bwh.harvard.edu){.email}\
**Institution**: Brigham and Women's Hospital | Harvard Medical School | Broad Institute of Harvard and MIT \
**Location**: Boston, MA, USA

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

Standard methods were used to cluster nuclei based on Uniform Manifold Approximation and Projection (UMAP), a non-linear dimensionality reduction technique that projects high-dimensional data into lower dimensions while preserving its local structure and relationships. It works by constructing a weighted graph representing the data’s local neighborhood and then optimizing its layout in a low-dimensional space to capture both local and global data patterns.


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
# # Save this instance to avoid needing to re-run:
# saveRDS(acta2_scaled, file = "ACTA2_clustered.rds")
```

## Annotated UMAP, Density Plots

Annotation of clusters was accomplished by comparing the expression patters of canonical cell type marker genes, shown below.


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
VSMC.1.markers <- FindMarkers(acta2_scaled, ident.1 = "VSMC_1", min.pct = 0.25)
VSMC.2.markers <- FindMarkers(acta2_scaled, ident.1 = "VSMC_2", min.pct = 0.25)
VSMC.3.markers <- FindMarkers(acta2_scaled, ident.1 = "VSMC_3", min.pct = 0.25)
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
  writeData(wb_DESeq, "VSMC1", VSMC.1.markers, startCo = 1, rowNames = T)
      addWorksheet(wb_DESeq, "VSMC2")
  writeData(wb_DESeq, "VSMC2", VSMC.2.markers, startCo = 1, rowNames = T)
      addWorksheet(wb_DESeq, "VSMC3")
  writeData(wb_DESeq, "VSMC3", VSMC.3.markers, startCo = 1, rowNames = T)
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
# Save file for downstream analysis
# saveRDS(acta2_scaled, file = "../1_Input/ACTA2_annotated.rds")
```


## VSMC - Proliferation scores
To explore whether differences in VSMC subpopulations reflects differential regulation of proliferative pathways, as previously reported by Kwartler et al. 2023, we used the Gene Ontology term GO0008283.


``` r
library(dplyr)
library(Seurat)
library(ggplot2)
# Load the necessary libraries
acta2_scaled<-readRDS(file = "../1_Input/ACTA2_annotated.rds")
vsmc_idents <- grep("VSMC", levels(acta2_scaled), value = TRUE)
acta2_scaled <- subset(acta2_scaled, idents = vsmc_idents)
acta2_scaled <- RunPCA(acta2_scaled) # re-cluster
acta2_scaled <- RunUMAP(acta2_scaled, assay = "RNA", dims = 1:40) # re-cluster
Idents(acta2_scaled) <- acta2_scaled$RNA_snn_res.0.3
library(clusterProfiler)
library(org.Hs.eg.db)  #gene database
# Define the Human GO term of interest, e.g., GO:0008150 (biological process)
go_term <- "GO:0008283" # used in the original 2023 paper
# Extract gene symbols associated with the specified GO term
genes_vector <- AnnotationDbi::select(
  org.Hs.eg.db,                 # gene database
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

# Contractility
go_term <- "GO:0006940" # used in the original 2023 paper
genes_vector <- AnnotationDbi::select(
  org.Hs.eg.db,                 # Human gene database
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
ann_colors = list(CellType = c(VSMC="coral2", 
                                 Fibroblast="wheat", 
                                 EC1="steelblue4",
                                 EC2="deepskyblue3",
                                 Macrophage="azure4",
                                 NKT="gray",
                                 Dendritic="tan2"))
acta2_scaled<-readRDS(file = "../1_Input/ACTA2_annotated.rds")
# Subsetting and re-normalization
acta2_scaled <- subset(acta2_scaled, idents = grep("VSMC", Idents(acta2_scaled), value = TRUE))
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
# saveRDS(TAAng_scaled, file = "../1_Input/TAA_nongenomic_clustered.rds")
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

# Contractility score
contractility_genes <- c("CNN1", "TAGLN", "ACTA2", "MYH11")
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
# saveRDS(TAAng_scaled, file = "../1_Input/TAAng_annotated.rds")
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
# cds <- order_cells(cds) # Used when setting the nodes; if already known, use the next line
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
Aortic_genes <- enrich_3 %>% filter(str_detect(Term,"Aort"))
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

Data from Chou et al. 2022 were incorporated to provide a more generalizeable understanding of cellular heterogeneity at single-nuclear resolution. Once the individual snRNA-seq datasets were interrogated for quality and normalized, they were fully integrated using standard protocol within Seurat. To accomplish this, shared cell populations were matched across datasets ("anchors") were used to correct for technical differences between datasets (i.e. batch effects), as well as to perform comparative sn-RNA-seq analysis across experimental conditions.


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
# saveRDS(ifnb, "../1_Input/FFPE_integrated.rds")
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
# saveRDS(TAA_integrated, file = "../1_Input/TAA_snRNA.rds")
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
  object = TAA.combined,
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
  dplyr::select(-count, -label.count.total.per.facet) %>% 
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


# Total number of cells
# Extract metadata
meta_df <- TAA.combined@meta.data
# Tabulate absolute counts of nuclei per cell type split by Disease
count_df <- meta_df %>%
  count(Disease, cell_type, name = "Nuclei_Count")
# Plot as barplot with absolute counts
ggplot(count_df, aes(x = Disease, y = Nuclei_Count, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  ggtitle(NULL) +
  scale_fill_manual(values = c(
    VSMC_1 = "coral2",
    VSMC_2 = "coral3",
    Fibroblast = "steelblue4",
    EC = "wheat",
    Macrophage = "azure4",
    NKT = "gray",
    Dendritic = "tan2",
    Fibromyocyte = "cadetblue3",
    `B Lymph` = "darkgray",
    `T Lymph` = "black"
  )) +
  ylab("Number of Nuclei") +
  xlab("Disease Condition") +
  theme_minimal()
```

## Trajectory - VSMC-specific Trajectory Analysis (Monocle 3)


``` r
library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(ggtrace)
library(ggplot2)
ann_colors = list(CellType = c(VSMC_1="coral2", 
                               VSMC_2="coral3",
                             Fibroblast="wheat", 
                             EC="steelblue4",
                             Macrophage="azure4",
                             NKT="gray",
                             Dendritic="tan2",
                             Fibromyocyte="darkcyan",
                             `B Lymph`="darkgray",
                             `T Lymph`="black"))
combined_scaled<-readRDS(file = "../1_Input/ACTA2_annotated.rds")
# Subsetting and re-normalization
combined_scaled <- subset(combined_scaled, idents = c("VSMC_1", "VSMC_2")) # select only the VSMCs
combined_scaled <- RunPCA(combined_scaled) # re-cluster
combined_scaled <- RunUMAP(combined_scaled, assay = "RNA", dims = 1:40) # re-cluster
# Idents(combined_scaled) <- combined_scaled$RNA_snn_res.0.4
# Clustering
pdf(file = "../2_Output/3_Merged/4_Trajectory/UMAP_VSMC.pdf", height = 4, width = 4)
DimPlot(combined_scaled, group.by = "RNA_snn_res.0.4", reduction = "umap", label = T, repel = T,label.size = 4, cols = c("coral2", "coral4", "coral3"))+NoLegend()+ xlim(-4,6) + ylim(-8,4)
dev.off()
#Contractility
library(ggplot2)
pdf(file = "../2_Output/3_Merged/4_Trajectory/Contractility.Score_UMAP.pdf", height = 4, width = 4)
FeaturePlot(combined_scaled, features = "contractility", reduction = "umap") + ggtitle("Contractility Score") +
scale_color_gradient2(low = "darkcyan", mid = "white", high = "darkred", midpoint = median(combined_scaled$contractility, na.rm = FALSE))+ xlim(-4,6) + ylim(-8,4)
dev.off()
#Osteogenic
pdf(file = "../2_Output/3_Merged/4_Trajectory/Osteogenicity.Score_UMAP.pdf", height = 4, width = 4)
FeaturePlot(combined_scaled, features = "osteogenic", reduction = "umap") + ggtitle("Osteogenicity Score") +
scale_colour_gradient2(low = "darkcyan", mid = "white", high = "darkred", midpoint = median(combined_scaled$osteogenic, na.rm = FALSE))+ xlim(-4,6) + ylim(-8,4)
dev.off()
#Synthetic
pdf(file = "../2_Output/3_Merged/4_Trajectory/Synthetic.Score_UMAP.pdf", height = 4, width = 4)
FeaturePlot(combined_scaled, features = "Synthetic", reduction = "umap", alpha = 0.6) + ggtitle("Synthetic Score") +
scale_colour_gradient2(low = "darkcyan", mid = "white", high = "darkred", midpoint = mean(combined_scaled$Synthetic, na.rm = FALSE))+ xlim(-4,6) + ylim(-8,4)
dev.off()
# 
DefaultAssay(combined_scaled) = "RNA" # This changed in seurat V5
# combined_scaled <- ScaleData(object = combined_scaled)
cds <- as.cell_data_set(combined_scaled)
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
root_cell <- "CATCTTGAGACTGAGT-1_2" # names(which(cds@principal_graph_aux$UMAP$pseudotime==0))
cds<-order_cells(cds, root_cells = root_cell)
# Plot the pseudotime on UMAP
plot_cells(cds, 
           color_cells_by = "pseudotime", 
           label_branch_points = FALSE, 
           label_leaves = FALSE) + xlim(-4,6) + ylim(-8,4)
pdf(file = "../2_Output/3_Merged/4_Trajectory/combined_UMAP_Trajectory_Partition.pdf", height = 3, width = 3)
plot_cells(cds, 
           color_cells_by = "partition", 
           label_branch_points = FALSE, 
           label_leaves = F,
           show_trajectory_graph = F,
           label_roots = F)+ xlim(-4,6) + ylim(-8,4)
dev.off()
pdf(file = "../2_Output/3_Merged/4_Trajectory/combined_UMAP_Trajectory_Pseudotime.pdf", height = 3, width = 4)
plot_cells(cds, 
           color_cells_by = "pseudotime", 
           label_branch_points = FALSE, 
           label_leaves = F,
           show_trajectory_graph = F,
           label_roots = F)+ xlim(-4,6) + ylim(-8,4)
dev.off()
# Examine specific genes
plot_cells(cds, 
           genes = c("ACTA2"),
           label_cell_groups = FALSE,
           show_trajectory_graph = FALSE,
           min_expr = 3)+ xlim(-4,6) + ylim(-8,4)
# Identify pseudotime
modulated_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=8) # Identify differentially-expressed genes with pseudotime
modulated_genes <- na.omit(modulated_genes) # remove NA's
modulated_genes <- modulated_genes %>% filter(modulated_genes$q_value < 0.05 & modulated_genes$status =="OK") # filter cds results down
modulated_genes <- modulated_genes[order(-modulated_genes$morans_test_statistic), ] # order by moran's test
modulated_genes <- top_n(modulated_genes, 1000, -q_value)
#### Create a heatmap of genes with similar pseudotime kinetics
genes <- row.names(subset(modulated_genes, q_value < 0.05))
openxlsx::write.xlsx(modulated_genes, "../2_Output/1_combined/Trajectory/Pseudotime_DEGs.xlsx")
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
Index$cell_type <- factor(Index$cell_type,
                         levels = c("VSMC_1",
                                    "VSMC_2",
                                    "Fibromyocyte",
                                    "Fibroblast", 
                                    "EC",
                                    "Macrophage",
                                    "NKT",
                                    "Dendritic",
                                    "B Lymph",
                                    "T Lymph"))
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
      # fontsize_col = 8,
      # color = myColor,
      # annotation_col = Index,
      # annotation_colors = ann_colors,
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
saveWorkbook(wb_DESeq, file = paste0("../2_Output/3_Merged/4_Trajectory/Trajectory_Pathway.xlsx"), overwrite = TRUE)
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
saveWorkbook(wb_DESeq, file = paste0("../2_Output/3_Merged/4_Trajectory/Trajectory_GWAS.Traits.xlsx"), overwrite = TRUE)
library(stringr)
Aortic_genes <- enrich_4 %>% filter(str_detect(Term,"Aort"))
Aortic_genes <- Aortic_genes$Genes
Aortic_genes <- unique(unlist(strsplit(Aortic_genes, ";", fixed = T)))
# proliferative score
gene_list <-Aortic_genes
gene_list <- gene_list[gene_list %in% rownames(combined_scaled)]
gene_expr <- GetAssayData(object = combined_scaled, assay = "RNA", slot = "data")[gene_list, ]
sum_expr <- colSums(gene_expr)
combined_scaled$Aortic <- scale(sum_expr)
Cluster1_Names<-unlist(strsplit(enrich_1$Genes[1], ";", fixed = T))
Cluster2_Names<-unlist(strsplit(enrich_2$Genes[1], ";", fixed = T))
Cluster3_Names<-unlist(strsplit(enrich_3$Genes[1], ";", fixed = T))
Cluster4_Names<-unlist(strsplit(enrich_4$Genes[1], ";", fixed = T))
pdf(file = "../2_Output/3_Merged/4_Trajectory/Aortic_genes.pdf")
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
############################################################################
GENES_HM<-c(Aortic_genes)
gene_list <- Aortic_genes[Aortic_genes %in% rownames(combined_scaled)]   # combine genes associated with aortic disease from the GWAS Map
gene_expr <- GetAssayData(object = combined_scaled, slot = "data.1")[gene_list, ]
sum_expr <- colSums(gene_expr)
combined_scaled$aortopathy <- scale(sum_expr)
pdf(file = "../2_Output/3_Merged/4_Trajectory/Aortopathy.Score_UMAP.pdf", height = 4, width = 4)
FeaturePlot(combined_scaled, features = "aortopathy", reduction = "umap") + ggtitle("Aortopathy Score") +
scale_color_gradient2(low = "darkcyan", mid = "white", high = "darkred", midpoint = median(combined_scaled$aortopathy, na.rm = FALSE))+ xlim(-4,6) + ylim(-8,4)
dev.off()

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
pdf(file = "../2_Output/3_Merged/4_Trajectory/Trajectory_Heatmap.pdf", height = 8, width = 5)
heatmap_combination
dev.off()

library(Seurat)
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
  pdf(paste0("../2_Output/3_Merged/3_Differential.Expression/Aortic_Genes.pdf"), height = 6, width = 10)
  print(wrap_plots(plots = plots_up, ncol = 3))
  dev.off()
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
  de_markers <- FindMarkers(subset_seurat, ident.1 = "TAA_ACTA2", ident.2 = "TAA_NG", min.pct = 0.25, test.use = "MAST")
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


# Trajectory Analysis

``` r
# --- QC SUPPLEMENT (ACTA2 + MERGE): populate pctRibo & DoubletScore ---

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(openxlsx)
  library(Matrix)
})

message("WD: ", getwd())
dir.create("../2_Output/QC", showWarnings = FALSE, recursive = TRUE)

fp_acta2  <- "../1_Input/ACTA2_annotated.rds"
fp_merged <- "../1_Input/TAA_snRNA.rds"
fp_dv200  <- "../1_Input/library_QC_DV200.csv"  # set to NULL to skip

# ---------- helpers ----------
safe_read_seurat <- function(fp, label){
  if (is.null(fp) || !nzchar(fp)) stop(sprintf("[%s] path is NULL/empty", label))
  if (!file.exists(fp)) stop(sprintf("[%s] file not found at '%s'", label, fp))
  obj <- tryCatch(readRDS(fp), error = function(e)
    stop(sprintf("[%s] readRDS failed: %s", label, conditionMessage(e))))
  if (!inherits(obj, "Seurat")) stop(sprintf("[%s] not a Seurat object (class: %s)", label, paste(class(obj), collapse=", ")))
  # Ensure RNA assay exists and is joined (Seurat v5 layers)
  if (!("RNA" %in% names(obj@assays))) stop(sprintf("[%s] RNA assay not present", label))
  if (!is.null(obj[["RNA"]]@layers) && length(obj[["RNA"]]@layers) > 0) {
    obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
  }
  DefaultAssay(obj) <- "RNA"
  obj
}

pct_feature <- function(seu, pattern, slot = "counts", ignore.case = TRUE){
  rn <- rownames(seu)
  feats <- rn[grepl(pattern, rn, ignore.case = ignore.case)]
  if (length(feats) == 0) return(list(values = rep(NA_real_, ncol(seu)), n_matches = 0L))
  m   <- GetAssayData(seu, assay = DefaultAssay(seu), slot = slot)
  tot <- Matrix::colSums(m)
  sub <- m[feats, , drop = FALSE]
  out <- rep(NA_real_, length(tot))
  ok  <- tot > 0
  out[ok] <- 100 * (Matrix::colSums(sub)[ok] / tot[ok])
  list(values = out, n_matches = length(feats))
}

qc_annotate <- function(seu, mito_pattern="^MT-", ribo_pattern="^RP[SL]"){
  mdn <- colnames(seu@meta.data)
  # mito %
  if (!"percent.mt" %in% mdn) {
    mt <- pct_feature(seu, mito_pattern)
    seu$percent.mt <- mt$values
    message("Added percent.mt (", mt$n_matches, " MT* genes matched)")
  }
  # ribo % (RPL/RPS)
  if (!"percent.ribo" %in% mdn) {
    rb <- pct_feature(seu, ribo_pattern)
    seu$percent.ribo <- rb$values
    message("Added percent.ribo (", rb$n_matches, " RP[SL]* genes matched)")
  }
  seu
}

qc_cell_table <- function(seu, dataset_label){
  md <- seu@meta.data
  tibble::tibble(
    Barcode    = colnames(seu),
    Dataset    = dataset_label,
    Sample     = if (!is.null(md$orig.ident)) as.character(md$orig.ident) else dataset_label,
    Disease    = if (!is.null(md$Disease))    as.character(md$Disease)    else NA_character_,
    CellType   = if (!is.null(md$cell_type))  as.character(md$cell_type)
                 else if (!is.null(Seurat::Idents(seu))) as.character(Seurat::Idents(seu)) else NA_character_,
    nUMI       = as.numeric(md$nCount_RNA),
    nGene      = as.numeric(md$nFeature_RNA),
    pctMT      = as.numeric(md$percent.mt),
    pctRibo    = as.numeric(md$percent.ribo)
  )
}

qc_sample_summary <- function(cell_tbl){
  cell_tbl %>%
    group_by(Dataset, Sample, Disease) %>%
    summarise(
      Nuclei         = n(),
      nUMI_median    = median(nUMI,  na.rm = TRUE),
      nUMI_IQR       = IQR(nUMI,     na.rm = TRUE),
      nGene_median   = median(nGene, na.rm = TRUE),
      nGene_IQR      = IQR(nGene,    na.rm = TRUE),
      pctMT_median   = median(pctMT, na.rm = TRUE),
      pctMT_IQR      = IQR(pctMT,    na.rm = TRUE),
      pctRibo_median = median(pctRibo, na.rm = TRUE),
      pctRibo_IQR    = IQR(pctRibo,    na.rm = TRUE),
      .groups = "drop"
    ) %>% arrange(Dataset, Sample, Disease)
}

qc_celltype_summary <- function(cell_tbl){
  cell_tbl %>%
    group_by(Dataset, Sample, Disease, CellType) %>%
    summarise(
      Nuclei       = n(),
      nUMI_median  = median(nUMI,  na.rm = TRUE),
      nGene_median = median(nGene, na.rm = TRUE),
      pctMT_median = median(pctMT, na.rm = TRUE),
      pctRibo_median = median(pctRibo, na.rm = TRUE),
      .groups = "drop"
    ) %>% arrange(Dataset, CellType, Sample, Disease)
}

# Doublets: scDblFinder if available; otherwise a proxy based on nUMI+nGene
add_doublets <- function(seu){
  if (requireNamespace("scDblFinder", quietly = TRUE) &&
      requireNamespace("SingleCellExperiment", quietly = TRUE) &&
      requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    m <- as.matrix(GetAssayData(seu, slot = "counts"))
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = m))
    sce <- scDblFinder::scDblFinder(sce, verbose = FALSE)
    seu$DoubletScore <- as.numeric(SummarizedExperiment::colData(sce)$scDblFinder.score)
    seu$DoubletCall  <- as.character(SummarizedExperiment::colData(sce)$scDblFinder.class)
    message("Doublets: scDblFinder populated (", sum(seu$DoubletCall == "doublet", na.rm = TRUE), " called)")
  } else {
    # proxy: rank-normalized sum of nUMI + nGene per sample
    md <- seu@meta.data
    nn <- scales::rescale(rank(md$nCount_RNA, ties.method = "average")) +
          scales::rescale(rank(md$nFeature_RNA, ties.method = "average"))
    # z-score the combined metric
    z  <- as.numeric(scale(nn))
    seu$DoubletScore <- z
    # call top 1.5% as doublets (tweakable)
    thr <- quantile(z, 0.985, na.rm = TRUE)
    seu$DoubletCall  <- ifelse(z >= thr, "doublet_proxy", "singlet")
    message("Doublets: proxy populated (threshold z >= ", round(thr, 2), "; ",
            sum(seu$DoubletCall != "singlet", na.rm = TRUE), " flagged)")
  }
  seu
}

# ---------- load & annotate ----------
ACTA2 <- safe_read_seurat(fp_acta2,  "ACTA2")  %>% qc_annotate("^MT-", "^RP[SL]") %>% add_doublets()
MERGE <- safe_read_seurat(fp_merged, "MERGE")  %>% qc_annotate("^MT-", "^RP[SL]") %>% add_doublets()

# ---------- tables ----------
tab_ACTA2 <- qc_cell_table(ACTA2, "TAA_ACTA2") %>%
  mutate(DoubletScore = ACTA2$DoubletScore, DoubletCall = ACTA2$DoubletCall)
tab_MERGE <- qc_cell_table(MERGE, "Merged") %>%
  mutate(DoubletScore = MERGE$DoubletScore, DoubletCall = MERGE$DoubletCall)

tab_all <- bind_rows(tab_ACTA2, tab_MERGE)

# ---------- summaries ----------
sum_sample    <- qc_sample_summary(tab_all)
sum_celltype  <- qc_celltype_summary(tab_all)

sum_disease_merge <- tab_MERGE %>%
  group_by(Disease) %>%
  summarise(
    Nuclei       = n(),
    nUMI_median  = median(nUMI,  na.rm = TRUE),
    nGene_median = median(nGene, na.rm = TRUE),
    pctMT_median = median(pctMT, na.rm = TRUE),
    pctRibo_median = median(pctRibo, na.rm = TRUE),
    .groups = "drop"
  )

# ---------- DV200 (optional) ----------
dv200_tbl <- NULL
if (!is.null(fp_dv200) && file.exists(fp_dv200)){
  dv200_raw <- read.csv(fp_dv200, stringsAsFactors = FALSE)
  cn <- tolower(colnames(dv200_raw))
  s_col <- which.max(cn %in% c("sample","sampleid","library","orig.ident"))
  d_col <- which.max(grepl("dv200", cn))
  if (s_col != 0 && d_col != 0) {
    dv200_tbl <- dv200_raw[, c(s_col, d_col)]
    colnames(dv200_tbl) <- c("Sample","DV200_percent")
    dv200_tbl$DV200_percent <- as.numeric(dv200_tbl$DV200_percent)
    sum_sample <- sum_sample %>% left_join(dv200_tbl, by = "Sample")
  } else {
    warning("DV200 file found but could not detect 'Sample' and 'DV200' columns; skipping merge")
  }
}

# ---------- plots ----------
plot_violin <- function(seu, features, title, out_pdf){
  p <- VlnPlot(seu, features = features, pt.size = 0, combine = TRUE) + ggtitle(title)
  ggsave(out_pdf, p, width = 8, height = 3)
}
plot_violin(ACTA2, c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo"), "ACTA2 QC",  "../2_Output/QC/ACTA2_QC_violin.pdf")
plot_violin(MERGE, c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo"), "Merged QC", "../2_Output/QC/MERGE_QC_violin.pdf")

# Doublet fraction by Sample
p_db <- tab_all %>%
  group_by(Dataset, Sample) %>%
  summarise(DoubletRate = mean(grepl("^doublet", DoubletCall), na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = Sample, y = DoubletRate, fill = Dataset)) +
  geom_col() + coord_flip() + ylab("Doublet fraction") + xlab(NULL) + theme_minimal()
ggsave("../2_Output/QC/Doublet_rates_by_sample.pdf", p_db, width = 6, height = 4)

# disease-split violins for MERGE
p_merge_split <- VlnPlot(
  MERGE, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo"),
  split.by = "Disease", pt.size = 0, combine = TRUE
) + ggtitle("MERGE QC by Disease")
ggsave("../2_Output/QC/MERGE_QC_violin_byDisease.pdf", p_merge_split, width = 10, height = 3.5)

# ---------- export Excel ----------
wb <- createWorkbook()
addWorksheet(wb, "Cell-level_QC");       writeData(wb, "Cell-level_QC", tab_all)
addWorksheet(wb, "Sample_Summary");      writeData(wb, "Sample_Summary", sum_sample)
addWorksheet(wb, "CellType_Summary");    writeData(wb, "CellType_Summary", sum_celltype)
addWorksheet(wb, "MERGE_by_Disease");    writeData(wb, "MERGE_by_Disease", sum_disease_merge)
if (!is.null(dv200_tbl)) { addWorksheet(wb, "DV200"); writeData(wb, "DV200", dv200_tbl) }
saveWorkbook(wb, file = "../2_Output/QC/Supplementary_QC_ACTA2+MERGE.xlsx", overwrite = TRUE)

message("QC summary written to ../2_Output/QC/Supplementary_QC_ACTA2+MERGE.xlsx")

# load tidyverse
library(dplyr)

# filter MERGE table for TAAng and ACTA2 samples
qc_summary <- tab_MERGE %>%
  group_by(Disease) %>%
  summarise(
    Cells = n(),
    Mean_nUMI = mean(nUMI, na.rm = TRUE), Median_nUMI = median(nUMI, na.rm = TRUE), SD_nUMI = sd(nUMI, na.rm = TRUE),
    Mean_nGene = mean(nGene, na.rm = TRUE), Median_nGene = median(nGene, na.rm = TRUE), SD_nGene = sd(nGene, na.rm = TRUE),
    Mean_pctMT = mean(pctMT, na.rm = TRUE), Median_pctMT = median(pctMT, na.rm = TRUE), SD_pctMT = sd(pctMT, na.rm = TRUE),
    Mean_DoubletScore = mean(DoubletScore, na.rm = TRUE), Median_DoubletScore = median(DoubletScore, na.rm = TRUE), SD_DoubletScore = sd(DoubletScore, na.rm = TRUE)
  )
qc_summary


# ---- QC VISUALS: nUMI, nGene, pctMT by Disease ----
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(openxlsx)
})

# 0) output dir ----------------------------------------------------------------
out_dir <- "../2_Output/QC/figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 1) reshape to long format ----------------------------------------------------
qc_long <- tab_MERGE %>%
  dplyr::select(Disease, nUMI, nGene, pctMT) %>%
  pivot_longer(cols = c(nUMI, nGene, pctMT), names_to = "Metric", values_to = "Value") %>%
  filter(is.finite(Value))

# Optional: set preferred Disease order if both exist
pref_levels <- c("TAA_NG", "TAA_ACTA2")
qc_long$Disease <- factor(qc_long$Disease, levels = intersect(pref_levels, unique(qc_long$Disease)))

# 2) summary stats (for caption/supplement) ------------------------------------
qc_stats <- qc_long %>%
  group_by(Metric, Disease) %>%
  summarise(
    n       = dplyr::n(),
    mean    = mean(Value, na.rm = TRUE),
    sd      = sd(Value, na.rm = TRUE),
    median  = median(Value, na.rm = TRUE),
    IQR     = IQR(Value, na.rm = TRUE),
    .groups = "drop"
  )

write.xlsx(qc_stats, file.path(out_dir, "QC_summary_stats.xlsx"), overwrite = TRUE)

# 3) mean ± SD helper for stat_summary ----------------------------------------
mean_sd <- function(x) {
  m <- mean(x, na.rm = TRUE); s <- sd(x, na.rm = TRUE)
  data.frame(y = m, ymin = m - s, ymax = m + s)
}

# 4) faceted violin + boxplot with mean ± SD ----------------------------------
p <- ggplot(qc_long, aes(x = Disease, y = Value, fill = Disease)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.6, linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9, linewidth = 0.3) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2,
               fill = "white", color = "black") +
  stat_summary(fun.data = mean_sd, geom = "errorbar", width = 0.1, size = 0.3) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 3) +
  labs(x = NULL, y = NULL, title = "QC metrics by Disease (nUMI, nGene, mito%)") +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey92", colour = NA),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(out_dir, "QC_violin_box_facet.pdf"), p,
       width = 10, height = 3.6, useDingbats = FALSE)
ggsave(file.path(out_dir, "QC_violin_box_facet.png"), p,
       width = 10, height = 3.6, dpi = 300)
```

# Heatmap - Proliferation


``` r
# ---- Reviewer 6: Heatmap of DEGs tied to VSMC proliferation -----------------
# This chunk builds a heatmap of differentially expressed genes (DEGs)
# associated with cell proliferation (GO:0008283) within VSMCs, comparing
# TAA_ACTA2 vs TAA_NG in your MERGE object.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)      # HUMAN gene annotations
  library(pheatmap)          # for pseudobulk heatmap
  library(openxlsx)
})

# Inputs/outputs
fp_merged <- "../1_Input/TAA_snRNA.rds"               # your merged object
out_dir   <- "../2_Output/3_Merged/3_Differential.Expression"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Load and prep MERGE
MERGE <- readRDS(fp_merged)
if (!is.null(MERGE[["RNA"]]@layers) && length(MERGE[["RNA"]]@layers) > 0) {
  MERGE[["RNA"]] <- JoinLayers(MERGE[["RNA"]])
}
DefaultAssay(MERGE) <- "RNA"

# ---- define VSMC subset (handles either VSMC, VSMC_1/2 labels) --------------
vsmc_levels <- c("VSMC", "VSMC_1", "VSMC_2")
if (!"cell_type" %in% colnames(MERGE@meta.data)) {
  stop("MERGE lacks 'cell_type' metadata required to subset VSMCs.")
}
vsmc <- subset(MERGE, subset = cell_type %in% vsmc_levels)
if (ncol(vsmc) == 0) stop("No VSMC cells found (checked labels: VSMC, VSMC_1, VSMC_2).")

# Set identities for DE test
if (!"Disease" %in% colnames(vsmc@meta.data)) {
  stop("MERGE lacks 'Disease' metadata required for group comparison.")
}
vsmc$Disease <- factor(vsmc$Disease, levels = c("TAA_NG","TAA_ACTA2"))
Idents(vsmc) <- vsmc$Disease

# ---- get GO:0008283 (cell proliferation) gene set (human) --------------------
go_term <- "GO:0008283"  # cell proliferation (as used in your earlier analysis)
go_genes_df <- suppressWarnings(
  AnnotationDbi::select(org.Hs.eg.db,
                        keys = go_term,
                        columns = "SYMBOL",
                        keytype = "GOALL")
)
go_genes <- unique(na.omit(go_genes_df$SYMBOL))
go_genes <- go_genes[go_genes %in% rownames(vsmc)]

if (length(go_genes) < 5) {
  warning("Few/no GO:0008283 symbols found in your dataset; heatmap may be sparse.")
}

# ---- DE analysis: TAA_ACTA2 vs TAA_NG within VSMC ---------------------------
de_vsmc <- FindMarkers(
  vsmc,
  ident.1 = "TAA_ACTA2",
  ident.2 = "TAA_NG",
  test.use = "MAST",
  min.pct = 0.1,
  logfc.threshold = 0,   # keep broad set, we filter below
  assay = "RNA"
)

de_vsmc$gene <- rownames(de_vsmc)

# filter to proliferation gene set + significance
de_vsmc_prolif <- de_vsmc %>%
  filter(gene %in% go_genes) %>%
  mutate(p_adj = ifelse(!is.na(p_val_adj), p_val_adj, p_val)) %>%
  arrange(desc(abs(avg_log2FC)))

# pick top genes for visualization
top_n_genes <- 30
sel_genes <- head(de_vsmc_prolif$gene[de_vsmc_prolif$p_adj < 0.05], top_n_genes)
if (length(sel_genes) < 5) {
  # fallback: take top by |log2FC| (even if not all pass adj p)
  sel_genes <- head(de_vsmc_prolif$gene, min(top_n_genes, nrow(de_vsmc_prolif)))
}
if (length(sel_genes) == 0) stop("No proliferation-associated DEGs found to plot.")

# ---- Option A: cell-level heatmap (Seurat::DoHeatmap) -----------------------
pdf(file.path(out_dir, "VSMC_Proliferation_DEGs_DoHeatmap.pdf"), width = 8, height = 6)
print(
  DoHeatmap(
    vsmc,
    features = sel_genes,
    group.by = "Disease",
    disp.min = -2, disp.max = 2
  ) + ggtitle("VSMC: Proliferation-associated DEGs (TAA_ACTA2 vs TAA_NG)")
)
dev.off()

# ---- Option B: pseudobulk (average by Disease) heatmap ----------------------
avg <- AverageExpression(
  vsmc,
  assays = "RNA",
  features = sel_genes,
  group.by = "Disease",
  slot = "data"  # log-normalized
)$RNA

# z-score by gene for visualization
mat <- t(scale(t(as.matrix(avg))))  # gene-wise z-score
mat[is.na(mat)] <- 0

pdf(file.path(out_dir, "VSMC_Proliferation_DEGs_pseudobulk.pdf"), width = 6, height = 8)
pheatmap::pheatmap(
  mat,
  color = colorRampPalette(c("dodgerblue4","white","coral2"))(101),
  border_color = NA,
  cluster_cols = FALSE,
  main = "VSMC: Proliferation-associated DEGs (pseudobulk by Disease)"
)
dev.off()

# ---- export the DEG table used for the heatmaps -----------------------------
openxlsx::write.xlsx(
  list(
    DE_full = de_vsmc,
    DE_proliferation = de_vsmc_prolif,
    Heatmap_genes = data.frame(gene = sel_genes)
  ),
  file = file.path(out_dir, "VSMC_Proliferation_DEGs.xlsx"),
  overwrite = TRUE
)

message("Saved heatmaps and tables to: ", out_dir)

# Define canonical proliferation-associated genes
prolif_genes <- c("KLF4","MYOCD","SRF","FOXO3",
                  "CCND1","CCNE1","CCNA2","CCNB1",
                  "CDK1","CDK2","CDK4","CDK6",
                  "E2F1","E2F2",
                  "PCNA","MKI67","MYC",
                  "PDGFRB","PDGFRA","PDGFA","PDGFB",
                  "FGF2","EGFR","STAT3",
                  "CDKN1A","CDKN1B","RB1")
Idents(MERGE) <- MERGE$Disease  # Set identities for differential expression
# Run differential expression test between groups (adjusted for your dataset)
deg_results <- FindMarkers(MERGE, ident.1 = "TAA_NG", ident.2 = "TAA_ACTA2", 
                           group.by = "Disease", features = prolif_genes, 
                           logfc.threshold = 0, min.pct = 0)

# Save to CSV for supplemental upload
write.csv(deg_results, "DEG_proliferation_genes.csv", row.names = TRUE)

# DotPlot separating groups by Disease (TAA_ACTA2 vs TAA_NG)
DotPlot(
  object = MERGE, 
  features = prolif_genes, 
  group.by = "Disease"
) + 
  RotatedAxis() + 
  ggtitle("Expression of proliferation-associated genes in TAA_ACTA2 vs TAA_NG") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Define a curated list of canonical apoptosis-related genes
apoptosis_genes <- c("BAX","BAK1","BCL2","BCL2L1","BCL2L11",
                     "CASP3","CASP7","CASP8","CASP9",
                     "FAS","FADD","TNFRSF10A","TNFRSF10B","TP53")

# Subset differentially expressed genes between TAA_ACTA2 and TAA_NG
deg_apoptosis <- subset(my_deg_results, gene %in% apoptosis_genes)

# Print table of apoptosis DEGs
deg_apoptosis

# Optionally, visualize apoptosis gene expression across TAA_ACTA2 vs TAA_NG
VlnPlot(MERGE, features = c("BAX","BCL2","CASP3","TP53"), 
        group.by = "cell_type", split.by = "Disease", 
        pt.size = 0, assay = "RNA") + 
  ggtitle("Apoptosis-related gene expression in TAA_ACTA2 vs TAA_NG")
```

# Supplemental Table: R Session Information

All packages and setting are acquired using the following command.


``` r
# Write 
options(kableExtra.latex.load_packages = FALSE)
library(kableExtra)
sinfo<-devtools::session_info()
sinfo$platform
```

```
##  setting  value
##  version  R version 4.5.0 (2025-04-11)
##  os       macOS Sequoia 15.6
##  system   aarch64, darwin20
##  ui       X11
##  language (EN)
##  collate  en_US.UTF-8
##  ctype    en_US.UTF-8
##  tz       America/Los_Angeles
##  date     2025-08-19
##  pandoc   3.4 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64/ (via rmarkdown)
##  quarto   1.6.42 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto
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
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/abind </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/abind </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-12 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AnnotationDbi </td>
   <td style="text-align:center;"> AnnotationDbi </td>
   <td style="text-align:center;"> 1.70.0 </td>
   <td style="text-align:center;"> 1.70.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/AnnotationDbi </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/AnnotationDbi </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-04-15 </td>
   <td style="text-align:center;"> Bioconductor 3.21 (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ape </td>
   <td style="text-align:center;"> ape </td>
   <td style="text-align:center;"> 5.8.1 </td>
   <td style="text-align:center;"> 5.8-1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/ape </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/ape </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-12-16 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> aplot </td>
   <td style="text-align:center;"> aplot </td>
   <td style="text-align:center;"> 0.2.8 </td>
   <td style="text-align:center;"> 0.2.8 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/aplot </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/aplot </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-07-02 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Biobase </td>
   <td style="text-align:center;"> Biobase </td>
   <td style="text-align:center;"> 2.68.0 </td>
   <td style="text-align:center;"> 2.68.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/Biobase </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/Biobase </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-04-15 </td>
   <td style="text-align:center;"> Bioconductor 3.21 (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> BiocGenerics </td>
   <td style="text-align:center;"> BiocGenerics </td>
   <td style="text-align:center;"> 0.54.0 </td>
   <td style="text-align:center;"> 0.54.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/BiocGenerics </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/BiocGenerics </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-04-15 </td>
   <td style="text-align:center;"> Bioconductor 3.21 (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> BiocParallel </td>
   <td style="text-align:center;"> BiocParallel </td>
   <td style="text-align:center;"> 1.42.1 </td>
   <td style="text-align:center;"> 1.42.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/BiocParallel </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/BiocParallel </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-05-29 </td>
   <td style="text-align:center;"> Bioconductor 3.21 (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Biostrings </td>
   <td style="text-align:center;"> Biostrings </td>
   <td style="text-align:center;"> 2.76.0 </td>
   <td style="text-align:center;"> 2.76.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/Biostrings </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/Biostrings </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-04-15 </td>
   <td style="text-align:center;"> Bioconductor 3.21 (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> bit </td>
   <td style="text-align:center;"> bit </td>
   <td style="text-align:center;"> 4.6.0 </td>
   <td style="text-align:center;"> 4.6.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/bit </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/bit </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-03-06 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> bit64 </td>
   <td style="text-align:center;"> bit64 </td>
   <td style="text-align:center;"> 4.6.0.1 </td>
   <td style="text-align:center;"> 4.6.0-1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/bit64 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/bit64 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-01-16 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> blob </td>
   <td style="text-align:center;"> blob </td>
   <td style="text-align:center;"> 1.2.4 </td>
   <td style="text-align:center;"> 1.2.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/blob </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/blob </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-03-17 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> bslib </td>
   <td style="text-align:center;"> bslib </td>
   <td style="text-align:center;"> 0.9.0 </td>
   <td style="text-align:center;"> 0.9.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/bslib </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/bslib </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-01-30 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cachem </td>
   <td style="text-align:center;"> cachem </td>
   <td style="text-align:center;"> 1.1.0 </td>
   <td style="text-align:center;"> 1.1.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/cachem </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/cachem </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-05-16 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cli </td>
   <td style="text-align:center;"> cli </td>
   <td style="text-align:center;"> 3.6.5 </td>
   <td style="text-align:center;"> 3.6.5 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/cli </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/cli </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-04-23 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cluster </td>
   <td style="text-align:center;"> cluster </td>
   <td style="text-align:center;"> 2.1.8.1 </td>
   <td style="text-align:center;"> 2.1.8.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/cluster </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/cluster </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-03-12 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> clusterProfiler </td>
   <td style="text-align:center;"> clusterProfiler </td>
   <td style="text-align:center;"> 4.16.0 </td>
   <td style="text-align:center;"> 4.16.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/clusterProfiler </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/clusterProfiler </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-04-15 </td>
   <td style="text-align:center;"> Bioconductor 3.21 (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> codetools </td>
   <td style="text-align:center;"> codetools </td>
   <td style="text-align:center;"> 0.2.20 </td>
   <td style="text-align:center;"> 0.2-20 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/codetools </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/codetools </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-03-31 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> colorspace </td>
   <td style="text-align:center;"> colorspace </td>
   <td style="text-align:center;"> 2.1.1 </td>
   <td style="text-align:center;"> 2.1-1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/colorspace </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/colorspace </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-26 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cowplot </td>
   <td style="text-align:center;"> cowplot </td>
   <td style="text-align:center;"> 1.2.0 </td>
   <td style="text-align:center;"> 1.2.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/cowplot </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/cowplot </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-07-07 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> crayon </td>
   <td style="text-align:center;"> crayon </td>
   <td style="text-align:center;"> 1.5.3 </td>
   <td style="text-align:center;"> 1.5.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/crayon </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/crayon </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-06-20 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:center;"> data.table </td>
   <td style="text-align:center;"> 1.17.8 </td>
   <td style="text-align:center;"> 1.17.8 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/data.table </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/data.table </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-07-10 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> DBI </td>
   <td style="text-align:center;"> DBI </td>
   <td style="text-align:center;"> 1.2.3 </td>
   <td style="text-align:center;"> 1.2.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/DBI </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/DBI </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-06-02 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> DelayedArray </td>
   <td style="text-align:center;"> DelayedArray </td>
   <td style="text-align:center;"> 0.34.1 </td>
   <td style="text-align:center;"> 0.34.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/DelayedArray </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/DelayedArray </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-04-17 </td>
   <td style="text-align:center;"> Bioconductor 3.21 (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> deldir </td>
   <td style="text-align:center;"> deldir </td>
   <td style="text-align:center;"> 2.0.4 </td>
   <td style="text-align:center;"> 2.0-4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/deldir </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/deldir </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-02-28 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> devtools </td>
   <td style="text-align:center;"> devtools </td>
   <td style="text-align:center;"> 2.4.5 </td>
   <td style="text-align:center;"> 2.4.5 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/devtools </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/devtools </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-10-11 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> dichromat </td>
   <td style="text-align:center;"> dichromat </td>
   <td style="text-align:center;"> 2.0.0.1 </td>
   <td style="text-align:center;"> 2.0-0.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/dichromat </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/dichromat </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-05-02 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> digest </td>
   <td style="text-align:center;"> digest </td>
   <td style="text-align:center;"> 0.6.37 </td>
   <td style="text-align:center;"> 0.6.37 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/digest </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/digest </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-08-19 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> dittoSeq </td>
   <td style="text-align:center;"> dittoSeq </td>
   <td style="text-align:center;"> 1.20.0 </td>
   <td style="text-align:center;"> 1.20.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/dittoSeq </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/dittoSeq </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-04-15 </td>
   <td style="text-align:center;"> Bioconductor 3.21 (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> DOSE </td>
   <td style="text-align:center;"> DOSE </td>
   <td style="text-align:center;"> 4.2.0 </td>
   <td style="text-align:center;"> 4.2.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/DOSE </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/DOSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-04-15 </td>
   <td style="text-align:center;"> Bioconductor 3.21 (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> dotCall64 </td>
   <td style="text-align:center;"> dotCall64 </td>
   <td style="text-align:center;"> 1.2 </td>
   <td style="text-align:center;"> 1.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/dotCall64 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/dotCall64 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-10-04 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> dplyr </td>
   <td style="text-align:center;"> dplyr </td>
   <td style="text-align:center;"> 1.1.4 </td>
   <td style="text-align:center;"> 1.1.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/dplyr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/dplyr </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-11-17 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ellipsis </td>
   <td style="text-align:center;"> ellipsis </td>
   <td style="text-align:center;"> 0.3.2 </td>
   <td style="text-align:center;"> 0.3.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/ellipsis </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/ellipsis </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-04-29 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> enrichplot </td>
   <td style="text-align:center;"> enrichplot </td>
   <td style="text-align:center;"> 1.28.4 </td>
   <td style="text-align:center;"> 1.28.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/enrichplot </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/enrichplot </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-07-14 </td>
   <td style="text-align:center;"> Bioconductor 3.21 (R 4.5.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> evaluate </td>
   <td style="text-align:center;"> evaluate </td>
   <td style="text-align:center;"> 1.0.4 </td>
   <td style="text-align:center;"> 1.0.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/evaluate </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/evaluate </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-06-18 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> farver </td>
   <td style="text-align:center;"> farver </td>
   <td style="text-align:center;"> 2.1.2 </td>
   <td style="text-align:center;"> 2.1.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/farver </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/farver </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-05-13 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fastDummies </td>
   <td style="text-align:center;"> fastDummies </td>
   <td style="text-align:center;"> 1.7.5 </td>
   <td style="text-align:center;"> 1.7.5 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/fastDummies </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/fastDummies </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-01-20 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fastmap </td>
   <td style="text-align:center;"> fastmap </td>
   <td style="text-align:center;"> 1.2.0 </td>
   <td style="text-align:center;"> 1.2.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/fastmap </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/fastmap </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-05-15 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fastmatch </td>
   <td style="text-align:center;"> fastmatch </td>
   <td style="text-align:center;"> 1.1.6 </td>
   <td style="text-align:center;"> 1.1-6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/fastmatch </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/fastmatch </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-12-23 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fgsea </td>
   <td style="text-align:center;"> fgsea </td>
   <td style="text-align:center;"> 1.34.2 </td>
   <td style="text-align:center;"> 1.34.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/fgsea </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/fgsea </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-07-10 </td>
   <td style="text-align:center;"> Bioconductor 3.21 (R 4.5.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fitdistrplus </td>
   <td style="text-align:center;"> fitdistrplus </td>
   <td style="text-align:center;"> 1.2.4 </td>
   <td style="text-align:center;"> 1.2-4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/fitdistrplus </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/fitdistrplus </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-07-03 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fs </td>
   <td style="text-align:center;"> fs </td>
   <td style="text-align:center;"> 1.6.6 </td>
   <td style="text-align:center;"> 1.6.6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/fs </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/fs </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-04-12 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> future </td>
   <td style="text-align:center;"> future </td>
   <td style="text-align:center;"> 1.67.0 </td>
   <td style="text-align:center;"> 1.67.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/future </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/future </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-07-29 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> future.apply </td>
   <td style="text-align:center;"> future.apply </td>
   <td style="text-align:center;"> 1.20.0 </td>
   <td style="text-align:center;"> 1.20.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/future.apply </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/future.apply </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-06-06 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> generics </td>
   <td style="text-align:center;"> generics </td>
   <td style="text-align:center;"> 0.1.4 </td>
   <td style="text-align:center;"> 0.1.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/generics </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/generics </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-05-09 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GenomeInfoDb </td>
   <td style="text-align:center;"> GenomeInfoDb </td>
   <td style="text-align:center;"> 1.44.1 </td>
   <td style="text-align:center;"> 1.44.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/GenomeInfoDb </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/GenomeInfoDb </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-07-21 </td>
   <td style="text-align:center;"> Bioconductor 3.21 (R 4.5.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GenomeInfoDbData </td>
   <td style="text-align:center;"> GenomeInfoDbData </td>
   <td style="text-align:center;"> 1.2.14 </td>
   <td style="text-align:center;"> 1.2.14 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/GenomeInfoDbData </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/GenomeInfoDbData </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-05-12 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GenomicRanges </td>
   <td style="text-align:center;"> GenomicRanges </td>
   <td style="text-align:center;"> 1.60.0 </td>
   <td style="text-align:center;"> 1.60.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/GenomicRanges </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/GenomicRanges </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-04-15 </td>
   <td style="text-align:center;"> Bioconductor 3.21 (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ggfun </td>
   <td style="text-align:center;"> ggfun </td>
   <td style="text-align:center;"> 0.2.0 </td>
   <td style="text-align:center;"> 0.2.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/ggfun </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/ggfun </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-07-15 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ggplot2 </td>
   <td style="text-align:center;"> ggplot2 </td>
   <td style="text-align:center;"> 3.5.2 </td>
   <td style="text-align:center;"> 3.5.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/ggplot2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/ggplot2 </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-04-09 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ggplotify </td>
   <td style="text-align:center;"> ggplotify </td>
   <td style="text-align:center;"> 0.1.2 </td>
   <td style="text-align:center;"> 0.1.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/ggplotify </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/ggplotify </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-08-09 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ggrepel </td>
   <td style="text-align:center;"> ggrepel </td>
   <td style="text-align:center;"> 0.9.6 </td>
   <td style="text-align:center;"> 0.9.6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/ggrepel </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/ggrepel </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-07 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ggridges </td>
   <td style="text-align:center;"> ggridges </td>
   <td style="text-align:center;"> 0.5.6 </td>
   <td style="text-align:center;"> 0.5.6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/ggridges </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/ggridges </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-01-23 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ggtangle </td>
   <td style="text-align:center;"> ggtangle </td>
   <td style="text-align:center;"> 0.0.7 </td>
   <td style="text-align:center;"> 0.0.7 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/ggtangle </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/ggtangle </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-06-30 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ggtree </td>
   <td style="text-align:center;"> ggtree </td>
   <td style="text-align:center;"> 3.16.3 </td>
   <td style="text-align:center;"> 3.16.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/ggtree </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/ggtree </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-07-14 </td>
   <td style="text-align:center;"> Bioconductor 3.21 (R 4.5.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> globals </td>
   <td style="text-align:center;"> globals </td>
   <td style="text-align:center;"> 0.18.0 </td>
   <td style="text-align:center;"> 0.18.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/globals </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/globals </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-05-08 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> glue </td>
   <td style="text-align:center;"> glue </td>
   <td style="text-align:center;"> 1.8.0 </td>
   <td style="text-align:center;"> 1.8.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/glue </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/glue </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-30 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO.db </td>
   <td style="text-align:center;"> GO.db </td>
   <td style="text-align:center;"> 3.21.0 </td>
   <td style="text-align:center;"> 3.21.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/GO.db </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/GO.db </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-05-12 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> goftest </td>
   <td style="text-align:center;"> goftest </td>
   <td style="text-align:center;"> 1.2.3 </td>
   <td style="text-align:center;"> 1.2-3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/goftest </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/goftest </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-10-07 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GOSemSim </td>
   <td style="text-align:center;"> GOSemSim </td>
   <td style="text-align:center;"> 2.34.0 </td>
   <td style="text-align:center;"> 2.34.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/GOSemSim </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/GOSemSim </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-04-15 </td>
   <td style="text-align:center;"> Bioconductor 3.21 (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> gridExtra </td>
   <td style="text-align:center;"> gridExtra </td>
   <td style="text-align:center;"> 2.3 </td>
   <td style="text-align:center;"> 2.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/gridExtra </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/gridExtra </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2017-09-09 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> gridGraphics </td>
   <td style="text-align:center;"> gridGraphics </td>
   <td style="text-align:center;"> 0.5.1 </td>
   <td style="text-align:center;"> 0.5-1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/gridGraphics </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/gridGraphics </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-12-13 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> gson </td>
   <td style="text-align:center;"> gson </td>
   <td style="text-align:center;"> 0.1.0 </td>
   <td style="text-align:center;"> 0.1.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/gson </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/gson </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-03-07 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> gtable </td>
   <td style="text-align:center;"> gtable </td>
   <td style="text-align:center;"> 0.3.6 </td>
   <td style="text-align:center;"> 0.3.6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/gtable </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/gtable </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-10-25 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> htmltools </td>
   <td style="text-align:center;"> htmltools </td>
   <td style="text-align:center;"> 0.5.8.1 </td>
   <td style="text-align:center;"> 0.5.8.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/htmltools </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/htmltools </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-04-04 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> htmlwidgets </td>
   <td style="text-align:center;"> htmlwidgets </td>
   <td style="text-align:center;"> 1.6.4 </td>
   <td style="text-align:center;"> 1.6.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/htmlwidgets </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/htmlwidgets </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-12-06 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> httpuv </td>
   <td style="text-align:center;"> httpuv </td>
   <td style="text-align:center;"> 1.6.16 </td>
   <td style="text-align:center;"> 1.6.16 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/httpuv </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/httpuv </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-04-16 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> httr </td>
   <td style="text-align:center;"> httr </td>
   <td style="text-align:center;"> 1.4.7 </td>
   <td style="text-align:center;"> 1.4.7 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/httr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/httr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-08-15 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ica </td>
   <td style="text-align:center;"> ica </td>
   <td style="text-align:center;"> 1.0.3 </td>
   <td style="text-align:center;"> 1.0-3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/ica </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/ica </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-07-08 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> igraph </td>
   <td style="text-align:center;"> igraph </td>
   <td style="text-align:center;"> 2.1.4 </td>
   <td style="text-align:center;"> 2.1.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/igraph </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/igraph </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-01-23 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IRanges </td>
   <td style="text-align:center;"> IRanges </td>
   <td style="text-align:center;"> 2.42.0 </td>
   <td style="text-align:center;"> 2.42.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/IRanges </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/IRanges </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-04-15 </td>
   <td style="text-align:center;"> Bioconductor 3.21 (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> irlba </td>
   <td style="text-align:center;"> irlba </td>
   <td style="text-align:center;"> 2.3.5.1 </td>
   <td style="text-align:center;"> 2.3.5.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/irlba </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/irlba </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-10-03 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> jquerylib </td>
   <td style="text-align:center;"> jquerylib </td>
   <td style="text-align:center;"> 0.1.4 </td>
   <td style="text-align:center;"> 0.1.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/jquerylib </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/jquerylib </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-04-26 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> jsonlite </td>
   <td style="text-align:center;"> jsonlite </td>
   <td style="text-align:center;"> 2.0.0 </td>
   <td style="text-align:center;"> 2.0.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/jsonlite </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/jsonlite </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-03-27 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> kableExtra </td>
   <td style="text-align:center;"> kableExtra </td>
   <td style="text-align:center;"> 1.4.0 </td>
   <td style="text-align:center;"> 1.4.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/kableExtra </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/kableExtra </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-01-24 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KEGGREST </td>
   <td style="text-align:center;"> KEGGREST </td>
   <td style="text-align:center;"> 1.48.1 </td>
   <td style="text-align:center;"> 1.48.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/KEGGREST </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/KEGGREST </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-06-19 </td>
   <td style="text-align:center;"> Bioconductor 3.21 (R 4.5.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KernSmooth </td>
   <td style="text-align:center;"> KernSmooth </td>
   <td style="text-align:center;"> 2.23.26 </td>
   <td style="text-align:center;"> 2.23-26 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/KernSmooth </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/KernSmooth </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-01-01 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> knitr </td>
   <td style="text-align:center;"> knitr </td>
   <td style="text-align:center;"> 1.50 </td>
   <td style="text-align:center;"> 1.50 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/knitr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/knitr </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-03-16 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> later </td>
   <td style="text-align:center;"> later </td>
   <td style="text-align:center;"> 1.4.2 </td>
   <td style="text-align:center;"> 1.4.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/later </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/later </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-04-08 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lattice </td>
   <td style="text-align:center;"> lattice </td>
   <td style="text-align:center;"> 0.22.7 </td>
   <td style="text-align:center;"> 0.22-7 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/lattice </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/lattice </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-04-02 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lazyeval </td>
   <td style="text-align:center;"> lazyeval </td>
   <td style="text-align:center;"> 0.2.2 </td>
   <td style="text-align:center;"> 0.2.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/lazyeval </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/lazyeval </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-03-15 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lifecycle </td>
   <td style="text-align:center;"> lifecycle </td>
   <td style="text-align:center;"> 1.0.4 </td>
   <td style="text-align:center;"> 1.0.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/lifecycle </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/lifecycle </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-11-07 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> listenv </td>
   <td style="text-align:center;"> listenv </td>
   <td style="text-align:center;"> 0.9.1 </td>
   <td style="text-align:center;"> 0.9.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/listenv </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/listenv </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-01-29 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lmtest </td>
   <td style="text-align:center;"> lmtest </td>
   <td style="text-align:center;"> 0.9.40 </td>
   <td style="text-align:center;"> 0.9-40 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/lmtest </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/lmtest </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-03-21 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> magrittr </td>
   <td style="text-align:center;"> magrittr </td>
   <td style="text-align:center;"> 2.0.3 </td>
   <td style="text-align:center;"> 2.0.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/magrittr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/magrittr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-03-30 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MASS </td>
   <td style="text-align:center;"> MASS </td>
   <td style="text-align:center;"> 7.3.65 </td>
   <td style="text-align:center;"> 7.3-65 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/MASS </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/MASS </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-02-28 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Matrix </td>
   <td style="text-align:center;"> Matrix </td>
   <td style="text-align:center;"> 1.7.3 </td>
   <td style="text-align:center;"> 1.7-3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/Matrix </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/Matrix </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-03-11 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MatrixGenerics </td>
   <td style="text-align:center;"> MatrixGenerics </td>
   <td style="text-align:center;"> 1.20.0 </td>
   <td style="text-align:center;"> 1.20.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/MatrixGenerics </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/MatrixGenerics </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-04-15 </td>
   <td style="text-align:center;"> Bioconductor 3.21 (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> matrixStats </td>
   <td style="text-align:center;"> matrixStats </td>
   <td style="text-align:center;"> 1.5.0 </td>
   <td style="text-align:center;"> 1.5.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/matrixStats </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/matrixStats </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-01-07 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> memoise </td>
   <td style="text-align:center;"> memoise </td>
   <td style="text-align:center;"> 2.0.1 </td>
   <td style="text-align:center;"> 2.0.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/memoise </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/memoise </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-11-26 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mime </td>
   <td style="text-align:center;"> mime </td>
   <td style="text-align:center;"> 0.13 </td>
   <td style="text-align:center;"> 0.13 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/mime </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/mime </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-03-17 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miniUI </td>
   <td style="text-align:center;"> miniUI </td>
   <td style="text-align:center;"> 0.1.2 </td>
   <td style="text-align:center;"> 0.1.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/miniUI </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/miniUI </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-04-17 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> nlme </td>
   <td style="text-align:center;"> nlme </td>
   <td style="text-align:center;"> 3.1.168 </td>
   <td style="text-align:center;"> 3.1-168 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/nlme </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/nlme </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-03-31 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> openxlsx </td>
   <td style="text-align:center;"> openxlsx </td>
   <td style="text-align:center;"> 4.2.8 </td>
   <td style="text-align:center;"> 4.2.8 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/openxlsx </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/openxlsx </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-01-25 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> org.Hs.eg.db </td>
   <td style="text-align:center;"> org.Hs.eg.db </td>
   <td style="text-align:center;"> 3.21.0 </td>
   <td style="text-align:center;"> 3.21.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/org.Hs.eg.db </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/org.Hs.eg.db </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-05-12 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> parallelly </td>
   <td style="text-align:center;"> parallelly </td>
   <td style="text-align:center;"> 1.45.1 </td>
   <td style="text-align:center;"> 1.45.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/parallelly </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/parallelly </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-07-24 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> patchwork </td>
   <td style="text-align:center;"> patchwork </td>
   <td style="text-align:center;"> 1.3.1 </td>
   <td style="text-align:center;"> 1.3.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/patchwork </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/patchwork </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-06-21 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pbapply </td>
   <td style="text-align:center;"> pbapply </td>
   <td style="text-align:center;"> 1.7.4 </td>
   <td style="text-align:center;"> 1.7-4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/pbapply </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/pbapply </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-07-20 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pheatmap </td>
   <td style="text-align:center;"> pheatmap </td>
   <td style="text-align:center;"> 1.0.13 </td>
   <td style="text-align:center;"> 1.0.13 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/pheatmap </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/pheatmap </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-06-05 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pillar </td>
   <td style="text-align:center;"> pillar </td>
   <td style="text-align:center;"> 1.11.0 </td>
   <td style="text-align:center;"> 1.11.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/pillar </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/pillar </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-07-04 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pkgbuild </td>
   <td style="text-align:center;"> pkgbuild </td>
   <td style="text-align:center;"> 1.4.8 </td>
   <td style="text-align:center;"> 1.4.8 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/pkgbuild </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/pkgbuild </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-05-26 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pkgconfig </td>
   <td style="text-align:center;"> pkgconfig </td>
   <td style="text-align:center;"> 2.0.3 </td>
   <td style="text-align:center;"> 2.0.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/pkgconfig </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/pkgconfig </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-09-22 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pkgload </td>
   <td style="text-align:center;"> pkgload </td>
   <td style="text-align:center;"> 1.4.0 </td>
   <td style="text-align:center;"> 1.4.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/pkgload </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/pkgload </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-06-28 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> plotly </td>
   <td style="text-align:center;"> plotly </td>
   <td style="text-align:center;"> 4.11.0 </td>
   <td style="text-align:center;"> 4.11.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/plotly </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/plotly </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-06-19 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> plyr </td>
   <td style="text-align:center;"> plyr </td>
   <td style="text-align:center;"> 1.8.9 </td>
   <td style="text-align:center;"> 1.8.9 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/plyr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/plyr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-10-02 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> png </td>
   <td style="text-align:center;"> png </td>
   <td style="text-align:center;"> 0.1.8 </td>
   <td style="text-align:center;"> 0.1-8 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/png </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/png </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-11-29 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polyclip </td>
   <td style="text-align:center;"> polyclip </td>
   <td style="text-align:center;"> 1.10.7 </td>
   <td style="text-align:center;"> 1.10-7 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/polyclip </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/polyclip </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-23 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> profvis </td>
   <td style="text-align:center;"> profvis </td>
   <td style="text-align:center;"> 0.4.0 </td>
   <td style="text-align:center;"> 0.4.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/profvis </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/profvis </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-20 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> progressr </td>
   <td style="text-align:center;"> progressr </td>
   <td style="text-align:center;"> 0.15.1 </td>
   <td style="text-align:center;"> 0.15.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/progressr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/progressr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-11-22 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> promises </td>
   <td style="text-align:center;"> promises </td>
   <td style="text-align:center;"> 1.3.3 </td>
   <td style="text-align:center;"> 1.3.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/promises </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/promises </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-05-29 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> purrr </td>
   <td style="text-align:center;"> purrr </td>
   <td style="text-align:center;"> 1.1.0 </td>
   <td style="text-align:center;"> 1.1.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/purrr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/purrr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-07-10 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> qvalue </td>
   <td style="text-align:center;"> qvalue </td>
   <td style="text-align:center;"> 2.40.0 </td>
   <td style="text-align:center;"> 2.40.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/qvalue </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/qvalue </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-04-15 </td>
   <td style="text-align:center;"> Bioconductor 3.21 (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> R.methodsS3 </td>
   <td style="text-align:center;"> R.methodsS3 </td>
   <td style="text-align:center;"> 1.8.2 </td>
   <td style="text-align:center;"> 1.8.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/R.methodsS3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/R.methodsS3 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-06-13 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> R.oo </td>
   <td style="text-align:center;"> R.oo </td>
   <td style="text-align:center;"> 1.27.1 </td>
   <td style="text-align:center;"> 1.27.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/R.oo </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/R.oo </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-05-02 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> R.utils </td>
   <td style="text-align:center;"> R.utils </td>
   <td style="text-align:center;"> 2.13.0 </td>
   <td style="text-align:center;"> 2.13.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/R.utils </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/R.utils </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-02-24 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> R6 </td>
   <td style="text-align:center;"> R6 </td>
   <td style="text-align:center;"> 2.6.1 </td>
   <td style="text-align:center;"> 2.6.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/R6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/R6 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-02-15 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RANN </td>
   <td style="text-align:center;"> RANN </td>
   <td style="text-align:center;"> 2.6.2 </td>
   <td style="text-align:center;"> 2.6.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/RANN </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/RANN </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-08-25 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RColorBrewer </td>
   <td style="text-align:center;"> RColorBrewer </td>
   <td style="text-align:center;"> 1.1.3 </td>
   <td style="text-align:center;"> 1.1-3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/RColorBrewer </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/RColorBrewer </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-04-03 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Rcpp </td>
   <td style="text-align:center;"> Rcpp </td>
   <td style="text-align:center;"> 1.1.0 </td>
   <td style="text-align:center;"> 1.1.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/Rcpp </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/Rcpp </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-07-02 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RcppAnnoy </td>
   <td style="text-align:center;"> RcppAnnoy </td>
   <td style="text-align:center;"> 0.0.22 </td>
   <td style="text-align:center;"> 0.0.22 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/RcppAnnoy </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/RcppAnnoy </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-01-23 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RcppHNSW </td>
   <td style="text-align:center;"> RcppHNSW </td>
   <td style="text-align:center;"> 0.6.0 </td>
   <td style="text-align:center;"> 0.6.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/RcppHNSW </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/RcppHNSW </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-02-04 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> remotes </td>
   <td style="text-align:center;"> remotes </td>
   <td style="text-align:center;"> 2.5.0 </td>
   <td style="text-align:center;"> 2.5.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/remotes </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/remotes </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-03-17 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> reshape2 </td>
   <td style="text-align:center;"> reshape2 </td>
   <td style="text-align:center;"> 1.4.4 </td>
   <td style="text-align:center;"> 1.4.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/reshape2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/reshape2 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-04-09 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> reticulate </td>
   <td style="text-align:center;"> reticulate </td>
   <td style="text-align:center;"> 1.43.0 </td>
   <td style="text-align:center;"> 1.43.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/reticulate </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/reticulate </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-07-21 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rlang </td>
   <td style="text-align:center;"> rlang </td>
   <td style="text-align:center;"> 1.1.6 </td>
   <td style="text-align:center;"> 1.1.6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/rlang </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/rlang </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-04-11 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rmarkdown </td>
   <td style="text-align:center;"> rmarkdown </td>
   <td style="text-align:center;"> 2.29 </td>
   <td style="text-align:center;"> 2.29 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/rmarkdown </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/rmarkdown </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-11-04 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ROCR </td>
   <td style="text-align:center;"> ROCR </td>
   <td style="text-align:center;"> 1.0.11 </td>
   <td style="text-align:center;"> 1.0-11 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/ROCR </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/ROCR </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-05-02 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RSpectra </td>
   <td style="text-align:center;"> RSpectra </td>
   <td style="text-align:center;"> 0.16.2 </td>
   <td style="text-align:center;"> 0.16-2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/RSpectra </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/RSpectra </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-18 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RSQLite </td>
   <td style="text-align:center;"> RSQLite </td>
   <td style="text-align:center;"> 2.4.2 </td>
   <td style="text-align:center;"> 2.4.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/RSQLite </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/RSQLite </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-07-18 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rstudioapi </td>
   <td style="text-align:center;"> rstudioapi </td>
   <td style="text-align:center;"> 0.17.1 </td>
   <td style="text-align:center;"> 0.17.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/rstudioapi </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/rstudioapi </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-10-22 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Rtsne </td>
   <td style="text-align:center;"> Rtsne </td>
   <td style="text-align:center;"> 0.17 </td>
   <td style="text-align:center;"> 0.17 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/Rtsne </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/Rtsne </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-12-07 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S4Arrays </td>
   <td style="text-align:center;"> S4Arrays </td>
   <td style="text-align:center;"> 1.8.1 </td>
   <td style="text-align:center;"> 1.8.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/S4Arrays </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/S4Arrays </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-05-29 </td>
   <td style="text-align:center;"> Bioconductor 3.21 (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S4Vectors </td>
   <td style="text-align:center;"> S4Vectors </td>
   <td style="text-align:center;"> 0.46.0 </td>
   <td style="text-align:center;"> 0.46.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/S4Vectors </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/S4Vectors </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-04-15 </td>
   <td style="text-align:center;"> Bioconductor 3.21 (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sass </td>
   <td style="text-align:center;"> sass </td>
   <td style="text-align:center;"> 0.4.10 </td>
   <td style="text-align:center;"> 0.4.10 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/sass </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/sass </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-04-11 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> scales </td>
   <td style="text-align:center;"> scales </td>
   <td style="text-align:center;"> 1.4.0 </td>
   <td style="text-align:center;"> 1.4.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/scales </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/scales </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-04-24 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> scattermore </td>
   <td style="text-align:center;"> scattermore </td>
   <td style="text-align:center;"> 1.2 </td>
   <td style="text-align:center;"> 1.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/scattermore </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/scattermore </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-06-12 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sctransform </td>
   <td style="text-align:center;"> sctransform </td>
   <td style="text-align:center;"> 0.4.2 </td>
   <td style="text-align:center;"> 0.4.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/sctransform </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/sctransform </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-04-30 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sessioninfo </td>
   <td style="text-align:center;"> sessioninfo </td>
   <td style="text-align:center;"> 1.2.3 </td>
   <td style="text-align:center;"> 1.2.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/sessioninfo </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/sessioninfo </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-02-05 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Seurat </td>
   <td style="text-align:center;"> Seurat </td>
   <td style="text-align:center;"> 5.3.0 </td>
   <td style="text-align:center;"> 5.3.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/Seurat </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/Seurat </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-04-23 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SeuratObject </td>
   <td style="text-align:center;"> SeuratObject </td>
   <td style="text-align:center;"> 5.1.0 </td>
   <td style="text-align:center;"> 5.1.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/SeuratObject </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/SeuratObject </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-04-22 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> shiny </td>
   <td style="text-align:center;"> shiny </td>
   <td style="text-align:center;"> 1.11.1.9000 </td>
   <td style="text-align:center;"> 1.11.1.9000 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/shiny </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/shiny </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-08-19 </td>
   <td style="text-align:center;"> Github (rstudio/shiny@0e355ed25cc1066d6894733f04f4b511a27acc53) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SingleCellExperiment </td>
   <td style="text-align:center;"> SingleCellExperiment </td>
   <td style="text-align:center;"> 1.30.1 </td>
   <td style="text-align:center;"> 1.30.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/SingleCellExperiment </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/SingleCellExperiment </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-05-05 </td>
   <td style="text-align:center;"> Bioconductor 3.21 (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sp </td>
   <td style="text-align:center;"> sp </td>
   <td style="text-align:center;"> 2.2.0 </td>
   <td style="text-align:center;"> 2.2-0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/sp </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/sp </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-02-01 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spam </td>
   <td style="text-align:center;"> spam </td>
   <td style="text-align:center;"> 2.11.1 </td>
   <td style="text-align:center;"> 2.11-1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/spam </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/spam </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-01-20 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SparseArray </td>
   <td style="text-align:center;"> SparseArray </td>
   <td style="text-align:center;"> 1.8.1 </td>
   <td style="text-align:center;"> 1.8.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/SparseArray </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/SparseArray </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-07-21 </td>
   <td style="text-align:center;"> Bioconductor 3.21 (R 4.5.1) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spatstat.data </td>
   <td style="text-align:center;"> spatstat.data </td>
   <td style="text-align:center;"> 3.1.6 </td>
   <td style="text-align:center;"> 3.1-6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/spatstat.data </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/spatstat.data </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-03-17 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spatstat.explore </td>
   <td style="text-align:center;"> spatstat.explore </td>
   <td style="text-align:center;"> 3.5.2 </td>
   <td style="text-align:center;"> 3.5-2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/spatstat.explore </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/spatstat.explore </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-07-22 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spatstat.geom </td>
   <td style="text-align:center;"> spatstat.geom </td>
   <td style="text-align:center;"> 3.5.0 </td>
   <td style="text-align:center;"> 3.5-0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/spatstat.geom </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/spatstat.geom </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-07-20 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spatstat.random </td>
   <td style="text-align:center;"> spatstat.random </td>
   <td style="text-align:center;"> 3.4.1 </td>
   <td style="text-align:center;"> 3.4-1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/spatstat.random </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/spatstat.random </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-05-20 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spatstat.sparse </td>
   <td style="text-align:center;"> spatstat.sparse </td>
   <td style="text-align:center;"> 3.1.0 </td>
   <td style="text-align:center;"> 3.1-0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/spatstat.sparse </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/spatstat.sparse </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-06-21 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spatstat.univar </td>
   <td style="text-align:center;"> spatstat.univar </td>
   <td style="text-align:center;"> 3.1.4 </td>
   <td style="text-align:center;"> 3.1-4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/spatstat.univar </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/spatstat.univar </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-07-13 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spatstat.utils </td>
   <td style="text-align:center;"> spatstat.utils </td>
   <td style="text-align:center;"> 3.1.5 </td>
   <td style="text-align:center;"> 3.1-5 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/spatstat.utils </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/spatstat.utils </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-07-17 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> stringi </td>
   <td style="text-align:center;"> stringi </td>
   <td style="text-align:center;"> 1.8.7 </td>
   <td style="text-align:center;"> 1.8.7 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/stringi </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/stringi </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-03-27 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> stringr </td>
   <td style="text-align:center;"> stringr </td>
   <td style="text-align:center;"> 1.5.1 </td>
   <td style="text-align:center;"> 1.5.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/stringr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/stringr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-11-14 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SummarizedExperiment </td>
   <td style="text-align:center;"> SummarizedExperiment </td>
   <td style="text-align:center;"> 1.38.1 </td>
   <td style="text-align:center;"> 1.38.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/SummarizedExperiment </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/SummarizedExperiment </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-04-28 </td>
   <td style="text-align:center;"> Bioconductor 3.21 (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> survival </td>
   <td style="text-align:center;"> survival </td>
   <td style="text-align:center;"> 3.8.3 </td>
   <td style="text-align:center;"> 3.8-3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/survival </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/survival </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-12-17 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> svglite </td>
   <td style="text-align:center;"> svglite </td>
   <td style="text-align:center;"> 2.2.1 </td>
   <td style="text-align:center;"> 2.2.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/svglite </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/svglite </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-05-12 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> systemfonts </td>
   <td style="text-align:center;"> systemfonts </td>
   <td style="text-align:center;"> 1.2.3 </td>
   <td style="text-align:center;"> 1.2.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/systemfonts </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/systemfonts </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-04-30 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> tensor </td>
   <td style="text-align:center;"> tensor </td>
   <td style="text-align:center;"> 1.5.1 </td>
   <td style="text-align:center;"> 1.5.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/tensor </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/tensor </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-06-17 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> textshaping </td>
   <td style="text-align:center;"> textshaping </td>
   <td style="text-align:center;"> 1.0.1 </td>
   <td style="text-align:center;"> 1.0.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/textshaping </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/textshaping </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-05-01 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> tibble </td>
   <td style="text-align:center;"> tibble </td>
   <td style="text-align:center;"> 3.3.0 </td>
   <td style="text-align:center;"> 3.3.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/tibble </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/tibble </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-06-08 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> tidyr </td>
   <td style="text-align:center;"> tidyr </td>
   <td style="text-align:center;"> 1.3.1 </td>
   <td style="text-align:center;"> 1.3.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/tidyr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/tidyr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-01-24 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> tidyselect </td>
   <td style="text-align:center;"> tidyselect </td>
   <td style="text-align:center;"> 1.2.1 </td>
   <td style="text-align:center;"> 1.2.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/tidyselect </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/tidyselect </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-03-11 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> tidytree </td>
   <td style="text-align:center;"> tidytree </td>
   <td style="text-align:center;"> 0.4.6 </td>
   <td style="text-align:center;"> 0.4.6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/tidytree </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/tidytree </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-12-12 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> treeio </td>
   <td style="text-align:center;"> treeio </td>
   <td style="text-align:center;"> 1.32.0 </td>
   <td style="text-align:center;"> 1.32.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/treeio </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/treeio </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-04-15 </td>
   <td style="text-align:center;"> Bioconductor 3.21 (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UCSC.utils </td>
   <td style="text-align:center;"> UCSC.utils </td>
   <td style="text-align:center;"> 1.4.0 </td>
   <td style="text-align:center;"> 1.4.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/UCSC.utils </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/UCSC.utils </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-04-17 </td>
   <td style="text-align:center;"> Bioconductor 3.21 (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> urlchecker </td>
   <td style="text-align:center;"> urlchecker </td>
   <td style="text-align:center;"> 1.0.1 </td>
   <td style="text-align:center;"> 1.0.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/urlchecker </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/urlchecker </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-11-30 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> usethis </td>
   <td style="text-align:center;"> usethis </td>
   <td style="text-align:center;"> 3.1.0 </td>
   <td style="text-align:center;"> 3.1.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/usethis </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/usethis </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-11-26 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> uwot </td>
   <td style="text-align:center;"> uwot </td>
   <td style="text-align:center;"> 0.2.3 </td>
   <td style="text-align:center;"> 0.2.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/uwot </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/uwot </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-02-24 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> vctrs </td>
   <td style="text-align:center;"> vctrs </td>
   <td style="text-align:center;"> 0.6.5 </td>
   <td style="text-align:center;"> 0.6.5 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/vctrs </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/vctrs </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-12-01 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> viridisLite </td>
   <td style="text-align:center;"> viridisLite </td>
   <td style="text-align:center;"> 0.4.2 </td>
   <td style="text-align:center;"> 0.4.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/viridisLite </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/viridisLite </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-05-02 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> withr </td>
   <td style="text-align:center;"> withr </td>
   <td style="text-align:center;"> 3.0.2 </td>
   <td style="text-align:center;"> 3.0.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/withr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/withr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-10-28 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> xfun </td>
   <td style="text-align:center;"> xfun </td>
   <td style="text-align:center;"> 0.52 </td>
   <td style="text-align:center;"> 0.52 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/xfun </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/xfun </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-04-02 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> xml2 </td>
   <td style="text-align:center;"> xml2 </td>
   <td style="text-align:center;"> 1.3.8 </td>
   <td style="text-align:center;"> 1.3.8 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/xml2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/xml2 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-03-14 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> xtable </td>
   <td style="text-align:center;"> xtable </td>
   <td style="text-align:center;"> 1.8.4 </td>
   <td style="text-align:center;"> 1.8-4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/xtable </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/xtable </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-04-21 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> XVector </td>
   <td style="text-align:center;"> XVector </td>
   <td style="text-align:center;"> 0.48.0 </td>
   <td style="text-align:center;"> 0.48.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/XVector </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/XVector </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-04-15 </td>
   <td style="text-align:center;"> Bioconductor 3.21 (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> yaml </td>
   <td style="text-align:center;"> yaml </td>
   <td style="text-align:center;"> 2.3.10 </td>
   <td style="text-align:center;"> 2.3.10 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/yaml </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/yaml </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-26 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> yulab.utils </td>
   <td style="text-align:center;"> yulab.utils </td>
   <td style="text-align:center;"> 0.2.0 </td>
   <td style="text-align:center;"> 0.2.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/yulab.utils </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/yulab.utils </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-01-29 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> zip </td>
   <td style="text-align:center;"> zip </td>
   <td style="text-align:center;"> 2.3.3 </td>
   <td style="text-align:center;"> 2.3.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/zip </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/zip </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-05-13 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> zoo </td>
   <td style="text-align:center;"> zoo </td>
   <td style="text-align:center;"> 1.8.14 </td>
   <td style="text-align:center;"> 1.8-14 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/zoo </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/zoo </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2025-04-10 </td>
   <td style="text-align:center;"> CRAN (R 4.5.0) </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library </td>
  </tr>
</tbody>
</table>
