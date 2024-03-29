library(Seurat)
library(tidyverse)
#setwd('/root/user/salma_tariq/Single_cell')
#rds_obj <- readRDS('913238fa-494f-4bfc-8a40-fb3005df448e.rds')
rds_obj <- readRDS('filt_rds_obj.rds')
cts <- rds_obj$RNA@data
rds_obj.seurat.obj <- CreateSeuratObject(counts = cts, project = 'BC', min.cells = 3, min.features = 200)
#str(rds_obj.seurat.obj)
#rds_obj.seurat.obj@meta.data
rds_obj.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(rds_obj.seurat.obj, pattern = "^MT-")
pdf(file = "SC_plots.pdf", width = 8, height = 6)
VlnPlot(rds_obj.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#dev.off()
#pdf(file = "SC_plots_2.pdf", width = 8, height = 6)
FeatureScatter(rds_obj.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')
#dev.off()
rds_obj.seurat.obj <- subset(rds_obj.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & 
                               percent.mt < 5)
rds_obj.seurat.obj <- NormalizeData(rds_obj.seurat.obj)
str(rds_obj.seurat.obj)
rds_obj.seurat.obj <- FindVariableFeatures(rds_obj.seurat.obj, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(rds_obj.seurat.obj), 10)

# plot variable features with and without labels
#pdf(file = "SC_plots_3.pdf", width = 8, height = 6)
plot1 <- VariableFeaturePlot(rds_obj.seurat.obj)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

# 5. Scaling -------------
#all.genes <- rownames(rds_obj.seurat.obj)
#rds_obj.seurat.obj <- ScaleData(rds_obj.seurat.obj, features = all.genes)
rds_obj.seurat.obj <- ScaleData(rds_obj.seurat.obj)
#Error: I got an error here due to huge size of the seurat object so i decided to subset the object 
#rds_obj@meta.data$subtype
#filt_rds_obj <- subset(rds_obj, subset = subtype %in% c("ER+", "TNBC"))
#subtype <- seurat_obj@meta.data$subtype
#unique(subtype)
#colnames(rds_obj@meta.data)
# Save the filtered Seurat object to a new RDS file
#saveRDS(filt_rds_obj, file = "filt_rds_obj.rds")
#str(rds_obj.seurat.obj)

# 6. Perform Linear dimensionality reduction --------------
rds_obj.seurat.obj<- RunPCA(rds_obj.seurat.obj, features = VariableFeatures(object = rds_obj.seurat.obj))

# visualize PCA results
print(rds_obj.seurat.obj[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(rds_obj.seurat.obj, dims = 1, cells = 500, balanced = TRUE)
# determine dimensionality of the data
ElbowPlot(rds_obj.seurat.obj)


# 7. Clustering ------------
rds_obj.seurat.obj <- FindNeighbors(rds_obj.seurat.obj, dims = 1:15)

# understanding resolution
rds_obj.seurat.obj <- FindClusters(rds_obj.seurat.obj, resolution = c(0.1,0.3, 0.5, 0.7, 1))
rds_obj.seurat.obj@meta.data

DimPlot(rds_obj.seurat.obj, group.by = "RNA_snn_res.0.5", label = TRUE)

# setting identity of clusters
Idents(rds_obj.seurat.obj)
Idents(rds_obj.seurat.obj) <- "RNA_snn_res.0.1"
Idents(rds_obj.seurat.obj)

# non-linear dimensionality reduction --------------
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
rds_obj.seurat.obj <- RunUMAP(rds_obj.seurat.obj, dims = 1:15)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(rds_obj.seurat.obj, reduction = "umap")
dev.off()