library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(data.table)
library(SingleCellExperiment)
sce <- readRDS('scRNA_sce.RDS')
#assays(sce)
#dim(counts(sce))
#counts(sce)[1:6, 1:6]
#dim(colData(sce))
#head(colData(sce))
cluster_names <- levels(colData(sce)$cluster_id)
#cluster_names
#length(cluster_names)
sample_names <- levels(colData(sce)$donor_id)
#sample_names
#length(sample_names)
groups <- colData(sce)[, c("cluster_id", "donor_id")]
#head(groups)
aggr_counts <- aggregate.Matrix(t(counts(sce)), 
                                groupings = groups, fun = "sum") 

#class(aggr_counts)
#dim(aggr_counts)
#aggr_counts[1:6, 1:6]
aggr_counts <- t(aggr_counts)
tstrsplit(colnames(aggr_counts), "_") %>% str()
head(colnames(aggr_counts), n = 10)
head(tstrsplit(colnames(aggr_counts), "_")[[1]], n = 10)
b_cell_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]] == "B-cells")
b_cell_idx
colnames(aggr_counts)[b_cell_idx]
aggr_counts[1:10, b_cell_idx]
#cluster_names
counts_ls <- list()

for (i in 1:length(cluster_names)) {
  column_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]] == cluster_names[i])
  counts_ls[[i]] <- aggr_counts[, column_idx]
  names(counts_ls)[i] <- cluster_names[i]
}

str(counts_ls)
'''
#head(colData(sce))
sce$sample_id <- paste(sce$subtype, sce$donor_id, sep = "_")
metadata <- colData(sce) %>% 
  as.data.frame() %>% 
  dplyr::select(subtype, donor_id, sample_id)
head(metadata)
metadata <- metadata[!duplicated(metadata), ]
dim(metadata)
head(metadata)
rownames(metadata) <- metadata$sample_id
head(metadata)
t <- table(colData(sce)$sample_id,
           colData(sce)$cluster_id)
t[1:6, 1:6]
# Creating metadata list

## Initiate empty list
metadata_ls <- list()

for (i in 1:length(counts_ls)) {
  
  ## Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
  df <- data.frame(cluster_sample_id = colnames(counts_ls[[i]]))
  df
  
  ## Use tstrsplit() to separate cluster (cell type) and sample IDs
  df$cluster_id <- tstrsplit(df$cluster_sample_id, "_")[[1]]
  df$sample_id  <- tstrsplit(df$cluster_sample_id, "_")[[2]]
  
  
  ## Retrieve cell count information for this cluster from global cell count table
  idx <- which(colnames(t) == unique(df$cluster_id))
  idx
  cell_counts <- t[, idx]
  
  ## Remove samples with zero cell contributing to the cluster
  cell_counts <- cell_counts[cell_counts > 0]
  
  ## Match order of cell_counts and sample_ids
  sample_order <- match(df$sample_id, names(cell_counts))
  cell_counts <- cell_counts[sample_order]
  
  ## Append cell_counts to data frame
  df$cell_count <- cell_counts
  
  
  ## Join data frame (capturing metadata specific to cluster) to generic metadata
  df <- plyr::join(df, metadata, 
                   by = intersect(names(df), names(metadata)))
  
  ## Update rownames of metadata to match colnames of count matrix, as needed later for DE
  rownames(df) <- df$cluster_sample_id
  
  ## Store complete metadata for cluster i in list
  metadata_ls[[i]] <- df
  names(metadata_ls)[i] <- unique(df$cluster_id)
  
}



str(metadata_ls)
# Select cell type of interest
cluster_names

# Double-check that both lists have same names
all(names(counts_ls) == names(metadata_ls))
idx <- which(names(counts_ls) == "B-cells")
cluster_counts <- counts_ls[[idx]]
cluster_metadata <- metadata_ls[[idx]]
# Check contents of extracted objects
cluster_counts[1:6, 1:6]
head(cluster_metadata)

# Check matching of matrix columns and metadata rows
all(colnames(cluster_counts) == rownames(cluster_metadata))
dds <- DESeqDataSetFromMatrix(cluster_counts, 
                              colData = cluster_metadata, 
                              design = ~ group_id)
rld <- rlog(dds, blind=TRUE)

# Plot PCA
DESeq2::plotPCA(rld, ntop = 500, intgroup = "group_id")
DESeq2::plotPCA(rld, ntop = 500, intgroup = "cell_count")
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

# Plot heatmap
pheatmap(rld_cor, annotation = cluster_metadata[, c("group_id"), drop=F])
dds <- DESeq(dds)
plotDispEsts(dds)
resultsNames(dds)

# Generate results object
res <- results(dds, 
               name = "ER+ vs TNBC",
               alpha = 0.05)

# Shrink the log2 fold changes to be more appropriate using the apeglm method - should cite [paper]() when using this method
res <- lfcShrink(dds, 
                 coef = "ER+ vs TNBC",
                 res=res,
                 type = "apeglm")
# Turn the DESeq2 results object into a tibble for use with tidyverse functions
res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  arrange(padj)

# Check results output
res_tbl 

# Write all results to file
write.csv(res_tbl,
          paste0("results/", unique(cluster_metadata$cluster_id), "_", 
                 levels(cluster_metadata$group_id)[2], "_vs_", levels(cluster_metadata$group_id)[1], "_all_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)
'''