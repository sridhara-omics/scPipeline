library(dplyr)
library(Seurat)
library(patchwork)


gse_counts <- read.table("C:/Users/vishs/Documents/R_projects/healthy_asthma_dropseq/h5ad/GSE130148_raw_counts.csv",sep=",",header=T,row.names=1)
gse_metadata <- read.table("C:/Users/vishs/Documents/R_projects/healthy_asthma_dropseq/h5ad/GSE130148_barcodes_cell_types.txt",sep="\t",header=T,row.names=1)
n_size <- dim(gse_counts)[2]
n_pcs <- round(n_size / 100)
n_dims <- round(n_pcs / 1.25)

seurat_object <- Seurat::CreateSeuratObject(
  counts = gse_counts,
  project = "project_title",
  min.cells = 3,
  min.features = 200,
  meta.data = gse_metadata
)

# Step 2: Normalize data
seurat_object <-
  NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = n_size)


# Step 3: Find variable features
seurat_object <-
  FindVariableFeatures(seurat_object)


# Step 4: Calculate mitochondrial percentage
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^mt-")

  # Step 5: Scale data
  seurat_object <- ScaleData(seurat_object, vars.to.regress = "percent.mt")

  # Step 6: Run PCA
  seurat_object <- RunPCA(seurat_object, npcs = n_pcs, ndims.print = 1:5, nfeatures.print = 5)

  # Step 7: Find neighbors for clustering
  seurat_object <- FindNeighbors(seurat_object, reduction = "pca", dims = 1:n_dims, nn.eps = 0.5)

  # Step 8: Find clusters
  seurat_object <- FindClusters(seurat_object, resolution = 3, n.start = 10)

  # Step 9:Run TSNE
  seurat_object <- RunTSNE(seurat_object, dims = 1:n_dims)

  # Step 10:Run UMAP
  seurat_object <- RunUMAP(seurat_object, dims = 1:n_dims, min.dist = 0.75)

  # step 11: UMAP plot by split.by
  UMAPPlot(seurat_object, split.by='ID')

  # step 12: Find All markers
  FindAllMarkers(seurat_object, min.pct = 0.3)

  # step 13: Reactome GSA analyse_sc_clusters
  gsva_result <- ReactomeGSA::analyse_sc_clusters(seurat_object)

  # step 14: Reactome GSA pathways
  pathway_expression <- ReactomeGSA::pathways(gsva_result)

  # step 15: Update pathway expression column names
  colnames(pathway_expression) <- gsub("\\.Seurat", "", colnames(pathway_expression))

  # Step 16: Calculate max differences in pathway expression
  max_difference <- do.call(rbind, apply(pathway_expression, 1, function(row) {
    values <- as.numeric(row[2:length(row)])
    return(data.frame(name = row[1], min = min(values), max = max(values)))
  }))

  # Step 17: Add difference column and sort
  max_difference$diff <- max_difference$max - max_difference$min
  max_difference <- max_difference[order(max_difference$diff, decreasing = TRUE), ]









