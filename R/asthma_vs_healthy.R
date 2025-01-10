
library(dplyr)
library(Seurat)
library(celldex)
library(SingleR)
library(patchwork)

# Load dropseq files
gse_counts <- read.table("C:/Users/vishs/Documents/R_projects/healthy_asthma_dropseq/h5ad/GSE130148_raw_counts.csv",sep=",",header=T,row.names=1)
gse_metadata <- read.table("C:/Users/vishs/Documents/R_projects/healthy_asthma_dropseq/h5ad/GSE130148_barcodes_cell_types.txt",sep="\t",header=T,row.names=1)
n_size <- dim(gse_counts)[2]
n_pcs <- round(n_size / 100)
n_dims <- round(n_pcs / 1.25)

# Step 1: Create Seurat object with dropseq counts
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

  
  # step 20: load celldex and singleR to run below:
  hpca.se <- celldex::HumanPrimaryCellAtlasData()
  
  # step 21: singleR annotation on summarized object
  sce <- Seurat::as.SingleCellExperiment(seurat_object)
  sce.pred <- SingleR(test = sce, ref = hpca.se, assay.type.test=1, labels = hpca.se$label.main)
  
  seurat_object@meta.data$SingleR_labels <- sce.pred$labels
  
  # step 22: Transfer annotations from both manually curated original annotations as well as SingleR annotations to map to seurat clusters
  seurat_object <- TransferAnnotations(seurat_object, "celltype", "seurat_clusters", output_col = "celltype_cluster_annotation")
  seurat_object <- TransferAnnotations(seurat_object, "SingleR_labels", "seurat_clusters", output_col = "SingleR_cluster_annotation")
  
  
  DimPlot(
         object = seurat_object,               # Your updated Seurat object
         #group.by = "seurat_clusters",
         group.by = "celltype_cluster_annotation",       # Column to use for labeling
         label = TRUE,                          # Display labels on the plot
         label.size = 4                         # Adjust label size
     ) + NoLegend() 

  
  # step 12: Find markers
  markers <- FindMarkers(seurat_object, ident.1 = "25", min.pct = 0.3)
  markers[1:5, ]
  #markers <- FindMarkers(seurat_object, ident.1 = "Dropseq_2", group.by = 'ID', subset.ident = "Ciliated")
  #FindAllMarkers(seurat_object, min.pct = 0.3)
  
  # subset seurat object
  #seurat_subset <- subset(seurat_object, idents = c("25", "26"))
  
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
  
  # step 18: Features for visualization
  features <- c("SCGB1A1", "LYZ", "ATP2B1", "IL4", "IL5", "IL13", "CCL21")
  #("AR", "NOS2")
  #, SCGB1A1", "IL1RL1", "POSTN", "IL12A")
  # IL1RL1, POSTN, SERPINB2, CLCA1, NOS2, and MUC5AC
  # IL12A and MUC5B
  
  # step 19: plots
  DotPlot(seurat_object, features = features, group.by = "celltype_cluster_annotation") + RotatedAxis()
  DotPlot(seurat_object, features = features, group.by = "SingleR_cluster_annotation") + RotatedAxis()
  
  VlnPlot(seurat_object, features = features, ncol = 2)
  FeaturePlot(seurat_object, features = features)
  DotPlot(seurat_object, features = features) + RotatedAxis()
  DoHeatmap(subset(seurat_object, downsample = 100), features = features, size = 3)
  
  VlnPlot(seurat_object, features = "percent.mt")
  
  DotPlot(seurat_object, features = features, group.by = "celltype_cluster_annotation") + RotatedAxis()
  
  
  



