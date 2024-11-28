#' Create a Low dimensional Seurat object from scaled seurat object
#'
#' This function converts the transformed data to low-dimensional data for downstream analysis
#' and provides detailed error/warning logs for each step.
#'
#' @export
#' @importFrom Seurat RunPCA
#' @importFrom Seurat FindNeighbors
#' @importFrom Seurat FindClusters
#' @importFrom Seurat RunTSNE
#' @importFrom Seurat RunUMAP
#' @param scaled_seurat_object A scaled Seurat object.
#' @param ... Additional arguments to be passed for downstream analyses.
#' @return A Seurat object.
SeuratLowDim <- function(scaled_seurat_object, ...) {
  log_messages <- c() # Initialize log for errors and warnings

  # Helper function to log messages
  log_step <- function(step_name, expr) {
    tryCatch(
      expr,
      error = function(e) {
        log_messages <<- c(log_messages, paste0("Error in ", step_name, ": ", e$message))
        stop(paste0("Step '", step_name, "' failed: ", e$message))
      },
      warning = function(w) {
        log_messages <<- c(log_messages, paste0("Warning in ", step_name, ": ", w$message))
        warning(paste0("Step '", step_name, "' raised a warning: ", w$message))
      }
    )
  }

  # Step 1: Compute number of PCs and dimensions
  n_size <- log_step("Compute n_size", {
    dim(scaled_seurat_object@meta.data)[1]
  })

  n_pcs <- log_step("Compute n_pcs", {
    round(n_size / 100)
  })

  n_dims <- log_step("Compute n_dims", {
    round(n_pcs / 1.25)
  })

  # Step 2: Run PCA
  seurat_object <- log_step("RunPCA", {
    RunPCA(scaled_seurat_object, npcs = n_pcs, ndims.print = 1:5, nfeatures.print = 5)
  })

  # Step 3: Find Neighbors
  seurat_object <- log_step("FindNeighbors", {
    FindNeighbors(seurat_object, reduction = "pca", dims = 1:n_dims, nn.eps = 0.5)
  })

  # Step 4: Find Clusters
  seurat_object <- log_step("FindClusters", {
    FindClusters(seurat_object, resolution = 3, n.start = 10)
  })

  # Step 5: Run t-SNE
  seurat_object <- log_step("RunTSNE", {
    RunTSNE(seurat_object, dims = 1:n_dims)
  })

  # Step 6: Run UMAP
  seurat_object <- log_step("RunUMAP", {
    RunUMAP(seurat_object, dims = 1:n_dims, min.dist = 0.75)
  })

  # Print log messages for user
  if (length(log_messages) > 0) {
    message("Summary of issues during execution:")
    message(paste(log_messages, collapse = "\n"))
  }

  return(seurat_object)
}
