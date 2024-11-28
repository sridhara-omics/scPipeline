#' Preprocess count data and create a Seurat object
#'
#' This function preprocesses count data, optionally applying batch correction using batchelor::fastMNN,
#' and creates a Seurat object. It reports errors or warnings at each step for better debugging.
#'
#' @export
#' @importFrom Seurat CreateSeuratObject
#' @importFrom Seurat NormalizeData
#' @importFrom Seurat FindVariableFeatures
#' @importFrom Seurat PercentageFeatureSet
#' @importFrom Seurat ScaleData
#' @importFrom Seurat as.Seurat
#' @importFrom batchelor fastMNN
#' @param counts_data A matrix or data frame of count data.
#' @param batch_column A vector or factor specifying batch assignments for each cell. Default is NULL.
#' @param use_fastMNN Logical. Whether to apply batch correction using fastMNN. Default is FALSE.
#' @param ... Additional arguments to be passed to Seurat::CreateSeuratObject.
#' @return A Seurat object.
SeuratPreprocess <- function(counts_data, batch_column = NULL, use_fastMNN = FALSE, ...) {
  n_size <- dim(counts_data)[2]
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

  # Step 1: Create Seurat object
  seurat_object <- log_step("CreateSeuratObject", {
    Seurat::CreateSeuratObject(
      counts = counts_data,
      project = "project_title",
      min.cells = 3,
      min.features = 200,
      ...
    )
  })

  # Step 2: Normalize data
  seurat_object <- log_step("NormalizeData", {
    NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = n_size)
  })

  # Step 3: Find variable features
  seurat_object <- log_step("FindVariableFeatures", {
    FindVariableFeatures(seurat_object)
  })

  # Step 4: Calculate mitochondrial percentage
  seurat_object <- log_step("PercentageFeatureSet", {
    seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^mt-")
    seurat_object
  })

  # Step 5: Apply batch correction (if requested)
  if (use_fastMNN) {
    seurat_object <- log_step("BatchCorrection", {
      if (is.null(batch_column)) {
        stop("Batch correction requires a 'batch_column' to specify batch assignments.")
      }

      # Add batch information to metadata
      seurat_object$batch <- batch_column

      # Split object by batch
      batches <- SplitObject(seurat_object, split.by = "batch")

      # Convert to SingleCellExperiment objects
      sce_list <- lapply(batches, Seurat::as.SingleCellExperiment)

      # Perform batch correction using fastMNN
      corrected <- batchelor::fastMNN(sce_list)

      # Convert back to Seurat object
      Seurat::as.Seurat(corrected)
    })
  }

  # Step 6: Scale data, regressing out mitochondrial percentage
  seurat_object <- log_step("ScaleData", {
    ScaleData(seurat_object, vars.to.regress = "percent.mt")
  })

  # Print log messages for user
  if (length(log_messages) > 0) {
    message("Summary of issues during execution:")
    message(paste(log_messages, collapse = "\n"))
  }

  return(seurat_object)
}
