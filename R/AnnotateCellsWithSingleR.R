#' Annotate cells in a Seurat object using SingleR
#'
#' This function annotates the cells in a Seurat object using the SingleR package
#' and provides detailed error/warning logs for each step.
#'
#' @export
#' @importFrom SingleR SingleR
#' @importFrom SingleR HumanPrimaryCellAtlasData
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom Seurat as.SingleCellExperiment AddMetaData
#' @param seurat_object A Seurat object to be annotated.
#' @param reference_data A reference dataset to use for annotation (e.g., HumanPrimaryCellAtlasData).
#'        If NULL, HumanPrimaryCellAtlasData is used by default.
#' @param assay The assay in the Seurat object to use for annotation. Default is "RNA".
#' @return The Seurat object with cell annotations added to the metadata.
AnnotateCellsWithSingleR <- function(seurat_object, reference_data = NULL, assay = "RNA") {
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

  # Step 1: Load default reference data if not provided
  reference_data <- log_step("LoadReferenceData", {
    if (is.null(reference_data)) {
      SingleR::HumanPrimaryCellAtlasData()
    } else {
      reference_data
    }
  })

  # Step 2: Convert Seurat object to SingleCellExperiment
  sce <- log_step("ConvertToSingleCellExperiment", {
    Seurat::as.SingleCellExperiment(seurat_object, assay = assay)
  })

  # Step 3: Run SingleR to annotate cells
  singleR_results <- log_step("RunSingleR", {
    SingleR::SingleR(
      test = SummarizedExperiment::SummarizedExperiment(list(counts = sce)),
      ref = reference_data,
      labels = reference_data$label.main
    )
  })

  # Step 4: Add SingleR annotations to Seurat metadata
  seurat_object <- log_step("AddSingleRAnnotations", {
    Seurat::AddMetaData(
      object = seurat_object,
      metadata = singleR_results$labels,
      col.name = "SingleR_Labels"
    )
  })

  # Print log messages for user
  if (length(log_messages) > 0) {
    message("Summary of issues during execution:")
    message(paste(log_messages, collapse = "\n"))
  }

  return(seurat_object)
}
