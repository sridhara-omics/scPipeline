#' Transfer annotations to Seurat clusters
#'
#' This function assigns cluster-level annotations in a Seurat object based on the majority 
#' annotation of cells within each cluster. It logs any errors or warnings encountered.
#'
#' @param seurat_object Seurat object containing cluster and annotation information.
#' @param annotation_col The name of the metadata column with annotations (character string).
#' @param cluster_col The name of the metadata column with cluster information (character string).
#' @importFrom dplyr group_by summarise
#' @importFrom Seurat DimPlot
#' @return The Seurat object with an additional column `cluster_annotation` in its metadata.
#' @export
TransferAnnotations <- function(seurat_object, annotation_col, cluster_col) {
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
  
  # Step 1: Validate input columns
  log_step("Validation of input columns", {
    if (!annotation_col %in% colnames(seurat_object@meta.data)) {
      stop(paste0("Annotation column '", annotation_col, "' not found in metadata."))
    }
    if (!cluster_col %in% colnames(seurat_object@meta.data)) {
      stop(paste0("Cluster column '", cluster_col, "' not found in metadata."))
    }
  })
  
  # Step 2: Calculate majority annotations
  majority_annotations <- log_step("Majority annotation calculation", {
    metadata <- seurat_object@meta.data
    metadata %>%
      group_by(!!sym(cluster_col)) %>%
      summarise(
        majority_annotation = names(sort(table(!!sym(annotation_col)), decreasing = TRUE)[1])
      )
  })
  
  # Step 3: Map annotations to all cells
  log_step("Mapping annotations to cells", {
    seurat_object@meta.data$cluster_annotation <- majority_annotations$majority_annotation[
      match(seurat_object@meta.data[[cluster_col]], majority_annotations[[cluster_col]])
    ]
  })
  
  # Print log messages for user
  if (length(log_messages) > 0) {
    message("Summary of issues during execution:")
    message(paste(log_messages, collapse = "\n"))
  }
  
  # Return the updated Seurat object
  return(seurat_object)
}
