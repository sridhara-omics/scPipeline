#' A thresholded markers list for better calculation of DE genes
#'
#' This function calculates differentially expressed genes using Seurat::FindAllMarkers
#' and logs any errors or warnings encountered during execution.
#'
#' @param lowdim_seurat_object Seurat object with cluster information
#' @importFrom Seurat FindAllMarkers
#' @return A list containing two marker lists:
#'         - Full markers list
#'         - Thresholded markers list with `min.pct = 0.1`
#' @export
SeuratMarkers <- function(lowdim_seurat_object) {
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

  # Step 1: Calculate the full markers list
  markers_list <- log_step("FindAllMarkers (full list)", {
    Seurat::FindAllMarkers(lowdim_seurat_object)
  })

  # Step 2: Calculate the thresholded markers list
  markers_list_threshold <- log_step("FindAllMarkers (thresholded)", {
    Seurat::FindAllMarkers(lowdim_seurat_object, min.pct = 0.1)
  })

  # Print log messages for user
  if (length(log_messages) > 0) {
    message("Summary of issues during execution:")
    message(paste(log_messages, collapse = "\n"))
  }

  # Return results as a list
  return(list(
    full_markers_list = markers_list,
    thresholded_markers_list = markers_list_threshold
  ))
}
