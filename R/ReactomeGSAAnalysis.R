#' Reactome Data Analysis for Seurat Object
#'
#' This function performs pathway analysis using ReactomeGSA on a Seurat object
#' with cluster information and logs any errors or warnings encountered.
#'
#' @param lowdim_seurat_object Seurat object that has clusters information
#' @return A list containing:
#'         - GSVA result (`gsva_result`)
#'         - Pathway expression data (`pathway_expression`)
#'         - Max difference between pathway expression values (`max_difference`)
#' @export
#' @importFrom ReactomeGSA analyse_sc_clusters pathways
ReactomeData <- function(lowdim_seurat_object) {
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

  # Step 1: Perform GSVA pathway analysis
  gsva_result <- log_step("analyse_sc_clusters", {
    ReactomeGSA::analyse_sc_clusters(lowdim_seurat_object)
  })

  # Step 2: Extract pathway expression data
  pathway_expression <- log_step("pathways", {
    ReactomeGSA::pathways(gsva_result)
  })

  # Step 3: Clean pathway expression column names
  pathway_expression <- log_step("CleanColumnNames", {
    colnames(pathway_expression) <- gsub("\\.Seurat", "", colnames(pathway_expression))
    pathway_expression
  })

  # Step 4: Calculate max differences in pathway expression
  max_difference <- log_step("CalculateMaxDifferences", {
    do.call(rbind, apply(pathway_expression, 1, function(row) {
      values <- as.numeric(row[2:length(row)])
      return(data.frame(name = row[1], min = min(values), max = max(values)))
    }))
  })

  # Step 5: Add difference column and sort
  max_difference <- log_step("SortMaxDifferences", {
    max_difference$diff <- max_difference$max - max_difference$min
    max_difference[order(max_difference$diff, decreasing = TRUE), ]
  })

  # Print log messages for user
  if (length(log_messages) > 0) {
    message("Summary of issues during execution:")
    message(paste(log_messages, collapse = "\n"))
  }

  # Return the results
  return(list(gsva_result = gsva_result, pathway_expression = pathway_expression, max_difference = max_difference))
}
