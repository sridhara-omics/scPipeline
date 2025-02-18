#' A thresholded markers list for better calculation of DE genes
#'
#' This function calculates differentially expressed genes using Seurat::FindAllMarkers.
#'
#' @param lowdim_seurat_object Seurat object with cluster information
#' @importFrom Seurat FindAllMarkers
#' @return A list containing two marker lists:
#'         - Full markers list
#'         - Thresholded markers list with `min.pct = 0.1`
#' @export
SeuratMarkers <- function(lowdim_seurat_object) {
  # Calculate the full markers list
  markers_list <- Seurat::FindAllMarkers(lowdim_seurat_object)

  # Calculate the thresholded markers list
  markers_list_threshold <- Seurat::FindAllMarkers(lowdim_seurat_object, min.pct = 0.1)

  # Return results as a list
  return(list(
    full_markers_list = markers_list,
    thresholded_markers_list = markers_list_threshold
  ))
}
