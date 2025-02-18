#' Annotate cells in a Seurat object using SingleR
#'
#' This function annotates the cells in a Seurat object using the SingleR package.
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
  # Step 1: Load default reference data if not provided
  if (is.null(reference_data)) {
    reference_data <- SingleR::HumanPrimaryCellAtlasData()
  }

  # Step 2: Convert Seurat object to SingleCellExperiment
  sce <- Seurat::as.SingleCellExperiment(seurat_object, assay = assay)

  # Step 3: Run SingleR to annotate cells
  singleR_results <- SingleR::SingleR(
    test = SummarizedExperiment::SummarizedExperiment(list(counts = sce)),
    ref = reference_data,
    labels = reference_data$label.main
  )

  # Step 4: Add SingleR annotations to Seurat metadata
  seurat_object <- Seurat::AddMetaData(
    object = seurat_object,
    metadata = singleR_results$labels,
    col.name = "SingleR_Labels"
  )

  return(seurat_object)
}
