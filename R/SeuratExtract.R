#' Extract Normalised Data From Seurat
#'
#' This function takes a Seurat object as an input, and returns
#' an expression matrix based on subsetting parameters. Either none, one, or
#' two metadata features can be selected for a given input.
#'
#' @param seu_obj Seurat object containing normalised counts
#'   stored in `seu_obj@assays$RNA@data`
#' @param meta1 Metadata column to subset
#' @param meta2 Metadata column to subset
#' @param value_meta1 Value to select within `meta1` column
#' @param value_meta2 Value to select within `meta2` column
#' @param pseudocount Pseudocount to add to data. Defaults to 0.001
#'
#' @examples \dontrun{
#' cd4 <- seurat_extract(
#'    Seurat_Object,
#'    meta1 = "Hour",
#'    value_meta1 = 12,
#'    meta2 = "Cell_Type",
#'    value_meta2 = "CD4"
#' )
#' }
#'
#' @return Matrix containing count values of selected populations
#' @export

seurat_extract <- function(seu_obj,
                           meta1 = NULL, value_meta1 = NULL,
                           meta2 = NULL, value_meta2 = NULL,
                           pseudocount = 0.001) {

  if (is.null(meta1) && is.null(meta2)) {
    message("No metadata selected. Converting whole Seurat object to matrix")
    seu_obj <- as.matrix(seu_obj@assays$RNA@data) + pseudocount
    return(seu_obj)
  }

  if (!is.null(meta1) && is.null(meta2)) {
    message(paste0("Extracting cells where ", meta1, " == ", value_meta1))
    met_1 <- Seurat::FetchData(object = seu_obj, vars = meta1)
    seu_obj <- seu_obj[, which(met_1 == value_meta1)]
    seu_obj <- as.matrix(seu_obj@assays$RNA@data) + pseudocount
    return(seu_obj)
  }

  if(!is.null(meta1) && !is.null(meta2)) {
    message(paste0("Extracting cells where ", meta1, " == ", value_meta1,
                   " AND ", meta2, " == ", value_meta2))
    met_1 <- Seurat::FetchData(object = seu_obj, vars = meta1)
    met_2 <- Seurat::FetchData(object = seu_obj, vars = meta2)
    seu_obj <- seu_obj[, which(met_1 == value_meta1 & met_2 == value_meta2)]
    seu_obj <- as.matrix(seu_obj@assays$RNA@data) + pseudocount
    return(seu_obj)
  }
}
