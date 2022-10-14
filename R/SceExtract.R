#' Extract data from a SingleCellExperiment object
#'
#' This function takes a SingleCellExperiment object as an input, and returns
#' an expression matrix based on subsetting parameters. Either none, one, or
#' two metadata features can be selected for a given input.
#'
#' @param sce_object SingleCellExperiment object containing expression data
#' @param assay_name Name of assay to pull from. Defaults to "logcounts"
#' @param meta1 Metadata column to subset
#' @param value_meta1 Value to select within `meta1` column
#' @param meta2 Metadata column to subset
#' @param value_meta2 Value to select within `meta2` column
#' @param pseudocount Pseudocount to add to data. Defaults to 0.001
#'
#' @examples \dontrun{
#' cd4 <- sce_extract(
#'    sce_object,
#'    meta1 = "Hour",
#'    value_meta1 = 12,
#'    meta2 = "Cell_Type",
#'    value_meta2 = "CD4"
#' )
#' }
#'
#' @return Matrix containing count values of selected populations
#' @export

sce_extract <- function(sce_object,
                        assay_name = "logcounts",
                        meta1 = NULL, value_meta1 = NULL,
                        meta2 = NULL, value_meta2 = NULL,
                        pseudocount = 0.001) {

  if (is.null(meta1) && is.null(meta2)) {
    message("No metadata selected. Converting whole SCE object to matrix")
    sub_sce <- assay(sce_object, assay_name)
    sub_sce <- as.matrix(sub_sce) + pseudocount
    return(sub_sce)
  }

  if (!is.null(meta1) && is.null(meta2)) {
    message(paste0("Extracting cells where ", meta1, " == ", value_meta1))
    sub_sce <- sce_object[, sce_object[[meta1]] == value_meta1]
    sub_sce <- assay(sub_sce, assay_name)
    sub_sce <- as.matrix(sub_sce) + pseudocount
    return(sub_sce)
  }

  if(!is.null(meta1) && !is.null(meta2)) {
    message(paste0("Extracting cells where ", meta1, " == ", value_meta1,
                   " AND ", meta2, " == ", value_meta2))
    sub_sce <- sce_object[, sce_object[[meta1]] == value_meta1 &
                            sce_object[[meta2]] == value_meta2]
    sub_sce <- assay(sub_sce, assay_name)
    sub_sce <- as.matrix(sub_sce) + pseudocount
    return(sub_sce)
  }
}




