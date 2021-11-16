#' Format the output of msigdbr to use with SCPA
#'
#' This function takes the output of msigdbr
#' and formats the pathways into something that
#' can be used within SCPA
#'
#' @param msigdbr_output output from msigdbr
#'
#' @examples \dontrun{
#' pathways <- format_pathways(
#'      pathways
#' )
#' }
#'
#' @return list of pathways with corresponding genes
#' @export

format_pathways <- function(msigdbr_output) {
  msigdbr_output <- msigdbr_output[, c("gs_name", "gene_symbol")]
  colnames(msigdbr_output) <- c("Pathway", "Genes")
  msigdbr_output <- dplyr::group_split(msigdbr_output, "Pathway")
  return(msigdbr_output)
}
