#' Convert gene sets to filtered pathway list
#'
#' This function takes one of three inputs:
#' 1) csv file where pathway title is in the first column
#' and subsequent columns contain genes of that pathway
#' 2) A gmt file
#' 3) A tidy formatted list where pathway name is contained in the first column
#' and the genes of that pathway are in the second column
#'
#' @param pathway_filepath filepath to csv or gmt file
#' @param min_genes Minimum number of genes required in a pathway for inclusion
#' @param max_genes Maximum number of genes required in a pathway for inclusion
#'
#' @examples \dontrun{
#' pathways <- get_paths(
#'      "Documents/gene_list.csv"
#' )
#' }
#'
#' @return list of pathways with corresponding genes
#' @export


get_paths <- function(pathway_filepath) {

  if (str_ends(pathway_filepath, "gmt")) {
    pathways <- utils::read.delim(pathway_filepath, row.names = 1, header = F)
    pathways <- as.data.frame(t(pathways))
    pathways <- tidyr::pivot_longer(pathways, cols = 1:length(pathways), names_to = "Pathway", values_to = "Genes")
    pathways <- dplyr::group_split(pathways, Pathway)
    pathways <- lapply(pathways, function(x) x[x$Genes != "",])
    return(pathways)

    } else if (str_ends(pathway_filepath, "csv")) {

      pathways <- utils::read.csv(pathway_filepath, row.names = 1, header = F)
      pathways <- as.data.frame(t(pathways))
      pathways <- tidyr::pivot_longer(pathways, cols = 1:length(pathways), names_to = "Pathway", values_to = "Genes")
      pathways <- dplyr::group_split(pathways, Pathway)
      pathways <- lapply(pathways, function(x) x[x$Genes != "",])
      return(pathways)

    } else {
      pathways <- pathway_filepath
      return(pathways)
    }
}




