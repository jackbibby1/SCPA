#' Create pathway matrices from gene sets
#'
#' This function takes a matrix and pathway csv/gmt files
#' and input and creates expression matrices for each
#' pathway. The resulting output is a nested list of expression
#' matrices for each pathway by each sample.
#'
#' @param samples Expression matrix for each population, with each population
#'   as a separate list. These matrices can be generated with e.g. `SCPA::seurat_extract()`
#' @param pathways Either a gene set gmt file, csv file, or output from `SCPA::format_pathways()`
#'   that contains your gene sets of interest. Information on formatting gene sets can be found
#'   at https://jackbibby1.github.io/SCPA/articles/using_gene_sets.html.
#' @param sample_names Names of the samples to be used in the output
#'
#'@examples \dontrun{
#' pathway_expression <- pathway_matrices(
#'      list(sample1, sample2),
#'      pathways = "your/path_to_gmt",
#'      sample_names = c("population1", "population2"))
#' }
#'
#' @export
#'
#'

pathway_matrices <- function(samples,
                             pathways,
                             sample_names = NULL) {

  if (class(pathways)[1] == "character") {
    pathways <- get_paths(pathways)
  }

  pop_paths <- vector(mode = "list", length = length(samples))
  for (i in 1:length(pop_paths)) {
    for (c in 1:length(pathways)) {
      pop_paths[[i]][[c]] <- samples[[i]][rownames(samples[[i]]) %in% pathways[[c]]$Genes, ]
    }
  }

  if (!is.null(sample_names)) {
    names(pop_paths) <- sample_names
  }

  path_names <- sapply(pathways, function(x) unique(x$Pathway))
  pop_paths <- lapply(pop_paths, function(x) magrittr::set_names(x, path_names))

  return(pop_paths)
}





