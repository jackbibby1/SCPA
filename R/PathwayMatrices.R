#' Create pathway matrices from gene sets
#'
#' This function takes a matrix and pathway list as
#' and input and creates expression matrices for each
#' pathway. The resulting output is a list of expression
#' matrices for each pathway
#'
#' @param samples List containing expression file of samples
#' @param pathways List of pathways to create pathway specific
#'   matrices. You can use the get_paths function to create
#'   this list if you are using a csv file generated from a gmt
#'   file.
#' @param sample_names Names of the samples to be used in output
#' @param pathway_names Names of pathways to be used in output
#'
#' @export
#'
#'

pathway_matrices <- function(samples,
                             pathways,
                             sample_names = NULL,
                             pathway_names = NULL) {

  pathways <- get_paths(pathways)

  pop_paths <- vector(mode = "list", length = length(samples))
  for (i in 1:length(pop_paths)) {
    for (c in 1:length(pathways)) {
      pop_paths[[i]][[c]] <- samples[[i]][rownames(samples[[i]]) %in% pathways[[c]]$Genes, ]
    }
  }

  if (!is.null(sample_names)) {
    names(pop_paths) <- sample_names
  }

  if (!is.null(pathway_names)) {
    path_names <- sapply(pathways, function(x) unique(x$Pathway))
  }

  return(pop_paths)
}

