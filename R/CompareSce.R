#' Use SCPA to compare pathways within a SingleCellExperiment object
#'
#' This function takes a SingleCellExperiment object as an input, and
#' compares gene sets over specified conditions/populations.
#'
#' @param sce_object SingleCellExperiment object with populations defined in the column data
#' @param assay_name Assay name to extract expression values from. Defaults to logcounts
#' @param group1 First comparison group as defined by `colData()` columns
#'   of SingleCellExperiment object e.g. cell_type
#' @param group1_population Populations within group1 to compare
#'   e.g. c("t_cell", "b_cell")
#' @param group2 Second comparison group as defined by `colData()` columns
#'   of SingleCellExperiment object e.g. hour
#' @param group2_population Population within group2 to compare
#'   e.g. 24
#' @param pathways Pathway gene sets with each pathway in a separate list. For formatting of
#'   gene lists, see documentation at https://jackbibby1.github.io/SCPA/articles/using_gene_sets.html
#' @param downsample Option to downsample cell numbers. Defaults to 500 cells per condition. If a population
#'   has < 500 cells, all cells from that condition are used.
#' @param min_genes Gene sets with fewer than this number of genes will be excluded
#' @param max_genes Gene sets with more than this number of genes will be excluded
#' @param parallel Should parallel processing be used?
#' @param cores The number of cores used for parallel processing
#'
#' @examples \dontrun{
#' scpa_out <- compare_sce(
#'      group1 = "cell",
#'      group1_population = c("t_cell", "b_cell"),
#'      group2 = "hour",
#'      group2_population = c("24"),
#'      pathways = pathways)
#' }
#'
#' @return Statistical results from the SCPA analysis. The qval should be the
#' primary metric that is used to interpret pathway differences i.e. a higher
#' qval translates to larger pathway differences between conditions.
#' If only two samples are provided, a fold change (FC) enrichment score will also be
#' calculated. The FC output is generated from a running sum of mean changes in gene
#' expression from all genes of the pathway. It's calculated from average pathway
#' expression in population1 - population2, so a negative FC means the pathway is
#' higher in population2.
#'
#' @export

compare_sce <- function(sce_object,
                        assay_name = "logcounts",
                        group1 = NULL,
                        group1_population = NULL,
                        group2 = NULL,
                        group2_population = NULL,
                        pathways,
                        min_genes = 15,
                        max_genes = 500,
                        downsample = 500,
                        parallel = FALSE,
                        cores = NULL) {

  ## Pathways
  if (class(pathways)[1] == "character") {
    pathways <- get_paths(pathways)
  }
  path_names <- sapply(pathways, function(x) unique(x$Pathway))

  ## Seurat extract
  if (is.null(group2)) {
    samples <- list()
    for (i in group1_population) {
      samples[[i]] <- sce_extract(sce_object,
                                  assay_name = assay_name,
                                  meta1 = group1,
                                  value_meta1 = i)
    }
  }

  if (!is.null(group2)) {
    samples <- list()
    for (i in group1_population) {
      samples[[i]] <- sce_extract(sce_object,
                                  assay_name = assay_name,
                                  meta1 = group1,
                                  value_meta1 = i,
                                  meta2 = group2,
                                  value_meta2 = group2_population)
    }
  }


  if (parallel == FALSE) {

    mcm_output <- compare_pathways(samples = samples,
                                   pathways = pathways,
                                   downsample = downsample,
                                   min_genes = min_genes,
                                   max_genes = max_genes)

  } else {

    mcm_output <- compare_pathways(samples = samples,
                                   pathways = pathways,
                                   downsample = downsample,
                                   min_genes = min_genes,
                                   max_genes = max_genes,
                                   parallel = TRUE,
                                   cores = cores)

  }


  return(mcm_output)

}




