#' Use SCPA to compare gene sets
#'
#' This function takes an input of samples and pathways
#' to compare gene set perturbations over different conditions with SCPA.
#'
#' @param samples List of samples, each supplied as an expression matrix with cells in columns
#'     and genes in rows.
#' @param pathways Pathways and their genes with each pathway in a separate list. For formatting of
#'     gene lists, see documentation at https://jackbibby1.github.io/SCPA/articles/using_gene_sets.html
#' @param downsample Option to downsample cell numbers. Defaults to 500 cells per condition. If a population
#'     has < 500 cells, all cells from that condition are used.
#' @param min_genes Gene sets with fewer than this number of genes will be excluded
#' @param max_genes Gene sets with more than this number of genes will be excluded
#'
#' @examples \dontrun{
#' scpa_result <- compare_pathways(
#'      list(sample1, sample2, sample3),
#'      pathways = pathways)
#' }
#'
#' @return Statistical results from the SCPA analysis. The qval should be the
#' primary metric that is used to interpret pathway differences i.e. a higher
#' qval translates to larger pathway differences between conditions.
#' If only two samples are provided, a fold change (FC) enrichment score will also be
#' calculated. The FC statistic is generated from a running sum of mean changes in gene
#' expression from all genes of the pathway. It's calculated from average pathway
#' expression in population1 - population2, so a negative FC means the pathway is
#' higher in population2.
#'
#' @export

compare_pathways <- function(samples,
                             pathways,
                             downsample = 500,
                             min_genes = 15,
                             max_genes = 500) {

  # get pathways for analysis
  if (class(pathways)[1] == "character") {
    pathways <- get_paths(pathways)
  }

  # define the number of cells in each condition
  cell_number <- sapply(samples, function(x) ncol(x))

  for (i in 1:length(cell_number)) {
    message(paste("Cell numbers in population", i, "=", cell_number[i]))
  }

  message("- If greater than ", downsample,
          " cells, these populations will be downsampled", "\n")

  # randomly sample cells
  for (i in 1:length(samples)) {

    samples[[i]] <- random_cells(samples[[i]], ifelse(cell_number[i] < 500, cell_number[i], downsample))

  }

  # only take shared genes
  genes <- lapply(samples, function(x) rownames(x))
  genes <- table(unlist(genes))
  genes <- genes[genes == length(samples)]
  genes <- names(genes)
  samples <- lapply(samples, function(x) x[rownames(x) %in% genes, ])

  # filter out pathways
  gene_numbers <- lapply(pathways, function(x) nrow(samples[[1]][rownames(samples[[1]]) %in% x$Genes, ]))
  keep_pathway <- sapply(gene_numbers, function(x) any(x >= min_genes & x <= max_genes))
  excluded_pathways <- sapply(pathways[!keep_pathway], function(x) unique(dplyr::pull(x, Pathway)))
  pathways_filtered <- pathways[keep_pathway]

  if (length(pathways_filtered) == 0) {

    stop(call. = F, "No pathways passed the min/max genes threshold")

  } else if (length(excluded_pathways) > 0) {

    message("Excluding ", length(excluded_pathways),
            " pathway(s) based on min/max genes parameter: ",
            paste(utils::head(excluded_pathways, 5), collapse = ", "), "...", "\n")

  } else {

    message("All ", length(pathways), " pathways passed the min/max genes threshold", "\n")

  }


  if (length(samples) > 2) {

    message("Performing a multisample analysis with SCPA...")

  } else {

    message("Calculating pathway fold changes...", "\n")
    message("Performing a two-sample analysis with SCPA...")

  }

  pb <- utils::txtProgressBar(min = 0, max = length(pathways_filtered),
                              style = 3, width = 50)

  scpa_result <- list()
  for (i in 1:length(pathways_filtered)) {

    utils::setTxtProgressBar(pb, i)

    # subset data to get one pathway
    path_subset <- lapply(samples, function(x) x[rownames(x) %in% pathways_filtered[[i]]$Genes, ])
    path_subset <- lapply(path_subset, function(x) t(x))
    path_subset <- lapply(path_subset, function(x) x[, sort(colnames(x))])

    if (length(path_subset) == 2) {

      avg_expression <- lapply(path_subset, function(x) data.frame(colMeans(x)))
      samp_combined <- cbind(avg_expression[[1]], avg_expression[[2]])
      samp_combined <- magrittr::set_colnames(samp_combined, c("Pop1", "Pop2"))
      samp_combined <- cbind(samp_combined, logFC = samp_combined[, "Pop1"]-samp_combined[, "Pop2"])
      path_fc <- sum(samp_combined[, "logFC"])

      scpa_result[[i]] <- multicross::mcm(path_subset, level = 0.05) %>%
        data.frame() %>%
        t() %>%
        data.frame() %>%
        dplyr::mutate(FC = path_fc) %>%
        dplyr::mutate(Pathway = pathways_filtered[[i]]$Pathway[1]) %>%
        dplyr::select(-X2) %>%
        dplyr::mutate(Pval = as.numeric(X1)) %>%
        dplyr::select(-X1) %>%
        dplyr::mutate(adjPval = stats::p.adjust(Pval , method = "bonferroni",
                                         n = length(pathways_filtered))) %>%
        dplyr::mutate(qval = sqrt(-log10(adjPval))) %>%
        dplyr::select(Pathway, Pval, adjPval, qval, FC)

    } else {

      scpa_result[[i]] <- multicross::mcm(path_subset, level = 0.05) %>%
        data.frame() %>%
        t() %>%
        data.frame() %>%
        dplyr::mutate(Pathway = pathways_filtered[[i]]$Pathway[1]) %>%
        dplyr::select(-X2) %>%
        dplyr::mutate(Pval = as.numeric(X1)) %>%
        dplyr::select(-X1) %>%
        dplyr::mutate(adjPval = stats::p.adjust(Pval , method = "bonferroni",
                                         n = length(pathways_filtered))) %>%
        dplyr::mutate(qval = sqrt(-log10(adjPval))) %>%
        dplyr::select(Pathway, Pval, adjPval, qval)

    }

  }

  close(pb)

  scpa_result <- scpa_result %>%
    dplyr::bind_rows() %>%
    tibble::remove_rownames() %>%
    dplyr::arrange(desc(qval))

  return(scpa_result)

}

#' Use SCPA to compare gene sets
#'
#' This function takes an input of samples and pathways
#' to compare gene set perturbations over different conditions with SCPA.
#'
#' @param samples List of samples, each supplied as an expression matrix with cells in columns
#'     and genes in rows.
#' @param pathways Pathways and their genes with each pathway in a separate list. For formatting of
#'     gene lists, see documentation at https://jackbibby1.github.io/SCPA/articles/using_gene_sets.html
#' @param downsample Option to downsample cell numbers. Defaults to 500 cells per condition. If a population
#'     has < 500 cells, all cells from that condition are used.
#' @param min_genes Gene sets with fewer than this number of genes will be excluded
#' @param max_genes Gene sets with more than this number of genes will be excluded
#' @param cores Number of cores to use for parallel processing
#'
#' @examples \dontrun{
#' scpa_result <- compare_pathways_parallel(
#'      list(sample1, sample2, sample3),
#'      pathways = pathways),
#'      cores = 2)
#' }
#'
#' @return Statistical results from the SCPA analysis. The qval should be the
#' primary metric that is used to interpret pathway differences i.e. a higher
#' qval translates to larger pathway differences between conditions.
#' If only two samples are provided, a fold change (FC) enrichment score will also be
#' calculated. The FC statistic is generated from a running sum of mean changes in gene
#' expression from all genes of the pathway. It's calculated from average pathway
#' expression in population1 - population2, so a negative FC means the pathway is
#' higher in population2.
#'
#' @export

compare_pathways_parallel <- function(samples,
                             pathways,
                             downsample = 500,
                             min_genes = 15,
                             max_genes = 500,
                             cores = 2) {

  # get pathways for analysis
  if (class(pathways)[1] == "character") {
    pathways <- get_paths(pathways)
  }

  # define the number of cells in each condition
  cell_number <- sapply(samples, function(x) ncol(x))

  for (i in 1:length(cell_number)) {
    message(paste("Cell numbers in population", i, "=", cell_number[i]))
  }

  message("- If greater than ", downsample,
          " cells, these populations will be downsampled", "\n")

  # randomly sample cells
  for (i in 1:length(samples)) {

    samples[[i]] <- random_cells(samples[[i]], ifelse(cell_number[i] < 500, cell_number[i], downsample))

  }

  # only take shared genes
  genes <- lapply(samples, function(x) rownames(x))
  genes <- table(unlist(genes))
  genes <- genes[genes == length(samples)]
  genes <- names(genes)
  samples <- lapply(samples, function(x) x[rownames(x) %in% genes, ])

  # filter out pathways
  gene_numbers <- lapply(pathways, function(x) nrow(samples[[1]][rownames(samples[[1]]) %in% x$Genes, ]))
  keep_pathway <- sapply(gene_numbers, function(x) any(x >= min_genes & x <= max_genes))
  excluded_pathways <- sapply(pathways[!keep_pathway], function(x) unique(dplyr::pull(x, Pathway)))
  pathways_filtered <- pathways[keep_pathway]

  if (length(pathways_filtered) == 0) {

    stop(call. = F, "No pathways passed the min/max genes threshold")

  } else if (length(excluded_pathways) > 0) {

    message("Excluding ", length(excluded_pathways),
            " pathway(s) based on min/max genes parameter: ",
            paste(utils::head(excluded_pathways, 5), collapse = ", "), "...", "\n")

  } else {

    message("All ", length(pathways), " pathways passed the min/max genes threshold", "\n")

  }


  if (length(samples) > 2) {

    message("Performing a multisample analysis with SCPA...")

  } else {

    message("Calculating pathway fold changes...", "\n")
    message("Performing a two-sample analysis with SCPA...")

  }

  if (!require(doParallel)) {
        stop('doParallel library not loaded. Please exeucte library("doParallel").')
  } else {
    cluster <- makeCluster(cores, type = "PSOCK")
    registerDoParallel(cluster)

    scpa_result <- foreach(pathway = pathways_filtered) %dopar% {
      res <- tryCatch(
        expr = {
          # subset data to get one pathway
          path_subset <- lapply(samples, function(x) x[rownames(x) %in% pathway$Genes, ])
          path_subset <- lapply(path_subset, function(x) t(x))
          path_subset <- lapply(path_subset, function(x) x[, sort(colnames(x))])

          if (length(path_subset) == 2) {

            avg_expression <- lapply(path_subset, function(x) data.frame(colMeans(x)))
            samp_combined <- cbind(avg_expression[[1]], avg_expression[[2]])
            samp_combined <- magrittr::set_colnames(samp_combined, c("Pop1", "Pop2"))
            samp_combined <- cbind(samp_combined, logFC = samp_combined[, "Pop1"]-samp_combined[, "Pop2"])
            path_fc <- sum(samp_combined[, "logFC"])

            multicross::mcm(path_subset, level = 0.05) %>%
              data.frame() %>%
              t() %>%
              data.frame() %>%
              dplyr::mutate(FC = path_fc) %>%
              dplyr::mutate(Pathway = pathway$Pathway[1]) %>%
              dplyr::select(-X2) %>%
              dplyr::mutate(Pval = as.numeric(X1)) %>%
              dplyr::select(-X1) %>%
              dplyr::mutate(adjPval = stats::p.adjust(Pval , method = "bonferroni",
                                              n = length(pathways_filtered))) %>%
              dplyr::mutate(qval = sqrt(-log10(adjPval))) %>%
              dplyr::select(Pathway, Pval, adjPval, qval, FC)

          } else {

            multicross::mcm(path_subset, level = 0.05) %>%
              data.frame() %>%
              t() %>%
              data.frame() %>%
              dplyr::mutate(Pathway = pathway$Pathway[1]) %>%
              dplyr::select(-X2) %>%
              dplyr::mutate(Pval = as.numeric(X1)) %>%
              dplyr::select(-X1) %>%
              dplyr::mutate(adjPval = stats::p.adjust(Pval , method = "bonferroni",
                                              n = length(pathways_filtered))) %>%
              dplyr::mutate(qval = sqrt(-log10(adjPval))) %>%
              dplyr::select(Pathway, Pval, adjPval, qval)
          }
        },
        error = function(e) {
          return("Error in pathway ", pathway$Pathway[1], ": ", e)
        })
    }

    stopCluster(cluster)

    scpa_result <- scpa_result %>%
      dplyr::bind_rows() %>%
      tibble::remove_rownames() %>%
      dplyr::arrange(desc(qval))

    return(scpa_result)

  }
}
