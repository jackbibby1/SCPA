#' Use SCPA to compare gene sets
#'
#' This function takes an input of samples and pathways
#' to compare gene set perturbations over different conditions.
#'
#' @param samples List of samples, each supplied as an expression matrix with cells in columns
#'     and genes in rows
#' @param pathways Pathway gene sets with each pathway in a separate list. For formatting of
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
#' calculated. The FC output is generated from a running sum of mean changes in gene
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
  path_names <- sapply(pathways, function(x) unique(x$Pathway))

  # define the number of cells in each condition
  cell_number <- lapply(samples, function(x) ncol(x))
  cell_number <- sapply(cell_number, function(x) x[1])

  for (i in 1:length(cell_number)) {
    message(paste("Cell numbers in population", i, "=", cell_number[i]))}
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

  # generate pathway matrices
  pop_paths <- vector(mode = "list", length = length(samples))
  for (i in 1:length(pop_paths)) {
    for (c in 1:length(pathways)) {
      pop_paths[[i]][[c]] <- samples[[i]][rownames(samples[[i]]) %in% pathways[[c]]$Genes, ]
    }
  }

  # set pathway names
  pop_paths <- lapply(pop_paths, function(x) purrr::set_names(x, path_names))

  # filter out pathways with < 15 genes or > 500
  filter_paths <- sapply(pop_paths[[1]], function(x) any(nrow(x) >= min_genes & nrow(x) <= max_genes))

  filtered_pathways <- names(filter_paths[filter_paths == "FALSE"])

  if (length(filtered_pathways) > 0) {
    message("Excluding ", length(filtered_pathways),
            " pathways based on min/max genes parameter: ",
            paste(utils::head(filtered_pathways, 5), collapse=", "), "...", "\n")
  } else {
    message("All ", length(pathways), " pathways passed the min/max genes threshold", "\n")
  }

  pop_paths <- lapply(pop_paths, function(x) x[unlist(filter_paths)])

  # transpose matrix
  pop_paths <- lapply(pop_paths, function(x) lapply(x, function(c) t(c)))

  # order columns
  pop_paths <- lapply(pop_paths, function(x) lapply(x, function(c) c[, sort(colnames(c))]))

  # fc calculation
  if (length(samples) == 2) {

    message("Calculating pathway fold changes...", "\n")

    avg_expression <- lapply(pop_paths, function(x) lapply(x, function(c) data.frame(colMeans(c))))
    samp_combined <- c()
    for (i in 1:length(pop_paths[[1]])) {
      samp_combined[[i]] <- cbind(avg_expression[[1]][[i]], avg_expression[[2]][[i]])
    }
    samp_combined <- lapply(samp_combined, function(x) magrittr::set_colnames(x, c("Pop1", "Pop2")))
    samp_combined <- lapply(samp_combined, function(x) cbind(x, logFC = x[, "Pop1"]-x[, "Pop2"]))
    path_fc <- sapply(samp_combined, function(x) sum(x[, "logFC"]))
  }

  # run scpa
  if (length(samples) > 2) {
    message("Performing a multisample analysis with SCPA...")
    pb <- utils::txtProgressBar(min = 0, max = length(pop_paths[[1]]),
                                style = 3, width = 50)
    mcm_result <- list()
    for (i in 1:length(pop_paths[[1]])) {
      Sys.sleep(0.1)
      utils::setTxtProgressBar(pb, i)
      mcm_result[[i]] <- multicross::mcm(lapply(pop_paths, function(x) x[[i]]), level = 0.05)
    }
    close(pb)

    mcm_output <- data.frame(t(sapply(mcm_result, c)), stringsAsFactors = F)
    mcm_output$Pathway <- names(pop_paths[[1]])
    mcm_output$X2 <- NULL
    colnames(mcm_output)[1] <- c("Pval")
    mcm_output$Pval <- replace(mcm_output$Pval, mcm_output$Pval == 0, 10^-300)
    mcm_output$Pval <- as.numeric(mcm_output$Pval)
    mcm_output$adjPval <- stats::p.adjust(mcm_output$Pval, method = "bonferroni",
                                          n = nrow(mcm_output))
    mcm_output$qval <- sqrt(-log10(mcm_output$adjPval))
    mcm_output <- mcm_output[, c(2, 1, 3, 4)]
    mcm_output <- mcm_output[order(-mcm_output$qval), ]
    return(mcm_output)

  } else {
    message("Performing a two-sample analysis with SCPA...")
    pb <- utils::txtProgressBar(min = 0, max = length(pop_paths[[1]]),
                                style = 3, width = 50)
    mcm_result <- list()
    for (i in 1:length(pop_paths[[1]])) {
      Sys.sleep(0.1)
      utils::setTxtProgressBar(pb, i)
      mcm_result[[i]] <- multicross::mcm(lapply(pop_paths, function(x) x[[i]]), level = 0.05)
    }
    close(pb)
    mcm_output <- data.frame(t(sapply(mcm_result, c)), stringsAsFactors = F)
    mcm_output$FC <- path_fc
    mcm_output$Pathway <- names(pop_paths[[1]])
    mcm_output$X2 <- NULL
    colnames(mcm_output)[1] <- c("Pval")
    mcm_output[mcm_output$Pval == 0] <- 10^-300
    mcm_output$Pval <- as.numeric(mcm_output$Pval)
    mcm_output$adjPval <- stats::p.adjust(mcm_output$Pval, method = "bonferroni",
                                          n = nrow(mcm_output))
    mcm_output$qval <- sqrt(-log10(mcm_output$adjPval))
    mcm_output <- mcm_output[, c(3, 1, 4, 5, 2)]
    mcm_output <- mcm_output[order(-mcm_output$qval), ]
    return(mcm_output)
  }
}




