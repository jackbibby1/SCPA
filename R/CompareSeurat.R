#' Use SCPA to compare pathways within a Seurat object
#'
#' This function takes a Seurat object as an input, and
#' compares gene sets over specified conditions or populations.
#'
#' @param seurat_object Seurat object with defined clusters
#' @param group1 First comparison group as defined by column names in
#'   Seurat object e.g. cell_type
#' @param group1_population Population within group1 to compare
#'   e.g. t_cell
#' @param group2 Second comparison group as defined by column names in
#'   Seurat object e.g. hour
#' @param group2_population Population within group2 to compare
#'   e.g. 24
#' @param pathways List of pathways and their genes
#' @param downsample Option to downsample cell numbers. Default is 500
#'
#' @examples \dontrun{
#' scpa_out <- compare_seurat(
#'      list(sample1, sample2, sample3),
#'      pathways = pathways)
#' }
#'
#' @return Statistical results from the multicross tool.
#' If only two samples are provided, an enrichment score will also be
#' calculated.
#'
#' @export

compare_seurat <- function(seurat_object,
                           group1 = NULL,
                           group1_population = NULL,
                           group2 = NULL,
                           group2_population = NULL,
                           pathways,
                           downsample = 500) {

  ## Pathways
  if (class(pathways)[1] == "character") {
    pathways <- get_paths(pathways)
  }
  path_names <- sapply(pathways, function(x) unique(x$Pathway))

  ## Seurat extract
  if (is.null(group2)) {
    samples <- list()
    for (i in group1_population) {
      samples[[i]] <- seurat_extract(seurat_object,
                                     meta1 = group1,
                                     value_meta1 = i)
    }
  }

  if (!is.null(group2)) {
    samples <- list()
    for (i in group1_population) {
      samples[[i]] <- seurat_extract(seurat_object,
                                     meta1 = group1,
                                     value_meta1 = i,
                                     meta2 = group2,
                                     value_meta2 = group2_population)
    }
  }

  ## Get cells in each sample
  cell_number <- lapply(samples, function(x) ncol(x))
  cell_number <- sapply(cell_number, function(x) x[1])

  ## Sample random cells to 500
  for (i in 1:length(samples)) {
    samples[[i]] <- random_cells(samples[[i]], ifelse(cell_number[i] < 500, cell_number[i], downsample))
  }

  pop_paths <- vector(mode = "list", length = length(samples))
  ## Pathway matrices
  for (i in 1:length(pop_paths)) {
    for (c in 1:length(pathways)) {
      pop_paths[[i]][[c]] <- samples[[i]][rownames(samples[[i]]) %in% pathways[[c]]$Genes, ]
    }
  }

  ##### Convert the matrices to get genes as cols and rows as cells #####
  pop_paths <- lapply(pop_paths, function(x) lapply(x, function(c) t(c)))

  ## Take genes present in all samples
  common_genes <- list()
  for (i in 1:length(pathways)) {
    common_genes[[i]] <- Reduce(intersect, lapply(pop_paths, function(x) colnames(x[[i]])))
  }
  pop_paths <- lapply(pop_paths, function(x) lapply(x, function(c) c[, colnames(c) %in% unlist(common_genes)]))

  ## Organise columns into same order
  pop_paths <- lapply(pop_paths, function(x) lapply(x, function(c) c[, sort(colnames(c))]))

  ## Calculate fold changes
  if (length(samples) == 2) {
    avg_expression <- lapply(pop_paths, function(x) lapply(x, function(c) data.frame(colMeans(c))))
    samp_combined <- c()
    for (i in 1:length(pathways)) {
      samp_combined[[i]] <- cbind(avg_expression[[1]][[i]], avg_expression[[2]][[i]])
    }
    samp_combined <- lapply(samp_combined, function(x) magrittr::set_colnames(x, c("Pop1", "Pop2")))
    samp_combined <- lapply(samp_combined, function(x) cbind(x, logFC = x[, "Pop1"]-x[, "Pop2"]))
    path_fc <- sapply(samp_combined, function(x) sum(x[, "logFC"]))
  }

  ## Test samples
  if (length(samples) > 2) {
    message("Conducting multivariate analysis")
    pb <- utils::txtProgressBar(min = 0, max = length(pathways),
                                style = 3, width = 50)
    mcm_result <- list()
    for (i in 1:length(pathways)) {
      Sys.sleep(0.1)
      utils::setTxtProgressBar(pb, i)
      mcm_result[[i]] <- multicross::mcm(lapply(pop_paths, function(x) x[[i]]), level = 0.05)
    }
    mcm_output <- data.frame(t(sapply(mcm_result, c)), stringsAsFactors = F)
    mcm_output$Pathway <- path_names
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
    message("Conducting 2-sample test")
    pb <- utils::txtProgressBar(min = 0, max = length(pathways),
                                style = 3, width = 50)
    mcm_result <- list()
    for (i in 1:length(pathways)) {
      Sys.sleep(0.1)
      utils::setTxtProgressBar(pb, i)
      mcm_result[[i]] <- multicross::mcm(lapply(pop_paths, function(x) x[[i]]), level = 0.05)
    }
    mcm_output <- data.frame(t(sapply(mcm_result, c)), stringsAsFactors = F)
    mcm_output$FC <- path_fc
    mcm_output$Pathway <- path_names
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
