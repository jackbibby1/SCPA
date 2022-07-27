#' Plot a heatmap of qvals from SCPA
#'
#' This function takes the output from SCPA and plots a
#' heatmap of the qvals
#'
#' @param scpa_out Data frame that contains a "Pathway" column, and a column or columns
#'     containing qvals. This can be the direct output from `compare_pathways()`,
#'     or a custom data frame.
#' @param highlight_pathways Pathway or pathways to annotate on the heatmap,
#'     supplied as character vector. If no argument is given, all pathways are
#'     shown
#' @param row_fontsize Font size of pathway names
#' @param column_fontsize Font size of the sample names
#' @param column_names Option to supply names of the heatmap
#'     columns if more than one populations is present
#' @param show_row_names Should row names be shown in the heatmap?
#' @param cluster_columns Should columns in the heatmap be clustered?
#' @param hm_colors Colors to be used in the heatmap
#' @param scale_breaks Breaks to be used in the colors of the heatmap.
#'     Length must be equal to the number of colors.
#'
#' @importFrom magrittr "%>%"
#'
#' @examples \dontrun{
#' plot_heatmap(
#'      scpa_out = scpa_result,
#'      pathway = "mtorc"
#' )
#' }
#'
#' @return list of pathways with corresponding genes
#' @export

plot_heatmap <- function(scpa_out,
                         highlight_pathways = NULL,
                         row_fontsize = 6,
                         column_fontsize = 10,
                         column_names = colnames(scpa_out),
                         show_row_names = TRUE,
                         cluster_columns = TRUE,
                         hm_colors = NULL,
                         scale_breaks = NULL) {

  if (!is.null(hm_colors)) {
    heatmap_colors <- hm_colors
  } else {
    heatmap_colors <- c("cornflowerblue", "white", "red")
  }

  scale <- scpa_out %>%
    dplyr::select(grep(pattern = "qval", x = colnames(scpa_out), ignore.case = T, value = T)) %>%
    as.matrix()
  scale_breaks <- c(min(scale), mean(scale), max(scale))

  hm_col <- circlize::colorRamp2(colors = heatmap_colors, breaks = scale_breaks)

  if (is.null(highlight_pathways) == F) {

    pathways <- scpa_out %>%
      dplyr::filter(grepl(paste(highlight_pathways, collapse = "|"), Pathway, ignore.case = T)) %>%
      dplyr::pull(Pathway)

    pathways <- gsub(pattern = "_", replacement = " ", x = pathways)

    scpa_out <- scpa_out %>%
      dplyr::mutate(Pathway = gsub(pattern = "_", replacement =  " ", x = Pathway))

    position <- scpa_out %>%
      dplyr::mutate(position = 1:nrow(.)) %>%
      dplyr::filter(grepl(paste(pathways, collapse = "|"), Pathway, ignore.case = T)) %>%
      dplyr::pull(position)

    scpa_out <- scpa_out %>%
      tibble::remove_rownames() %>%
      tibble::column_to_rownames("Pathway") %>%
      dplyr::select(grep(pattern = "qval", x = colnames(scpa_out), ignore.case = T, value = T))

    row_an <- ComplexHeatmap::rowAnnotation(Genes = ComplexHeatmap::anno_mark(at =  which(rownames(scpa_out) %in% pathways),
                                              labels = rownames(scpa_out)[position],
                                              labels_gp = grid::gpar(fontsize = 7),
                                              link_width = grid::unit(2.5, "mm"),
                                              padding = grid::unit(1, "mm"),
                                              link_gp = grid::gpar(lwd = 0.5)))

    ComplexHeatmap::Heatmap(scpa_out,
            name = "Qval",
            border = T,
            rect_gp = grid::gpar(col = "white", lwd = 0.1),
            row_names_gp = grid::gpar(fontsize = row_fontsize),
            column_names_gp = grid::gpar(fontsize = column_fontsize),
            right_annotation = row_an,
            column_labels = column_names,
            col = hm_col,
            show_row_names = show_row_names,
            row_dend_width = grid::unit(6, "mm"),
            cluster_columns = cluster_columns)

  } else {

    scpa_out <- scpa_out %>%
      dplyr::mutate(Pathway = gsub(pattern = "_", replacement = " ", x = Pathway)) %>%
      tibble::remove_rownames() %>%
      tibble::column_to_rownames("Pathway") %>%
      dplyr::select(grep(pattern = "qval", x = colnames(scpa_out), ignore.case = T, value = T))

    ComplexHeatmap::Heatmap(scpa_out,
            name = "Qval",
            border = T,
            rect_gp = grid::gpar(col = "white", lwd = 0.1),
            row_names_gp = grid::gpar(fontsize = row_fontsize),
            column_names_gp = grid::gpar(fontsize = column_fontsize),
            show_row_names = show_row_names,
            column_labels = column_names,
            col = hm_col,
            row_dend_width = grid::unit(6, "mm"),
            cluster_columns = cluster_columns)

  }

}
