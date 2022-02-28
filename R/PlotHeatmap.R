#' Plot a heatmap of qvals from SCPA
#'
#' This function takes the output from SCPA and plots a
#' heatmap of the qvals
#'
#' @param scpa_out Data frame containing Pathways and qvals generated from
#'     compare_pathways
#' @param highlight_pathways Chosen pathway or pathways to highlight on the heatmap,
#'     supplied as character vector. If no argument is given, all pathways are
#'     shown
#' @param row_fontsize Font size of pathway names
#' @param column_fontsize Font size of the sample names
#' @param column_names Option to supply names of the heatmap
#'     columns if more than one populations is present
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
                         show_row_names = FALSE,
                         cluster_columns = TRUE,
                         hm_colors = NULL,
                         scale_breaks = NULL) {

  hm_colors <- c("cornflowerblue", "white", "red")
  scale <- scpa_out %>%
    select(grep(pattern = "qval", x = colnames(.), ignore.case = T, value = T)) %>%
    as.matrix()
  scale_breaks <- c(min(scale), mean(scale), max(scale))
  hm_col <- colorRamp2(colors = hm_colors, breaks = scale_breaks)

  if (is.null(highlight_pathways) == F) {

    pathways <- scpa_out %>%
      filter(grepl(paste(highlight_pathways, collapse = "|"), Pathway, ignore.case = T)) %>%
      pull(Pathway)

    pathways <- gsub(pattern = "_", replacement = " ", x = pathways)

    scpa_out <- scpa_out %>%
      mutate(Pathway = gsub(pattern = "_", replacement =  " ", x = Pathway))

    position <- scpa_out %>%
      mutate(position = 1:nrow(.)) %>%
      filter(grepl(paste(pathways, collapse = "|"), Pathway, ignore.case = T)) %>%
      pull(position)

    scpa_out <- scpa_out %>%
      remove_rownames() %>%
      column_to_rownames("Pathway") %>%
      select(grep(pattern = "qval", x = colnames(.), ignore.case = T, value = T))

    row_an <- rowAnnotation(Genes = anno_mark(at =  which(rownames(scpa_out) %in% pathways),
                                              labels = rownames(scpa_out)[position],
                                              labels_gp = gpar(fontsize = 7),
                                              link_width = unit(2.5, "mm"),
                                              padding = unit(1, "mm"),
                                              link_gp = gpar(lwd = 0.5)))

    ht_opt$message = FALSE

    Heatmap(scpa_out,
            name = "Qval",
            border = T,
            rect_gp = gpar(col = "white", lwd = 0.1),
            row_names_gp = gpar(fontsize = row_fontsize),
            column_names_gp = gpar(fontsize = column_fontsize),
            right_annotation = row_an,
            column_labels = column_names,
            col = hm_col,
            show_row_names = show_row_names,
            row_dend_width = unit(6, "mm"),
            cluster_columns = cluster_columns)


  } else {

    scpa_out <- scpa_out %>%
      mutate(Pathway = gsub(pattern = "_", replacement = " ", x = Pathway)) %>%
      remove_rownames() %>%
      column_to_rownames("Pathway") %>%
      select(grep(pattern = "qval", x = colnames(.), ignore.case = T, value = T))

    ht_opt$message = FALSE

    Heatmap(scpa_out,
            name = "Qval",
            border = T,
            rect_gp = gpar(col = "white", lwd = 0.1),
            row_names_gp = gpar(fontsize = row_fontsize),
            column_names_gp = gpar(fontsize = column_fontsize),
            show_row_names = show_row_names,
            column_labels = column_names,
            col = hm_col,
            row_dend_width = unit(6, "mm"),
            cluster_columns = cluster_columns)

  }

}

plot_heatmap(df,
             highlight_pathways = c("myc", "mtorc"),
             cluster_columns = F,
             column_names = c("t cell", "b cell"))






