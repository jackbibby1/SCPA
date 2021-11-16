#' Plot the rank of specific pathways from SCPA output
#'
#' This function takes the output from SCPA and plots the
#' rank of a user defined pathway.
#'
#' @param scpa_out Data frame containing Pathways and qvals generated from
#'     compare_pathways
#' @param pathway Chosen pathway or pathways to plot the rank of. This can
#'     be specific e.g. HALLMARK_GLYCOLYSIS to plot a specific result, or
#'     generic e.g. glycolysis to plot all glycolysis pathways
#' @param population_name Column name of the population to plot, if more
#'     than one population present in the data frame
#' @param base_point_size Size of the base points to plot on the graph
#' @param base_point_color Color of the base points to plot on the graph
#' @param highlight_point_size Size of the highlighted points to plot on the graph
#' @param highlight_point_color Color of the highlighted points to plot on the graph
#'
#' @examples \dontrun{
#' plot_rank(
#'      scpa_out = scpa_result,
#'      pathway = "interferon",
#'      population_name = cd4_qval,
#' )
#' }
#'
#' @return list of pathways with corresponding genes
#' @export

plot_rank <- function(scpa_out,
                      pathway,
                      population_name = "qval",
                      base_point_size = 2,
                      base_point_color = "royalblue2",
                      highlight_point_size = 3,
                      highlight_point_color = "orangered2",
                      label_pathway = NULL,
                      label_size = 4) {

  selected_paths <- grep(pattern = paste(pathway, collapse = "|"), x = scpa_out$Pathway, ignore.case = T,  value = T)

  path_ranking <- arrange(scpa_out, desc(population_name))
  path_ranking$path_rank <- percent_rank(path_ranking[[population_name]])*100

  df_sub <- subset(path_ranking, path_ranking$Pathway %in% selected_paths)

  if (label_pathway == T) {
    path_lab <- grep(pathway, path_ranking$Pathway, value = T, ignore.case = T)
    path_lab <- path_ranking[path_ranking$Pathway %in% path_lab, ]

    ggplot(path_ranking, aes(.data[[population_name]], path_rank)) +
      geom_hline(yintercept = c(0, 25, 50, 75, 100), linetype = 'dotted', lwd = 0.3, color = 'gray40') +
      geom_point(shape = 21, cex = base_point_size, color = 'black', fill = base_point_color, stroke = 0.05) +
      ggrepel::geom_label_repel(data = path_lab, label = path_lab$Pathway,
                                size = label_size, label.padding = unit(0.7, "mm"),
                                label.r = unit(0.3, "mm"), nudge_x = -30) +
      geom_point(data = df_sub, shape = 21, cex = highlight_point_size, color = 'black', fill = highlight_point_color) +
      xlab("Qval") +
      ylab("Pathway rank") +
      scale_y_continuous(expand = c(0.03, 0.03), breaks = c(0, 25, 50, 75, 100)) +
      scale_x_continuous(expand = c(0.2, 0.2)) +
      theme(panel.border = element_rect(fill = NA),
            panel.background = element_blank(),
            title = element_text(size = 9),
            axis.title = element_text(size = 11))

  } else {

    ggplot(path_ranking, aes(.data[[population_name]], path_rank)) +
      geom_hline(yintercept = c(0, 25, 50, 75, 100), linetype = 'dotted', lwd = 0.3, color = 'gray40') +
      geom_point(shape = 21, cex = base_point_size, color = 'black', fill = base_point_color, stroke = 0.05) +
      geom_point(data = df_sub, shape = 21, cex = highlight_point_size, color = 'black', fill = highlight_point_color) +
      xlab("Qval") +
      ylab("Pathway rank") +
      scale_y_continuous(expand = c(0.03, 0.03), breaks = c(0, 25, 50, 75, 100)) +
      scale_x_continuous(expand = c(0.2, 0.2)) +
      theme(panel.border = element_rect(fill = NA),
            panel.background = element_blank(),
            title = element_text(size = 9),
            axis.title = element_text(size = 11))

  }
}
