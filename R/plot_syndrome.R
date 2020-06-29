# plot_agebreakdown(pne6, plot_type = "pie")
# plot_agebreakdown(pne6, plot_type = "waffle")
# plot_agebreakdown(pne6, plot_type = "bar", include_legend = FALSE)
# plot_agebreakdown(pne6, plot_type = "tree", include_legend = FALSE)


#' Plot the syndrome distribution from a computed CHAMPS table
#'
#' @param x A computed CHAMPS object.
#' @param legend_location Can be one of "none", "bottom", "right", "left". Defaults to none.
#' @param plot_type One of four optins. Pie Chart 'pie', Bar Chart 'bar', Treemap 'tree', Waffle 'waffle'
#' @param include_pie_numbers Whether to include the numbers within the slices of the pie chart. Defaults to FALSE.
#' @param include_bar_text Whether to include the category text within the bar. Defaults to TRUE
#' @param waffle_cols The number of columns of the waffle chart.  
#' @importFrom ggplot2 ggplot aes geom_text theme_bw theme element_text labs ylim
#' @importFrom forcats fct_reorder
# @import treemapify
#' @import ggfittext
#' @export
plot_syndrome <- function(x, legend_location = "none", plot_type = "pie", 
                              include_pie_numbers = FALSE, include_bar_text = TRUE,
                              waffle_cols = 12) {
  
  check_champs_object(x)
  
  d1 <- dplyr::ungroup(x$table) %>%
    dplyr::mutate(percent = round(100*n/sum(n),0),
                  percent_label = dplyr::case_when(percent <= 2 ~ "" ,
                                            percent <= 5 ~ paste0(percent),
                                            TRUE ~ paste0(percent, "%")))
  
  if(nrow(d1) > 6) stop("Including more than 6 groups is not supported")
  
  
  piecolors <- c("#4c98a4","#938595","#405b9f", "#5494c8", "#1e73cb", "#73849e")

  
  if (!any(plot_type %in% c("pie", "bar", "tree", "waffle"))) stop("Expecting one of 'pie', 'bar', 'treemap' or 'waffle'")
  
  if (plot_type == "pie") {
    p_out <- d1  %>%
      ggplot2::ggplot(ggplot2::aes(x = "", y = n, fill = syndrome)) +
      ggplot2::geom_col(position = "fill") +
      ggplot2::coord_polar("y", start = 0) +
      ggplot2::scale_fill_manual(values = piecolors) +
      ggplot2::guides(fill = ggplot2::guide_legend(nrow=2,byrow=TRUE)) +
      ggplot2::labs(fill = "") +
      ggplot2::theme_minimal()+
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(size=14, face="bold"),
        legend.position = legend_location)
    
    if (include_pie_numbers) {
      p_out <- p_out + 
        ggplot2::geom_text(ggplot2::aes(label = percent_label), color = "white", 
                           position = ggplot2::position_fill(vjust = 0.5)) 
    }
    
  } else if (plot_type == "bar" ) {
    # https://github.com/wilkox/ggfittext
   p_out <- d1  %>%
      dplyr::mutate(xuse = forcats::fct_reorder(syndrome, n, .desc = TRUE)) %>%
      ggplot2::ggplot(ggplot2::aes(x = xuse, 
                                   y = percent, fill = xuse, label = xuse)) +
      ggplot2::geom_col() +
      ggplot2::scale_fill_manual(values = piecolors) +
      ggplot2::guides(fill = ggplot2::guide_legend(nrow=2,byrow=TRUE)) +
      ggplot2::labs(fill = "", y = "Percentage") +
      ggplot2::theme_minimal()+
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(),
        panel.grid.major.x =  ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(size=14, face="bold"),
        legend.position = legend_location) 
   
     if (include_bar_text) {
       p_out <- p_out + 
         ggfittext::geom_bar_text(reflow = TRUE, contrast = TRUE)

     }
   
  } else if (plot_type == "tree") {
    # # https://github.com/wilkox/treemapify
    # # https://github.com/wilkox/ggfittext
    # p_out <- d1 %>%
    #   dplyr::mutate(xuse = forcats::fct_reorder(syndrome, n, .desc = TRUE)) %>%
    #   ggplot2::ggplot(ggplot2::aes(area = n, fill = factor(n), label = paste0(xuse,"\n", percent, "%"))) +
    #   treemapify::geom_treemap(color = "black") +
    #   treemapify::geom_treemap_text(fontface = "italic", place = "centre", reflow = TRUE) +
    #   ggplot2::scale_fill_manual(values = piecolors) +
    #   ggplot2::guides(fill = ggplot2::guide_legend(nrow=2,byrow=TRUE)) +
    #   ggplot2::labs(fill = "Count") +
    #   ggplot2::theme(plot.title = ggplot2::element_text(size=14, face="bold"),
    #     legend.position = legend_location) +
    #   ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1))
    
  } else if (plot_type == "waffle") {
    # # https://github.com/liamgilbey/ggwaffle
    # total <- sum(d1$n)
    # col_use <- waffle_cols
    # rows_use <- ceiling(total / col_use)
    
    # p_out <- d1 %>%
    #   dplyr::group_by(syndrome) %>% 
    #   tidyr::expand(count = seq(1:n)) %>%
    #   dplyr::mutate(n = c(max(count), rep(NA, dplyr::n()-1)),
    #                 sort = max(count)) %>%
    #   dplyr::arrange(desc(sort)) %>%
    #   dplyr::ungroup() %>%
    #   dplyr::mutate(x = rep(1:col_use, each = rows_use)[1:dplyr::n()],
    #                 y = rep(1:rows_use, col_use)[1:dplyr::n()],
    #                 syndrome = forcats::fct_reorder(syndrome, sort, .desc = TRUE)) %>%
    #   ggplot2::ggplot(ggplot2::aes(x, y, fill = syndrome)) + 
    #   ggwaffle::geom_waffle() + 
    #   ggplot2::coord_equal() + 
    #   ggplot2::scale_fill_manual(values = piecolors) +
    #   ggplot2::theme_minimal() +
    #   ggplot2::theme(axis.title.x = ggplot2::element_blank(),
    #                  axis.text = ggplot2::element_blank(),
    #                  panel.border = ggplot2::element_blank(),
    #                  panel.grid =  ggplot2::element_blank(),
    #                  axis.ticks = ggplot2::element_blank(),
    #                  plot.title = ggplot2::element_text(size=14, face="bold"),
    #                  legend.position = legend_location) +
    #   ggplot2::labs(fill = "", x = "", y = "" )
  }

  p_out
}
