
#' Plot the etiologies of a computed CHAMPS table
#'
#' @param x A computed CHAMPS object.
#' @param etiologies The etiologies to include in the graphic with their shortened names as a named vector. Defaults to NULL.  If NULL than the top variable is used and the names are shortened automatically. 
#' @param top The number of etiologies to include in the graphic
#' @param bar_width defines the size of the bar the chart. Defaults to 0.5. 
# @examples ekeep <- c(`Acinetobacter baumannii` = "A. baumannii", `Salmonella spp.` = "Salmonella spp.", Parechovirus = "Parechovirus")
#' @importFrom ggplot2 ggplot aes geom_text theme_bw theme element_text labs ylim
#' @export
plot_etiologies <- function(x, etiologies = NULL, top = 6, bar_width = .5) {
  
  check_champs_object(x)
  
  d1 <- dplyr::ungroup(x$etiol_counts)
  
  # total deaths
  d1_deaths <- sum(d1$n)
  
  if (is.null(etiologies)) {
    common_etiol <- d1 %>%
      dplyr::group_by(etiol) %>%
      dplyr::summarise(count = sum(n)) %>%
      dplyr::arrange(desc(count)) %>%
      dplyr::filter(!is.na(etiol)) %>%
      dplyr::slice(1:top) %>%
      dplyr::pull(etiol)
    
    names_table <- stringr::str_split_fixed(common_etiol, pattern = " ", n = 2)
    
    # find the cases that don't need to be shortened.  
    one_word_name <- names_table[,2] == ""
    spp_name <- stringr::str_detect(names_table[,2], "spp")
    
    # spp names shouldn't be shortened.
    names_table[,2][spp_name] <- paste(names_table[,1][spp_name], names_table[,2][spp_name])
    
    # one word names shouldn't be shortened.  Put word in second column then make first word NA
    names_table[,2][one_word_name] <- names_table[,1][one_word_name]
    names_table[,1] <- paste0(stringr::str_trunc(names_table[,1], 1, ellipsis = ""), ".")
    
    
    names_table[,1][one_word_name] <- ""
    names_table[,1][spp_name] <-  ""
      
    short_names <- paste(names_table[,1], names_table[,2]) %>% stringr::str_trim()
    
    # Now build etiolgies object to be like non-NULL value potentially entered by user.    
    etiologies <- short_names
    names(etiologies) <- common_etiol
    
  } else {
    
    common_etiol <- names(etiolgies)
    short_names <- etiologies
    
  }
    
  group_colors <- c("#405b9f", "#5494c8", "#1e73cb", "#73849e")
    
  p1 <- d1 %>%
    dplyr::filter(etiol %in% common_etiol) %>%
    dplyr::mutate(eshort = etiologies[etiol]  %>% forcats::fct_reorder(n, .fun = sum, .desc = TRUE)) %>%
    ggplot2::ggplot(ggplot2::aes(x = eshort, y = n, fill = case_type_desc2)) +
    ggplot2::geom_col(width = bar_width) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom", panel.grid.major.x = ggplot2::element_blank()) +
    ggplot2::scale_fill_manual(values = group_colors) +
    ggplot2::scale_y_continuous(breaks = seq(0, d1_deaths, by = 20)) +
    ggplot2::labs(y = "", x = "", fill = "",
                  title = paste0("Common Etiologies of CHAMPS\nNeonatal Infectious Deaths: All Sites (N=", d1_deaths ,")"))

p1


}

# plot_etiologies(pne3) # pne3$etiol_counts for slide 4
