

plotMotif <- function(x,
                      ylim,
                      x_title = "Motif",
                      plot_type = c("Fraction"),
                      bar_color = c("A" = "cornflowerblue",
                                     "C" = "darksalmon",
                                     "G" = "palevioletred2",
                                     "T" = "forestgreen"),
                      motif_levels = c("C", "G", "A", "T"),
                      ...){
  
  y_use <- NULL
  x <- tibble::as_tibble(x)
  plot_type <-  stringr::str_to_lower(plot_type)
  max_count <- max(dplyr::pull(x, .data$n))
  max_fraction <- max(dplyr::pull(x, .data$fraction ))
  
  
  if(plot_type == c("fraction")) {
    
    y_use <- "fraction"
    
    if(missing(ylim)) {
      ylim <- c(0, max_fraction * 1.1)
      message("ylim", " was set as: ", ylim[1]  ," - ", ylim[2] )
    }
    
    
  } else if (plot_type == c("n")) {
    y_use <- "n"
    
    if(missing(ylim)) {
      ylim <- c(0, max_count * 1.2)
      
      message("ylim", " was set as: ", ylim[1]  ," - ", ylim[2] )
      
    }
    
  }
  
  x <- x %>%
    dplyr::mutate(group  = stringr::str_extract(.data$motif, pattern = "^[ATCG]")) %>%
    dplyr::filter(.data$group %in% motif_levels)
  
  x$group <- factor(x$group, levels = motif_levels)
  
  x <- x %>%
    dplyr::arrange(.data$group)
  
  x$motif <-  factor(x$motif, levels = x$motif)
  
  
  p <- ggplot(data = x, 
              mapping = aes(x = .data$motif, 
                            y =  .data[[y_use]], 
                            fill = .data$group)) +
    geom_col() +
    labs(x = x_title, y = stringr::str_to_title(plot_type)) +
    scale_fill_manual(values = bar_color) +
    scale_y_continuous(limits = ylim) +
    theme_motif_plot()
  return(p)
  
}


