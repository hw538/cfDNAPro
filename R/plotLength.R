
plotLength <- function(x,
                       plot_type = c("fraction"),
                       xlim = c(30, 500),
                       vline = c(167),
                       ylim = NA_real_ ,
                       x_breaks = c(50, 100, 167, 200, 300, 400 , 500),
                       x_labels = c("50", "100", "167", "200", "300", "400", "500"),
                       line_color = "grey1",
                       line_size = 0.8,
                       line_alpha = 0.8,
                       vline_color = "grey",
                       vline_type = "dashed",
                       vline_size = 0.7,
                       ...
                       ){
  
  
  
  y_use <- NULL
  x <- tibble::as_tibble(x)
  plot_type <-  stringr::str_to_lower(plot_type)
  max_count <- max(dplyr::pull(x, All_Reads.fr_count))
  max_fraction <- max(dplyr::pull(x, prop ))
  

  
  if(plot_type == c("fraction")) {
    ylim <- c(0, max_fraction + 0.0015)
    y_use <- "prop"
  } else if (plot_type == c("count")) {
    ylim <- c(0, max_count + 1000)
    y_use <- "All_Reads.fr_count"
  }
    
  x$insert_size <- as.numeric(x$insert_size)
  plot_type <- stringr::str_to_title(plot_type)
  
      

      p <-  ggplot(x) +
        geom_line(aes(.data$insert_size, .data[[y_use]]), 
                  color = line_color, 
                  size = line_size, 
                  alpha = line_alpha) +
        geom_vline(xintercept = vline, linetype = vline_type, color = vline_color, size = vline_size) + 
        coord_cartesian(xlim = xlim) +
        labs(x = "cfDNA Fragment Length (bp)", y = plot_type) +
        scale_x_continuous(limits = xlim, 
                           breaks = x_breaks, 
                           labels = x_labels) + 
        scale_y_continuous(expand = c(0, 0)) +
        theme_length_plot()

      
      
      return(p)
  
}
