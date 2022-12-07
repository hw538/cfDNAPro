
plotLength <- function(x,
                       plot_type = c("fraction"),
                       xlim = c(30, 500),
                       vline = c(167),
                       ylim = NA_real_ ,
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
  x_breaks <- c(xlim[[1]], 50, 100, 167, 200, 300, 400 , 500, 600, 700, 800, 900, 1000)
  x_labels <- c(xlim[[1]], "50", "100", "167", "200", "300", "400", "500", "600", "700", "800", "900", "1000")
      

      p <-  ggplot(x) +
        geom_line(aes(.data$insert_size, .data[[y_use]]), color = "grey2", size = 0.8, alpha = 0.8) +
        geom_vline(xintercept = vline, linetype = "dashed", color = "grey") + 
        coord_cartesian(xlim = xlim) +
        labs(x = "cfDNA Fragment Length (bp)", y = plot_type) +
        scale_x_continuous(limits = xlim, 
                           breaks = x_breaks, 
                           labels = x_labels) + 
        scale_y_continuous(expand = c(0, 0)) +
        theme_length_plot()

      
      
      return(p)
  
}
