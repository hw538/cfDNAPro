
plotLength <- function(x,
                       plot_type = c("Fraction"),
                       xlim = c(30, 500),
                       ylim  ,
                       x_breaks,
                       x_labels,
                       vline,
                       add_vline = TRUE,
                       line_color = "grey1",
                       line_size = 0.6,
                       line_alpha = 0.8,
                       vline_color = "grey",
                       vline_type = "dashed",
                       vline_size = 0.7,
                       area_highlight = NULL,
                       ...
                       ){
  
  
  
  y_use <- NULL
  x <- tibble::as_tibble(x)
  plot_type <-  stringr::str_to_lower(plot_type)
  max_count <- max(dplyr::pull(x, All_Reads.fr_count))
  max_fraction <- max(dplyr::pull(x, prop ))

  
  if(plot_type == c("fraction")) {
    
    y_use <- "prop"
    
    if(missing(ylim)) {
    ylim <- c(0, max_fraction * 1.1)
    message("ylim", " was set as: ", ylim[1]  ," - ", ylim[2] )
    }
    
    
  } else if (plot_type == c("count")) {
    y_use <- "All_Reads.fr_count"
    
    if(missing(ylim)) {
    ylim <- c(0, max_count * 1.2)
    
    message("ylim", " was set as: ", ylim[1]  ," - ", ylim[2] )
    
    }
    
  }
    
  x$insert_size <- as.numeric(x$insert_size)
  plot_type <- stringr::str_to_title(plot_type)
  
  # find the modal length 
    peak_isize_breaks <- dplyr::filter(x, !!rlang::sym(y_use) == max(x[[y_use]])) %>%
      dplyr::pull(.data$insert_size)
    
  # define the vline as modal length
  if(missing(vline)){
    vline <- as.vector(peak_isize_breaks)
  }
  
  # set the x breaks 
  if(missing(x_breaks) & missing(x_labels)) {
    ref_breaks <- seq(100, 1000, by = 100)
    
    i <- which(ref_breaks > xlim[1] & ref_breaks < xlim[2])
    middle_isize_breaks <- ref_breaks[i]
    
    x_breaks <- c(xlim[1], middle_isize_breaks, peak_isize_breaks, xlim[2]) %>% 
      sort()
    x_labels <- as.character(x_breaks)
  }
  
  if(!missing(x_breaks) & missing(x_labels)) {
    x_labels <- as.character(x_breaks)
  }
  
  if(missing(x_breaks) & !missing(x_labels)) {
    x_breaks <- as.numeric(x_labels)
  }
  
  # plot
  p <-  ggplot(data = x, aes(.data$insert_size, .data[[y_use]])) +
    geom_line(
              color = line_color, 
              size = line_size, 
              alpha = line_alpha)
  
  if(add_vline) {
    p <- p +
    geom_vline(xintercept = vline, linetype = vline_type, color = vline_color, size = vline_size)
  }
  
  
  if(is.list(area_highlight)){
    
   plot_length_add_hightlight_area <- function(params, plot) {
     
    if (rlang::has_name(params, "fill")) {
     fill <- params[["fill"]] 
    } else {
      fill <- "lightgrey"
    }
     
    if (rlang::has_name(params, "alpha")) {
     alpha <- params[["alpha"]]
    } else {
     alpha <- 0.7 
    }
     
     left_boundary <- params[["range"]][[1]] %>% as.numeric()
     right_boundary <- params[["range"]][[2]] %>% as.numeric()
    
     p <- plot +
       geom_area(
         mapping = aes(x = ifelse( .data$insert_size >= left_boundary & .data$insert_size <= right_boundary , insert_size, 0)), 
         fill = fill,
         alpha  = alpha)
     
     return(p)
   }
   
   for (area in area_highlight){
     p <- plot_length_add_hightlight_area(params = area, plot = p)
   }
   
    
  }
  
  p <- p + 
    coord_cartesian(xlim = xlim) +
    labs(x = "cfDNA Fragment Length (bp)", y = stringr::str_to_title(plot_type)) +
    scale_x_continuous(limits = xlim, 
                       breaks = x_breaks, 
                       labels = x_labels) + 
    scale_y_continuous(limits = ylim) +
    theme_length_plot()
  
      
  return(p)
  
}


