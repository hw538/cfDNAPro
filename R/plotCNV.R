
#' plot copy number profile 
#' @import ggplot2
#' @import tibble
#' @import dplyr
#' @importFrom rlang has_name
#' @param x 
#' @param ylim 
#' @param chromosome 
#' @param point_color 
#' @param x_title 
#' @param y_title 
#' @param point_size 
#' @param point_alpha 
#' @param chr_edge_color 
#' @param chr_edge_line_size 
#' @param chr_edge_alpha 
#' @param chr_edge_type 
#' @param segment_color 
#' @param segment_alpha 
#' @param segment_line_end 
#' @param segment_line_size 
#' @param legend_position 
#' @param x_axis_expand 
#' @param y_axis_expand 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
plotCNV <- function(x, 
                    ylim = c(-2, 2), 
                    chromosome = c(seq(1, 22, 1), "X"), 
                    point_color = c("Loss" = "royalblue", 
                                    "Deletion" = "darkblue", 
                                    "Neutral" = "darkgrey", 
                                    "Amplification" = "darkorange", 
                                    "Gain" = "orange3"),
                    x_title = "Chromosome",
                    y_title = "Log2 Ratio",
                    point_size = 0.3, 
                    point_alpha = 0.9,
                    chr_edge_color = "black",
                    chr_edge_line_size = 0.2, 
                    chr_edge_alpha = 0.8,
                    chr_edge_type = "dotted",
                    segment_color = "red",
                    segment_alpha = 1, 
                    segment_line_end = "round",
                    segment_line_size = 0.75,
                    legend_position = "none",
                    x_axis_expand = c(0.1, 0.1),
                    y_axis_expand = c(0, 0),
                    ...) {
  
  names(point_color) <- stringr::str_to_title(names(point_color))
  
  if(!rlang::has_name(point_color, "Loss")){
    point_color["Loss"] <- "royalblue"
  }
  if(!rlang::has_name(point_color, "Deletion")){
    point_color["Deletion"] <- "darkblue"
  }
  if(!rlang::has_name(point_color, "Neutral")){
    point_color["Neutral"] <- "darkgrey"
  }
  if(!rlang::has_name(point_color, "Amplification")){
    point_color["Amplification"] <- "darkorange"
  }
  if(!rlang::has_name(point_color, "Gain")){
    point_color["Gain"] <- "orange3"
  }
   
  sample_name <- x@phenoData@data$name
  
  chr_levels <- as.character(c(seq(1, 22, 1), "X", "Y"))
  
  cnv_tibble <- tibble::tibble(
    chr = x@featureData@data$chromosome,
    start = x@featureData@data$start,
    end = x@featureData@data$end,
    use = x@featureData@data$use,
    copynumber = x@assayData$copynumber[, sample_name],
    seg = x@assayData$segmented[, sample_name],
    call = x@assayData$calls[, sample_name]
  )
  
  cnv_tibble$chr <- factor(cnv_tibble$chr, levels = chr_levels)
  
  
  cnv_tibble2 <- cnv_tibble %>%
    dplyr::filter(.data$chr %in% chromosome) %>%
    dplyr::filter(.data$use == TRUE) %>%
    dplyr::filter(is.finite(.data$copynumber)) %>%
    dplyr::mutate(call = dplyr::case_when(
      .data$call == -2 ~ "Deletion",
      .data$call == -1 ~ "Loss",
      .data$call == 0 ~ "Neutral",
      .data$call == 1 ~ "Gain",
      .data$call == 2 ~ "Amplification"
    )) %>%
    dplyr::mutate(x_index = dplyr::row_number()) %>%
    dplyr::mutate(logratio = log2(.data$copynumber)) %>%
    dplyr::mutate(logseg = log2(.data$seg)) %>%
    dplyr::mutate(seg_id = rleid_cfdnapro(.data$logseg))
  
  cnv_tibble2$call <- factor(cnv_tibble2$call, 
                             levels = c("Amplification", "Gain", "Neutral", "Loss", "Deletion"))

  
  # create chromosome boundries
  chr_edge <- cnv_tibble2 %>% 
    dplyr::group_by(chr) %>%
    dplyr::filter(row_number() == n()) %>%
    dplyr::pull(x_index)
  
  chr_edge <- c(1, chr_edge)
  
  x_breaks <- cnv_tibble2 %>% 
    dplyr::group_by(chr) %>%
    dplyr::filter(row_number() == ceiling(n()/2)) %>%
    dplyr::pull(x_index)
  
  x_labels <- as.character(unique(cnv_tibble2$chr))
    
  # ggplot2 plot 
  p <- ggplot(data = cnv_tibble2) +
    # copy number points
    geom_point(aes(x = .data$x_index, y = .data$logratio, color = .data$call),
               size = point_size, 
               alpha = point_alpha) +
    scale_color_manual(values = point_color) +
    # segmentation line
    geom_line(aes(x = .data$x_index, y = .data$logseg, group = .data$seg_id), 
              colour = segment_color, 
              alpha = segment_alpha, 
              size = segment_line_size,
              lineend = segment_line_end) +
    # chr edges 
    geom_vline(xintercept = chr_edge, 
               linetype = chr_edge_type, 
               color = chr_edge_color, 
               size = chr_edge_line_size) +
    # chr names 
    scale_x_continuous(
                       breaks = x_breaks, 
                       labels = x_labels,
                       expand = x_axis_expand) + 
    scale_y_continuous(limits = ylim,
                       expand = y_axis_expand) +
    labs(x = x_title, y = y_title) +
    theme_cnv_plot()  +
    theme(legend.position = legend_position)
  
  return(p)
  
}
