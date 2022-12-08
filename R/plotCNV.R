
plotCNV <- function(x, ylim = c(-2, 2)) {
  
  
  sample_name <- x@phenoData@data$name
  
  cnv_tibble <- tibble::tibble(
    chr = x@featureData@data$chromosome,
    start = x@featureData@data$start,
    end = x@featureData@data$end,
    use = x@featureData@data$use,
    copynumber = x@assayData$copynumber[, sample_name],
    seg = x@assayData$segmented[, sample_name],
    call = x@assayData$calls[, sample_name]
  )
  
  
  
  cnv_tibble2 <- cnv_tibble %>%
    dplyr::filter(chr %in% seq(1, 22, by = 1)) %>%
    dplyr::mutate(x_index = dplyr::row_number()) %>%
    dplyr::mutate(logratio = log2(copynumber)) %>%
    dplyr::mutate(logseg = log2(seg)) %>%
    dplyr::mutate(color = case_when(
      call == -1 ~ "lightred",
      call == 0 ~ "grey2",
      call == 1 ~ "lightorange",
      call == -2 ~ "red",
      call == 2 ~ "orange"
    ))
  
  cnv_tibble2$call <- forcats::as_factor(cnv_tibble2$call)
  
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
    
  
  p <- ggplot(data = cnv_tibble2) +
    geom_point(
               aes(x = .data[["x_index"]], 
                   y = .data[["logratio"]], 
                   color = .data[["call"]]),
               size = 0.35) +
    scale_color_manual(values = c("-1" = "royalblue", 
                                 "-2" = "darkblue", 
                                 "0" = "grey", 
                                 "2" = "darkorange", 
                                 "1" = "orange3")) +
    geom_line(aes(x = .data[["x_index"]], y = .data[["logseg"]], group = 1)) +
    geom_vline(xintercept = chr_edge, linetype = "dotted", color = "black", size = 0.2) +
    labs(x = "Chromosome", y = "Log2Ratio") +
    scale_x_continuous( 
                       breaks = x_breaks, 
                       labels = x_labels) + 
    theme_classic() +
    theme(axis.ticks.x = element_blank(),
          axis.line.x = element_blank())
  
  p
  return(p)
  
}
