#' Plot Motif Count or Fraction as Bar Plot
#'
#' Creates a bar plot to visualize motif counts or fractions from genomic data.
#' Uses `ggplot2` for plotting and `tidyr` and `dplyr` for data manipulation.
#'
#' @import ggplot2
#' @import tidyr
#' @import dplyr
#' @importFrom tibble as_tibble
#' @importFrom stringr str_to_lower str_extract str_to_title
#' @importFrom ggplot2 geom_col aes labs scale_fill_manual scale_y_continuous
#'
#' @param x Object containing the data to plot.
#' @param ylim Numeric vector of length 2 providing y axis limits.
#' @param x_title Character string, title for the x axis.
#' @param plot_type Character string specifying the type of plot: 'count' or 
#' 'fraction'.
#' @param bar_color Color to be used for bars in the plot.
#' @param motif_levels Character vector specifying the order and inclusion of
#' motif levels.
#' @param output_file Optional string specifying the path and filename 
#' where the plot should be saved. Include the file extension, such as '.png' or 
#' '.pdf'.
#' @param ggsave_params List of parameters to be passed to `ggplot2::ggsave()`.
#' This list can include any of the arguments that `ggsave()` accepts.
#' Example: `list(width = 10, height = 8, dpi = 300)`
#' @param ... Additional arguments passed on to `geom_bar()`.
#'
#' @return Returns a `ggplot2` object.
#' @export
#'
#' @examples
#' \dontrun{
#'   plot <- plotMotif(data, ylim = c(0, 20), x_title = "Motif",
#'                     plot_type = "count", bar_color = "blue")
#' }

plotMotif <- function(x,
                      ylim = NULL,
                      x_title = "Motif",
                      plot_type = c("Fraction"),
                      bar_color = c("A" = "cornflowerblue",
                                    "C" = "darksalmon",
                                    "G" = "palevioletred2",
                                    "T" = "forestgreen"),
                      motif_levels = c("C", "G", "A", "T"),
                      output_file = NULL,
                      ggsave_params = list(width = 16,
                                           height = 4,
                                           unit = "cm",
                                           bg = "white",
                                           dpi = 500),
                      ...) {

  # Conditional check for mutational data
  has_mutational_info <- attr(x, "mutational_info")

  if (is.null(has_mutational_info) || !has_mutational_info) {
    y_use <- NULL
    x <- tibble::as_tibble(x)
    plot_type <-  stringr::str_to_lower(plot_type)
    max_count <- max(dplyr::pull(x, .data$n))
    max_fraction <- max(dplyr::pull(x, .data$fraction))

    if (plot_type == c("fraction")) {

      y_use <- "fraction"

      if (missing(ylim)) {
        ylim <- c(0, max_fraction * 1.1)
        message("ylim", " was set as: ", ylim[1], " - ", ylim[2])
      }

    } else if (plot_type == c("n")) {

      y_use <- "n"

      if (missing(ylim)) {
        ylim <- c(0, max_count * 1.2)
        message("ylim", " was set as: ", ylim[1], " - ", ylim[2])
      }

    }

    x <- x %>%
      dplyr::mutate(group = stringr::str_extract(.data$motif,
                                                 pattern = "^[ATCG]")) %>%
      dplyr::filter(.data$group %in% motif_levels)

    x$group <- factor(x$group, levels = motif_levels)

    x <- x %>%
      dplyr::arrange(.data$group)

    x$motif <-  factor(x$motif, levels = x$motif)

    p <- ggplot(data = x,
                mapping = aes(x = .data$motif,
                              y = .data[[y_use]],
                              fill = .data$group)) +
      geom_col() +
      labs(x = x_title, y = stringr::str_to_title(plot_type)) +
      scale_fill_manual(values = bar_color) +
      scale_y_continuous(limits = ylim) +
      theme_motif_plot()

    # Check if 'output_file' is not NULL and save the plot to the specified path
    if (!is.null(output_file)) {
      do.call("ggsave", c(list(plot = p, filename = output_file), ggsave_params))
      message("Plot saved to ", output_file)
    }

    return(p)

  } else if (has_mutational_info) {
    plot_motif_mut(x, ylim, x_title, plot_type,
                   bar_color, motif_levels,
                   output_file, ggsave_params)
  }
}

theme_motif_plot <- function(){

  theme_classic() %+replace%
    theme(
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 8),
      axis.title = element_text(size = 14, face = "bold"),
      legend.position = "none"
    )
}

theme_motif_plot_mut <- function() {
  theme_classic() %+replace%
    theme(
      # X-axis text
      axis.text.x = element_text(angle = 45, size = 5),
      
      # Y-axis text
      axis.text.y = element_text(size = 5),
      
      # X-axis title
      axis.title.x = element_text(
        size = 5, 
        face = "bold",
        margin = margin(t = 3, r = 0, b = 0, l = 0)
      ),
      
      # Y-axis title
      axis.title.y = element_text(
        size = 5,
        face = "bold",
        angle = 90,
        margin = margin(t = 0, r = 3, b = 0, l = 0)
      ),
      
      # Legend title and text
      legend.title = element_text(size = 4),
      legend.text = element_text(size = 4),
      
      # Panel and plot background
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      
      # Legend key and margin
      legend.key.size = unit(0.15, "cm"),
      legend.margin = margin(t = 0, r = 0.1, b = 0, l = -0.3, unit = "cm")
    )
}
