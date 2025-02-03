#' Plot fragment length profile
#' @import ggplot2
#' @import stringr
#' @import tibble
#' @importFrom rlang has_name
#' @importFrom purrr map
#'
#' @param x Data to be plotted.
#' @param plot_type Type of plot to generate (e.g., 'count' or 'fraction').
#' @param xlim Limits for the x-axis.
#' @param ylim Limits for the y-axis.
#' @param x_breaks Breaks for the x-axis.
#' @param x_labels Labels for the x-axis.
#' @param vline Position of vertical lines.
#' @param add_vline Logical, whether to add vertical lines.
#' @param line_color Color of the plot lines.
#' @param line_size Size of the plot lines.
#' @param line_alpha Transparency of the plot lines.
#' @param vline_color Color of the vertical lines.
#' @param vline_type Type of the vertical lines (e.g., 'dashed').
#' @param vline_size Size of the vertical lines.
#' @param area_highlight List specifying areas to highlight; each list element
#'   is a sublist specifying 'range', 'fill', and 'alpha'.
#' @param output_file File path and name for saving the output plot. Ensure the
#'   directory exists and is writable.
#' @param ggsave_params A list of parameters for ggplot2::ggsave(). May include
#'   'width', 'height', 'dpi', etc. Use named elements in the list to specify.
#' @param ... Additional parameters for the plotting function.
#'
#' @return A ggplot2 object that can be further modified or saved.
#' @export
#'
#' @examples
#' # Assuming 'data' is a dataframe with 'insert_size' and 'value' columns:
#' plotLength(data, plot_type = "fraction", xlim = c(100, 500), ylim = c(0, 1),
#'            x_breaks = seq(100, 500, by = 50), add_vline = TRUE, vline = 150)
plotLength <- function(x,
                       plot_type = c("Fraction"),
                       xlim = c(30, 500),
                       ylim = NULL,
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
                       output_file = NULL,
                       ggsave_params = list(width = 10,
                                            height = 11,
                                            dpi = 500,
                                            bg = "white",
                                            unit = "cm"),
                       ...
                       ) {

  # Conditional check for mutational data
  has_mutational_info <- attr(x, "mutational_info")

  if (is.null(has_mutational_info) || !has_mutational_info) {
    y_use <- NULL
    x <- tibble::as_tibble(x)
    plot_type <-  stringr::str_to_lower(plot_type)
    max_count <- max(dplyr::pull(x, All_Reads.fr_count))
    max_fraction <- max(dplyr::pull(x, prop))


    if (plot_type == c("fraction")) {

      y_use <- "prop"

      if (missing(ylim)) {
        ylim <- c(0, max_fraction * 1.1)
        message("ylim", " was set as: ", ylim[1], " - ", ylim[2])
      }

    } else if (plot_type == c("count")) {
      y_use <- "All_Reads.fr_count"

      if (missing(ylim)) {
        ylim <- c(0, max_count * 1.2)

        message("ylim", " was set as: ", ylim[1], " - ", ylim[2])

      }

    }

    x$insert_size <- as.numeric(x$insert_size)
    plot_type <- stringr::str_to_title(plot_type)

    # find the modal length
    peak_isize_breaks <- dplyr::filter(x, !!rlang::sym(y_use) == max(x[[y_use]])) %>%
      dplyr::pull(.data$insert_size)

    # define the vline as modal length
    if (missing(vline)) {
      vline <- as.vector(peak_isize_breaks)
    }

    # set the x breaks
    if (missing(x_breaks) & missing(x_labels)) {
      ref_breaks <- seq(100, 1000, by = 100)

      i <- which(ref_breaks > xlim[1] & ref_breaks < xlim[2])
      middle_isize_breaks <- ref_breaks[i]

      x_breaks <- c(xlim[1],
                    middle_isize_breaks,
                    peak_isize_breaks,
                    xlim[2]) %>% sort()
      x_labels <- as.character(x_breaks)
    }

    if (!missing(x_breaks) && missing(x_labels)) {
      x_labels <- as.character(x_breaks)
    }

    if (missing(x_breaks) && !missing(x_labels)) {
      x_breaks <- as.numeric(x_labels)
    }

    # plot
    p <- ggplot(data = x, aes(.data$insert_size, .data[[y_use]])) +
      geom_line(color = line_color,
                size = line_size,
                alpha = line_alpha)

    if (add_vline) {
      p <- p +
        geom_vline(xintercept = vline,
                   linetype = vline_type,
                   color = vline_color,
                   size = vline_size)
    }

    if (is.list(area_highlight)) {

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
          geom_area(mapping = aes(x = ifelse( .data$insert_size >= left_boundary & .data$insert_size <= right_boundary,
                                             insert_size,
                                             0)),
                    fill = fill,
                    alpha  = alpha)

        return(p)
      }

      for (area in area_highlight) {
        p <- plot_length_add_hightlight_area(params = area, plot = p)
      }


      sum_hightlight_area <- function(params, length_tibble) {

        left_boundary <- params[["range"]][[1]] %>% as.numeric()
        right_boundary <- params[["range"]][[2]] %>% as.numeric()

        x_filtered <- length_tibble %>%
          dplyr::filter(.data$insert_size >= left_boundary) %>%
          dplyr::filter(.data$insert_size <= right_boundary)

        if (plot_type == c("fraction")) {
          y_use <- "prop"
        } else if (plot_type == c("count")) {
          y_use <- "All_Reads.fr_count"
        }

        area_sum <- sum(x_filtered[[y_use]])

        params["area_sum"] <- area_sum
        return(params)

      }

      area_highlight <- purrr::map(area_highlight,
                                   sum_hightlight_area,
                                   length_tibble = x)

    }

    p <- p +
      coord_cartesian(xlim = xlim) +
      labs(x = "cfDNA Fragment Length (bp)",
           y = stringr::str_to_title(plot_type)) +
      scale_x_continuous(limits = xlim,
                         breaks = x_breaks,
                         labels = x_labels) +
      scale_y_continuous(limits = ylim) +
      theme_length_plot()
    
    # Check if 'output_file' is not NULL and save the plot to the specified path
    if (!is.null(output_file)) {
      do.call("ggsave", c(list(plot = p, filename = output_file), ggsave_params))
      message("Plot saved to ", output_file)

    }
    return(p)

  } else if (has_mutational_info) {
    plot_length_mut(x, ylim, output_file, ggsave_params)
  }
}

