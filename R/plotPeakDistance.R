#' Plot the inter-peak distance of fragment size distance distribution
#'
#' @importFrom magrittr %>%
#' @import ggplot2
#' @import stringr
#' @import dplyr
#' @importFrom rlang .data
#' @param x A long-format dataframe contains the inter-peak distance,
#'    a template please refer to the result of 'callPeakDistance' function.
#' @param order The groups show in the final plot,
#'    the input value should be vector,
#'    e.g. `groups = c('group1','group2')`,
#'  default is all folders in the folder path.
#' @param summarized Logical value,
#'    describe whether the x is summarzied already.
#'    summarized means the count and proportion of each interpeak_dist.
#' @param type The plot type, default is line plot, now only support line plot.
#'    Don't change this parameter in this version, keep it as default.
#' @param mincount Minimum count value of inter peak distance,
#'    count number less than this value will be removed first,
#'    then proportion of each count value will be calculated.
#'    Default value is 0.
#' @param xlim The x axis range shown in the plot. Default is `c(8,13)`.
#' @param ... Further arguments passed to or from other methods.
#' @return  The function returns the line plot of inter peak distance.
#'
#' @examples
#' # Get the path to example data.
#' path <- examplePath("groups_picard")
#'
#' # Calculate the inter-peak distance.
#' df <- callPeakDistance(path = path)
#'
#' # Plot the inter-peak distance.
#' plot <- plotPeakDistance(df,
#'     xlim = c(8, 13),
#'     mincount = 2
#' )
#' @author Haichao Wang
#'
#' @export

plotPeakDistance <- function(x,
                            summarized,
                            order,
                            type,
                            mincount,
                            xlim,
                            ...) {
    group <- prop <- interpeak_dist <- NULL

    if (missing(summarized)) {
          summarized <- FALSE
      }

    if (missing(order)) {
          order <- as.vector(unique(x$group))
      }

    if (missing(type)) {
          type <- "line"
      }

    if (missing(mincount)) {
        mincount <- 0
        cat("setting the mincount to 0.\n ")
    }

    if (missing(xlim)) {
        xlim <- c(7, 13)
        cat("setting the xlim to c(7,13). \n ")
    }
    # Summarize interpeak dist if summarized = FALSE.
    if (!summarized) {
        summary <-
            x %>%
            group_by(.data$group, .data$interpeak_dist) %>%
            summarise(count = n()) %>%
            na.omit() %>%
            filter(count >= mincount) %>%
            mutate(prop = count / sum(count)) %>%
            ungroup() %>%
            filter(.data$group %in% as.vector(order))
        # Give x to summary if summarized = TRUE.
    } else if (summarized) {
        summary <- x
    }
    # Set the factor level of group.
    summary$group <- factor(summary$group, levels = order)
    # Interpeak line plot.
    interpeak_dist_plot <- ggplot() +
        geom_point(
            data = summary,
            aes(x = interpeak_dist, y = prop),
            color = "gray30",
            size = 3,
            alpha = 0.8
        ) +
        geom_line(
            data = summary,
            aes(
                x = interpeak_dist,
                y = prop,
                color = group
            ),
            size = 2
        ) +
        labs(
            x = "Interpeak distance",
            y = "Percent"
        ) +
        theme_classic() +
        geom_vline(
            xintercept = 10,
            linetype = "dashed",
            color = "gray30",
            size = 0.6
        ) +
        geom_vline(
            xintercept = 11,
            linetype = "dashed",
            color = "gray30",
            size = 0.6
        ) +
        scale_x_continuous(
            limits = xlim,
            breaks = c(xlim[1]:xlim[2]),
            labels = c(xlim[1]:xlim[2])
        )

    return(interpeak_dist_plot)
}
