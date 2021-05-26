#' Plot the inter-valley distance of fragment size distance distribution
#'
#' @importFrom magrittr %>%
#' @import ggplot2
#' @import stringr
#' @import dplyr
#' @param x A long-format dataframe contains the inter-valley distance,
#'   a template please refer to the result of 'callValleyDistance' function.
#' @param order The groups show in the final plot,
#'   the input value should be vector,
#'   e.g. `groups=c('group1','group2')`,
#'   default is all folders in the folder path.
#' @param type The plot type, default is line plot, now only support line plot.
#'   Don't change this parameter in this version, keep it as default.
#' @param mincount Minimum count value of inter valley distance,
#'   count number less than this value will be removed first,
#'   then proportion of each count value will be calculated. Default value is 0.
#' @param xlim The x axis range shown in the plot. Default is c(8,13).
#' @param ... Further arguments passed to or from other methods.
#' @return  The function returns the line plot of inter valley distance.
#'
#' @examples
#' # Get the path to example data.
#' path <- examplePath("groups_picard")
#'
#' # Calculate the inter-valley distance.
#' df <- callValleyDistance(path = path)
#'
#' # Plot the inter-valley distance.
#' plot <- plotValleyDistance(df,
#'     xlim = c(8, 13),
#'     mincount = 2
#' )
#' @author Haichao Wang
#' @export

plotValleyDistance <-
    function(x, order, type, mincount, xlim, ...) {
        group <- intervalley_dist <- prop <- NULL

        if (missing(order)) {
              order <- as.vector(unique(x$group))
          }

        if (missing(type)) {
              type <- "line"
          }

        if (missing(mincount)) {
            mincount <- 0
            cat("setting the mincount to 0. \n ")
        }

        if (missing(xlim)) {
            xlim <- c(7, 13)
            cat("setting the xlim to c(7,13). \n ")
        }
        # Summarize intervalley dist.
        summary <-
            x %>%
            group_by(group, intervalley_dist) %>%
            summarise(count = n()) %>%
            na.omit() %>%
            filter(count >= mincount) %>%
            mutate(prop = count / sum(count)) %>%
            ungroup() %>%
            filter(group %in% as.vector(order))

        # Set the factor level of group.
        summary$group <- factor(summary$group, levels = order)

        # Intervalley line plot.
        intervalley_dist_plot <-
            ggplot() +
            geom_point(
                data = summary,
                aes(x = intervalley_dist, y = prop),
                color = "gray30",
                size = 3,
                alpha = 0.8
            ) +
            geom_line(
                data = summary,
                aes(x = intervalley_dist, y = prop, color = group),
                size = 2
            ) +
            labs(x = "Intervalley distance", y = "Percent") +
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

        return(intervalley_dist_plot)
    }
