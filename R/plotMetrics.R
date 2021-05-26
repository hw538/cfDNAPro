#' Plot the fragment size metrics (i.e. proportion, cdf and 1-cdf)
#'
#' @importFrom magrittr %>%
#' @import dplyr
#' @import ggplot2
#' @param x A long-format dataframe contains the metrics of different cohort.
#' @param order The groups show in the final plot,
#'   the input value should be vector, e.g. `groups = c('group1','group2')``,
#'   default is all folders in the folder path
#' @param plot The plot type,
#'     default is 'all': both median and mean metrics will be shown.
#'     They will include: mean_prop, mean_cdf, mean_1-cdf,
#'     median_prop, median_cdf, median_1-cdf.
#'     Could also specify as "median" or "mean".
#' @param vline Vertical dashed lines, default value is c(81,167).
#' @param xlim The x axis range shown in the plot. Default is c(0,500).
#' @param ylim The y axis range shown
#'   in the fraction of fragment size plots. Default is c(0,0.0125).
#' @param ... Further arguments passed to or from other methods.
#' @return  The function returns a list plots.
#'
#' @examples
#' # Get the path to example data.
#' path <- examplePath("groups_picard")
#' # Calculate the metrics.
#' df <- callMetrics(path = path)
#' # Plot metrics.
#' plot <- plotMetrics(df,
#'     plot = "median",
#'     order = c("cohort_1", "cohort_2")
#' )
#' @author Haichao Wang
#'
#' @export


plotMetrics <- function(x, order, plot, vline, xlim, ylim, ...) {
    group <- insert_size <- NULL
    prop_median <- cdf_median <- one_minus_cdf_median <- NULL
    prop_mean <- cdf_mean <- one_minus_cdf_mean <- NULL

    # Setting vline.
    if (missing(vline)) {
          vline <- c(81, 167)
      }

    # Setting the plot.
    if (missing(plot)) {
          plot <- "all"
      }

    # Setting the orders of groups in the plot.
    if (missing(order)) {
          order <- as.vector(unique(x$group))
      }

    # Setting xlim.
    if (missing(xlim)) {
          xlim <- c(0, 500)
      }

    # setting ylim for proportion plots.
    if (missing(ylim)) {
          ylim <- c(0, 0.0125)
      }

    # Filtering the groups.
    x <- filter(x, group %in% as.vector(order))

    # Setting the factor level of the groups.
    x$group <- factor(x$group, levels = order)

    # median ------------------------------------------------------------------


    # Generate the median prop plots.
    median_prop_plot <-
        ggplot2::ggplot(
            data = x,
            aes(x = insert_size, y = prop_median, color = group)
        ) +
        geom_line(size = 1) +
        xlim(xlim) +
        ylim(ylim) +
        xlab("Fragment size (bp)") +
        ylab("Median proportion of reads") +
        theme_classic()
    # Generate the median cdf plots.
    median_cdf_plot <-
        ggplot2::ggplot(
            data = x,
            aes(
                x = insert_size,
                y = cdf_median,
                color = group
            )
        ) +
        geom_line(size = 1) +
        xlim(xlim) +
        xlab("Fragment size (bp)") +
        ylab("Median cumulative fraction of reads") +
        theme_classic()
    # Generate the median 1-cdf plots.
    median_one_minus_cdf_plot <-
        ggplot2::ggplot(
            data = x,
            aes(
                x = insert_size,
                y = one_minus_cdf_median,
                color = group
            )
        ) +
        geom_line(size = 1) +
        xlim(xlim) +
        xlab("Fragment size (bp)") +
        ylab("Median cumulative fraction of reads > Fragment size") +
        theme_classic()

    # mean --------------------------------------------------------------------

    # Generate the mean prop plots.
    mean_prop_plot <-
        ggplot(
            data = x,
            aes(x = insert_size, y = prop_mean, color = group)
        ) +
        geom_line(size = 1) +
        xlim(xlim) +
        ylim(ylim) +
        xlab("Fragment size (bp)") +
        ylab("Mean proportion of reads") +
        theme_classic()
    # Generate the mean cdf plots.
    mean_cdf_plot <-
        ggplot(
            data = x,
            aes(x = insert_size, y = cdf_mean, color = group)
        ) +
        geom_line(size = 1) +
        xlim(xlim) +
        xlab("Fragment size (bp)") +
        ylab("Mean cumulative fraction of reads") +
        theme_classic()
    # Generate the mean 1-cdf plots.
    mean_one_minus_cdf_plot <-
        ggplot(
            data = x,
            aes(x = insert_size, y = one_minus_cdf_mean, color = group)
        ) +
        geom_line(size = 1) +
        xlim(xlim) +
        xlab("Fragment size (bp)") +
        ylab("Mean cumulative fraction of reads > Fragment size") +
        theme_classic()

    # Add vertical lines to the plot.
    for (i in c(seq_len(length(vline)))) {
        median_prop_plot <-
            median_prop_plot + geom_vline(
                xintercept = vline[i],
                linetype = "dashed",
                color = "red",
                size = 0.6
            )

        median_cdf_plot <-
            median_cdf_plot + geom_vline(
                xintercept = vline[i],
                linetype = "dashed",
                color = "red",
                size = 0.6
            )

        median_one_minus_cdf_plot <-
            median_one_minus_cdf_plot + geom_vline(
                xintercept = vline[i],
                linetype = "dashed",
                color = "red",
                size = 0.6
            )

        mean_prop_plot <-
            mean_prop_plot + geom_vline(
                xintercept = vline[i],
                linetype = "dashed",
                color = "red",
                size = 0.6
            )

        mean_cdf_plot <-
            mean_cdf_plot + geom_vline(
                xintercept = vline[i],
                linetype = "dashed",
                color = "red",
                size = 0.6
            )

        mean_one_minus_cdf_plot <-
            mean_one_minus_cdf_plot + geom_vline(
                xintercept = vline[i],
                linetype = "dashed",
                color = "red",
                size = 0.6
            )
    }

    # Return the plots in a list.
    output_all <-
        list(
            median_prop_plot = median_prop_plot,
            median_cdf_plot = median_cdf_plot,
            median_one_minus_cdf_plot = median_one_minus_cdf_plot,
            mean_prop_plot = mean_prop_plot,
            mean_cdf_plot = mean_cdf_plot,
            mean_one_minus_cdf_plot = mean_one_minus_cdf_plot
        )

    output_mean <-
        list(
            mean_prop_plot = mean_prop_plot,
            mean_cdf_plot = mean_cdf_plot,
            mean_one_minus_cdf_plot = mean_one_minus_cdf_plot
        )

    output_median <-
        list(
            median_prop_plot = median_prop_plot,
            median_cdf_plot = median_cdf_plot,
            median_one_minus_cdf_plot = median_one_minus_cdf_plot
        )

    if (plot == "all") {
        return(output_all)
    } else if (plot == "median") {
        return(output_median)
    } else if (plot == "mean") {
        return(output_mean)
    }
}
