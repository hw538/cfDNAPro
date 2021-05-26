#' Plot the raw fragment size metrics (e.g. proportion, cdf and 1-cdf) of
#' all groups with different colors in a single plot
#'
#' @importFrom magrittr %>%
#' @import dplyr
#' @import ggplot2
#' @param x A long-format dataframe contains the metrics of different cohort.
#' @param order The groups show in the final plot,
#'   the input value should be vector,
#'   e.g. `groups=c('group1','group2')`,
#'   default is all folders in the folder path.
#' @param plot The plot type, default is 'all' which means
#'   all of proportion, cdf and 1-cdf plots will be shown.
#' @param vline Vertical dashed lines, default value is `c(81,167)`.
#' @param xlim The x axis range shown in the plot. Default is `c(0,500)`.
#' @param ylim The y axis range shown in the fraction of fragment size plots.
#'   Default is `c(0,0.035)`.
#' @param ... Further arguments passed to or from other methods.
#' @return  The function returns a list plots.
#'
#' @examples
#' # Get the path to example data.
#' path <- examplePath("groups_picard")
#' # Calculate the sizes.
#' df <- callSize(path = path)
#' # Plot all samples from multiple groups into one figure.
#' plot <- plotAllToOne(df)
#' @author Haichao Wang
#' @export
plotAllToOne <- function(x, order, plot, vline, xlim, ylim, ...) {
    group <-
        insert_size <- prop <- file_name <- cdf <- one_minus_cdf <- NULL

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

    # Setting ylim for proportion plots.
    if (missing(ylim)) {
          ylim <- c(0, 0.035)
      }

    # Filtering the groups.
    x <- filter(x, group %in% as.vector(order))

    # Setting the factor level of the groups.
    x$group <- factor(x$group, levels = order)

    # Prop plot.
    prop_plot <-
        ggplot(data = x, aes(x = insert_size, y = prop, color = group)) +
        geom_line(size = 1, aes(group = file_name)) +
        xlim(xlim) +
        ylim(ylim) +
        xlab("Fragment size (bp)") +
        ylab("Proportion of reads") +
        theme_classic()

    # Cdf plot.
    cdf_plot <-
        ggplot(data = x, aes(x = insert_size, y = cdf, color = group)) +
        geom_line(size = 1, aes(group = file_name)) +
        xlim(xlim) +
        xlab("Fragment size (bp)") +
        ylab("Cumulative proportion of reads") +
        theme_classic()
    # 1-cdf plot.
    one_minus_cdf_plot <-
        ggplot(data = x, aes(x = insert_size, 
                             y = one_minus_cdf, 
                             color = group)) +
        geom_line(size = 1, aes(group = file_name)) +
        xlim(xlim) +
        xlab("Fragment size (bp)") +
        ylab("Cumulative proportion of reads > fragment size") +
        theme_classic()

    # Add vertical lines to the plot.
    for (i in c(seq_len(length(vline)))) {
        prop_plot <-
            prop_plot + geom_vline(
                xintercept = vline[i],
                linetype = "dashed",
                color = "red",
                size = 0.6
            )

        cdf_plot <-
            cdf_plot + geom_vline(
                xintercept = vline[i],
                linetype = "dashed",
                color = "red",
                size = 0.6
            )

        one_minus_cdf_plot <-
            one_minus_cdf_plot + geom_vline(
                xintercept = vline[i],
                linetype = "dashed",
                color = "red",
                size = 0.6
            )
    }

    output_all <-
        list(
            prop_plot = prop_plot,
            cdf_plot = cdf_plot,
            one_minus_cdf_plot = one_minus_cdf_plot
        )

    if (plot == "all") {
        return(output_all)
    } else if (plot == "prop") {
        return(prop_plot)
    } else if (plot == "cdf") {
        return(cdf_plot)
    } else if (plot == "1-cdf") {
        return(one_minus_cdf_plot)
    }
}
