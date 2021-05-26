#' Plot the raw fragment size metrics of single group in a single plot,
#' colored by samples.
#'
#' @importFrom magrittr %>%
#' @import ggplot2
#' @import stringr
#' @import dplyr
#' @param x A long-format dataframe contains the metrics of different cohort.
#' @param order The groups show in the final plot,
#'     the input value should be vector,
#'     e.g. `order = c('group1')``,
#'     default is all groups/cohorts in the folder path.
#' @param plot The plot type, default is 'all' which means both proportion,
#'     cdf and 1-cdf plots will be shown.
#' @param vline Vertical dashed lines, default value is `c(81,167)`.
#' @param xlim The x axis range shown in the plot. Default is `c(0,500)`.
#' @param ylim The y axis range shown in the fraction of fragment size plots.
#'     Default is `c(0,0.035)`.
#' @param ... Further arguments passed to or from other methods.
#' @return  The function returns a list plots.
#' @examples
#' # Get the path to example data.
#' path <- examplePath("groups_picard")
#'
#' # Calculate the metrics.
#' df <- callMetrics(path = path)
#'
#' # Plot the only the group specified..
#' plot <- plotSingleGroup(x = df, order = c("cohort_1"))
#' @author Haichao Wang
#'
#' @export

plotSingleGroup <-
    function(x, xlim, ylim, vline, order, plot, ...) {
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
              ylim <- c(0, 0.03)
          }

        # Filtering the groups.
        x <- filter(x, group %in% as.vector(order))

        # Setting the factor level of the groups.
        x$group <- factor(x$group, levels = order)
        # Prop plot.
        prop_plot <-
            ggplot(
                data = x,
                aes(x = insert_size, y = prop, color = file_name)
            ) +
            geom_line(size = 1) +
            xlim(xlim) +
            ylim(ylim) +
            xlab("Fragment size (bp)") +
            ylab("Proportion of reads") +
            theme_classic() +
            theme(
                axis.text = element_text(size = 10),
                axis.title = element_text(size = 14, face = "bold"),
                legend.position = "none"
            )
        # Cdf plot
        cdf_plot <-
            ggplot(
                data = x,
                aes(x = insert_size, y = cdf, color = file_name)
            ) +
            geom_line(size = 1) +
            xlim(xlim) +
            xlab("Fragment size (bp)") +
            ylab("Cumulative proportion of reads") +
            theme_classic() +
            theme(
                axis.text = element_text(size = 10),
                axis.title = element_text(size = 14, face = "bold"),
                legend.position = "none"
            )
        # 1-cdf plot
        one_minus_cdf_plot <-
            ggplot(
                data = x,
                aes(x = insert_size, y = one_minus_cdf, color = file_name)
            ) +
            geom_line(size = 1) +
            xlim(xlim) +
            xlab("Fragment size (bp)") +
            ylab("Cumulative proportion of reads > fragment size") +
            theme_classic() +
            theme(
                axis.text = element_text(size = 10),
                axis.title = element_text(size = 14, face = "bold"),
                legend.position = "none"
            )

        # Add vertical lines to the plot
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

        return(
            list(
                prop_plot = prop_plot,
                cdf_plot = cdf_plot,
                one_minus_cdf_plot = one_minus_cdf_plot
            )
        )
    }
