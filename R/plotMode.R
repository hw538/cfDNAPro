#' Plot mode fragment size
#'
#' @importFrom magrittr %>%
#' @import ggplot2
#' @import stringr
#' @import dplyr
#' @param x A long-format dataframe contains the interpeak distance,
#'   a template please refer to the result of "callPeakdist" function.
#' @param order The groups show in the final plot,
#'  the input value should be vector, e.g. `groups = c("group1","group2")`,
#'  default is all folders in the folder path.
#' @param type The plot type, could choose "bin" or "stacked" chart.
#'   Default is bin plot.
#' @param mincount Minimum count of mode fragment size that will be included.
#'   Count number smaller than this value will be removed first,
#'   then proportion of each count value will be calculated.
#'   Default value is 0.
#' @param hline The horizontal lines added to the bin plot.
#'   Default lines will be `c(81,112,170)`.
#' @param ... Further arguments passed to or from other methods.
#' @return  The function returns the plot.
#'
#' @examples
#' # Get the path to example data.
#' path <- examplePath("groups_picard")
#' # Calculate the modes.
#' df <- callMode(path = path)
#' # Plot modes.
#' plot <- plotMode(df, hline = c(80, 111, 170))
#' @author Haichao Wang
#' @export


plotMode <- function(x, order, type, mincount, hline, ...) {
    group <- insert_size <- x_index <- file_name <- prop <- NULL

    if (missing(order)) order <- as.vector(unique(x$group))
    if (missing(type)) type <- "bin"

    if (missing(mincount)) {
        mincount <- 0
        cat("setting default mincount as 0.\n")
    }

    if (missing(hline)) {
        hline <- c(81, 112, 170)
        cat("setting default horizontal lines: y = 81, 112, 170. \n")
    }

    # Generate bin plot.
    if (type == "bin") {
        if (missing(order)) {
            order <- as.vector(unique(x$group))
        }

        # Sort the dataframe based on the order vector.
        df <- x %>%
            filter(group %in% as.vector(order)) %>%
            # Group by "group" for sort
            group_by(group) %>%
            # Sort by insert_size
            arrange(insert_size, .by_group = TRUE) %>%
            # Change group order which will alfter the order of groups.
            arrange(factor(group, levels = order)) %>%
            ungroup() %>%
            mutate(x_index = row_number())
        # Set the levels of "group".
        df$group <- factor(df$group, levels = order)

        # Start the plot.
        p1 <- ggplot(df, aes(x = x_index, y = insert_size, fill = group)) +
            geom_bar(stat = "identity", color = "gray25", show.legend = FALSE) +
            theme(
                plot.background = element_blank(),
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
                axis.ticks = element_blank(),
                # axis.text = element_blank(),
                axis.title = element_blank()
            ) +
            geom_vline(xintercept = 0, 
                       colour = "gray15", 
                       linetype = 1, 
                       size = 2) +
            geom_hline(yintercept = 0, 
                       colour = "gray15", 
                       linetype = 1, 
                       size = 2)

        # Add horizontal lines to the plot.
        for (i in c(seq_len(length(hline)))) {
            p1 <- p1 +
                geom_hline(
                    yintercept = hline[i],
                    linetype = "dashed",
                    color = "gray30",
                    size = 0.6
                )
        }

        # Calculate breaks.
        breaks <- vector()
        # Set the first value in breaks.
        if (length(order) == 1) {
            breaks[1] <- nrow(filter(df, group == order[1])) / 2
        } else if (length(order) > 1) {
            breaks[1] <- nrow(filter(df, group == order[1])) / 2
            for (i in c(2:length(order))) {
                breaks[i] <- nrow(filter(df, group %in%
                    as.vector(order[seq_len(i - 1)]))) +
                    nrow(filter(df, group == order[i])) / 2
            }
        } else {
            print("No group(s) was included in the analysis!")
        }

        # Calculate labels.
        labels <- vector()
        for (i in c(seq_len(length(order)))) {
            labels[i] <- paste0(
                stringr::str_to_title(order[i]),
                "\n",
                "(",
                "n=",
                filter(df, group == order[i]) %>%
                    group_by(file_name) %>%
                    count() %>%
                    nrow(),
                ")"
            )
        }

        # Set the x-lab.
        result <- p1 +
            scale_y_continuous(
                expand = c(0, 0),
                breaks = append(hline, 0, after = 0),
                labels = as.character(append(hline, 0, after = 0))
            ) +
            scale_x_continuous(
                expand = c(0, 0),
                breaks = breaks,
                labels = labels
            ) +
            theme(
                axis.text = element_text(size = 10),
                axis.title = element_text(size = 14, face = "bold")
            ) +
            labs(x = "Group", y = "Mode Fragment size (bp)", fill = "Group")
    } else if (type == "stacked") {

        # Generate summary.
        summary <-
            x %>%
            filter(group %in% as.vector(order)) %>%
            group_by(group, insert_size) %>%
            summarise(count = n()) %>%
            filter(count >= mincount) %>%
            mutate(prop = count / sum(count))

        # Generate stacked bar chart plot.
        result <- ggplot(
            summary,
            aes(
                x = factor(group, levels = order),
                y = prop * 100,
                fill = factor(insert_size)
            )
        ) +
            geom_bar(stat = "identity", width = 0.7, color = "gray24") +
            labs(x = "Group", y = "Percent", fill = "Mode Fragment Size") +
            theme_classic(base_size = 15) +
            theme(
                axis.text = element_text(size = 10),
                axis.title = element_text(size = 14, face = "bold")
            )
    }

    return(result)
}
