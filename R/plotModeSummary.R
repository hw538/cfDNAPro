#' Summarize and plot mode fragment size in a stacked bar chart
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr n
#' @import ggplot2
#' @import stringr
#' @import dplyr
#' @importFrom rlang .data
#' @param x A long-format dataframe contains mode fragment size,
#'   a template please refer to the result of `callMode` function.
#' @param order The groups show in the final plot,
#'   the input value should be vector,
#'   e.g. `groups = c('group1','group2')`,
#'  default is all folders in the folder path.
#' @param summarized Logical value, default is False.
#' @param mode_partition This should be a list.
#'   This decides how the modes are
#'   partitioned in each stacked bar. Default value is
#'   `list(c(80, 81), c(111, 112), c(167))`.
#'   Also this function will automatically calculate an 'Others' group
#'   which includes
#'   the modes not mentioned by users.
#' @param ... Further arguments passed to or from other methods.
#' @return  The function returns the plot.
#'
#' @examples
#' # Get the path to example data.
#' path <- examplePath("groups_picard")
#' # Calculate the modes.
#' df <- callMode(path = path)
#' # Plot mode summary.
#' plot <- plotModeSummary(df,
#'     mode_partition = list(c(80, 81), c(111, 112), c(167))
#' )
#' @author Haichao Wang
#' @export
#'
#'
plotModeSummary <- function(x,
                            order,
                            summarized,
                            mode_partition,
                            ...) {
    insert_size <- group <- prop <- NULL


    if (missing(summarized)) {
          summarized <- FALSE
      }
    if (missing(order)) {
          order <- as.vector(unique(x$group))
      }
    if (missing(mode_partition)) {
          mode_partition <- list(c(80, 81), c(111, 112), c(167))
      }

    # This step is to calculated the proportion of each modes.
    if (summarized == FALSE) {
        x <-
            x %>%
            dplyr::filter(.data$group %in% order) %>%
            dplyr::group_by(.data$group, .data$insert_size) %>%
            dplyr::summarise(count = n()) %>%
            dplyr::mutate(prop = count / sum(count))
    }



    extract1 <- function(partition) {
        result <- x %>%
            filter(.data$insert_size %in% partition) %>%
            group_by(.data$group) %>%
            summarise(
                insert_size = paste(partition, collapse = "/"),
                count = sum(count)
            )
    }
    # Extract the mode partitions set by users.
    r1 <- lapply(mode_partition, extract1)

    # Extract the Others.

    Others <-
        x %>%
        filter(!insert_size %in% unlist(mode_partition)) %>%
        group_by(group) %>%
        summarise(
            insert_size = "Others",
            count = sum(count)
        )

    # Append the 'Others' group to the r1 list if it is not blank.

    if (length(Others) == 0) {
        message("No other mode groups detected except those indicated by you 
              in the parameter mode_partition.")
    } else {
        r1[["Others"]] <- Others
    }


    # Bind the rows of dataframes in r1 list and calculate the prop.
    mode_final <-
        bind_rows(r1) %>%
        filter(group %in% order) %>%
        group_by(group) %>%
        mutate(prop = count / sum(count))



    # Plot the stacked bar chart.

    plot <- ggplot(
        mode_final,
        aes(
            x = factor(group, levels = order),
            y = prop * 100,
            fill = factor(insert_size)
        )
    ) +
        geom_bar(
            stat = "identity",
            width = 0.7,
            color = "gray24"
        ) +
        labs(
            x = "Group",
            y = "Percent (%)",
            fill = "Mode Fragment Size"
        ) +
        theme_classic(base_size = 15) +
        theme(
            axis.text = element_text(size = 10),
            axis.title = element_text(size = 14, face = "bold")
        )

    return(plot)
}
