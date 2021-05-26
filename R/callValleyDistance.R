#' Calculate the inter-valley distance of insert size
#'
#' @importFrom magrittr %>%
#' @importFrom  dplyr group_by ungroup arrange mutate filter summarise bind_rows
#' @importFrom rlang .data
#' @param path The root folder containing all groups folders,
#'   default is the present working folder.
#' @param groups The name of the groups, the input value should be vector,
#'   e.g. groups = c('group1','group2'), default is all folders
#'   in the folder path.
#' @param limit The insert size range that will be focused on,
#'   default value is `limit = c(35,135)`.
#' @param outfmt The output format,
#'   could specify as 'list' or 'dataframe' or 'df',
#'   default is dataframe.
#' @param summary If TRUE, summarize the output.
#' @param mincount The minimum count value of inter-valley distance.
#' @param input_type Character. The input file format
#'   should be 'picard' or 'bam'. The bam files has to be 
#'   marked duplicates.
#' @param ... Further arguments passed to or from other methods.
#' @return  The inter-valley distance in a list or dataframe.
#'
#' @examples
#' # Get the path to example data.
#' path <- examplePath("groups_picard")
#' # Calculate the inter-valley distance.
#' df <- callValleyDistance(path = path)
#' @author Haichao Wang
#'
#' @export

callValleyDistance <- function(path = getwd(),
                              groups,
                              limit,
                              outfmt,
                              summary,
                              mincount,
                              input_type,
                              ...) {
    # Define missing parameters.
    if (missing(groups)) {
        groups <- list.dirs(
            path = path,
            recursive = FALSE,
            full.names = FALSE
        )
    }

    if (missing(limit)) {
        limit <- c(35, 135)
        cat("setting default limit to c(35,135).\n")
    }

    if (missing(outfmt)) {
        outfmt <- "df"
        cat("setting default outfmt to df.\n")
    }

    if (missing(summary)) {
          summary <- FALSE
      }

    if (missing(mincount)) {
        mincount <- 0
        cat("setting the mincount to 0.\n")
    }

    if (missing(input_type)) {
        input_type <- "picard"
        cat("setting default input_type to picard.\n")
    }

    # Output format == dataframe.
    if (outfmt == "df" || outfmt == "dataframe") {
        dist <-
            GroupValleyFocus(
                path = path,
                groups = groups,
                limit = limit,
                input_type = input_type,
                ...
            ) %>%
            lapply(lapply, calcualte_intervalley_dist) %>%
            lapply(bind_rows, .id = "file_name") %>%
            bind_rows(.id = "group")

        # If summary = T, summarize the data.
        if (summary) {
            dist <-
                dist %>%
                group_by(.data$group, .data$intervalley_dist) %>%
                summarise(count = n()) %>%
                na.omit() %>%
                filter(count >= mincount) %>%
                mutate(prop = count / sum(count))
        }

        return(dist)

        # Output format == list.
    } else if (outfmt == "list") {
        dist <-
            GroupValleyFocus(
                path = path,
                groups = groups,
                limit = limit,
                input_type = input_type,
                ...
            ) %>%
            lapply(lapply, calcualte_intervalley_dist)

        # If summary = T, summarize the data.
        if (summary) {
            dist <- lapply(dist, function(x) {
                x %>%
                    group_by(.data$group, .data$intervalley_dist) %>%
                    summarise(count = n()) %>%
                    na.omit() %>%
                    filter(count >= mincount) %>%
                    mutate(prop = count / sum(count))
            })
        }
        return(dist)
    }
}
