#' Calculate the inter-peak distance of insert size
#'
#' @importFrom magrittr %>%
#' @importFrom  dplyr group_by ungroup arrange mutate filter summarise bind_rows
#' @importFrom rlang .data
#' @param path The root folder containing all groups folders.
#'   Default is the present working folder.
#' @param groups The name of the groups, the input value should be vector,
#'   e.g. groups=c('group1','group2').
#'   Default is all folders in the folder path.
#' @param limit The insert size range that will be focused on.
#'   Default value is `limit = c(35,135)`.
#' @param outfmt The output format, a 'list' or 'dataframe'.
#'   Default is dataframe.
#' @param summary If TRUE, summarize the output.
#' @param mincount The minimum count value of inter-peak distance
#'     in the summary.
#' @param input_type Character. The input file format, 'picard' or 'bam'.
#'     The bam files has to be marked duplicates.
#' @param ... Further arguments passed to or from other methods.
#' @return  The function returns the inter peak distance in list
#'   or dataframe format.
#' @examples
#' # Get the path to example data.
#' path <- examplePath("groups_picard")
#' # Calculate the inter-peak distance.
#' df <- callPeakDistance(path = path)
#' @author Haichao Wang
#' @export

callPeakDistance <-
    function(path = getwd(),
            groups,
            limit,
            outfmt,
            summary,
            mincount,
            input_type,
            ...) {
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
            cat("Setting default mincount to 0.\n")
        }
        if (missing(input_type)) {
            input_type <- "picard"
            cat("setting default input_type to picard.\n")
        }
        if (outfmt == "df" || outfmt == "dataframe") {
            dist1 <-
                GroupPeakFocus(
                    path = path,
                    groups = groups,
                    limit = limit,
                    input_type = input_type,
                    ...
                ) %>%
                lapply(lapply, calcualte_interpeak_dist) %>%
                lapply(bind_rows, .id = "file_name") %>%
                bind_rows(.id = "group")

            # If summary = T, summarize the data.
            if (summary) {
                dist1 <-
                    dist1 %>%
                    group_by(.data$group, .data$interpeak_dist) %>%
                    summarise(count = n()) %>%
                    na.omit() %>%
                    filter(count >= mincount) %>%
                    mutate(prop = count / sum(count))
            }
            return(dist1)
        } else if (outfmt == "list") {
            dist2 <-
                GroupPeakFocus(
                    path = path,
                    groups = groups,
                    limit = limit,
                    input_type = input_type,
                    ...
                ) %>%
                lapply(lapply, calcualte_interpeak_dist)
            # If summary = T, summarize the data.
            if (summary) {
                dist2 <- lapply(dist2, function(x) {
                    x %>%
                        group_by(.data$group, .data$interpeak_dist) %>%
                        summarise(count = n()) %>%
                        na.omit() %>%
                        filter(.data$count >= mincount) %>%
                        mutate(prop = count / sum(count))
                })
            }
            return(dist2)
        }
    }
