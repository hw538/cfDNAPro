#' Calculate the mode fragment size of each sample
#'
#' @importFrom magrittr %>%
#' @importFrom  dplyr group_by ungroup arrange mutate filter summarise
#' @importFrom rlang .data
#'
#' @param path The root folder containing all groups folders,
#'     default is the present working folder.
#' @param groups The name of the groups, the input value should be vector,
#'     e.g. groups=c('group1','group2'),
#'     default is all folders in the folder path.
#' @param outfmt The output format, 'list' or 'dataframe' or 'df',
#'     default is dataframe.
#' @param order The order in the sorted output,
#'     default value equals to 'groups' parameter.
#' @param summary Summarize the dataframe result
#'     by calculating each mode size and its count number.
#'     Default value is False.
#' @param mincount Minimum count number of each mode size
#'     in the summarized output.
#'     Only significant when 'summary = TRUE'.
#' @param input_type Character. The input file format,
#'     should be one of these: 'picard', 'bam'. The bam files has to be 
#'     marked duplicates.
#' @param ... Further arguments passed to or from other methods.
#' @return  The function returns the inter valley distance
#'     in list or dataframe format.
#' @examples
#' # Get the path to example data.
#' path <- examplePath("groups_picard")
#' # Calculate the mode.
#' df <- callMode(path = path)
#' @author Haichao Wang
#' @export


callMode <-
    function(path,
            groups,
            outfmt = "df",
            order = groups,
            summary,
            mincount,
            input_type,
            ...) {
        # Check the path.
        if (missing(path)) {
            path <- getwd()
        }

        # Check the groups.
        if (missing(groups)) {
            groups <-
                list.dirs(
                    path = path,
                    recursive = FALSE,
                    full.names = FALSE
                )
        }
        if (missing(summary)) {
            summary <- FALSE
        }
        if (missing(mincount)) {
            mincount <- 0
        }
        if (missing(input_type)) {
            input_type <- "picard"
            cat("setting default input_type to picard.\n")
        }

        # Generating sorted df.
        if (outfmt == "df" || outfmt == "dataframe") {
            result <-
                AllGroupMode_Dataframe(
                    path = path,
                    groups = groups,
                    input_type = input_type,
                    ...
                ) %>%
                # Group by 'group' for sort.
                group_by(.data$group) %>%
                # Sort by insert_size.
                arrange(.data$insert_size, .by_group = TRUE) %>%
                # Change group order which will alter the order of groups.
                arrange(factor(.data$group, levels = order)) %>%
                ungroup() %>%
                mutate(x_index = row_number()) %>%
                mutate(duplication = duplicated(.data$file_name))

            # Checking the existence of dual modes in each samples.
            check_dup_mode(result)

            if (summary) {
                result <-
                    result %>%
                    filter(.data$group %in% as.vector(order)) %>%
                    group_by(.data$group, .data$insert_size) %>%
                    summarise(count = n()) %>%
                    filter(count >= mincount) %>%
                    mutate(prop = count / sum(count))
            }

            # Generate list format result.
        } else if (outfmt == "list") {
            result <-
                AllGroupMode_List(
                    path = path,
                    groups = groups,
                    input_type = input_type,
                    ...
                ) %>% lapply(function(x) {
                    mutate(x, duplication = duplicated(.data$file_name))
                })

            # Check the existence of dual modes.
            lapply(result, check_dup_mode)

            # Summarize the df.
            if (summary) {
                result <- lapply(result, function(x) {
                    x %>%
                        group_by(.data$group, .data$insert_size) %>%
                        summarise(count = n()) %>%
                        filter(count >= mincount) %>%
                        mutate(prop = count / sum(count))
                })
            }
        }
        return(result)
    }
