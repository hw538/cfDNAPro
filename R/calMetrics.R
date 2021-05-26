#' Calculate the metrics of insert size
#'
#' @importFrom magrittr %>%
#' @import dplyr
#' @param path The root folder containing all groups folders,
#'     default is the present working folder.
#' @param groups The name of the groups, the input value should be vector,
#'     e.g. groups=c('group1','group2'),
#'     default is all sub-folders in the 'path'.
#' @param fun String value, the types of metrics to be calculated.
#'    Default is 'all', which means both median and mean values
#'    will be returned.
#' @param outfmt The output format, a 'list' or 'dataframe' or 'df',
#'     default is dataframe.
#' @param input_type Character. The input file format,
#'    should be one of these: 'picard', 'bam'. The bam files has to be 
#'    marked duplicates.
#' @param ... Further arguments passed to or from other methods.
#' @return  The inter valley distance in list or dataframe format.
#' @examples
#' # Get the path to example data.
#' path <- examplePath("groups_picard")
#' # Calculate the metrics.
#' df <- callMetrics(path = path)
#' @author Haichao Wang
#' @export

callMetrics <-
    function(path = getwd(),
            groups,
            fun = "all",
            outfmt = "df",
            input_type,
            ...) {
        if (missing(groups)) {
            groups <- list.dirs(
                path = path,
                recursive = FALSE,
                full.names = FALSE
            )
        }

        if (missing(input_type)) {
            input_type <- "picard"
            cat("setting default input_type to picard.\n")
        }


        if (fun == "all" & outfmt == "df") {
            result <-
                AllGroupPropCdfMeanMedian_Dataframe(
                    path = path,
                    groups = groups,
                    input_type = input_type,
                    ...
                )
        } else if (fun == "all" & outfmt == "list") {
            result <-
                AllGroupPropCdfMeanMedian_List(
                    path = path,
                    groups = groups,
                    input_type = input_type,
                    ...
                )
        } else if (fun == "median" & outfmt == "df") {
            result <-
                AllGroupPropCdfMeanMedian_Dataframe(
                    path = path,
                    groups = groups,
                    input_type = input_type,
                    ...
                ) %>% select(-ends_with("mean"))
        } else if (fun == "median" & outfmt == "list") {
            result <-
                AllGroupPropCdfMeanMedian_List(
                    path = path,
                    groups = groups,
                    input_type = input_type,
                    ...
                ) %>%
                lapply(select, -ends_with("mean"))
        } else if (fun == "mean" & outfmt == "df") {
            result <-
                AllGroupPropCdfMeanMedian_Dataframe(
                    path = path,
                    groups = groups,
                    input_type = input_type,
                    ...
                ) %>%
                select(-ends_with("median"))
        } else if (fun == "mean" & outfmt == "list") {
            result <-
                AllGroupPropCdfMeanMedian_List(
                    path = path,
                    groups = groups,
                    input_type = input_type,
                    ...
                ) %>%
                lapply(select, -ends_with("median"))
        }

        return(result)
    }
