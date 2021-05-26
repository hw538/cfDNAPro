#' Calculate the insert size metrics (i.e. prop, cdf, 1-cdf) or each group
#'
#' @import dplyr
#' @param path The root folder containing all groups folders,
#'    default is the present working folder.
#' @param groups The name of the groups, the input value should be vector,
#'    e.g. `groups=c('group1','group2')`,
#'    default is all folders in the folder path.
#' @param outfmt The output format,
#'    could specify as 'list' or 'dataframe' or 'df',
#'    default is dataframe.
#' @param input_type Character. The input file format,
#'    should be one of these: 'picard', 'bam'. The bam files has to be 
#'    marked duplicates.
#' @param ... Further arguments passed to or from other methods.
#' @return  The function returns the insert size metrics
#'   of each group in list or dataframe format.
#' @examples
#' # Get the path to example data.
#' path <- examplePath("groups_picard")
#' # Calculate the size.
#' df <- callSize(path = path)
#' @author Haichao Wang
#'
#' @export

callSize <- function(path, groups, outfmt, input_type, ...) {
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

    if (missing(outfmt)) {
        outfmt <- "df"
        cat("setting default outfmt to df.\n")
    }

    if (missing(input_type)) {
        input_type <- "picard"
        cat("setting default input_type to picard.\n")
    }

    if (outfmt == "df" || outfmt == "dataframe") {
        result <-
            AllGroupPropCdf_Dataframe(
                path = path,
                groups = groups,
                input_type = input_type,
                ...
            )
    } else if (outfmt == "list") {
        result <-
            AllGroupPropCdf_List(
                path = path,
                groups = groups,
                input_type = input_type,
                ...
            )
    }
    return(result)
}
