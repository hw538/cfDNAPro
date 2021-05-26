#' Get path to cfDNAPro example folder.
#'
#' cfDNAPro package has sample files in `inst/extdata`
#' directory. This function helps get the path to the data.
#'
#' @param data Name of data set. Such as "groups_picard" or "step6".
#' If `NULL`, the path of extdata folder will be returned.
#' @export
#' @examples
#' examplePath()
#' examplePath("groups_picard")
#' examplePath("step6")
#' @return A string. (i.e. the path.)
#' @export



examplePath <- function(data = NULL) {
    if (is.null(data)) {
        system.file("extdata", package = "cfDNAPro")
    } else {
        system.file("extdata", data, package = "cfDNAPro")
    }
}
