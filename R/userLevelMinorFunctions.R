# These functions might be helpful for users,
# but not important enough to export them.

#' @noRd
#' @importFrom stats median na.omit
#' @importFrom utils read.table
#' @importFrom rlang .data
#' @importFrom stringr str_detect
#' @importFrom dplyr filter select rename group_by mutate summarise summarize
#' @importFrom dplyr lag
#' @importFrom magrittr %>%

# InterPeakDist -----------------------------------------------------------

InterPeakDist <-
    function(path = getwd(),
             groups,
             limit,
             outfmt,
             input_type,
             ...) {
        if (missing(groups)) {
            groups <- list.dirs(path = path, 
                                recursive = FALSE, 
                                full.names = FALSE)
        }

        if (missing(limit)) {
            limit <- c(35, 135)
            cat("setting default limit to c(35,135).\n")
        }
        if (missing(outfmt)) {
            outfmt <- "df"
            cat("setting default outfmt to df.\n")
        }


        if (outfmt == "df" || outfmt == "dataframe") {
            dist <-
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
            return(dist)
        } else if (outfmt == "list") {
            dist <-
                GroupPeakFocus(
                    path = path,
                    groups = groups,
                    limit = limit,
                    input_type = input_type,
                    ...
                ) %>%
                lapply(lapply, calcualte_interpeak_dist)
            return(dist)
        }
    }



InterValleyDist <-
    function(path = getwd(),
             groups,
             limit,
             outfmt,
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
            cat("setting default limit to c(35,135)")
        }
        if (missing(outfmt)) {
            outfmt <- "df"
            cat("setting default outfmt to df")
        }

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
            return(dist)
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
            return(dist)
        }
    }



# InsertSizePropMode ------------------------------------------------------

InsertSizePropMode <- function(path, groups, outfmt, ...) {
    if (missing(path)) {
        path <- getwd()
    }
    if (missing(groups)) {
        groups <- list.dirs(
            path = path,
            recursive = FALSE,
            full.names = FALSE
        )
    }


    if (missing(outfmt) ||
        outfmt == "df" || outfmt == "dataframe") {
        result <- AllGroupMode_Dataframe(path = path, groups = groups, ...)
        return(result)
    } else if (outfmt == "list") {
        result <- AllGroupMode_List(path = path, groups = groups, ...)
        return(result)
    }
}


# InsertSizePropMedian ----------------------------------------------------


InsertSizePropMedian <-
    function(path = getwd(), groups, outfmt, ...) {
        if (missing(groups)) {
            groups <- list.dirs(
                path = path,
                recursive = FALSE,
                full.names = FALSE
            )
        }
        if (outfmt == "list") {
            result <- AllGroupPropMedian_List(path = path, groups = groups)
            return(result)
        } else if (missing(outfmt) ||
            outfmt == "df" || outfmt == "dataframe") {
            result <- AllGroupPropMedian_Dataframe(path = path, groups = groups)
            return(result)
        }
    }

# InsertSizePropMean ------------------------------------------------------


InsertSizePropMean <-
    function(path = getwd(), groups, outfmt, ...) {
        if (missing(groups)) {
            groups <- list.dirs(
                path = path,
                recursive = FALSE,
                full.names = FALSE
            )
        }

        if (outfmt == "list") {
            result <- AllGroupPropMean_List(path = path, groups = groups, ...)
            return(result)
        } else if (missing(outfmt) ||
            outfmt == "df" || outfmt == "dataframe") {
            result <-
                AllGroupPropMean_Dataframe(path = path, groups = groups, ...)
            return(result)
        }
    }


## Prop AllGroupPropMedian --------------------------------------------


AllGroupPropMedian_Dataframe <-
    function(path = getwd(), groups, ...) {
        if (missing(groups)) {
            groups <- list.dirs(
                path = path,
                recursive = FALSE,
                full.names = FALSE
            )
        }

        result <-
            AllGroupProp_Dataframe(path = path, ...) %>%
            calculate_median_for_df()
        return(result)
    }


AllGroupPropMedian_List <- function(path = getwd(), groups, ...) {
    if (missing(groups)) {
        groups <- list.dirs(
            path = path,
            recursive = FALSE,
            full.names = FALSE
        )
    }

    result <-
        AllGroupProp_List(path = path, ...) %>%
        lapply(bind_rows, .id = "file_name") %>%
        lapply(calculate_median_for_list)
    return(result)
}


# AllGroupPropMean -------------------------------------------------------


AllGroupPropMean_Dataframe <-
    function(path = getwd(), groups, ...) {
        if (missing(groups)) {
            groups <- list.dirs(
                path = path,
                recursive = FALSE,
                full.names = FALSE
            )
        }

        result <-
            AllGroupProp_Dataframe(path = path, ...) %>%
            calculate_mean_for_df()

        return(result)
    }


AllGroupPropMean_List <- function(path = getwd(), groups, ...) {
    if (missing(groups)) {
        groups <- list.dirs(
            path = path,
            recursive = FALSE,
            full.names = FALSE
        )
    }

    result <-
        AllGroupProp_List(path = path, ...) %>%
        lapply(bind_rows, .id = "file_name") %>%
        lapply(calculate_mean_for_list)

    return(result)
}



# AllGroupCdfMedian_Dataframe ----------------------------------------------


# AllGroupCdfMedian_List ---------------------------------------------------


# AllGroupCdfMean_Dataframe ------------------------------------------------



# AllGroupCdfMean_List -----------------------------------------------------



## PropCdf

## PropCdf


# AllGroupPropCdfMedian ----------------------------------------------------

AllGroupPropCdfMedian_Dataframe <-
    function(path = getwd(), groups, ...) {
        group <- insert_size <- prop <- cdf <- one_minus_cdf <- NULL

        if (missing(groups)) {
            groups <- list.dirs(
                path = path,
                recursive = FALSE,
                full.names = FALSE
            )
        }

        result <-
            AllGroupPropCdf_Dataframe(path = path, ...) %>%
            group_by(group, insert_size) %>%
            summarise(
                prop_median = median(prop,
                    na.rm = TRUE
                ),
                cdf_median = median(cdf, na.rm = TRUE),
                one_minus_cdf_median = median(one_minus_cdf, na.rm = TRUE)
            )

        return(result)
    }

AllGroupPropCdfMedian_List <-
    function(path = getwd(), groups, ...) {
        if (missing(groups)) {
            groups <- list.dirs(
                path = path,
                recursive = FALSE,
                full.names = FALSE
            )
        }

        result <-
            AllGroupPropCdf_List(path = path, ...) %>%
            lapply(bind_rows, .id = "file_name") %>%
            lapply(calculate_propcdf_median_for_list)

        return(result)
    }


# AllGroupPropCdfMean ------------------------------------------------------

AllGroupPropCdfMean_Dataframe <-
    function(path = getwd(), groups, ...) {
        group <- insert_size <- prop <- cdf <- one_minus_cdf <- NULL

        if (missing(groups)) {
            groups <- list.dirs(
                path = path,
                recursive = FALSE,
                full.names = FALSE
            )
        }

        result <-
            AllGroupPropCdf_Dataframe(path = path, ...) %>%
            group_by(group, insert_size) %>%
            summarise(
                prop_mean = mean(prop,
                    na.rm = TRUE
                ),
                cdf_mean = mean(cdf, na.rm = TRUE),
                one_minus_cdf_mean = mean(one_minus_cdf, na.rm = TRUE)
            )

        return(result)
    }

AllGroupPropCdfMean_List <- function(path = getwd(), groups, ...) {
    if (missing(groups)) {
        groups <- list.dirs(
            path = path,
            recursive = FALSE,
            full.names = FALSE
        )
    }

    result <-
        AllGroupPropCdf_List(path = path, ...) %>%
        lapply(bind_rows, .id = "file_name") %>%
        lapply(calculate_propcdf_mean_for_list)

    return(result)
}

# AllGroupPropCdfMeanMedian -------------------------------------------------

AllGroupPropCdfMeanMedian_Dataframe <-
    function(path = getwd(), groups, ...) {
        group <- insert_size <- prop <- cdf <- one_minus_cdf <- NULL

        if (missing(groups)) {
            groups <- list.dirs(
                path = path,
                recursive = FALSE,
                full.names = FALSE
            )
        }

        result <-
            AllGroupPropCdf_Dataframe(path = path, groups = groups, ...) %>%
            group_by(group, insert_size) %>%
            summarise(
                prop_mean = mean(prop, na.rm = TRUE),
                prop_median = median(prop, na.rm = TRUE),
                cdf_mean = mean(cdf, na.rm = TRUE),
                cdf_median = median(cdf, na.rm = TRUE),
                one_minus_cdf_mean = mean(one_minus_cdf, na.rm = TRUE),
                one_minus_cdf_median = median(one_minus_cdf, na.rm = TRUE)
            )

        return(result)
    }

AllGroupPropCdfMeanMedian_List <-
    function(path = getwd(), groups, ...) {
        if (missing(groups)) {
            groups <- list.dirs(
                path = path,
                recursive = FALSE,
                full.names = FALSE
            )
        }

        result <-
            AllGroupPropCdf_List(path = path, groups = groups, ...) %>%
            lapply(bind_rows, .id = "file_name") %>%
            lapply(calculate_propcdf_meanmedian_for_list)

        return(result)
    }




# Others ------------------------------------------------------------------

InsertSizeProp <- function(path = getwd(), groups, outfmt, ...) {
    if (missing(groups)) {
        groups <- list.dirs(
            path = path,
            recursive = FALSE,
            full.names = FALSE
        )
    }


    if (missing(outfmt) ||
        outfmt == "df" ||
        outfmt == "dataframe") {
        result <- AllGroupProp_Dataframe(path = path, groups = groups, ...)
        return(result)
    } else if (outfmt == "list") {
        result <- AllGroupProp_List(path = path, groups = groups, ...)
        return(result)
    }
}

GroupValleyFocus <- function(path, groups, limit, input_type, ...) {
    if (missing(groups)) {
        groups <- list.dirs(
            path = path,
            recursive = FALSE,
            full.names = FALSE
        )
    }

    if (missing(limit)) {
        limit <- c(35, 135)
        cat("setting default limit to c(35,135)")
    }
    group_valleys_focus <-
        AllGroupProp_List(
            path = path,
            groups = groups,
            input_type = input_type,
            ...
        ) %>%
        lapply(lapply, calculate_valleys) %>%
        lapply(lapply, filter_insert_size, limit = limit)

    return(group_valleys_focus)
}

GroupPeakFocus <- function(path, groups, limit, input_type, ...) {
    if (missing(groups)) {
        groups <- list.dirs(
            path = path,
            recursive = FALSE,
            full.names = FALSE
        )
    }

    if (missing(limit)) {
        limit <- c(35, 135)
        cat("setting default limit to c(35,135)")
    }
    group_peak_focus <-
        AllGroupProp_List(
            path = path,
            groups = groups,
            input_type = input_type,
            ...
        ) %>%
        lapply(lapply, calculate_peaks) %>%
        lapply(lapply, filter_insert_size, limit = limit)

    return(group_peak_focus)
}

AllGroupMetrics <- function(path, groups, input_type, ...) {
    if (missing(groups)) {
        groups <- list.dirs(
            path = path,
            recursive = FALSE,
            full.names = FALSE
        )
    }

    result1 <- read_raw_data(
        path = path,
        groups = groups,
        input_type = input_type,
        ...
    )

    return(result1)
}

# AllGroupMode_List----------------------------------------------------------

AllGroupMode_List <- function(path, groups, input_type, ...) {
    if (missing(groups)) {
        groups <- list.dirs(
            path = path,
            recursive = FALSE,
            full.names = FALSE
        )
    }


    result1 <- read_raw_data(
        path = path,
        groups = groups,
        input_type = input_type,
        ...
    )

    result2 <- lapply(
        result1,
        function(x) {
            lapply(x, calculate_sample_max_prop)
        }
    ) %>%
        lapply(bind_rows, .id = "file_name")

    return(result2)
}

# AllGroupMode_Dataframe-----------------------------------------------------

AllGroupMode_Dataframe <- function(path, groups, input_type, ...) {
    if (missing(groups)) {
        groups <- list.dirs(
            path = path,
            recursive = FALSE,
            full.names = FALSE
        )
    }

    result1 <- read_raw_data(
        path = path,
        groups = groups,
        input_type = input_type,
        ...
    )

    result2 <- lapply(
        result1,
        function(x) {
            lapply(x, calculate_sample_max_prop)
        }
    ) %>%
        lapply(bind_rows, .id = "file_name") %>%
        bind_rows(.id = "group")

    return(result2)
}
