
#' Call Fragment Length
#'
#' Computes distribution of fragment sizes from BAM/GRanges data. Optionally
#' integrates mutational data, normalizes fragment counts, and adjusts for
#' specific reference types.
#'
#' @import ggplot2
#' @import plyranges
#' @importFrom dplyr n sample_n
#' @importFrom magrittr '%>%'
#'
#' @param frag_obj GRanges object containing DNA fragment information.
#' @param isize_min Minimum insert size for filtering (default 1).
#' @param isize_max Maximum insert size for filtering (default 1000).
#' @param integrate_mut Logical, TRUE integrates mutational data.
#' @param ref_type Specifies the type of reference fragments for comparison
#' with MUT fragments: 'locus_fragment' refers to fragments/read-pairs
#' overlapping specific mutation loci, 'outer_fragment' refers to
#' fragments/read-pairs that do not overlap indicated mutation loci.
#' @param downsample_ref Logical, TRUE to match REF and MUT fragment counts.
#' @param ... Additional arguments for downstream functions.
#'
#' @return Returns a tibble with summarized fragment data, optionally
#' including mutational analysis.
#' @export
#'
#' @examples
#' callLength(frag_obj = myFragData, isize_min = 1L, isize_max = 1000L,
#'            integrate_mut = TRUE)
callLength <- function(frag_obj,
                       isize_min = 1L,
                       isize_max = 1000L,
                       integrate_mut = FALSE,
                       ref_type = "locus_fragment",
                       downsample_ref = FALSE,
                       ...) {

  if (integrate_mut == FALSE) {

    frag <- frag_obj
    # Calculating insert sizes
    message("Calculating insert sizes...")
    frag$insert_size <- BiocGenerics::width(frag)

    # Size analysis
    frag <- plyranges::filter(frag,
                              insert_size >= isize_min &
                                insert_size <= isize_max)
    isize <- frag$insert_size

    isize_tibble <- tibble("insert_size" = isize, "count" = 1) %>%
      dplyr::filter(!is.na(insert_size))

    result <- isize_tibble %>%
      dplyr::group_by(.data$insert_size) %>%
      dplyr::summarise("All_Reads.fr_count" = sum(count))

    # Quality control results
    # Create a vector of elements
    isize_ref <- seq.int(isize_min, isize_max, by = 1L) %>%
      as_tibble()

    colnames(isize_ref) <- c("insert_size")

    # Report abnormal isizes
    missing_isize <- dplyr::anti_join(isize_ref, result, by = "insert_size")

    # Handle missing isize(s)
    if (nrow(missing_isize) != 0) {
      message("Missing isize detected: ")
      print(dplyr::pull(missing_isize, insert_size))

      result <- dplyr::right_join(result, isize_ref, by = "insert_size") %>%
        tidyr::replace_na(replace = list(All_Reads.fr_count = 0)) %>%
        dplyr::arrange(insert_size) %>%
        dplyr::mutate(prop = All_Reads.fr_count / sum(All_Reads.fr_count))

      message(paste(" Missing isize(s) added back",
                    "to the final result with count of 0!"))

      message("Job completed successfully. ")
    }

    return(result)

  } else if (integrate_mut == TRUE) {

    return(process_length_mut(frag_obj, ref_type, downsample_ref))

  }

}



