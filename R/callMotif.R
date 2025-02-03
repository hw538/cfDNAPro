#' Call Motif
#'
#' This function analyzes motifs within DNA sequences
#' based on specified criteria.
#' It supports integration with mutation data and allows for downsampling of
#' reference based (i.e., non-mutant) motifs to balance the dataset.
#'
#' @import tibble
#' @import purrr
#' @import tidyr
#' @importFrom magrittr '%>%'
#' @importFrom dplyr all_of
#' @importFrom tidyr expand_grid unite
#'
#' @param frag_obj GRanges object containing DNA fragment information.
#' @param motif_type Character string describing the type of motifs to analyze.
#' @param motif_length An integer specifying the length of the motifs.
#' @param integrate_mut Logical, indicates whether to integrate mutational data.
#' @param ref_type Specifies the type of reference fragments for comparison
#' with MUT fragments: 'locus_fragment' refers to fragments/read-pairs
#' overlapping specific mutation loci, 'outer_fragment' refers to
#' fragments/read-pairs that do not overlap indicated mutation loci.
#' @param downsample_ref Logical, if TRUE, downsamples reference motifs to match
#'        the count of mutant motifs.
#' @param ... Additional parameters passed to underlying functions.
#'
#' @return Returns a tibble containing the results of the motif analysis.
#' @export
#'
#' @examples
#' callMotif(frag_obj = myData, genome_label = "hg38", motif_type = "s",
#'           motif_length = 2, integrate_mut = TRUE, downsample_ref = FALSE)
callMotif <- function(frag_obj,
                      motif_type = "s",
                      motif_length = 3L,
                      integrate_mut = FALSE,
                      ref_type = "locus_fragment",
                      downsample_ref = FALSE,
                      ...) {

  motif_length <- as.numeric(motif_length)

  message("Started to extract ", paste0(motif_type, motif_length), " motif...")

  # Select the appropriate genome reference
  bsgenome_obj <- get_genome_reference(frag_obj_mut)

  genome_label <- genome@metadata$genome

  message("You reference genome was set to ", genome_label, "!")

  if (integrate_mut == TRUE) {

    # Integrating mutational data
    integrate_motif_mut(
      frag_obj, downsample_ref, ref_type,
      bsgenome_obj, motif_type, motif_length)

  } else {

    # Default processing without mutational data
    processMotif(
      frag_obj, bsgenome_obj,
      motif_type, motif_length)

  }
}

#' Helper function to process the motif data
processMotif <- function(frag_obj, bsgenome_obj, motif_type, motif_length) {

  # Extract the motif to a vector
  motif <- get_motif(obj = frag_obj,
                     genome = bsgenome_obj,
                     motif_type = motif_type,
                     motif_length = motif_length)

  # Summarize the vector and convert to tibble format
  result <- table(motif) %>%
    tibble::as_tibble()

  # Create a vector of elements
  pos <- c("C", "G", "A", "T")
  base_index <- seq.int(1, motif_length, by = 1)

  letter_list <- lapply(base_index, function(x, y) return(y), y = pos) %>%
    setNames(paste0("base", base_index))

  motif_ref <- expand_grid(!!!letter_list) %>%
    unite(col = "motif", all_of(base_index), sep = "")

  # Report abnormal motifs.
  # Handle ambiguous motif(s)
  ambiguous_motif <- dplyr::anti_join(result, motif_ref, by = "motif")
  missing_motif <- dplyr::anti_join(motif_ref, result, by = "motif")

  if (nrow(ambiguous_motif) != 0) {
    message("ambiguous motif (i.e., 'N' base exists) detected: ")
    print(ambiguous_motif)

    result <- dplyr::filter(result, !stringr::str_detect(motif, "N"))
    message("Ambiguous motif(s) removed from result!")
  }

  # Handle missing motif(s)
  if (nrow(missing_motif) != 0) {
    message("Missing motif detected: ")
    print(missing_motif)

    result <- dplyr::right_join(result, motif_ref, by = "motif") %>%
      tidyr::replace_na(list(n = 0))

    message("Missing motif added back to the final result with count of 0!")
  }

  # Calculate the fraction
  result_frac <- result %>%
    dplyr::mutate(fraction = .data$n / sum(.data$n))

  # Set the factor level of motifs
  result_frac$motif <- factor(result_frac$motif, levels = sort(motif_ref$motif))

  message("Job completed successfully. ")

  return(result_frac)
}

