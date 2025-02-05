#' Call trinucleotides and summarise the cfDNA information
#' for each target mutation locus.
#'
#' This function processes a GRanges object, summarizing cfDNA fragment
#' information for each target mutation locus.
#' 
#' It annotates each mutation locus with the number
#' and type of supporting fragments according
#' to their read-pair overlap status.
#'
#' The median fragment length is also annotated
#' for each read-pair overlap type.
#'
#' Locus-based consensus mutation is determined
#' by selecting the most frequent mismatch type,
#' with priority given to concordant read-pair mutations (CO_MUT)
#' followed by single-read mutations (SO_MUT).
#' The consensus mutations are then used to derive
#' the trinucleotide substitution types (SBS96).
#'
#' @import GenomicRanges
#' @import IRanges
#' @import stringr
#' @importFrom dplyr select
#' @importFrom stats median
#' @importFrom magrittr '%>%'
#'
#' @param frag_obj_mut GRanges object containing fragment and mutation data.
#'
#' @return dataframe with summarised mutational and trinucleotide data.
#' @export
#'
#' @examples
#' trinuc_df <-  callTrinucleotide(frag_obj_mut)
callTrinucleotide <- function(frag_obj_mut) {

  # Check if GRanges has mutational information
  check_mutation_in_metadata(frag_obj = frag_obj_mut)

  # Get the appropriate genome reference
  genome <- get_genome_reference(frag_obj_mut)

  # Convert relevant data to dataframe
  gr_df <- prepare_data(frag_obj_mut)
  # Check if any fragments are present within dataframe
  check_mutation_status(gr_df)
  # Stratify mutational data into separate columns
  gr_df <- stratify_mutation_status(gr_df)
  # Update the format of mutational columns
  gr_df <- update_locus_status(gr_df)

  # Summarise mutation counts and fragment lengths stratified by overlap type
  df1 <- summarize_mutational_data(gr_df)
  df2 <- summarize_fragment_lengths(gr_df)
  # Merge the two summaries
  merged_table_df <- merge(df1, df2, by = 'target_mutation')

  # Define mapping for calculation of consensus mismatch type per locus
  mapping <- list(
    CO_MUT = "MUT",
    SO_MUT = "MUT",
    DO = "discordant",
    SO_OTHER = "other_base_single_read",
    CO_OTHER = "other_base_read_pair"
  )

  # Find the consensus mismatch type per locus based on counts
  merged_table_df <- update_consensus_mismatch(merged_table_df, gr_df, mapping)

  # Add the SBS 96 annotation for consensus locus-based mutations
  merged_table_df <- get_trinucleotide(merged_table_df, genome)

  return(merged_table_df)
}
