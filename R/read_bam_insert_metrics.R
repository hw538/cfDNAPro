#' Calculate insert sizes from a curated GRanges object
#' 
#' @param genome A BSgenome object, required when
#'   \code{genome_label = "user_define"}. The BAM and reference must use
#'   matching sequence names and lengths. For built-in labels a supplied
#'   BSgenome is ignored with a warning. If \code{genome_label} is omitted,
#'   the older shorthand \code{genome = "hg38"} still selects that built-in
#'   reference. References and dependencies must be installed before offline use.
#' @param bamfile The bam file name.
#' @param genome_label Reference label: "hg19", "hg38", "hg38-NCBI", or
#'   "user_define". Defaults to "hg19". For "user_define", supply \code{genome}.
#' @param outdir The path for saving rds file. Default is FALSE, i.e. not saving.
#' @param strand_mode Usually the strand_mode  = 1 means the First read is 
#'    aligned to positive strand. Details please see GenomicAlignments docs.
#' @param chromosome_to_keep Should be a character vector containing the 
#'    seqnames to be kept in the GRanges object. 
#'    Default is paste0("chr", 1:22).
#' @param isize_min min fragment length to keep, default is 1L.
#' @param isize_max max fragment length to keep, default is 1000L.
#' @param ... Further arguments passed to or from other methods.

#'
#' @importFrom Rsamtools ScanBamParam
#' @importFrom Rsamtools scanBam
#' @importFrom dplyr filter rename select group_by summarise 
#' @importFrom BiocGenerics start end strand width
#' @importFrom rlang .data
#' @importFrom tibble tibble
#' 
#' @return This function returns a dataframe with two columns: "insert_size" 
#'    and "All_Reads.fr_count".
#' @export
#' @author Haichao Wang
#'
#' @examples 
#' \dontrun{
#' 
#' object <- read_bam_insert_metrics(bamfile = "/path/to/bamfile.bam")
#' }
#' 
read_bam_insert_metrics <- function(bamfile = NULL,
                                    fragment_obj = NULL,
                                    chromosome_to_keep =paste0("chr", 1:22),
                                    strand_mode = 1,
                                    genome_label = "hg19",
                                    outdir = FALSE,
                                    isize_min = 1L,
                                    isize_max = 1000L,
                                    genome = NULL,
                                    ...) {
 
  insert_size <- NULL
  
  # users are required to supply only one of 'bamfile' or 'fragment_obj" from
  # readBam() function
  switch(
    rlang::check_exclusive(bamfile, fragment_obj),
    bamfile = message("`bamfile` was supplied."),
    fragment_obj = message("`fragment_obj` was supplied.")
  )
  
  if(!is.null(bamfile)) {
    reference <- .cfDNAPro_resolve_genome(genome_label, genome, missing(genome_label))
    frag <- readBam(bamfile = bamfile, 
                  chromosome_to_keep = chromosome_to_keep,
                  strand_mode = strand_mode,
                  genome_label = reference$label,
                  genome = if (reference$label == "user_define") reference$genome else NULL,
                  outdir = outdir)
  } else if(!is.null(fragment_obj)) {
    
    frag <- fragment_obj
    
  }
  
  
  # calculating insert sizes
  message("Calculating insert sizes...")
  frag$insert_size <- BiocGenerics::width(frag)
  
  # size analysis
  frag <- frag[frag$insert_size >= isize_min & frag$insert_size <= isize_max]
  isize <- frag$insert_size 
  
  isize_tibble <- tibble("insert_size" = isize, "count" = 1 ) %>%
    dplyr::filter(!is.na(insert_size))
  
  result <- isize_tibble %>%
    dplyr::group_by(.data$insert_size) %>%
    dplyr::summarise("All_Reads.fr_count" = sum(count))
  
  
  # quality control results
  # Create a vector of elements
  isize_ref <- seq.int(isize_min, isize_max, by = 1L) %>% 
    as_tibble()
  
  colnames(isize_ref) <- c("insert_size")
  
  # report abnormal isizes
  
  missing_isize <- dplyr::anti_join(isize_ref, result, by = "insert_size")

  # handle missing isize(s)
  
  if(nrow(missing_isize) != 0) {
    message("Missing isize detected: ")
    print(dplyr::pull(missing_isize, insert_size))
    
    result <- dplyr::right_join(result, isize_ref, by = "insert_size") %>%
      tidyr::replace_na(replace = list(All_Reads.fr_count = 0)) %>% 
      dplyr::arrange(insert_size)
    
    message("Missing isize(s) added back to the final result with count of 0!")
    
    message("Job completed successfully. ")
    
  }
  

  
  return(result)
}
