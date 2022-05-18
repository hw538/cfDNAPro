#' Calculate insert sizes from a curated GRanges object
#' 
#' @param bamfile The bam file name.
#' @param genome_label The Genome you used in the alignment. 
#'    Should be "hg19" or "hg38". Default is "hg19".
#' @param outdir The path for saving rds file. Default is NA, i.e. not saving.
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
#' @importFrom plyranges filter
#' @importFrom dplyr filter rename select group_by summarise n
#' @importFrom tibble as_tibble

read_bam_insert_metrics <- function(bamfile,
                                    chromosome_to_keep =paste0("chr", 1:22),
                                    strand_mode = 1,
                                    genome_label = "hg19",
                                    outdir = NA,
                                    isize_min = 1L,
                                    isize_max = 1000L,
                                    ...) {
 
  
  frag <- readBam(bamfile = bamfile, 
                  chromosome_to_keep = chromosome_to_keep,
                  strand_mode = strand_mode,
                  genome_label = genome_label,
                  outdir = outdir)
  
  # calculating insert sizes
  message("Calculating insert sizes...")
  frag$insert_size <- width(frag)
  
  # size analysis
  frag <- plyranges::filter(frag, 
                          insert_size >= isize_min & insert_size <= isize_max)
  isize <- frag$insert_size 
  
  isize_tibble <- tibble::tibble("insert_size" = isize, "count" = 1 )
  
  result <- isize_tibble %>%
    dplyr::group_by(.data$insert_size) %>%
    dplyr::summarise("All_Reads.fr_count" = sum(count))
  
  return(result)
}
