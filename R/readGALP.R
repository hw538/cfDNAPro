
#' Read bam file into GAlignmentPairs object
#' @import magrittr
#' @import GenomeInfoDb
#' @import GenomicAlignments
#' @import S4Vectors
#' @importFrom methods is
#' @importFrom Rsamtools scanBamFlag ScanBamParam
#' @importFrom GenomicRanges GRanges seqnames
#' @importFrom IRanges IRanges
#'
#' @param genome_label The genome you used for alignment: `"hg19"`,
#'   `"hg38"`, `"hg38-NCBI"`, or `"user_define"`. Default is `"hg19"`.
#'   For `"user_define"`, pass a matching `BSgenome` via `genome`.
#' @param genome Used only when `genome_label = "user_define"`: a `BSgenome`
#'   object. Ignored with a warning for built-in `genome_label` values if
#'   a `BSgenome` is supplied.
#' @param bamfile The bam file name.
#' @param outdir The path for saving rds file. Default is NA, i.e. not saving.
#' @param strand_mode Usually the strand_mode  = 1 means the First read is 
#'    aligned to positive strand. Details please see GenomicAlignments docs.
#' @param chromosome_to_keep Should be a character vector containing the 
#'    seqnames to be kept in the GRanges object or boolean FALSE. FALSE means not filtering.
#'    Default is paste0("chr", 1:22).
#' @param use_names 
#' @param galp_flag 
#' @param galp_what A character vector naming the fields to return 
#' Rsamtools::scanBamWhat() returns a vector of available fields. Fields are 
#' described on the Rsamtools::scanBam help page.
#' @param galp_tag 
#' @param galp_mapqFilter 
#' @param ... Further arguments passed to or from other methods.
#'
#' @return This function returns curated GRanges object.
#' @export
#' @author Haichao Wang
#'
#' @examples 
#' \dontrun{
#' 
#' object <- readGALP(bamfile = "/path/to/bamfile.bam", 
#'                    outdir = "./",
#'                    chromosome_to_keep = c("chr1", "chr2", "chr3"))
#' }
#' 

readGALP <- function(
    bamfile,
    use_names = TRUE,
    chromosome_to_keep = paste("chr", 1:22, sep = ""),
    strand_mode = 1,
    genome_label = "hg19",
    genome = NULL,
    outdir = FALSE,
    galp_flag =  Rsamtools::scanBamFlag(
      isPaired = TRUE,
      isDuplicate = FALSE,
      isSecondaryAlignment = FALSE,
      isUnmappedQuery = FALSE,
      isSupplementaryAlignment = FALSE),
    galp_what =  c("cigar", "mapq", "isize"),
    galp_tag = c("NM", "MD"),
    galp_mapqFilter = 30,
    ...){
  
  ############################################################################
  # Check parameters
  ############################################################################
  
  # check if the input bam is paired-end
  if (!Rsamtools::testPairedEndBam(bamfile)) {
    stop("Input is not paired-end Bam file.")
  }
  
  .valid_genome_labels <- c("hg19", "hg38", "hg38-NCBI", "user_define")
  if (!genome_label %in% .valid_genome_labels || length(genome_label) != 1L ||
        is.na(genome_label)) {
    stop(
      "genome_label must be one of: ",
      paste0("'", .valid_genome_labels, "'", collapse = ", "),
      "."
    )
  }

  user_genome_arg <- genome

  if (identical(genome_label, "user_define")) {
    if (is.null(user_genome_arg)) {
      stop("When genome_label is 'user_define', argument `genome` must be a BSgenome object.")
    }
    if (!methods::is(user_genome_arg, "BSgenome")) {
      stop("When genome_label is 'user_define', `genome` must inherit from BSgenome.")
    }
    genome <- user_genome_arg
  } else {
    if (!is.null(user_genome_arg) && methods::is(user_genome_arg, "BSgenome")) {
      warning(
        "Argument `genome` is ignored unless genome_label is \"user_define\".",
        call. = FALSE
      )
    }
    if (genome_label == "hg19") {
      genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    } else if (genome_label == "hg38") {
      genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    } else if (genome_label == "hg38-NCBI") {
      genome <- BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38
    }
  }

  genome_name <- unique(GenomeInfoDb::seqinfo(genome)@genome)
  if (length(genome_name) != 1L || is.na(genome_name)) {
    stop(
      "BSgenome seqinfo must have a single non-NA genome assembly name ",
      "(GenomeInfoDb::seqinfo(genome)@genome)."
    )
  }
  
  ############################################################################
  # Read bam into galp
  ############################################################################
  
  galp_param <- Rsamtools::ScanBamParam(flag = galp_flag, 
                                        what = galp_what,
                                        tag = galp_tag,
                                        mapqFilter = galp_mapqFilter, ...)
  
  galp <- bam_to_galp2(bamfile = bamfile, 
                       use_names = use_names,
                       chromosome_to_keep = chromosome_to_keep,
                       strand_mode = strand_mode, 
                       genome = genome_name,
                       param = galp_param
  ) 
  
  ############################################################################
  # Remove outward facing pairs
  ############################################################################
  
  galp <- remove_outward_facing_readpairs(galp)
  
  
  ############################################################################
  # Save galp object if given outdir 
  ############################################################################
  if(!isFALSE(outdir)) {
    
    bamfile_no_suffix <- gsub(bamfile, ".bam", "")
    out_rds_file_name <- paste0(bamfile_no_suffix, "_galp.rds")
    saveRDS(object = galp, file = file.path(outdir, out_rds_file_name))
    message("Saved ")
    message(out_rds_file_name)
  }
  
  return(galp)
  
  

}
