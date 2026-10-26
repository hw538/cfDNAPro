
#' Read bam file into GAlignmentPairs object
#' @import magrittr
#' @import GenomeInfoDb
#' @import GenomicAlignments
#' @import S4Vectors
#' @importFrom Rsamtools scanBamFlag ScanBamParam
#' @importFrom GenomicRanges GRanges seqnames
#' @importFrom IRanges IRanges
#'
#' @param genome_label The Genome you used in the alignment. 
#'    Should be "hg19", "hg38", "hg38-NCBI", or "user_define". Default is "hg19".
#'    Note: "hg19" will load BSgenome.Hsapiens.UCSC.hg19 package, which is 
#'    Full genome sequences for Homo sapiens (Human) as provided by 
#'    UCSC (hg19, based on GRCh37.p13) and stored in Biostrings objects; 
#'    "hg38" will load BSgenome.Hsapiens.UCSC.hg38 package, which is 
#'    Full genome sequences for Homo sapiens (Human) as provided by 
#'    UCSC (hg38, based on GRCh38.p13) and stored in Biostrings objects.
#'    "hg38-NCBI" will load BSgenome.Hsapiens.NCBI.GRCh38 package, which is 
#'    full genome sequences for Homo sapiens (Human) as provided by 
#'    NCBI (GRCh38, 2013-12-17) and stored in Biostrings objects.
#' @param genome A BSgenome object, required when
#'   \code{genome_label = "user_define"}. The BAM and reference must use
#'   matching sequence names and lengths. For built-in labels a supplied
#'   BSgenome is ignored with a warning. If \code{genome_label} is omitted,
#'   the older shorthand \code{genome = "hg38"} still selects that built-in
#'   reference. References and dependencies must be installed before offline use.
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
    genome = NULL,
    ...){
  
  ############################################################################
  # Check parameters
  ############################################################################
  
  # check if the input bam is paired-end
  if (!Rsamtools::testPairedEndBam(bamfile)) {
    stop("Input is not paired-end Bam file.")
  }
  
  reference <- .cfDNAPro_resolve_genome(genome_label, genome, missing(genome_label))
  genome_label <- reference$label
  genome <- reference$genome
  
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
                       genome = seqinfo(genome),
                       param = galp_param
  ) 
  
  ############################################################################
  # Remove outward facing pairs
  ############################################################################
  
  galp <- remove_outward_facing_readpairs(galp)
  galp <- .cfDNAPro_attach_reference(galp, genome_label, genome)
  
  
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
