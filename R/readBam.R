
#' Read bam file into a curated GRanges object
#' @importFrom  Rsamtools scanBamFlag ScanBamParam
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @import GenomicAlignments 
#' @import magrittr
#' @import BiocGenerics
#'
#' @param genome_label The Genome you used in the alignment. 
#'    Should be "hg19" or "hg38". Default is "hg19".
#' @param bamfile The bam file name.
#' @param outdir The path for saving rds file. Default is NA, i.e. not saving.
#' @param strand_mode Usually the strand_mode  = 1 means the First read is 
#'    aligned to positive strand. Details please see GenomicAlignments docs.
#' @param chromosome_to_keep Should be a character vector containing the 
#'    seqnames to be kept in the GRanges object. Default is paste0("chr", 1:22).
#' @param ... Further arguments passed to or from other methods.
#'
#' @return This function returns curated GRanges object.
#' @export
#' @author Haichao Wang
#'
#' @examples 
#' \dontrun{
#' 
#' object <- read_bam(bamfile = "/path/to/bamfile.bam", 
#'                    outdir = "./",
#'                    chromosome_to_keep = c("chr1", "chr2", "chr3"))
#' }
#' 

readBam <- function(
                     bamfile,
                     chromosome_to_keep = paste0("chr", 1:22),
                     strand_mode = 1,
                     genome_label = "hg19",
                     outdir = NA,
                     ...){
  
  ############################################################################
  # Check parameters
  ############################################################################
  
  if (is.na(outdir)) {message("GRanges object won't be saved as an RDS file.")}
  
  if (!genome_label %in% c("hg19", "hg38")) {
    stop("Only accept hg19 or hg38 as the genome_label...")
  }
  
  if (genome_label == "hg19") {
    
  genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  
  } else if (genome_label == "hg38") {
    
  genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  
  }
  
  ############################################################################
  # Read bam into galp
  ############################################################################
  
  galp <- bam_to_galp(bamfile = bamfile, 
                       chromosome_to_keep = chromosome_to_keep,
                       strand_mode = strand_mode) 
  
  ############################################################################
  # Remove outward facing pairs
  ############################################################################
  
  galp <- remove_outward_facing_readpairs(galp)
  
  
  ############################################################################
  # Curate starts and ends 
  ############################################################################
  
  fragmentwise <- curate_start_and_end(galp = galp) 
  
  names(fragmentwise) <- names(galp)
  seqlengths(fragmentwise) <- seqlengths(genome)[1:22]
  genome(fragmentwise) <- seqinfo(genome)@genome %>% unique()
  
  #############################################################################
  # remove out-of-bound reads
  #############################################################################
  
  # remove out-of-bound fragments and sort the galp
  frag <- remove_out_of_bound_reads(fragmentwise) %>% GenomicRanges::sort() 

  #############################################################################
  # Saving RDS file 
  #############################################################################
  
  if(!is.na(outdir)) {
    
    bamfile_no_suffix <- gsub(bamfile, ".bam", "")
    out_rds_file_name <- paste0(bamfile_no_suffix, "_GRanges_clean.rds")
    saveRDS(object = frag, file = file.path(outdir, out_rds_file_name))
    message(paste0("Saved", " ", out_rds_file_name, "."))
  }
  
  return(frag)
}
