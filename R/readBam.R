
#' Read bam file into a curated GRanges object
#' @import magrittr
#' @import GenomeInfoDb
#' @import GenomicAlignments
#' @import S4Vectors
#' @importFrom Rsamtools scanBamFlag ScanBamParam
#' @importFrom GenomicRanges GRanges seqnames
#' @importFrom IRanges IRanges
#'
#' @param genome_label The Genome you used in the alignment. 
#'    Should be "hg19" or "hg38" or "hg38-NCBI". Default is "hg19". 
#'    Note: "hg19" will load BSgenome.Hsapiens.UCSC.hg19 package, which is 
#'    Full genome sequences for Homo sapiens (Human) as provided by 
#'    UCSC (hg19, based on GRCh37.p13) and stored in Biostrings objects; 
#'    "hg38" will load BSgenome.Hsapiens.UCSC.hg38 package, which is 
#'    Full genome sequences for Homo sapiens (Human) as provided by 
#'    UCSC (hg38, based on GRCh38.p13) and stored in Biostrings objects.
#'    "hg38-NCBI" will load BSgenome.Hsapiens.NCBI.GRCh38 package, which is 
#'    full genome sequences for Homo sapiens (Human) as provided by 
#'    NCBI (GRCh38, 2013-12-17) and stored in Biostrings objects.
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

readBam <- function(
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
                     ...){
  
  ############################################################################
  # Check parameters
  ############################################################################
  
  # check if the input bam is paired-end
  if (!Rsamtools::testPairedEndBam(bamfile)) {
    stop("Input is not paired-end Bam file.")
  }
  
  if (!genome_label %in% c("hg19", "hg38", "hg38-NCBI") | is.na(genome_label)) {
    stop("Only accept hg19 or hg38 or hg38-NCBI as the genome_label...")
  }
  
  if (genome_label == "hg19") {
    
  genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  
  } else if (genome_label == "hg38") {
    
  genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  
  } else if (genome_label == "hg38-NCBI") {
    
  genome <- BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38
  
  }
  
  genome_name <- seqinfo(genome)@genome %>% unique()
  
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
  # Curate starts and ends 
  ############################################################################
  
  fragmentwise <- curate_start_and_end(galp = galp) 
  
  names(fragmentwise) <- names(galp)
  seqlengths(fragmentwise) <- seqlengths(genome)[1:22]
  genome(fragmentwise) <- genome_name
  
  #############################################################################
  # remove out-of-bound reads
  #############################################################################
  
  # remove out-of-bound fragments and sort the galp
  frag <- remove_out_of_bound_reads(fragmentwise) %>% GenomicRanges::sort() 

  #############################################################################
  # Saving RDS file 
  #############################################################################
  
  if(!isFALSE(outdir) | !is.na(outdir)) {
    
    bamfile_no_suffix <- gsub(bamfile, ".bam", "")
    out_rds_file_name <- paste0(bamfile_no_suffix, "_GRanges_clean.rds")
    saveRDS(object = frag, file = file.path(outdir, out_rds_file_name))
    message("Saved ")
    message(out_rds_file_name)
  }
  
  return(frag)
}




#-----------------------------------------------------------------------------


#' @import magrittr
#' @import GenomeInfoDb
#' @import GenomicAlignments
#' @import S4Vectors
#' @import Rsamtools 

bam_to_galp2 <- function(bamfile,
                         use_names  = TRUE, 
                         param = Rsamtools::ScanBamParam(...), 
                         chromosome_to_keep = FALSE,
                         strand_mode = 1,
                         genome = NA_character_,
                         ...) {
  # Check parameters
  stopifnot(file.exists(bamfile))
  stopifnot(isSingleStringOrNA(genome) || is(genome, "Seqinfo"))
  
  # Read bam into galp
  message("Reading bam into galp...")
  
  galp <- readGAlignmentPairs(file = bamfile, 
                              use.names = use_names, 
                              strandMode = strand_mode, 
                              param = param)
  # add genome information
  if (isSingleStringOrNA(genome)) {
    genome <- Seqinfo(genome=genome)
  }
  seqinfo(galp) <- merge(seqinfo(galp), genome)
  
  # strandMode should be one for downstream operations
  stopifnot(GenomicAlignments::strandMode(galp) == 1)
  
  # only keep needed seqnames
  if (!isFALSE(chromosome_to_keep)) {
  galp <- keepSeqlevels(galp, chromosome_to_keep, pruning.mode = "coarse")
    
  }
  
  message("Curating seqnames and strand information...")
  # remove read pairs without correct seqnames and strand information
  galp2 <- galp[!is.na(GenomicAlignments::seqnames(galp))]
  galp3 <- galp2[GenomicAlignments::strand(galp2) != '*']
  
  
  return(galp3)
  
}
