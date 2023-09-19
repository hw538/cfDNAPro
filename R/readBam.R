#' Read bam file into a curated GRanges object
#' @import magrittr
#' @import GenomeInfoDb
#' @import GenomicAlignments
#' @import S4Vectors
#' @importFrom Rsamtools scanBamFlag ScanBamParam
#' @importFrom GenomicRanges GRanges seqnames
#' @importFrom IRanges IRanges
#' @import(dplyr)
#' @import(ggplot2)
#' @import(magrittr)
#' @import(plyranges)
#' @import(purrr)
#' @import(stringr)
#' @import(tibble)
#' @import(tidyr)
#' @import(data.table)
#' @import(IRanges)
#' @import(gsubfn)
#' @import(stringi)
#' @import(reshape2)
#' @import(data.table)
#' @import(DT)
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
#' @param strand_mode The strand_mode = 1 means
#'    the strand of the pair is strand of its first alignment.
#'    For More Details please see GenomicAlignments docs.
#' @param chromosome_to_keep Should be a character vector containing the
#'    seqnames to be kept in the GRanges object or boolean FALSE.
#'    FALSE means not filtering.
#'    Default is paste0("chr", 1:22).
#' @param use_names
#' @param galp_flag
#' @param galp_what A character vector naming the fields to return
#' Rsamtools::scanBamWhat() returns a vector of available fields.
#' Fields are described on the Rsamtools::scanBam help page.
#' @param galp_tag
#' @param galp_mapqFilter
#' @param mutation_file An optional file containing a list of mutations.
#' The fragments that overlap the mutation loci are mutationally annotated.
#' @param denovo If TRUE:
#' will use pileup to find mismatches and use them as mutations.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return This function returns curated GRanges object.
#' @export
#' @author Haichao Wang
#' @author Paulius D. Mennea
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
                     outdir = NULL,
                     galp_flag =  Rsamtools::scanBamFlag(
                       isPaired = TRUE,
                       isDuplicate = FALSE,
                       isSecondaryAlignment = FALSE,
                       isUnmappedQuery = FALSE,
                       isSupplementaryAlignment = FALSE),
                     galp_what =  c("cigar", "mapq", "isize"),
                     galp_tag = c("NM", "MD"),
                     galp_mapqFilter = 30,
                     mutation_file = NULL,
                     denovo = FALSE,
                     denovo_max_depth = 1000,
                     denovo_min_base_quality = 10,
                     denovo_min_nucleotide_depth = 3,
                     denovo_min_minor_allele_depth = 3,
                     outer_downstream_lim = 1000,
                     outer_upstream_lim = 1000,
                     ...) {

  ############################################################################
  # Check parameters
  ############################################################################

  # Check if the input bam is paired-end
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
  # Fragmentomic Annotation
  ############################################################################
  galp_param <- Rsamtools::ScanBamParam(
                                  flag = galp_flag,
                                  what = galp_what,
                                  tag = galp_tag,
                                  mapqFilter = galp_mapqFilter)

  galp <- bam_to_galp2(bamfile = bamfile,
                        use_names = use_names,
                        chromosome_to_keep = chromosome_to_keep,
                        strand_mode = strand_mode,
                        genome = genome_name,
                        param = galp_param)

  # Remove indels for mutational annotation
  if (!is.null(mutation_file) || denovo == "TRUE") {

    # Identify indel fragments
      has_indels <- grepl("[ID]", cigar(galp@first)) |
                    grepl("[ID]", cigar(galp@last))

    # Remove rows with indels
      galp <- galp[!has_indels]
  }

  ###########################################################################
  # Remove outward facing pairs
  ###########################################################################
  galp <- remove_outward_facing_readpairs(galp)

  ###########################################################################
  # Curate starts and ends
  ###########################################################################
  fragmentwise <- curate_start_and_end(galp = galp)

  names(fragmentwise) <- names(galp)

  seqlengths(fragmentwise) <- seqlengths(genome)[1:22]

  genome(fragmentwise) <- genome_name

  ###########################################################################
  # Remove out-of-bound reads
  ###########################################################################
  # Remove out-of-bound fragments and sort the galp
  frag <- remove_out_of_bound_reads(fragmentwise) %>% GenomicRanges::sort()

  ###########################################################################
  # Optional Mutational Annotation
  ###########################################################################

  if (!is.null(mutation_file)) {

  # Mutation File used to obtain fragment-level mutational information
  frag <- mutationFile(
                      bamfile = bamfile,
                      use_names = use_names,
                      chromosome_to_keep = chromosome_to_keep,
                      strand_mode = strand_mode,
                      genome = genome,
                      galp_flag = galp_flag,
                      galp_what = galp_what,
                      galp_tag = galp_tag,
                      galp_mapqFilter = galp_mapqFilter,
                      mutation_file = mutation_file,
                      outer_downstream_lim = outer_downstream_lim,
                      outer_upstream_lim = outer_upstream_lim,
                      frag = frag
                      )

  } else if ((denovo == "TRUE")) {

  # If no mutation file is provided, generate it by calling mismatches
  # Rsamtools::pileup() is used to call mismatches
  frag <- denovoPileup(
                      bamfile = bamfile,
                      use_names = use_names,
                      strand_mode = strand_mode,
                      galp_what = galp_what,
                      galp_tag = galp_tag,
                      genome = genome,
                      chromosome_to_keep = chromosome_to_keep,
                      galp_flag = galp_flag,
                      max_depth = denovo_max_depth,
                      min_base_quality = denovo_min_base_quality,
                      min_mapq = galp_mapqFilter,
                      min_nucleotide_depth = denovo_min_nucleotide_depth,
                      min_minor_allele_depth = denovo_min_minor_allele_depth,
                      outer_downstream_lim = outer_downstream_lim,
                      outer_upstream_lim = outer_upstream_lim,
                      frag
                      )

  }

  ###########################################################################
  # Saving RDS file
  ###########################################################################
  bamfile_base <- basename(bamfile)
  bamfile_no_suffix <- stringr::str_remove(bamfile_base, "\\.bam$")
  out_rds_file_name <- paste0(bamfile_no_suffix, "_GRanges_clean.rds")

  if (!is.null(outdir)) {
    saveRDS(object = frag, file = file.path(outdir, out_rds_file_name))
    message("Saved in specified directory: ")
    message(out_rds_file_name)
  } else {
  # Save in the current directory if outdir is not provided
    saveRDS(object = frag, file = out_rds_file_name)
    message("Saved in current directory: ")
    message(out_rds_file_name)
  }

# Return the data
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
                         param = galp_param,
                         chromosome_to_keep = FALSE,
                         strand_mode = 1,
                         genome = NA_character_) {
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
  galp3 <- galp2[GenomicAlignments::strand(galp2) != "*"]


  return(galp3)

}

#supporting function which is internal
bam_to_galp3 <- function(bamfile,
                         use_names = TRUE,
                         param = galp_param,
                         chromosome_to_keep = FALSE,
                         strand_mode = 1,
                         genome = NA_character_
                         ) {
  # Check parameters
  stopifnot(file.exists(bamfile))
  stopifnot(isSingleStringOrNA(genome) || is(genome, "Seqinfo"))

  # Read bam into galp
  message("Reading bam into galp...")

  galp <- readGAlignmentPairs(file = bamfile,
                              use.names = use_names,
                              strandMode = strand_mode,
                              param = param,
                              with.which_label = TRUE)
  # add genome information
  if (isSingleStringOrNA(genome)) {
    genome <- Seqinfo(genome = genome)
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
  galp3 <- galp2[GenomicAlignments::strand(galp2) != "*"]

  return(galp3)

}

# mage granges for mismatch and  mutation positions
make_granges <- function(muts) {
  which_processed <- GRanges(seqnames = muts[, 1],
          ranges = IRanges(start = muts[, 2], end = muts[, 3]))

  return(which_processed)
}

