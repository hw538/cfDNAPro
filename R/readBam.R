#' Read bam file into a curated GRanges object
#' @import magrittr
#' @import GenomeInfoDb
#' @import GenomicAlignments
#' @import S4Vectors
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb seqlengths genome Seqinfo keepSeqlevels
#' @importFrom GenomicAlignments readGAlignmentPairs
#' @importFrom Rsamtools scanBamFlag ScanBamParam testPairedEndBam
#' @importFrom GenomicRanges GRanges seqnames sort granges
#' @importFrom stringr str_remove
#' @importFrom methods is
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
#' @param curate_start_and_end A logical (TRUE or FALSE) parameter
#'   for curating alignment start and end coordinates.
#'   The default value is TRUE.
#' @param outdir The path for saving rds file. Default is NA, i.e. not saving.
#' @param strand_mode The strand_mode = 1 means
#'    the strand of the pair is strand of its first alignment.
#'    For More Details please see GenomicAlignments docs.
#' @param chromosome_to_keep Should be a character vector containing the
#'    seqnames to be kept in the GRanges object or boolean FALSE.
#'    FALSE means not filtering.
#'    Default is paste0("chr", 1:22).
#' @param use_names Logical, assigns read names to GRanges, default TRUE.
#' @param galp_flag ScanBamFlag object for BAM file scanning flags.
#' @param galp_what A character vector of fields to return from the BAM file.
#'   The specific fields returned depend on the analysis context:
#'   - For analyses with mutational files or when using `call_mutations`,
#'     `galp_what` should include "cigar", "mapq", "isize", "seq", and "qual".
#'   - For fragmentomic analysis without mismatch mutational data,
#'     `galp_what` should include just "cigar", "mapq", and "isize".
#'   This parameter's values and their descriptions are available on the
#'   Rsamtools::scanBamWhat() function and detailed on the Rsamtools::scanBam
#'   help page.
#' @param galp_tag Vector specifying optional fields (tags) to retrieve.
#' @param galp_mapqFilter Minimum mapping quality to include a read.
#' @param galp_bqFilter Base quality threshold for filtering read pairs.
#'   When a mutation_file is provided, this parameter specifies the minimum
#'   base quality required at the mutation locus. Read pairs not meeting
#'   this threshold at the mutation site are discarded, irrespective of
#'   the base.
#' @param mutation_file An optional file containing a list of mutations.
#'    Fragments that overlap the mutation loci are mutationally annotated.
#' @param mut_fragments_only A logical (TRUE or FALSE) parameter
#'    that determines the type of genomic alignments
#'    to retrieve from the BAM file.
#'    If set to TRUE, the function will only retrieve alignments
#'    that overlap with specified mutation locations.
#'    If set to FALSE, it will retrieve all alignments
#'    that pass the selected filters.
#'    The default value is FALSE.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return This function returns a curated GRanges object.
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
readBam <- function(bamfile,
                    curate_start_and_end = TRUE,
                    use_names = TRUE,
                    chromosome_to_keep = paste0("chr", 1:22),
                    strand_mode = 1,
                    genome_label = "hg19",
                    outdir = NA,
                    galp_flag = Rsamtools::scanBamFlag(
                      isPaired = TRUE,
                      isDuplicate = FALSE,
                      isSecondaryAlignment = FALSE,
                      isUnmappedQuery = FALSE,
                      isSupplementaryAlignment = FALSE),
                    galp_what =  c("cigar", "mapq", "isize"),
                    galp_tag = c("NM", "MD"),
                    galp_mapqFilter = 40,
                    galp_bqFilter = 20,
                    mutation_file = NULL,
                    mut_fragments_only = FALSE,
                    ...) {

  ############################################################################
  # Check parameters
  ############################################################################

  # Check if the input bam is paired-end
  if (!Rsamtools::testPairedEndBam(bamfile)) {
    stop("Input is not paired-end Bam file.")
  }

  if (!genome_label %in% c("hg19", "hg38", "hg38-NCBI") || is.na(genome_label)) {
      stop("Only 'hg19', 'hg38', or 'hg38-NCBI' are accepted as genome labels.")
  }

  if (genome_label == "hg19") {

    genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

  } else if (genome_label == "hg38") {

    genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

  } else if (genome_label == "hg38-NCBI") {

    genome <- BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38

  }

  genome_name <- seqinfo(genome)@genome %>% unique()

  ###########################################################################
  # Read BAM into galp based on user-defined parameters
  ###########################################################################
  galp <- process_bam_file(bamfile = bamfile,
                           mutation_file = mutation_file,
                           mut_fragments_only = mut_fragments_only,
                           use_names = use_names,
                           chromosome_to_keep = chromosome_to_keep,
                           strand_mode = strand_mode,
                           genome_name = genome_name,
                           galp_flag = galp_flag,
                           galp_what = galp_what,
                           galp_tag = galp_tag,
                           galp_mapqFilter = galp_mapqFilter)

  # Verify if the GALP object contains no data
  if (length(galp) == 0) {
      stop("The GALP object is empty. Please revise the 'galp' parameter",
              " or other input arguments and re-execute the function.")
  } else {

  ###########################################################################
  # Remove outward facing pairs
  ###########################################################################
    galp <- remove_outward_facing_readpairs(galp)

  ###########################################################################
  # Curate starts and ends
  ###########################################################################
    if (curate_start_and_end) {

        fragmentwise <- curate_start_and_end(galp = galp)

        names(fragmentwise) <- names(galp)
        # Get sequence names from genome (adjust 'genome' as needed)
        genome_seqnames <- names(seqlengths(genome))

        # Create named vector with chromosome names as indices
        genome_indices <- setNames(seq_along(genome_seqnames), genome_seqnames)

        # Extract indices for required chromosomes
        chromosome_indices <- genome_indices[chromosome_to_keep]

        # Update sequence lengths and genome name
        seqlengths(fragmentwise) <- seqlengths(genome)[chromosome_indices]
        genome(fragmentwise) <- genome_name

    } else {

    message("Skipped curating start and end coordinates.")
    fragmentwise <- granges(galp)

    }

    ###########################################################################
    # Remove out-of-bound fragments and sort the galp
    ###########################################################################
    frag <- remove_out_of_bound_reads(fragmentwise) %>% GenomicRanges::sort()

  }

  ###########################################################################
  # Mutational Annotation of read-pairs/fragments (Optional)
  ###########################################################################
  # Check if mutation file is provided and use it for annotation
  if (!is.null(mutation_file)) {
      frag <- process_mutation_fragments(
        bamfile = bamfile,
        genome = genome,
        galp = galp,
        mutation_file = mutation_file,
        frag = frag,
        galp_bqFilter = galp_bqFilter,
        chromosome_to_keep = chromosome_to_keep,
        mut_fragments_only = mut_fragments_only
      )
  }

  ###########################################################################
  # Saving RDS file
  ###########################################################################
  if (is.na(outdir) || isFALSE(outdir)) {

    return(frag)

  } else {

    bamfile_no_suffix <- stringr::str_remove(bamfile, "\\.bam$")
    out_rds_file_name <- paste0(basename(bamfile_no_suffix), "_GRanges.rds")
    saveRDS(object = frag, file = file.path(outdir, out_rds_file_name))
    message("Saved GRanges object to")
    message(file.path(outdir, out_rds_file_name))

    return(frag)

  }

}

#-----------------------------------------------------------------------------

#' Function to read bam file
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
                         genome = NA_character_,
                         ...) {
  # Check parameters
  stopifnot(file.exists(bamfile))
  stopifnot(isSingleStringOrNA(genome) || is(genome, "Seqinfo"))

  # Read bam into galp
  message("Reading bam into galp...")

  galp <- GenomicAlignments::readGAlignmentPairs(file = bamfile,
                                                 use.names = use_names,
                                                 strandMode = strand_mode,
                                                 param = param,
                                                 with.which_label = FALSE)

  # Message indicating successful read and possibly some details about 'galp'
  message("BAM file has been read into galp.")

  # add genome information
  if (isSingleStringOrNA(genome)) {
    genome <- Seqinfo(genome = genome)[chromosome_to_keep]
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
