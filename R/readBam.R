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
#'    Default is paste0("chr", 1:22). Select "chrX", "chrY", or "chrM"
#'    explicitly to include sex or mitochondrial chromosomes. Mitochondrial
#'    aliases ("chrM", "chrMT", "MT", "M") are resolved against the BAM header
#'    and returned using the selected reference genome's name. Nuclear names
#'    must match the BAM header. BAM and reference sequence lengths must agree.
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
#' object <- readBam(bamfile = "/path/to/bamfile.bam",
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
                    genome = NULL,
                    ...) {

  ############################################################################
  # Check parameters
  ############################################################################

  # Check if the input bam is paired-end
  if (!Rsamtools::testPairedEndBam(bamfile)) {
    stop("Input is not paired-end Bam file.")
  }

  reference <- .cfDNAPro_resolve_genome(genome_label, genome, missing(genome_label))
  genome_label <- reference$label
  genome <- reference$genome

  genome_info <- seqinfo(genome)

  ###########################################################################
  # Read BAM into galp based on user-defined parameters
  ###########################################################################
  galp <- process_bam_file(bamfile = bamfile,
                           mutation_file = mutation_file,
                           mut_fragments_only = mut_fragments_only,
                           use_names = use_names,
                           chromosome_to_keep = chromosome_to_keep,
                           strand_mode = strand_mode,
                           genome_name = genome_info,
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
        # Preserve the validated BAM/reference lengths and circularity.
        seqinfo(fragmentwise) <- seqinfo(galp)[seqlevels(fragmentwise)]

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
  frag <- .cfDNAPro_attach_reference(frag, genome_label, genome)
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
#' @param bamfile Path to a paired-end BAM file.
#' @param use_names Retain read-pair names in the returned object.
#' @param param A ScanBamParam object controlling fields and read filters.
#' @param chromosome_to_keep Sequence names to retain, or FALSE for all sequences.
#' @param strand_mode Strand mode passed to readGAlignmentPairs().
#' @param genome Reference Seqinfo, an assembly label, or NA for no reference.
#' @param ... Additional arguments reserved for compatibility.
#' @keywords internal
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

  galp <- curate_bam_seqinfo(galp, genome, chromosome_to_keep)

  # strandMode should be one for downstream operations
  stopifnot(GenomicAlignments::strandMode(galp) == 1)

  message("Curating seqnames and strand information...")

  # remove read pairs without correct seqnames and strand information
  galp2 <- galp[!is.na(GenomicAlignments::seqnames(galp))]
  galp3 <- galp2[GenomicAlignments::strand(galp2) != "*"]


  return(galp3)

}

# Resolve mitochondrial aliases against the names in a BAM or mutation table.
# Nuclear chromosome names must match the input naming style exactly.
resolve_chromosome_selection <- function(chromosomes, available) {
  if (isFALSE(chromosomes)) return(FALSE)
  if (!is.character(chromosomes) || length(chromosomes) == 0L ||
      anyNA(chromosomes) || any(!nzchar(chromosomes))) {
    stop("'chromosome_to_keep' must be a non-empty character vector or FALSE.")
  }
  mt_names <- c("chrM", "chrMT", "MT", "M")
  available_mt <- intersect(available, mt_names)
  needs_alias <- chromosomes %in% mt_names & !chromosomes %in% available
  if (any(needs_alias) && length(available_mt) == 1L) {
    chromosomes[needs_alias] <- available_mt
  }
  unique(chromosomes)
}

# Filter using BAM names before reconciling the selected reference metadata.
curate_bam_seqinfo <- function(galp, genome, chromosome_to_keep) {
  selected <- resolve_chromosome_selection(chromosome_to_keep, seqlevels(galp))
  if (!isFALSE(selected)) {
    missing <- setdiff(selected, seqlevels(galp))
    present <- intersect(selected, seqlevels(galp))
    # A regional BAM may omit reference chromosomes from its header. Retain
    # the available selection, while rejecting names absent from the reference.
    absent_reference <- if (is(genome, "Seqinfo")) setdiff(missing, seqlevels(genome)) else missing
    if (!length(present) || length(absent_reference)) {
      stop("Chromosomes not found in the BAM header: ",
           paste(missing, collapse = ", "), ". Check the chromosome naming style.")
    }
    galp <- keepSeqlevels(galp, present, pruning.mode = "coarse")
  }

  # readBam supplies the already-loaded BSgenome Seqinfo, avoiding a second
  # reference lookup over the network. Retain string support for direct callers.
  if (isSingleStringOrNA(genome)) {
    if (is.na(genome)) return(galp)
    genome <- Seqinfo(genome = genome)
  }

  mt_names <- c("chrM", "chrMT", "MT", "M")
  bam_mt <- intersect(seqlevels(galp), mt_names)
  genome_mt <- intersect(seqlevels(genome), mt_names)
  if (length(bam_mt) && length(genome_mt)) {
    # hg19 has both chrM and chrMT with different sequences and lengths.
    # Preserve exact reference names before considering an alias.
    mt_targets <- vapply(bam_mt, function(chr) {
      if (chr %in% genome_mt) return(chr)
      if ("chrM" %in% genome_mt) return("chrM")
      if (length(genome_mt) == 1L) return(genome_mt)
      stop("The reference contains multiple mitochondrial sequences; ",
           "use an exact reference chromosome name.")
    }, character(1))
    if (anyDuplicated(unname(mt_targets))) {
      stop("Multiple mitochondrial aliases would map to the same reference ",
           "sequence; select one explicitly.")
    }
    galp <- GenomeInfoDb::renameSeqlevels(galp, stats::setNames(mt_targets, bam_mt))
  }

  common <- intersect(seqlevels(galp), seqlevels(genome))
  # merge rejects conflicting lengths, including hg19 chrM (16571) versus
  # a 16569-base mitochondrial reference. An alias is not a reference conversion.
  seqinfo(galp) <- merge(seqinfo(galp), genome[common])
  galp
}
