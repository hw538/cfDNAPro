
#' Call Copy Number Variation Using QDNAseq Utilities
#'
#' This function leverages QDNAseq to process BAM files for CNV analysis.
#' It applies a sequence of operations including binning, filtering,
#' correcting, normalizing, smoothing, segmenting, and calling to
#' generate a QDNAseqCopyNumbers object.
#'
#' @import QDNAseq
#' @importFrom dplyr filter
#' @importFrom magrittr `%>%`
#' @importFrom QDNAseq getBinAnnotations binReadCounts applyFilters 
#' @importFrom QDNAseq estimateCorrection correctBins normalizeBins 
#' @importFrom QDNAseq smoothOutlierBins segmentBins normalizeSegmentedBins callBins
#' @param bamfile A string specifying the path to a single BAM file.
#' @param bin_size A numerical value indicating the size of the bins (default is 1000).
#' @param genome_label A string specifying the genome label (default is "hg19").
#'
#' @return A `QDNAseqCopyNumbers` object containing the copy number variation results.
#' @export
#'
#' @examples
#' callCNV(bamfile = "path/to/your/file.bam", bin_size = 1000, genome_label = "hg19")
#'
callCNV <- function(
  bamfile,
  bin_size = 1000,
  genome_label = "hg19"
) {

  bins <- QDNAseq::getBinAnnotations(binSize = bin_size, genome = genome_label)

  count <- QDNAseq::binReadCounts(bins, bamfiles = bamfile) %>%
    QDNAseq::applyFilters(residual = TRUE, blacklist = TRUE) %>%
    QDNAseq::estimateCorrection()

  seg <- QDNAseq::correctBins(count) %>%
    QDNAseq::normalizeBins() %>%
    QDNAseq::smoothOutlierBins() %>%
    QDNAseq::segmentBins(transformFun = "sqrt") %>%
    QDNAseq::normalizeSegmentedBins() %>%
    QDNAseq::callBins()

  return(seg)
}