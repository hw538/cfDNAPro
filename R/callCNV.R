

#' call copy number variation using QDNAseq utilities
#' @import QDNAseq
#' @param bamfile string. path to a single bam file.
#' @param bin_size numerical. available options are 
#' @param genome_label 
#'
#' @return QDNAseqCopyNumbers object
#' @export
#'
#' @examples

callCNV <- function(bamfile,
                    bin_size = 1000,
                    genome_label = "hg19"
                    ){
  
  bins <- QDNAseq::getBinAnnotations(binSize = bin_size, genome = genome_label)
  
  count <- QDNAseq::binReadCounts(bins, bamfiles = bamfile) %>%
    QDNAseq::applyFilters(residual=TRUE, blacklist=TRUE) %>%
    QDNAseq::estimateCorrection()
  
  
  seg <- QDNAseq::correctBins(count) %>%
    QDNAseq::normalizeBins() %>%
    QDNAseq::smoothOutlierBins() %>%
    QDNAseq::segmentBins(transformFun="sqrt") %>%
    QDNAseq::normalizeSegmentedBins() %>%
    QDNAseq::callBins()
  
  return(seg)
}
  