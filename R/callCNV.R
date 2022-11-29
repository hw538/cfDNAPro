
bam <- "/Users/wang04/Documents/phd_projects/cfdnapro_debug/SLX-11379_D703_D504.mrkdup.BL.filtered.bam.markduplicates.bam.downsamp_0.1x.mrkdup.bam"

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
    QDNAseq::normalizeSegmentedBins()
  
  return(seg)
}
  