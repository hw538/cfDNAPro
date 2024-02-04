
#' Summarise descriptive Bam stats
#'
#' @import Biostrings
#' @import Rsamtools
#' @import dplyr
#' @param bamfile  Bam file
#' @param total_count Boolean. default = TRUE, which means calculating the total 
#' number of reads. 
#' @param total_mapped_count Boolean. Default = TRUE, which mean calculating the
#' number of mapped reads, and these reads must have mate reads.
#' @param chrM_count 
#' @param duplicate_count 
#' @param coverage_by_mapped 
#' @param genome_length_bp 
#' @param gc_metrics
#' @param loci_coverage_metrics
#' @param customized_count 
#' @param customized_count_mapqFilter 
#' @param customized_count_isPaired 
#' @param customized_count_isProperPair 
#' @param customized_count_isUnmappedQuery 
#' @param customized_count_hasUnmappedMate 
#' @param customized_count_isMinusStrand 
#' @param customized_count_isMateMinusStrand 
#' @param customized_count_isFirstMateRead 
#' @param customized_count_isSecondMateRead 
#' @param customized_count_isSecondaryAlignment 
#' @param customized_count_isNotPassingQualityControls 
#' @param customized_count_isDuplicate 
#' @param customized_count_isSupplementaryAlignment 
#' @param ... 
#'
#' @return tibble object (i.e. a dataframe)
#' @export
#' @author Haichao Wang
#'
#' @examples
#' \dontrun{
#' 
#'  summarizeBam(bamfile = "/path/to/bamfile.bam")
#' }
#' 

summarizeBam <- function(bamfile = NULL,
                         total_count = TRUE,
                         total_mapped_count = TRUE,
                         chrM_count = TRUE,
                         duplicate_count = TRUE,
                         coverage_by_mapped = TRUE,
                         genome_length_bp = 3200000000,
                         gc_metrics = FALSE,
                         loci_coverage_metrics = FALSE,
                         customized_count = FALSE,
                         customized_count_mapqFilter=NA_integer_,
                         customized_count_isPaired = NA, 
                         customized_count_isProperPair = NA, 
                         customized_count_isUnmappedQuery = NA, 
                         customized_count_hasUnmappedMate = NA, 
                         customized_count_isMinusStrand = NA, 
                         customized_count_isMateMinusStrand = NA,
                         customized_count_isFirstMateRead = NA, 
                         customized_count_isSecondMateRead = NA, 
                         customized_count_isSecondaryAlignment = NA, 
                         customized_count_isNotPassingQualityControls = NA,
                         customized_count_isDuplicate = NA, 
                         customized_count_isSupplementaryAlignment = NA,
                         ...) {
  
  
  #if (is.null(bamfile))
  #  bamfiles <- list.files(ifelse(is.null(path), '.', path),
  #                         pattern=sprintf('%s$', ext), full.names=TRUE)
  #if (length(bamfile) == 0L)
  #  stop('No files to process.')
  
  summary_metrics <- tibble::tibble(file = bamfile) 
  
  
  if(total_count) {
    
    message("Count total reads (mapped + unmapped)...")
    
    total_count <- bam_count(bamfile = bamfile)
    
    summary_metrics <- dplyr::full_join(summary_metrics, total_count, 
                                        by = "file")
  }
  
  
  if(total_mapped_count) {
    
    message("Count total mapped reads...")
    param <- Rsamtools::ScanBamParam(flag = 
                                       Rsamtools::scanBamFlag(isPaired = TRUE, 
                                                              isProperPair = NA, 
                                                              isUnmappedQuery = FALSE, 
                                                              hasUnmappedMate = FALSE, 
                                                              isMinusStrand = NA, 
                                                              isMateMinusStrand = NA,
                                                              isFirstMateRead = NA, 
                                                              isSecondMateRead = NA, 
                                                              isSecondaryAlignment = NA, 
                                                              isNotPassingQualityControls = NA,
                                                              isDuplicate = NA, 
                                                              isSupplementaryAlignment = NA
                                                              ))
    
    total_mapped_count <- bam_count(bamfile = bamfile, param = param) %>%
      dplyr::rename(n_read_mapped = n_read, 
                    n_nucleotide_mapped = n_nucleotide)
    
    summary_metrics <- dplyr::full_join(summary_metrics, total_mapped_count, 
                                        by = "file")
    
  }
  
  if(duplicate_count) {
    
    message("Count duplicate reads...")
    
    param <- Rsamtools::ScanBamParam(flag = scanBamFlag(isDuplicate = TRUE))
    
    duplicate_count <- bam_count(bamfile = bamfile, param = param) %>%
      dplyr::rename(n_read_duplicate = n_read, 
                    n_nucleotide_duplicate = n_nucleotide)
    
    summary_metrics <- dplyr::full_join(summary_metrics, duplicate_count, 
                                        by = "file")
  }
  
  
  if(chrM_count) {
    
    message("Count chrM reads...")
    
    chrM_count <- chr_count(bamfile = bamfile, chr = c("M", "chrM", "MT") ) %>%
      dplyr::rename(n_read_mapped_chrM = chr_count_mapped)
    summary_metrics <- dplyr::full_join(summary_metrics, chrM_count, 
                                        by = "file")
  }
  
  
  
  if(coverage_by_mapped) {
    
    if (isFALSE(total_mapped_count)) {
      stop("Please set total_mapped_count = TRUE in order to calculate coverage.")
    }
    
    
    message("Calculate coverage...")
    summary_metrics <- dplyr::mutate(summary_metrics, 
                                     coverage_by_mapped_reads = n_nucleotide_mapped / !!genome_length_bp )
  }
  
  if(gc_metrics) {
    
    message("Calculate GC metrics...")
    
    gc_metrics <- gc_count(bamfile = bamfile)
    
    summary_metrics <- dplyr::full_join(summary_metrics, gc_metrics, 
                                        by = "file")
  }
  
  
  if(loci_coverage_metrics) {
    
    message("Calculate loci coverage metrics...")
    
    ans <- get_loci_cov(bamfile = bamfile)
    
    summary_metrics <- dplyr::full_join(summary_metrics, ans, 
                                        by = "file")
  }
  
  
  
  
  if(customized_count) {
    
    message("Count reads with customized filtering criteria...")
    
    flag_customized <- scanBamFlag(isPaired = customized_count_isPaired, 
                        isProperPair = customized_count_isProperPair, 
                        isUnmappedQuery = customized_count_isUnmappedQuery, 
                        hasUnmappedMate = customized_count_hasUnmappedMate, 
                        isMinusStrand = customized_count_isMinusStrand, 
                        isMateMinusStrand = customized_count_isMateMinusStrand,
                        isFirstMateRead = customized_count_isFirstMateRead, 
                        isSecondMateRead = customized_count_isSecondMateRead, 
                        isSecondaryAlignment = customized_count_isSecondaryAlignment, 
                        isNotPassingQualityControls = customized_count_isNotPassingQualityControls,
                        isDuplicate = customized_count_isDuplicate, 
                        isSupplementaryAlignment = customized_count_isSupplementaryAlignment)
    
    mapqFilter_customized <- customized_count_mapqFilter
    
    param <- Rsamtools::ScanBamParam(flag = flag_customized,
                                     mapqFilter = mapqFilter_customized)
    
    customized_count <- bam_count(bamfile = bamfile, param = param) %>%
      dplyr::rename(n_read_customized = n_read, 
                    n_nucleotide_customized = n_nucleotide)
    
    summary_metrics <- dplyr::full_join(summary_metrics, customized_count, 
                                        by = "file")
  }
  
  message("Done!")
  return(summary_metrics)
  
  
}


# add synonym function

#' @rdname summariseBam
#' @export

summariseBam <- summarizeBam



###############################################################################
# helpers for summariseBam function
###############################################################################

chr_count <- function(bamfile, chr) {
  library(Rsamtools)
  library(tidyverse)
  a <- Rsamtools::idxstatsBam(bamfile)
  b <- a %>% 
    dplyr::filter(seqnames %in% as.vector(chr)) %>% 
    dplyr::select(mapped) %>%
    dplyr::rename(chr_count_mapped = mapped) %>%
    dplyr::mutate(file = bamfile)
  return(b)
}


bam_count <- function(bamfile, param , ...){
  
  if(missing(param)) {
    param <- Rsamtools::ScanBamParam()
  } 
  
  suppressWarnings(
    result <- Rsamtools::countBam(file = bamfile, param = param)  %>%
      tibble::as_tibble() %>%
      dplyr::select(.data$file, .data$records, .data$nucleotides)  %>%
      dplyr::mutate(file = bamfile)  %>%
      dplyr::rename(
        n_read = records,
        n_nucleotide = nucleotides
      )
  )
  
  result$n_read <- as.numeric(result$n_read)
  result$n_nucleotide = as.numeric(result$n_nucleotide)
  return(result)
}


gc_count <- function(bamfile){
  
  gcFunction <- function(x){
    alf <- alphabetFrequency(x, as.prob=FALSE)
    rowSums(alf[, c("C", "G")])
  }
  
  atgcFunction <- function(x){
    alf <- Biostrings::alphabetFrequency(x, as.prob=FALSE)
    rowSums(alf[, c("A", "T", "G", "C")])
  }
  param <- ScanBamParam(what="seq")
  seqs <- scanBam(bamfile, param=param)
  
  readGC <- gcFunction(seqs[[1]][["seq"]])
  readATGC <- atgcFunction(seqs[[1]][["seq"]])
  
  per_read_gc <- readGC / readATGC
  per_bam_gc <- sum(readGC) / sum(readATGC)
  mean_read_gc <- mean(per_read_gc, na.rm = TRUE)
  sd_read_gc <- sd(per_read_gc, na.rm = TRUE)
  median_read_gc <- median(per_read_gc, na.rm = TRUE)
  
  result <- tibble::tibble(
    file = bamfile,
    per_bam_gc = per_bam_gc,
    mean_read_gc = mean_read_gc,
    sd_read_gc = sd_read_gc,
    median_read_gc = median_read_gc
  )
  
  return(result)
}


# whole genome level metrics
get_loci_cov <- function(
    bamfile,
    chr = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
            "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", 
            "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
            "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"),
    quantile_positions = c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99)){
  
  res_cov <- GenomicAlignments::coverage(bamfile)
  res_cov_filtered <- res_cov[names(res_cov) %in% chr]
  combined_rle <- Reduce(c, res_cov_filtered) %>% as.numeric()
  whole_genome_mean_cov <- mean(combined_rle, na.rm = TRUE)
  whole_genome_median_cov <- median(combined_rle, na.rm = TRUE)
  #whole_genome_sd_cov <- sd(combined_rle, na.rm = TRUE)
  cov_metrics <- quantile(combined_rle, probs = quantile_positions) %>%
    as.data.frame() %>%
    t() %>% 
    tibble::as_tibble()
  cov_metrics <- cov_metrics %>%
    add_column(all_loci_mean_cov = whole_genome_mean_cov, 
               all_loci_median_cov = whole_genome_median_cov
    ) %>%
    add_column(file = bamfile)
  return(cov_metrics)
}


# loci level metrics
get_loci_cov_for_each_chr <- function(
    bamfile,
    chr = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
            "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", 
            "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
            "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"),
    quantile_positions = c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99)){
  
  res_cov <- GenomicAlignments::coverage(bamfile)
  res_cov_filtered <- res_cov[names(res_cov) %in% chr]
  #combined_rle <- Reduce(c, res_cov_filtered) %>% as.numeric()
  mean_cov <- mean(res_cov_filtered, na.rm = TRUE)
  median_cov <- median(res_cov_filtered, na.rm = TRUE)
  sd_cov <- sd(res_cov_filtered, na.rm = TRUE)
  cov_metrics <- quantile(res_cov_filtered, probs = quantile_positions) %>%
    tibble::as_tibble(rownames = "chr")
  cov_metrics <- cov_metrics %>%
    add_column(all_loci_mean_cov = mean_cov, 
               all_loci_median_cov = median_cov,
               all_loci_sd_cov = sd_cov
    ) %>%
    add_column(file = bamfile)
  return(cov_metrics)
}



