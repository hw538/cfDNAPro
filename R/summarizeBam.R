
#' Summarise descriptive Bam stats
#'
#' @param bamfile A Bam file
#' @param total_count Boolean. default = TRUE, which means calculating the total 
#' number of reads. 
#' @param total_mapped_count Boolean. Default = TRUE, which mean calculating the
#' number of mapped reads.
#' @param chrM_count 
#' @param duplicate_count 
#' @param coverage_by_mapped 
#' @param genome_length_bp 
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
                                       Rsamtools::scanBamFlag(isUnmappedQuery = FALSE))
    
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
    
    chrM_count <- chr_count(bamfile = bamfile, chr = c("M", "chrM", "ChrM") ) %>%
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


