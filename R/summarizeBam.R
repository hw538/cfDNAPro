
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
  
  
  result <- Rsamtools::countBam(file = bamfile, param = param)  %>%
    tibble::as_tibble() %>%
    dplyr::select(.data$file, .data$records, .data$nucleotides)  %>%
    dplyr::mutate(file = bamfile)  %>%
    dplyr::rename(
      n_read = records,
      n_nucleotide = nucleotides
    )
  
  
  result$n_read <- as.integer(result$n_read)
  result$n_nucleotide = as.integer(result$n_nucleotide)
  return(result)
}




#' Summarise descriptive Bam stats
#'
#' @param bamfile 
#' @param total_count 
#' @param total_mapped_count 
#' @param chrM_count 
#' @param duplicate_count 
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

summarizeBam <- function(bamfile,
                         total_count = TRUE,
                         total_mapped_count = TRUE,
                         chrM_count = TRUE,
                         duplicate_count = TRUE,
                         ...) {
  
  
  summary_metrics <- tibble::tibble(file = bamfile) 
  
  
  if(total_count) {
    
    total_count <- bam_count(bamfile = bamfile)
    
    summary_metrics <- dplyr::full_join(summary_metrics, total_count, 
                                        by = "file")
  }
  
  
  if(total_mapped_count) {
    
    param <- Rsamtools::ScanBamParam(flag = scanBamFlag(isUnmappedQuery = FALSE))
    
    total_mapped_count <- bam_count(bamfile = bamfile, param = param) %>%
      dplyr::rename(n_read_mapped = n_read, 
                    n_nucleotide_mapped = n_nucleotide)
    
    summary_metrics <- dplyr::full_join(summary_metrics, total_mapped_count, 
                                        by = "file")
    
  }
  
  if(duplicate_count) {
    
    param <- Rsamtools::ScanBamParam(flag = scanBamFlag(isDuplicate = TRUE))
    
    duplicate_count <- bam_count(bamfile = bamfile, param = param) %>%
      dplyr::rename(n_read_duplicate = n_read, 
                    n_nucleotide_duplicate = n_nucleotide)
    
    summary_metrics <- dplyr::full_join(summary_metrics, duplicate_count, 
                                        by = "file")
  }
  
  
  if(chrM_count) {
    
    chrM_count <- chr_count(bamfile = bamfile, chr = c("M", "chrM", "ChrM") ) %>%
      dplyr::rename(chrM_count_mapped = chr_count_mapped)
    summary_metrics <- dplyr::full_join(summary_metrics, chrM_count, 
                                        by = "file")
  }
  
  return(summary_metrics)
  
  
}

#' @rdname summariseBam
#' @export

summariseBam <- summarizeBam




