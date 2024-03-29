
#' Randomly downsample bam file to a target depth
#' 
#' @importFrom rtracklayer export
#' @importFrom Rsamtools BamFile 
#'
#' @param bamfile String. A bam file. 
#' @param genome Single String or NA. Supported values: "hg19", "hg38", or "GRCh38".
#' Used for setting the seqinfo in the resulting GALP object.
#' @param input_depth Numerical. If not supplied, infer the depth based on the 
#' total mapped paired reads in the input bam file.
#' @param input_mapped_reads_count Numerical. If not supplied, infer the depth based on the 
#' total mapped reads in the input bam file.
#' @param output_depth Numerical. Default is 0.1, which mean 0.1x target depth 
#' after downsampling. 
#' @param output_type Vector. Default is c("bam", "rds"), which means saving 
#' the downsampled bam as bam file AND rds file (i.e. GAlignmentPairs object).
#' @param output_dir String. If not supplied, it will be set as the same folder 
#' with input bamfile. 
#' @param return_result Boolean. TRUE(default) means the result will be returned
#' as an R object to console or the designated object.
#' @param ... further parameters for bam_to_galp2 function.
#'
#' @return GAlignmentPairs object containg the downsampled reads.
#' @export 
#'
#' @examples

downsampleBam <- function(bamfile, 
                          genome = NA_character_,
                          input_depth = NULL, 
                          input_mapped_reads_count = NULL, 
                          output_depth = 0.1, 
                          output_type = c("bam", "rds"),
                          output_dir = NULL,
                          return_result = TRUE,
                          ...) {
  
  #----------------------------------------------------------------------------
  # check params
  #----------------------------------------------------------------------------
  if(is.null(input_depth) & is.null(input_mapped_reads_count)){
    
    bam_stats <- summarizeBam(bamfile = bamfile,
                              total_count = FALSE,
                              total_mapped_count = TRUE, 
                              coverage_by_mapped = TRUE,
                              chrM_count = FALSE, 
                              duplicate_count = FALSE, 
                              customized_count = FALSE) 
    
    input_depth <- bam_stats$coverage_by_mapped_reads[[1]]
    message("Input Bam depth: ", input_depth, " x.")
    
    
    input_mapped_reads_count <- bam_stats$n_read_mapped[[1]]
    message("Input Bam N mapped reads (excluding singletons ): ", input_mapped_reads_count, ".")
    
    
  }
  
  if(!is.null(input_mapped_reads_count) & is.null(input_depth)){
    
    bam_stats <- summarizeBam(bamfile = bamfile,
                              total_count = FALSE,
                              total_mapped_count = TRUE, 
                              coverage_by_mapped = FALSE,
                              chrM_count = FALSE, 
                              duplicate_count = FALSE, 
                              customized_count = FALSE) 
    
    input_depth <- bam_stats$coverage_by_mapped_reads[[1]]
    message("Input Bam depth: ", input_depth, " x.")
    
  }
  
  if(is.null(input_mapped_reads_count) & !is.null(input_depth)){
    
    bam_stats <- summarizeBam(bamfile = bamfile,
                              total_count = FALSE,
                              total_mapped_count = TRUE, 
                              coverage_by_mapped = FALSE,
                              chrM_count = FALSE, 
                              duplicate_count = FALSE, 
                              customized_count = FALSE) 
    
    input_mapped_reads_count <- bam_stats$n_read_mapped[[1]]
    message("Input Bam N mapped reads (must be paired): ", input_mapped_reads_count, ".")
    
  }
  
  if(is.null(output_dir)){
    output_dir <- dirname(bamfile)
  }
  
  if(input_depth < output_depth){
    stop("Not enought reads for targeted output depth/reads.")
  }
  
  #----------------------------------------------------------------------------
  # downsample 
  #----------------------------------------------------------------------------
  
  keep <- round((output_depth/input_depth) * input_mapped_reads_count * 0.5)
  message("Randomly selecting ", keep, " read-pairs...")
  
  
  #important to set 'what' param here because it needed to be save as bam file. 
  bam_galp_object <- bam_to_galp2(bamfile, 
                                  genome = genome,
                                  ...) %>% 
    remove_outward_facing_readpairs()
  bam_sample <- sample(bam_galp_object,keep,replace=FALSE)
  
  #----------------------------------------------------------------------------
  # save outputfile 
  #----------------------------------------------------------------------------
  if(!isFALSE(output_type)) {
    output_type <- stringr::str_to_lower(output_type)
    
    if("bam" %in% output_type){
      output_file_name <- paste(basename(bamfile), "_samp_", output_depth, "x.bam" ,sep = "")   
      output_file <- file.path(output_dir, output_file_name )
      message("Saving bam file...")
      rtracklayer::export(bam_sample, Rsamtools::BamFile(output_file))
      message(output_file, "     Saved.")
    }
    
    if("rds" %in% output_type){
      output_file_name <- paste(basename(bamfile), "_samp_", output_depth, "x.rds" ,sep = "")   
      output_file <- file.path(output_dir, output_file_name )
      message("Saving rds file...")
      saveRDS(object = bam_sample, file = output_file)
      message(output_file, "     Saved.")
    }
    
  }
  #----------------------------------------------------------------------------
  # return result 
  #----------------------------------------------------------------------------
  
  if(return_result) {
    return(bam_sample)
  }
  
  message("Done.")
  
}




