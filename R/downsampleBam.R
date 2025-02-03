#' Randomly Downsample BAM File to a Target Depth
#'
#' Downsampling a BAM file to achieve a specified sequencing depth. Exports the
#' data as a BAM file and/or an RDS file containing a `GAlignmentPairs` object.
#'
#' @importFrom rtracklayer export
#' @importFrom Rsamtools BamFile 
#'
#' @param bamfile Path to a BAM file.
#' @param genome Supported values: "hg19", "hg38", "GRCh38", or NA. Sets seqinfo.
#' @param input_depth Target depth, inferred if not supplied.
#' @param input_mapped_reads_count Mapped reads count to infer depth if not supplied.
#' @param output_depth Target output depth after downsampling, default is 0.1.
#' @param output_type Output file formats, default is c("bam", "rds").
#' @param output_dir Directory for output files, defaults to input file directory.
#' @param return_result If TRUE (default), returns result as an R object.
#' @param ... Additional parameters for `bam_to_galp2` function.
#'
#' @return `GAlignmentPairs` object with the downsampled reads if return_result.
#' @export 
#'
#' @examples
#' \dontrun{
#'   # Downsample a BAM file to 10% of its original depth
#'   result <- downsampleBam("/path/to/your.bam", genome = "hg38",
#'                           output_depth = 0.1, output_type = c("bam", "rds"),
#'                           return_result = TRUE)
#' }

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




