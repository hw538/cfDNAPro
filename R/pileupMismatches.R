#' Pileup Mismatches from BAM Files
#'
#' This function conducts a genomic pileup on specified chromosomes or regions
#' from a BAM file. It leverages multiple cores for parallel processing using
#' Rsamtools pileup methods. User-defined parameters customize the operation.
#' Outputs a .tsv file with columns: chr, start, end, ref, alt. This file can
#' serve as input to the readBAM() function's mutation_file parameter.
#' @importFrom parallel mclapply detectCores makeCluster clusterEvalQ
#' @importFrom parallel clusterExport parLapply stopCluster
#' @importFrom Rsamtools pileup BamFile ScanBamParam PileupParam
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom dplyr select
#' @importFrom BSgenome getSeq
#' @importFrom utils write.table
#' @importFrom GenomeInfoDb seqlengths genome seqinfo
#'
#' @param bamfile The path to the BAM file or a BamFile object.
#' @param genome The BSgenome object corresponding to the BAM file's genome.
#' @param chromosome_to_keep Should be a character vector containing the
#'    seqnames to be kept in the GRanges object or boolean FALSE.
#'    FALSE means not filtering.
#'    Default is paste0("chr", 1:22).
#' @param galp_flag ScanBamFlag object for BAM file scanning flags.
#' @param pileup_params A list containing parameters for the pileup operation,
#'        including 'max_depth', 'min_base_quality', 'min_mapq',
#'        'min_nucleotide_depth', 'min_minor_allele_depth',
#'        'distinguish_strands', 'distinguish_nucleotides',
#'        'ignore_query_Ns', 'include_deletions', 'include_insertions',
#'        and optional binning parameters 'left_bins', 'query_bins',
#'        'cycle_bins'. For parameter descriptions, refer to
#'        `Rsamtools::pileup` documentation.
#' @param num_cores Number of processor cores to use for parallel processing.
#'        Defaults to all available cores minus one.
#' @param outfile The path and filename for the output .tsv file.
#' @return A dataframe of processed mutation data from specified chromosomes.

#' @examples
#' bamfile <- "/path/to/your.bam"
#' genome <- BSgenome.Hsapiens.UCSC.hg19
#' chromosomes <- c("chr1", "chr2", "chrX")
#' flags <- ScanBamFlag(isPaired=TRUE, isDuplicate=FALSE)
#' pileup_params <- list(max_depth=250, min_base_quality=20)
#' outfile <- "path/to/output.tsv"
#'
#' result <- pileupMismatches(
#'     bamfile=bamfile,
#'     genome=genome,
#'     chromosome_to_keep=chromosomes,
#'     galp_flag=flags,
#'     pileup_params=pileup_params,
#'     num_cores=3
#'     outfile=outfile
#' )
#' @export

pileupMismatches <- function(
    bamfile,
    genome,
    chromosome_to_keep,
    galp_flag,
    pileup_params,
    num_cores = detectCores() - 1,
    outfile
) {
  
  # Extract parameters from the list with default values if not specified
  params_list <- list(
    yieldSize = pileup_params$yieldSize %||% 10000,
    max_depth = pileup_params$max_depth %||% 250,
    min_base_quality = pileup_params$min_base_quality %||% 13,
    min_mapq = pileup_params$min_mapq %||% 30,
    min_nucleotide_depth = pileup_params$min_nucleotide_depth %||% 1,
    min_minor_allele_depth = pileup_params$min_minor_allele_depth %||% 0,
    distinguish_strands = pileup_params$distinguish_strands %||% FALSE,
    distinguish_nucleotides = pileup_params$distinguish_nucleotides %||% TRUE,
    ignore_query_Ns = pileup_params$ignore_query_Ns %||% FALSE,
    include_deletions = pileup_params$include_deletions %||% FALSE,
    include_insertions = pileup_params$include_insertions %||% FALSE
  )
  
  # Setup parallel processing
  cl <- makeCluster(num_cores)
  clusterEvalQ(cl, {
    library(Rsamtools)
    library(GenomicRanges)
    library(BSgenome)
    library(dplyr)
    library(IRanges)
  })
  
  clusterExport(cl,
                varlist = c("bamfile", "chromosome_to_keep",
                            "galp_flag", "genome", "params_list"),
                envir = environment())
  
  # Function to process each chromosome
  process_chromosome <- function(chr) {
    pileupFreq <- function(pileupres) {
      # Step 1: Get nucleotides
      nucleotides <- levels(pileupres$nucleotide)
      
      # Step 2: Split by seqnames
      res <- split(pileupres, pileupres$seqnames)
      
      # Step 3: Split by pos
      res <- lapply(res, function(x) {
        split(x, x$pos)
      })
      
      # Step 4: Process each positionsplit
      res <- lapply(res, function(positionsplit) {
        # Step 4.1: Calculate tablecounts for each nucleotide
        nuctab <- lapply(positionsplit, function(each) {
          chr <- as.character(unique(each$seqnames))
          pos <- as.character(unique(each$pos))
          tablecounts <- sapply(nucleotides, function(n) {
            sum(each$count[each$nucleotide == n])
          })
          c(chr, pos, tablecounts)
        })
        
        # Step 4.2: Combine nuctab into a data frame
        nuctab <- data.frame(do.call("rbind", nuctab),
                             stringsAsFactors = FALSE)
        rownames(nuctab) <- NULL
        nuctab
      })
      
      # Step 5: Combine res into a data frame
      res <- data.frame(do.call("rbind", res), stringsAsFactors = FALSE)
      rownames(res) <- NULL
      colnames(res) <- c("seqnames", "start", levels(pileupres$nucleotide))
      
      # Step 6: Convert columns 3 to ncol(res) to numeric
      res[3:ncol(res)] <- apply(res[3:ncol(res)], 2, as.numeric)
      
      # Step 7: Return the final result
      return(res)
    }
    
    # Extract chromosome length from the BSgenome object
    chr_length <- GenomeInfoDb::seqlengths(genome)[chr]
    
    # Check if chromosome exists in the BSgenome object
    if (is.na(chr_length)) {
      stop("Chromosome ", chr, " not found in genome data.")
    }
    
    bf <- BamFile(file = bamfile, yieldSize = params_list$yieldSize)
    sbp <- ScanBamParam(
      flag = galp_flag,
      mapqFilter = params_list$min_mapq,
      tagFilter = list(NM = c(1, 2, 3, 4, 5, 6, 7)),
      which = GRanges(seqnames = chr, ranges = IRanges(start = 1,
                                                       end = chr_length))
    )
    p_param <- PileupParam(
      max_depth = params_list$max_depth,
      min_base_quality = params_list$min_base_quality,
      min_mapq = params_list$min_mapq,
      min_nucleotide_depth = params_list$min_nucleotide_depth,
      min_minor_allele_depth = params_list$min_minor_allele_depth,
      distinguish_strands = params_list$distinguish_strands,
      distinguish_nucleotides = params_list$distinguish_nucleotides,
      ignore_query_Ns = params_list$ignore_query_Ns,
      include_deletions = params_list$include_deletions,
      include_insertions = params_list$include_insertions)

    message("Running pileup...")
    res <- pileup(bf, scanBamParam = sbp, pileupParam = p_param)

    message("Processing frequencies...")
    res <- pileupFreq(res)

    message("Processing mismatches...")
    res <- res[rowSums(res[, c("A", "C", "G", "T")] != 0) >= 2, ]

    # Process mismatch information
    res <- res[rowSums(res[, c("A", "C", "G", "T")] != 0) >= 2, ]
    
    loci <- res %>% select(seqnames, start)
    loci$end <- loci$start
    loci$end <- as.numeric(loci$end)
    loci$start <- as.numeric(loci$start)
    
    which_loci <- GRanges(
      seqnames = loci[, 1],
      ranges = IRanges(start = loci[, 2], end = loci[, 3])
    )
    
    ref_bases <- as.vector(BSgenome::getSeq(genome, names = which_loci))
    
    ref_bases_df <- data.frame(
      seqnames = loci$seqnames,
      start = loci$start,
      ref_base = ref_bases
    )
    
    res <- merge(res, ref_bases_df, by = c("seqnames", "start"))
    
    res$max1 <- apply(res[, c("A", "C", "T", "G")],
                     1,
                     function(x) names(sort(x, decreasing = TRUE)[1]))
    
    res$max2 <- apply(res[, c("A", "C", "T", "G")],
                     1,
                     function(x) names(sort(x, decreasing = TRUE)[2]))
    
    res$alt <- ifelse(res$max1 != res$ref_base & res$max2 != res$ref_base,
                     res$max1,
                     ifelse(res$max1 != res$ref_base,
                            res$max1,
                            ifelse(res$max2 != res$ref_base, res$max2, "NA")))
    
    res <- res %>% select(seqnames, start, ref_base, alt)
    
    colnames(res) <- c("chr", "pos", "ref", "alt")
    
    res$pos <- as.numeric(res$pos)    
    res <- res %>% select(chr, pos, ref, alt)
    
    # Filter out unnecessary chromosomes
    res <- res[
      res$chr %in% chromosome_to_keep,
    ]
    
    return(res)
  }
  
  # Process each chromosome in parallel
  results <- parLapply(cl, chromosome_to_keep, process_chromosome)
  
  # Stop the cluster
  stopCluster(cl)
  
  # Combine results from all chromosomes
  combined_results <- do.call(rbind, results)

  # Write to outfile as a .tsv file
  write.table(combined_results,
              file = outfile,
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE,
              quote = FALSE)
  
  # Return the combined results
  return(combined_results)
}