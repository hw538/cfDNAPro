
callMotif <- function(fragment_obj, 
                      genome_label  = "hg19", 
                      motif_type = "s",
                      motif_length = 3L, 
                      ...) {
  
  motif_length <- as.numeric(motif_length)
  
  message("Started to extract ", paste0(motif_type, motif_length)," motif..." )
  message("You reference genome was set to ", genome_label, "!")
  
  # select the correct reference genome.
  
  if (genome_label == "hg19") {
    
    bsgenome_obj <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    
  }else if (genome_label == "hg38") {
    
    bsgenome_obj <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    
  }else if (genome_labe == "hg38-NCBI") {
    
    bsgenome_obj <-BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38
    
  }
  
  # extract the motif to a vector
  
  motif <- get_motif(obj = fragment_obj, 
                      genome = bsgenome_obj, 
                      motif_type = motif_type,
                      motif_length = motif_length
                      )
  
  # summarize the vector and convert to tibble format 
  
  result <- table(motif) %>% 
    tibble::as_tibble() 
  
  # validate results
  
  # Create a vector of elements
  pos <- c("C", "G", "A", "T")
  base_index <- seq.int(1, motif_length, by = 1)
  
  letter_list <- lapply(base_index, function(x, y) return(y), y = pos) %>%
    setNames(paste0("base", base_index))
  
  motif_ref <- purrr::cross_df(letter_list)  %>%
    tidyr::unite(col = "motif", all_of(base_index), sep = "")
  
  # report abnormal motifs. 
  # i.e. ambiguous base in motif and/or missing motif
  
  #handle ambiguous motif(s)
  
  ambiguous_motif <- dplyr::anti_join(result, motif_ref, by = "motif")
  missing_motif <- dplyr::anti_join(motif_ref, result, by = "motif")
  
  if(nrow(ambiguous_motif) != 0) {
    print("ambiguous motif (i.e. 'N' base exists) detected: ")
    print(ambiguous_motif)
    
    
    result <- dplyr::filter(result, 
                            !stringr::str_detect(motif, "N"))
    
    print("Ambiguous motif(s) removed from result!")
  }
  
  # handle missing motif(s)
  
  if(nrow(missing_motif) != 0) {
    print("Missing motif detected: ")
    print(missing_motif)
    
    result <- dplyr::right_join(result, motif_ref, by = "motif") %>%
      tidyr::replace_na(replace = list(n = 0)) 
    
    print("Missing motif added back to the final result with count of 0!")
    
  }

  
  # calculate the fraction
  result_frac <- result %>%
    dplyr::mutate(fraction = n / sum(n))
  
  # set the factor level of motifs
  
  result_frac$motif <- factor(result_frac$motif,  levels = sort(motif_ref$motif))
  
  return(result_frac)
  
}










