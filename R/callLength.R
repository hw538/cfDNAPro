
#' call fragment length
#' @import ggplot2
#' @import dplyr
#' @import plyranges
#'
#' @param fragment_obj 
#' @param isize_min 
#' @param isize_max 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
callLength <- function(fragment_obj, 
                       isize_min = 1L,
                       isize_max = 1000L,
                       ...){
  
  frag <- fragment_obj
  # calculating insert sizes
  message("Calculating insert sizes...")
  frag$insert_size <- BiocGenerics::width(frag)
  
  # size analysis
  frag <- plyranges::filter(frag, 
                            insert_size >= isize_min & insert_size <= isize_max)
  isize <- frag$insert_size 
  
  isize_tibble <- tibble("insert_size" = isize, "count" = 1 ) %>%
    dplyr::filter(!is.na(insert_size))
  
  result <- isize_tibble %>%
    dplyr::group_by(.data$insert_size) %>%
    dplyr::summarise("All_Reads.fr_count" = sum(count))
  
  
  # quality control results
  # Create a vector of elements
  isize_ref <- seq.int(isize_min, isize_max, by = 1L) %>% 
    as_tibble()
  
  colnames(isize_ref) <- c("insert_size")
  
  # report abnormal isizes
  
  missing_isize <- dplyr::anti_join(isize_ref, result, by = "insert_size")
  
  # handle missing isize(s)
  
  if(nrow(missing_isize) != 0) {
    message("Missing isize detected: ")
    print(dplyr::pull(missing_isize, insert_size))
    
    result <- dplyr::right_join(result, isize_ref, by = "insert_size") %>%
      tidyr::replace_na(replace = list(All_Reads.fr_count = 0)) %>% 
      dplyr::arrange(insert_size) %>%
      dplyr::mutate(prop = All_Reads.fr_count / sum(All_Reads.fr_count))
    
    message("Missing isize(s) added back to the final result with count of 0!")
    
    message("Job completed successfully. ")
    
  }
  
  
  return(result)
  
  
  
}
