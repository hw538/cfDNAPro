

plotMotif <- function(x){
  
  
  
  
  
  p <- ggplot(x, aes(.data$motif, .data$fraction)) +
    geom_col() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5))
  
  return(p)
  
}


