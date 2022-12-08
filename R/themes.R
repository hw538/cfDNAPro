
theme_length_plot <- function(){

  theme_classic() %+replace%
    theme(
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
      axis.title = element_text(size = 14, face = "bold"),
      legend.position = "none"
    )
}


theme_cnv_plot <- function(){
  theme_classic() %+replace%
    theme(axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title = element_text(face = "bold"))
}
