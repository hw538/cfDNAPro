
theme_length_plot <- function(){

  theme_classic() %+replace%
    theme(
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
      axis.title = element_text(size = 14, face = "bold"),
      legend.position = "none"
    )
}