
library(cfDNAPro)
library(tidyverse)
bam <- "/Users/wang04/Documents/phd_projects/cfdnapro_debug/SLX-11379_D703_D504.mrkdup.BL.filtered.bam.markduplicates.bam.downsamp_0.1x.mrkdup.bam"

frag <-  readBam(bamfile = bam)

# callSize

isize1 <- read_bam_insert_metrics(fragment_obj = frag)

isize2 <- callLength(fragment_obj = frag)

plotLength(isize2, plot_type = "fraction", xlim = c(20, 290), line_color = "grey4", vline_size =  0.51 )
plotLength(isize2, plot_type = "count", xlim = c(1, 500), line_color = "grey4", add_vline = FALSE)

# Motif

motifs1 <- callMotif(frag, motif_type = "s", motif_length = 1L)
motifs1_plot <- plotMotif(test2)

# CNV
cnv1 <- callCNV(bamfile = bam) 
cnv1_plot <- plotCNV(cnv1)


###############################################################################
#QCit
###############################################################################

cfqc < -function(
    bamfile, 
    
    plot_file) {
  
  frag <- readBam(bamfile = bamfile)
  
  p_isize <- callLength(frag) %>% plotLength()
    
  p_motif <- callMotif(frag) %>% plotMotif()
    
  p_cnv <- callCNV(bamfile = bamfile) %>% plotCNV()
  
  cnv_obj <- callCNV(bamfile = bam)
  
  
}



