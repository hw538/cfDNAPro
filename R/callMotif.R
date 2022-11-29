
callMotif <- function(bamfile, 
                      motif_type,
                      motif_length,
                      genome_label  = "hg19", 
                      chromosome_to_keep = paste0("chr", 1:22),
                      strand_mode = 1,
                      outdir = NA,
                      ...) {
  
  
  fragments <- cfDNAPro::readBam(bamfile = bamfile, 
                    genome_label = genome_label,
                    chromosome_to_keep = chromosome_to_keep,
                    strand_mode = strand_mode,
                    outdir = outdir)
  
  
  
  
  
  
}




plotMotif <- function(){
  
  
  
  
  
}