
# cfDNAPro 1.5.4
* In addition to "bam" and "picard" files as the input, now we accept 
"cfdnapro" as input_type to various functions, this 'cfdnapro' input is exactly 
the output of `read_bam_insert_metrics` function in cfDNAPro package. It is a 
tsv file containing two columns, i.e., "insert_size" (fragment length) and 
"All_Reads.fr_count" (the count of the fragment length).

# cfDNAPro 1.5.3
* added support for hg38-NCBI version, i.e. GRCh38
# cfDNAPro 0.99.3 (19 July 2021)
* Modified vignette.
# cfDNAPro 0.99.2 (5 July 2021)
* Modified vignette.
# cfDNAPro 0.99.1 (27 May 2021)
* Added 'cfDNAPro' into the "watched tag".

# cfDNAPro 0.99.0 (26 May 2021)
* Now cfDNAPro supports bam file as input for data characterisation.
* Coding style improvements.
* Documentation improvements.
* Submitted to Bioconductor.


