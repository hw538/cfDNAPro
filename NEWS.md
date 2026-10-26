# cfDNAPro 1.7.4 (2026-10-26)
* Regenerated roxygen help pages and declared direct namespace dependencies.
* Added versioned GitHub Actions container builds with installed, offline tests.
* Added explicit custom BSgenome references with offline readers, saved-reference
  provenance, and safe motif/mutation reference recovery. Inspired by PR #11.
* Fixed mitochondrial alias resolution and `chromosome_to_keep = FALSE` in
  `readBam()`. Preserve reference lengths and circularity after curation, use
  local reference metadata, and retain length checks for incompatible references.
* Added hg19/hg38 regression tests for autosomal, sex-chromosome, mitochondrial,
  and mutation-aware filtering.

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
