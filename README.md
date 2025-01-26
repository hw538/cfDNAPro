# cfDNAPro  <img src="vignettes/logo.png" width="150" align="right">
[![Anaconda-Server Badge](https://anaconda.org/bioconda/bioconductor-cfdnapro/badges/downloads.svg)](https://anaconda.org/bioconda/bioconductor-cfdnapro)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/bioconductor-cfdnapro/badges/latest_release_date.svg)](https://anaconda.org/bioconda/bioconductor-cfdnapro)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/bioconductor-cfdnapro/badges/license.svg)](https://anaconda.org/bioconda/bioconductor-cfdnapro)
[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2Fhw538%2FcfDNAPro&count_bg=%2379C83D&title_bg=%23555555&icon=github.svg&icon_color=%23E7E7E7&title=GitHub+view&edge_flat=true)](https://hits.seeyoufarm.com)

## Official tutorials

For detailed documentation, please visit: https://cfdnapro.readthedocs.io/en/latest/ 

About R version: As of 11 Jan 2025, some R package dependencies have not been able to support R version 4.4 yet. Please try R version 4.2 or 4.3 instead.

## Declaration  
cfDNAPro is designed for research only.

## Challenges in the cfDNA fragment length calculation
Unlike genomic DNA, cfDNA has specific fragmentation patterns. Ambiguous definition of "fragment length" by various alignment software is raising concerns: see page 9 footnote in SAM file format spec:  https://samtools.github.io/hts-specs/SAMv1.pdf   
Cell-free DNA data fragmentomic analysis requires single-molecule level resolution, which further emphasizes the importance of accurate/un-biased feature extraction.

`cfDNAPro` is designed to resolve this issue and standardize the cfDNA fragmentomic analysis using the bioconductor R ecosystem.

## Input
A paired-end sequencing bam file, with duplicates marked. (e.g., using the MarkDuplicates function from Picard).  
Please do not impose any filtering on the bam files; For example, do not filter by the proper-pairs flag.  
`cfDNAPro` filters the reads by following default criteria (You can toggle those criteria using parameters built in `readBam()` function):
(1) Reads mapping qualities less than 30 were discarded;  
(2) Reads must be paired. Of note, by default, cfDNAPro doesn’t impose filtration by “proper pair”;   
(3) No duplicate;  
(4) No secondary alignment;  
(5) No supplementary alignment;  
(6) No unmapped reads.  

Note: remember to choose the correct `genome_label`, a parameter in `readBam()` function, based on the ref genome you used for alignment.
At the moment, it supports three different ref genomes, hg19, hg38 and hg38-NCBI, For details see  readBam() R documentation by typing `?readBam` in the R console or 
see source code:https://github.com/hw538/cfDNAPro/blob/master/R/readBam.R

## Output
`cfDNAPro` can extract (i.e., "quantify in a standandised and robust way") these features/bio-markers:
- fragment length
- fragment start/end/upstream/downstream motifs
- copy number variation
- single nucleotide mutation
- more...

Feature extraction depends on essential data objects/R packages in the _Bioconductor_ ecosystem, such as `Rsamtools`, `plyranges`, `GenomicAlignments`, `GenomeInfoDb` and `Biostrings`.  
Data engineering depends on packges in the _tidyverse_ ecosystem, such as `dplyr`, and `stringr`.  
All plots depend on `ggplot2` R packge.  

For issues/feature request etc., please contact:   
__Nitzan Rosenfeld Lab admin mailbox__:  
Rosenfeld.LabAdmin@cruk.cam.ac.uk  
__Author__: Haichao Wang  
wanghaichao2014@gmail.com  
__Author__: Paulius D. Mennea  
paulius.mennea@cruk.cam.ac.uk   


## Quick Start 1
Read bam file, return the fragment name (i.e. read name in bam file) and alignment coordinates in GRanges object in R.
If needed, you can convert the GRanges into a dataframe and the fragment length is stored in the "width" column.

```R

library(cfDNAPro)

# read bam file, do alignment curation
 frags <- readBam(bamfile = "/path/to/bamfile.bam")
# convert GRanges object to a dataframe in R
 frag_df <- as.data.frame(frags)

# You can calculate fragment length and motifs from the frags object (i.e., the output of readBam() function)

frag_length <- callLength(frags)
frag_motif <- callMotif(frags)
```
A screenshot of the output:  

<img width="545" alt="image" src="https://github.com/hw538/cfDNAPro/assets/15274940/49f9cc93-d6af-4503-9b65-bbfea7b5ba87">

## Quick Start 2
Read in bam file, return the fragment length counts.
A straightforward and frequent user case: calculate the fragment size of a bam file, use the following code:

```R

# install cfDNAPro newest version 

if (!require(devtools)) install.packages("devtools")
devtools::install_github("hw538/cfDNAPro", build_vignettes = FALSE)

# calculate insert size of a bam file

library(cfDNAPro)
 frag_lengths <- read_bam_insert_metrics(bamfile = "/path/to/bamfile.bam")
```
The returned dataframe contains two columns, i.e., "insert_size" (fragment length) and "All_Reads.fr_count" (the count of the fragment length). A screenshot of the output:  
<img width="298" alt="image" src="https://github.com/hw538/cfDNAPro/assets/15274940/cba6709d-c49c-4c0d-8ae3-4ee7e82884f0">




## News
### cfDNAAPro 1.7.1 (Aug 2024)
* multiple updates

### cfDNAPro 1.7.1 (May 2023)
* Resolved issues when building vignette
* Various updates
* Added/Updated readBam() functions
### cfDNAPro 1.5.4 (Nov 2022)
* In addition to "bam" and "picard" files as the input, now we accept 
"cfdnapro" as input_type to various functions, this 'cfdnapro' input is exactly 
the output of `read_bam_insert_metrics` function in cfDNAPro package. It is a 
tsv file containing two columns, i.e., "insert_size" (fragment length) and 
"All_Reads.fr_count" (the count of the fragment length).
### cfDNAPro 1.5.3 (Oct 2022)
* added support for hg38-NCBI version, i.e. GRCh38
### cfDNAPro 0.99.3 (July 2021)
* Modified vignette.
### cfDNAPro 0.99.2 (July 2021)
* Modified vignette.
### cfDNAPro 0.99.1 (May 2021)
* Added 'cfDNAPro' into the "watched tag".
### cfDNAPro 0.99.0 (May 2021)
* Now cfDNAPro supports bam file as input for data characterisation.
* Coding style improvements.
* Documentation improvements.
* Submitted to Bioconductor.






## Installation

Please install our latest version(highly recommended):

```R
if (!require(devtools)) install.packages("devtools")
library(devtools)
devtools::install_github("hw538/cfDNAPro", build_vignettes = TRUE, dependencies = TRUE)
# run below instead if you don't want to build vignettes inside R
# devtools::install_github("hw538/cfDNAPro", build_vignettes = FALSE, dependencies = FALSE)

```


Or install the released/steady version (i.e., not newest version, some functions might be missing in comparison to functions shown in this webpage) 
via Bioconductor:
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("cfDNAPro")
```

## Vignettes/tutorials
visit: https://cfdnapro.readthedocs.io/en/latest/ 


## Citation

Please cite package ‘cfDNAPro’ in publications:

Haichao Wang, Paulius Mennea, Elkie Chan, Hui Zhao, Christopher G. Smith, Tomer Kaplan, Florian Markowetz, Nitzan Rosenfeld(2024). cfDNAPro:An R/Bioconductor package to extract and visualise cell-free DNA biological features. R package version 1.7.1 <https://github.com/hw538/cfDNAPro>
