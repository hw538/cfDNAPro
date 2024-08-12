# cfDNAPro  <img src="vignettes/logo.png" width="150" align="right">
[![Anaconda-Server Badge](https://anaconda.org/bioconda/bioconductor-cfdnapro/badges/downloads.svg)](https://anaconda.org/bioconda/bioconductor-cfdnapro)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/bioconductor-cfdnapro/badges/latest_release_date.svg)](https://anaconda.org/bioconda/bioconductor-cfdnapro)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/bioconductor-cfdnapro/badges/license.svg)](https://anaconda.org/bioconda/bioconductor-cfdnapro)
[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2Fhw538%2FcfDNAPro&count_bg=%2379C83D&title_bg=%23555555&icon=github.svg&icon_color=%23E7E7E7&title=GitHub+view&edge_flat=true)](https://hits.seeyoufarm.com)

## Official tutorials

For detailed documentation please visit: https://cfdnapro.readthedocs.io/en/latest/ 

## Declaration  
cfDNAPro is designed for research only.

## Abstract
**An R/Bioconductor package to extract and visualise cell-free DNA biological features in an open, standardized, robust and reproducible way.**

Cell-free DNA (cfDNA) enters human blood circulation by various biological processes, and includes tumour-derived circulating tumour DNA (ctDNA). There is increasing evidence that differences in biological features between cfDNA and ctDNA could be exploited to improve cancer detection, treatment selection and minimal residual disease detection. However, there are currently no R packages that support analysis of cfDNA biological features such as fragment length, nucleotide frequency, nucleosome occupancy etc.  

Here we present a Bioconductor R package, cfDNAPro, which provides an easy-to-use framework for automated data characterisation and visualisation of cfDNA sequencing data. The cfDNAPro R package implements functions for calculating overall, median and modal fragment size distributions, calculating the peaks and troughs, as well as the periodicity of oscillations in the fragment size profile and it includes functions for data visualisation. 

As the first R package for the analysis of cfDNA fragmentation profiles, we anticipate that cfDNAPro will improve the efficiency and reproducibility of cfDNA fragmentation analyses. We plan to regularly add support for other analyses and visualisations such as copy number variation, nucleosome position calling, GC content, fragment end motif analysis of fragments and others. cfDNAPro provides a foundation for follow-up analyses of fragmentation patterns by more advanced machine learning and data science methods. The package has been accepted by Bioconductor: https://bioconductor.org/packages/release/bioc/html/cfDNAPro.html 

_Challenges in the cfDNA fragment length calculation_  
Ambiguous definition of "fragment length" by various alignment software is raising concerns: see page 9 footnote in SAM file format spec:  https://samtools.github.io/hts-specs/SAMv1.pdf   
Cell-free DNA data fragmentomic analysis requires single-molecule level resolution, which further emphasizes the importance of accurate/un-biased feature extraction.

_cfDNAPro is designed to resolve this issue and standardize the cfDNA fragmentomic analysis._

## Input
cfDNAPro is specifically written for cell-free DNA paire-ed whole-genome sequencing data. 
Its ensures accurate (i.e. up-to-standard) calculation of fragmentomic features (e.g., fragment lengths and motif)

## Output
cfDNApro supports extration and visualisation of features/markers:
- fragment length
- fragment start/end/upstream/downstream motifs
- copy number variation
- single nucleotide mutation
- more...

All plots depend on ggplot2 R packge.

For issues/feature requests/comments, please raise an issue or email:   
__Author__: Haichao Wang  
wanghaichao2014@gmail.com  
__Author__:Paulius D. Mennea  
paulius.mennea@cruk.cam.ac.uk 
__Rosenfeld Lab admin mailbox__:  
Rosenfeld.LabAdmin@cruk.cam.ac.uk  

## Quick Start 1
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


## Quick Start 2
Read bam file, return the fragment name (i.e. read name in bam file) and alignment coordinates in GRanges object in R.
If needed, you can convert the GRanges into a dataframe and the fragment length is stored in the "width" column.

```R

library(cfDNAPro)

# read bam file, do alignment curation
 frags <- readBam(bamfile = "/path/to/bamfile.bam")
# convert GRanges object to a dataframe in R
 frag_df <- as.data.frame(frags)

```
A screenshot of the output:  

<img width="545" alt="image" src="https://github.com/hw538/cfDNAPro/assets/15274940/49f9cc93-d6af-4503-9b65-bbfea7b5ba87">



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

Haichao Wang, Elkie Chan, Hui Zhao, Christopher G. Smith, Tomer Kaplan, Florian Markowetz, Nitzan Rosenfeld(2020). cfDNAPro:An R/Bioconductor package to extract and visualise cell-free DNA biological features. R package version 1.7 <https://github.com/hw538/cfDNAPro>
