# cfDNAPro  <img src="vignettes/logo.png" width="150" align="right">
[![Anaconda-Server Badge](https://anaconda.org/bioconda/bioconductor-cfdnapro/badges/downloads.svg)](https://anaconda.org/bioconda/bioconductor-cfdnapro)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/bioconductor-cfdnapro/badges/latest_release_date.svg)](https://anaconda.org/bioconda/bioconductor-cfdnapro)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/bioconductor-cfdnapro/badges/license.svg)](https://anaconda.org/bioconda/bioconductor-cfdnapro)
[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2Fhw538%2FcfDNAPro&count_bg=%2379C83D&title_bg=%23555555&icon=github.svg&icon_color=%23E7E7E7&title=GitHub+view&edge_flat=true)](https://hits.seeyoufarm.com)


**An R/Bioconductor package to extract and visualise cell-free DNA biological features in an open, standardized, robust and reproducible way.**

Cell-free DNA (cfDNA) enters human blood circulation by various biological processes, and includes tumour-derived circulating tumour DNA (ctDNA). There is increasing evidence that differences in biological features between cfDNA and ctDNA could be exploited to improve cancer detection, treatment selection and minimal residual disease detection. However, there are currently no R packages that support analysis of cfDNA biological features such as fragment length, nucleotide frequency, nucleosome occupancy etc. Here we present a Bioconductor R package, cfDNAPro, which provides an easy-to-use framework for automated data characterisation and visualisation of cfDNA sequencing data. The cfDNAPro R package implements functions for calculating overall, median and modal fragment size distributions, calculating the peaks and troughs, as well as the periodicity of oscillations in the fragment size profile and it includes functions for data visualisation. 

As the first R package for the analysis of cfDNA fragmentation profiles, we anticipate that cfDNAPro will improve the efficiency and reproducibility of cfDNA fragmentation analyses. We plan to regularly add support for other analyses and visualisations such as copy number variation, nucleosome position calling, GC content, fragment end motif analysis of fragments and others. cfDNAPro provides a foundation for follow-up analyses of fragmentation patterns by more advanced machine learning and data science methods. The package has been accepted by Bioconductor: https://bioconductor.org/packages/release/bioc/html/cfDNAPro.html 

## Highlights

cfDNAPro is under active development, its internal quality control steps ensures accurate calculation of fragment lengths. 
More feature extraction utilities will be added. For issues/feature requests/comments, please raise an issue or email me: haichao.wang@cruk.cam.ac.uk


## News

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



## Quick Start

A straightforward and frequent user case: calculate the fragment size of a bam file, use the following code:

```R

# install cfDNAPro newest version 

if (!require(devtools)) install.packages("devtools")
devtools::install_github("hw538/cfDNAPro", build_vignettes = TRUE)

# calculate insert size of a bam file

library(cfDNAPro)
 dataframe <- read_bam_insert_metrics(bamfile = "/path/to/bamfile.bam")
```
The returned dataframe contains two columns, i.e., "insert_size" (fragment length) and "All_Reads.fr_count" (the count of the fragment length).


## Installation

Please install our latest version(recommended):

```R
if (!require(devtools)) install.packages("devtools")
library(devtools)
devtools::install_github("hw538/cfDNAPro", build_vignettes = TRUE)
```

Or install the released/steady version (i.e., not newest version) 
via Bioconductor:
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("cfDNAPro")
```

## Vignettes

See Bioconductor official documentation:  
https://bioconductor.org/packages/release/bioc/vignettes/cfDNAPro/inst/doc/cfDNAPro.html

## Citation

Please cite package ‘cfDNAPro’ in publications:

Haichao Wang, Elkie Chan, Hui Zhao, Christopher G. Smith, Tomer Kaplan, Florian Markowetz, Nitzan Rosenfeld(2020). cfDNAPro:An R/Bioconductor package to extract and visualise cell-free DNA biological features. R package version 1.5 <https://github.com/hw538/cfDNAPro>
