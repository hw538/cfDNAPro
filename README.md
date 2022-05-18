# cfDNAPro  <img src="vignettes/logo.png" width="150" align="right">

 An R/Bioconductor package to extract and visualise cell-free DNA biological features.

> Cell-free DNA (cfDNA) enters human blood circulation by various biological processes, and includes tumour-derived circulating tumour DNA (ctDNA). There is increasing evidence that differences in biological features between cfDNA and ctDNA could be exploited to improve cancer detection, treatment selection and minimal residual disease detection. However, there are currently no R packages that support analysis of cfDNA biological features such as fragment length, nucleotide frequency, nucleosome occupancy etc. Here we present a Bioconductor R package, cfDNAPro, which provides an easy-to-use framework for automated data characterisation and visualisation of cfDNA sequencing data. The cfDNAPro R package implements functions for calculating overall, median and modal fragment size distributions, calculating the peaks and troughs, as well as the periodicity of oscillations in the fragment size profile and it includes functions for data visualisation. 

> As the first R package for the analysis of cfDNA fragmentation profiles, we anticipate that cfDNAPro will improve the efficiency and reproducibility of cfDNA fragmentation analyses. We plan to regularly add support for other analyses and visualisations such as copy number variation, nucleosome position calling, GC content, fragment end motif analysis of fragments and others. cfDNAPro provides a foundation for follow-up analyses of fragmentation patterns by more advanced machine learning and data science methods. The package has been accepted by Bioconductor: https://bioconductor.org/packages/release/bioc/html/cfDNAPro.html 

## Highlights

cfDNAPro is under active development, its internal quality control steps ensures accurate calculation of fragment lengths. 
More feature extraction utilities will be added.

## Quick Start

A straight forward user case: calculate the fragment size of a bam file, use the following code:

```R
library(cfDNAPro)
 dataframe <- read_bam_insert_metrics(bamfile = "/path/to/bamfile.bam")
```
The dataframe contains two columns, i.e., "insert_size" (fragment length) and "All_Reads.fr_count" (the count of the fragment length).

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

See bioconductor documentation for cfDNAPro here: https://bioconductor.org/packages/release/bioc/vignettes/cfDNAPro/inst/doc/cfDNAPro.html

Or to see the vignettes in Rstudio (you have to indicate `build_vignettes = TRUE` during aforementioned installation step), use the command:

```R
browseVignettes("cfDNAPro")
```
## Citation

Please cite package ‘cfDNAPro’ in publications:

Haichao Wang, Hui Zhao, Elkie Chan, Christopher G. Smith, Tomer Kaplan, Florian Markowetz, Nitzan Rosenfeld(2020). cfDNAPro:An R/Bioconductor package to extract and visualise cell-free DNA biological features. R package version 1.3 <https://github.com/hw538/cfDNAPro>
