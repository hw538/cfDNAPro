% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/callLength.R
\name{callLength}
\alias{callLength}
\title{Call Fragment Length}
\usage{
callLength(
  frag_obj,
  isize_min = 1L,
  isize_max = 1000L,
  integrate_mut = FALSE,
  ref_type = "locus_fragment",
  downsample_ref = FALSE,
  ...
)
}
\arguments{
\item{frag_obj}{GRanges object containing DNA fragment information.}

\item{isize_min}{Minimum insert size for filtering (default 1).}

\item{isize_max}{Maximum insert size for filtering (default 1000).}

\item{integrate_mut}{Logical, TRUE integrates mutational data.}

\item{ref_type}{Specifies the type of reference fragments for comparison
with MUT fragments: 'locus_fragment' refers to fragments/read-pairs
overlapping specific mutation loci, 'outer_fragment' refers to
fragments/read-pairs that do not overlap indicated mutation loci.}

\item{downsample_ref}{Logical, TRUE to match REF and MUT fragment counts.}

\item{...}{Additional arguments for downstream functions.}
}
\value{
Returns a tibble with summarized fragment data, optionally
including mutational analysis.
}
\description{
Computes distribution of fragment sizes from BAM/GRanges data. Optionally
integrates mutational data, normalizes fragment counts, and adjusts for
specific reference types.
}
\examples{
callLength(frag_obj = myFragData, isize_min = 1L, isize_max = 1000L,
           integrate_mut = TRUE)
}
