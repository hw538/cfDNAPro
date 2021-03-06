% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotPeakDistance.R
\name{plotPeakDistance}
\alias{plotPeakDistance}
\title{Plot the inter-peak distance of fragment size distance distribution}
\usage{
plotPeakDistance(x, summarized, order, type, mincount, xlim, ...)
}
\arguments{
\item{x}{A long-format dataframe contains the inter-peak distance,
a template please refer to the result of 'callPeakDistance' function.}

\item{summarized}{Logical value,
describe whether the x is summarzied already.
summarized means the count and proportion of each interpeak_dist.}

\item{order}{The groups show in the final plot,
  the input value should be vector,
  e.g. `groups = c('group1','group2')`,
default is all folders in the folder path.}

\item{type}{The plot type, default is line plot, now only support line plot.
Don't change this parameter in this version, keep it as default.}

\item{mincount}{Minimum count value of inter peak distance,
count number less than this value will be removed first,
then proportion of each count value will be calculated.
Default value is 0.}

\item{xlim}{The x axis range shown in the plot. Default is `c(8,13)`.}

\item{...}{Further arguments passed to or from other methods.}
}
\value{
The function returns the line plot of inter peak distance.
}
\description{
Plot the inter-peak distance of fragment size distance distribution
}
\examples{
# Get the path to example data.
path <- examplePath("groups_picard")

# Calculate the inter-peak distance.
df <- callPeakDistance(path = path)

# Plot the inter-peak distance.
plot <- plotPeakDistance(df,
    xlim = c(8, 13),
    mincount = 2
)
}
\author{
Haichao Wang
}
