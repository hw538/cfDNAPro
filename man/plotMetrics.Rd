% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotMetrics.R
\name{plotMetrics}
\alias{plotMetrics}
\title{Plot the fragment size metrics (i.e. proportion, cdf and 1-cdf)}
\usage{
plotMetrics(x, order, plot, vline, xlim, ylim, ...)
}
\arguments{
\item{x}{A long-format dataframe contains the metrics of different cohort.}

\item{order}{The groups show in the final plot,
the input value should be vector, e.g. `groups = c('group1','group2')``,
default is all folders in the folder path}

\item{plot}{The plot type,
default is 'all': both median and mean metrics will be shown.
They will include: mean_prop, mean_cdf, mean_1-cdf,
median_prop, median_cdf, median_1-cdf.
Could also specify as "median" or "mean".}

\item{vline}{Vertical dashed lines, default value is c(81,167).}

\item{xlim}{The x axis range shown in the plot. Default is c(0,500).}

\item{ylim}{The y axis range shown
in the fraction of fragment size plots. Default is c(0,0.0125).}

\item{...}{Further arguments passed to or from other methods.}
}
\value{
The function returns a list plots.
}
\description{
Plot the fragment size metrics (i.e. proportion, cdf and 1-cdf)
}
\examples{
# Get the path to example data.
path <- examplePath("groups_picard")
# Calculate the metrics.
df <- callMetrics(path = path)
# Plot metrics.
plot <- plotMetrics(df,
    plot = "median",
    order = c("cohort_1", "cohort_2")
)
}
\author{
Haichao Wang
}
