% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotValleyDistance.R
\name{plotValleyDistance}
\alias{plotValleyDistance}
\title{Plot the inter-valley distance of fragment size distance distribution}
\usage{
plotValleyDistance(x, order, type, mincount, xlim, ...)
}
\arguments{
\item{x}{A long-format dataframe contains the inter-valley distance,
a template please refer to the result of 'callValleyDistance' function.}

\item{order}{The groups show in the final plot,
the input value should be vector,
e.g. `groups=c('group1','group2')`,
default is all folders in the folder path.}

\item{type}{The plot type, default is line plot, now only support line plot.
Don't change this parameter in this version, keep it as default.}

\item{mincount}{Minimum count value of inter valley distance,
count number less than this value will be removed first,
then proportion of each count value will be calculated. Default value is 0.}

\item{xlim}{The x axis range shown in the plot. Default is c(8,13).}

\item{...}{Further arguments passed to or from other methods.}
}
\value{
The function returns the line plot of inter valley distance.
}
\description{
Plot the inter-valley distance of fragment size distance distribution
}
\examples{
# Get the path to example data.
path <- examplePath("groups_picard")

# Calculate the inter-valley distance.
df <- callValleyDistance(path = path)

# Plot the inter-valley distance.
plot <- plotValleyDistance(df,
    xlim = c(8, 13),
    mincount = 2
)
}
\author{
Haichao Wang
}
