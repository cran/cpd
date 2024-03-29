% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitctp.R
\name{plot.fitCBP}
\alias{plot.fitCBP}
\title{Plot of observed and theoretical frequencies for a CBP fit}
\usage{
\method{plot}{fitCBP}(x, plty = "FREQ", maxValue = NULL, ...)
}
\arguments{
\item{x}{An object of class \code{'fitCBP'}}

\item{plty}{Plot type to be shown. Default is \code{"FREQ"} which shows the observed and theoretical frequencies for each value of the variable; \code{"CDF"} and \code{"PP"} are also available for plotting the empirical and theoretical cumulative distribution functions or the theoretical cumulative probabilities against the empirical cumulative probabilities, respectively.}

\item{maxValue}{maxValue you want to appear in the plot}

\item{...}{Additional parameters.}
}
\description{
Plot of observed and theoretical frequencies for a CBP fit
}
\examples{
set.seed(123)
x <- rcbp(500, 1.75, 3.5)
fit <- fitcbp(x)
plot(fit)
plot(fit, plty = "CDF")
plot(fit, plty = "PP")
 
}
