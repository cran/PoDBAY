% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FunctionsSupplementary.R
\name{numToBool}
\alias{numToBool}
\title{Numeric to boolean}
\usage{
numToBool(x)
}
\arguments{
\item{x}{numeric value (0, 1)}
}
\value{
boolean value (T, F)
}
\description{
Converts numeric format to boolean format.
}
\details{
If the function is supposed to be used on a vector, the form \code{sapply("vector", numToBool)} needs to be applied.
}
\examples{
dStatus <- c(0,0,1,1,0,1)
sapply(dStatus, numToBool)

}
