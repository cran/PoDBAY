% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RefClassPopulation.R
\name{getDiseasedTiters}
\alias{getDiseasedTiters}
\title{Diseased titers}
\value{
numeric vector: titers of diseased subjects in the \code{Population-class} object
}
\description{
Function returns titers of diseased subjects (disease status = TRUE) in the \code{Population-class} object.
}
\details{
Input into the function, "diseaseStatus", is taken from the \code{Population-class} object attribute. Information about disease status is written into the \code{Population-class} object by the \code{ClinicalTrial()} function.
}
