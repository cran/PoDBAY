#' Dataset containing the information for vaccinated subjects
#'
#' A dataset containing the N, mean, stdDev, titers of vaccinated subjects. The dataset is provided in the form of population class object (see the \code{Population-class} function for more details).
#'
#' @format Population class object:
#' \describe{
#'   \item{N}{number of subjects}
#'   \item{mean}{mean of titers}
#'   \item{stdDev}{standard deviation of titers}
#'   \item{titers}{subject level titers}
#' }
"vaccinated"

#' Dataset containing the information for control subjects
#'
#' A dataset containing the N, mean, stdDev, titers of control subjects. The dataset is provided in the form of population class object (see the \code{Population-class} function for more details).
#'
#' @format Population class object:
#' \describe{
#'   \item{N}{number of subjects}
#'   \item{mean}{mean of titers}
#'   \item{stdDev}{standard deviation of titers}
#'   \item{titers}{subject level titers}
#' }
"control"

#' Dataset containing the information for diseased subjects
#'
#' A dataset containing the N, mean, stdDev, titers of diseased subjects. The dataset is provided in the form of population class object (see the \code{Population-class} function for more details).
#'
#' @format Population class object:
#' \describe{
#'   \item{N}{number of subjects}
#'   \item{mean}{mean of titers}
#'   \item{stdDev}{standard deviation of titers}
#'   \item{titers}{subject level titers}
#' }
"diseased"

#' Dataset containing the information for non-diseased subjects
#'
#' A dataset containing the N, mean, stdDev, titers of non-diseased subjects. The dataset is provided in the form of population class object (see the \code{Population-class} function for more details).
#'
#' @format Population class object:
#' \describe{
#'   \item{N}{number of subjects}
#'   \item{mean}{mean of titers}
#'   \item{stdDev}{standard deviation of titers}
#'   \item{titers}{subject level titers}
#' }
"nondiseased"

#' PoD curve parameters
#'
#' A dataset containing PoD curve parameters
#'
#' @format data frame
#' \describe{
#'   \item{pmax}{pmax: maximum PoD}
#'   \item{et50}{et50: titer value corresponding to the pmax/2}
#'   \item{slope}{slope: slope of the PoD curve}
#' }
"PoDParams"

#' Estimated PoD curve parameters
#'
#' A dataset containing estimated set of PoD curve parameters. (Set of PoD curve parameters is a vector obtained by number of replications specified by repeatCount. These replications are performed for calculation of a confidence interval. For more details, see the supplementary material of the article).
#'
#' @format data frame
#' \describe{
#'   \item{pmax}{pmax: maximum PoD}
#'   \item{et50}{et50: titer value corresponding to the pmax/2}
#'   \item{slope}{slope: slope of the PoD curve}
#' }
"estimatedParameters"

#' Estimated PoDBAY efficacies
#'
#' A dataset containing estimated set of PoDBAY efficacies. (Set of efficacies is a vector obtained by number of replications specified by repeatCount. These replications are performed for calculation of a confidence interval. For more details, see the supplementary material of the article).
#'
#' @format vector
#' \describe{
#'   \item{numeric vector}{PoDBAY efficacies}
#' }
"efficacySet"
