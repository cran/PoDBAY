#' @title Clinical trial: estimation of case-count efficacy
#'
#' @description Function assigns disease status (DS) to vaccinated and control groups and based on that calculates the case-count efficacy. Vaccinated and control groups are provided in the form of population class objects (see the \code{Population-class} function for more details).
#'
#' Input populations need to contain information about Probability of disease (PoD) for each subject - calculated using \code{population$assignPoD(PoD(x))}. See \code{PoD} function for further details.
#'
#' @param vaccinated \code{Population-class} object: vaccinated subjects with assigned PoD
#' @param control \code{Population-class} object: control subjects with assigned PoD
#' @param CI numeric: value from (0, 1) interval, confidence level of interest
#'
#' @return
#' \itemize{
#'   \item vaccinated: vaccinated subjects with assigned DS, \code{Population-class} object
#'
#'   \item control: control subjects with assigned DS, \code{Population-class} object
#'
#'   \item efficacy: case-count efficacy
#'
#'   \item confidenceInterval: case-count efficacy confidence interval calculated with \code{waldCI()} function
#' }
#'
#' @usage
#' ClinicalTrial(vaccinated, control, CI = 0.95)
#'
#' @examples
#' # Loading vaccinated, control population data with PoD information
#' data(vaccinated)
#' data(control)
#'
#' # Estimating the disease status and case-count efficacy with 95\% confidence interval
#' CT <- ClinicalTrial(vaccinated, control, CI = 0.95)
#'
#' CT$efficacy
#' CT$confidenceInterval
#'
#' CT$vaccinated
#'
#' @export
ClinicalTrial <- function(vaccinated,
                          control,
                          CI = 0.95){

  if (!is(vaccinated, "Population")) {
    incorrectPopulationInput("vaccinated")
  }
  if (!is(control, "Population")) {
    incorrectPopulationInput("control")
  }

  dsVacc <- mapply(rbinom, n = 1, size = 1, prob = vaccinated$PoDs)
  vaccinated$diseaseStatus <- sapply(dsVacc, numToBool)
  dsCont <- mapply(rbinom, n = 1, size = 1, prob = control$PoDs)
  control$diseaseStatus <- sapply(dsCont, numToBool)

  casecountEfficacy <- 1 - (vaccinated$getDiseasedCount() / vaccinated$N) /
    (control$getDiseasedCount() / control$N)

  confidenceInterval <- waldCI(vaccinated, control, CI)

  return(list(
    vaccinated = vaccinated,
    control = control,
    efficacy = casecountEfficacy,
    confidenceInterval = confidenceInterval
  ))
}

#' @title Clinical trial function expanded for usage in simulations when the calculation of coverage probability is needed for three confidence intervals: 80\%, 90\%, and user-defined  
#' 
#' 
#' @description Function works the same way as \code{ClinicalTrial} function but it also calculates 80\% and 90\% confidence intervals.
#'
#' @param vaccinated \code{Population-class} object: vaccinated subjects with assigned PoD
#' @param control \code{Population-class} object: control subjects with assigned PoD
#' @param CI numeric: value from (0, 1) interval, confidence level of interest 
#'
#' @return
#' \itemize{
#'   \item vaccinated: vaccinated subjects with assigned DS, \code{Population-class} object
#'
#'   \item control: control subjects with assigned DS, \code{Population-class} object
#'
#'   \item efficacy: case-count efficacy
#'
#'   \item confidenceInterval: confidence interval calculated with \code{waldCI} function
#'
#'   \item confidenceInterval90: 90\% confidence interval calculated with \code{waldCI} function
#'
#'   \item confidenceInterval80: 80\% confidence interval calculated with \code{waldCI} function
#' }
#'
#' @usage
#' ClinicalTrialCoverage(vaccinated, control, CI = 0.95)
#'
#' @export
ClinicalTrialCoverage <- function(vaccinated,
                                  control,
                                  CI = 0.95){

  if (!is(vaccinated, "Population")) {
    incorrectPopulationInput("vaccinated")
  }
  if (!is(control, "Population")) {
    incorrectPopulationInput("control")
  }

  dsVacc <- mapply(rbinom, n = 1, size = 1, prob = vaccinated$PoDs)
  vaccinated$diseaseStatus <- sapply(dsVacc, numToBool)
  dsCont <- mapply(rbinom, n = 1, size = 1, prob = control$PoDs)
  control$diseaseStatus <- sapply(dsCont, numToBool)

  casecountEfficacy <- 1 - (vaccinated$getDiseasedCount() / vaccinated$N) /
    (control$getDiseasedCount() / control$N)

  confidenceInterval95 <- waldCI(vaccinated, control, CI)
  confidenceInterval90 <- waldCI(vaccinated, control, 0.90)
  confidenceInterval80 <- waldCI(vaccinated, control, 0.80)

  return(list(
    vaccinated = vaccinated,
    control = control,
    efficacy = casecountEfficacy,
    confidenceInterval95 = confidenceInterval95,
    confidenceInterval90 = confidenceInterval90,
    confidenceInterval80 = confidenceInterval80
  ))
}

#' @title Wald confidence interval estimation
#'
#' @description Function calculates and returns case-count efficacy confidence intervals estimated using Wald's method.
#'
#' Input data need to contain information about disease status on individual level.
#'
#' @param vaccinated \code{Population-class} object: vaccinated subjects, containing information about disease status
#' @param control \code{Population-class} object: control subjects, containing information about disease status
#' @param confLevel numeric: value from (0, 1) interval, confidence level of interest
#'
#' @return Named list of lower and upper confidence interval bound
#'
#' @details  Confidence interval of the relative risk is calculated using the Wald method. (Wald, A. Tests of statistical hypotheses concerning several parameters when the number of observations is large. Transactions of the American Mathematical Society 54, 426-482 (1943)).
#'
#' @examples
#' # Loading vaccinated and control populations data with PoD information
#' data(vaccinated)
#' data(control)
#'
#' # Estimating the disease status and case-count efficacy with 95\% confidence interval
#' set.seed(1)
#' CT <- ClinicalTrial(vaccinated, control, CI = 0.95)
#'
#' waldCI(vaccinated, control)
#'
#' @export
waldCI <- function(vaccinated, 
                   control, 
                   confLevel = 0.95){
  
  if (!is(vaccinated, "Population")) {
    incorrectPopulationInput("vaccinated")
  }
  if (!is(control, "Population")) {
    incorrectPopulationInput("control")
  }
  
  # significance level z value
  N. <- 1 - ((1 - confLevel)/2)
  z <- stats::qnorm(N., mean = 0, sd = 1)
  
  # Get matrix of vaccine/control + diseased/non-diseased 
  A <- vaccinated$getDiseasedCount()
  B <- vaccinated$getNondiseasedCount()
  C <- control$getDiseasedCount()
  D <- control$getNondiseasedCount()
  
  # Get size of vaccinated and control
  N_VACC <- A + B # N VACC
  N_CONTROL <- C + D # N CONT
  
  # calculate RR point estimate
  RR_point_t <- (A/N_VACC) / (C/N_CONTROL)
  ln_RR_point <- log(RR_point_t)
  
  # Wald RR se  
  ln_RR_point_se <- sqrt((1/A) - (1/N_VACC) + (1/C) - (1/N_CONTROL))
  RR_point_se <- exp(ln_RR_point_se)
  
  # Wald RR confidence intervals
  RR_lower <- 1 - exp(ln_RR_point + (z * ln_RR_point_se))
  RR_upper <- 1 - exp(ln_RR_point - (z * ln_RR_point_se))
  RR_point <- 1 - RR_point_t
  
  # summary
  wald_ci_point <- RR_point
  wald_ci_lower <- RR_lower
  wald_ci_upper <- RR_upper
  
  return(list(
    lowerBound = wald_ci_lower,
    upperBound = wald_ci_upper
  )
  )
}
