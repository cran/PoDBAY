#' @title PoDBAY efficacy estimation
#'
#' @description
#' Function calculates the PoDBAY efficacy based on the set of PoD curve parameters calculated in \code{PoDParamEstimation} function, vaccinated and control immunogenicity subset means and standard deviations.
#'
#' @param estimatedParameters named data frame ("pmax", "slope", "et50"): set of estimated PoD curve parameters
#' @param blindVaccinated \code{Population-class} object: vaccinated subjects from immunogenicity subset, containing N, mean, standard deviation information
#' @param blindControl \code{Population-class} object: control subjects from immunogenicity subset, containing N, mean, standard deviation information
#' @param adjustTiters boolean: set to TRUE if titer values should be adjusted, for details see \code{PoD} function
#' @param adjustFrom numeric: value specifying the detection limit, all values below the detection limit will be adjusted to adjustTo value
#' @param adjustTo numeric: value to which titers below the detection limit will be adjusted
#'
#' @return
#' efficacySet, set of PoDBAY effficacies corresponding to estimated set of PoD curve parameters
#'
#' @usage
#' PoDBAYEfficacy(estimatedParameters,
#'               blindVaccinated,
#'               blindControl,
#'               adjustTiters = FALSE,
#'               adjustFrom = log2(10),
#'               adjustTo = log2(5))
#'
#' @examples
#' ## Data preparation
#' data(diseased)
#' data(nondiseased)
#' data(estimatedParameters)
#'
#' ## Example 1
#' # Creating imunogenicity subset, method = "Ratio", value = 4
#' ImmunogenicitySubset <- 
#'   BlindSampling(diseased, 
#'                 nondiseased, 
#'                 method = list(name = "Ratio", 
#'                               value = 4))
#'                               
#' # Estimating PoD curve parameters
#' nondiseasedGenerationCount <- nondiseased$N
#'
#' estimatedParameters <- PoDParamEstimation(diseased$titers,
#'                        ImmunogenicitySubset$ImmunogenicityNondiseased$titers,
#'                        nondiseasedGenerationCount,
#'                        repeatCount = 10)
#'                        
#' # Estimating PoDBAY efficacy  
#' PoDBAYEfficacy(estimatedParameters$results,
#'               ImmunogenicitySubset$ImmunogenicityVaccinated,
#'               ImmunogenicitySubset$ImmunogenicityControl)
#'
#' @details
#' Application of \code{efficacyComputation} function to the all PoD curves (each characterized by three PoD parameters) estimated by \code{PoDParamEstimation} function.
#'
#' Inputs into the \code{efficacyComputation} are:
#' \itemize{
#'   \item PoDParameters: i'th estimated PoD parameters from \code{PoDParamEstimation}. i = 1, ..., N, where N = number of estimations in which MLE converges. See \code{PoDMLE} for details.
#'
#'   \item means: jittered means of immunogenicity subset. See \code{JitterMeans} for details.
#'
#'   \item standardDeviations: standard deviations of the vaccinated and control subjects from the immunogenicity subset.
#' }
#'
#' @export
PoDBAYEfficacy <- function(estimatedParameters,
                           blindVaccinated,
                           blindControl,
                           adjustTiters = FALSE,
                           adjustFrom = log2(10),
                           adjustTo = log2(5)) {
  
  if (!is(blindVaccinated, "Population")) {
    incorrectPopulationInput("blindVaccinated")
  }
  if (!is(blindControl, "Population")) {
    incorrectPopulationInput("blindControl")
  }
  if(blindVaccinated$mean < blindControl$mean) {
    warning(paste("Mean of vaccinated population is lower than mean of control population"))
  }
  
  efficacySet <- numeric()
  
  for (i in 1:nrow(estimatedParameters)) {
    jitterVaccinated <- JitterMean(blindVaccinated)
    jitterControl <- JitterMean(blindControl)
    means <- list(
      vaccinated = jitterVaccinated,
      control = jitterControl
    )
    stdDevs <- list(
      vaccinated = blindVaccinated$stdDev,
      control = blindControl$stdDev
    )
    
    efficacy <- efficacyComputation(
      estimatedParameters[i, ],
      means,
      stdDevs,
      adjustTiters = adjustTiters,
      adjustFrom = adjustFrom,
      adjustTo = adjustTo
    )
    efficacySet <- c(efficacySet, efficacy)
  }
  
  return(efficacySet)
}

#' @title PoDBAY efficacy summary: mean, median, confidence intervals
#'
#' @description
#' Function summarizes PoDBAY efficacy statistics (mean, median, confidence intervals) based on the set of estimated efficacies and chosen condfidence level. (Set of efficacies is a vector obtained by number of replications specified by repeatCount. These replications are performed for calculation of a confidence interval. For more details, see the supplementary material of the article).
#'
#' @param efficacySet numeric vector: estimated PoDBAY efficacies from \code{PoDBAYEfficacy} function.
#' @param ci numeric: required confidence level
#'
#' @return
#' named list: mean, median, CILow, CIHigh
#'
#' @usage
#' EfficacyCI(efficacySet, ci = 0.95)
#'
#' @examples
#' ## Data preparation
#' data(efficacySet)
#'
#' ## Example 1
#' EfficacyCI(efficacySet, ci = 0.95)
#'
#' @details
#' Confidence intervals are calculated using quantiles of estimated efficacies.
#'
#' @export
EfficacyCI <- function(efficacySet, ci = 0.95) {
  CILow <- quantile(efficacySet, (1 - ci) / 2, names = F)
  CIHigh <- quantile(efficacySet, ci + (1 - ci) / 2, names = F)
  return(
    list(
      mean = mean(efficacySet),
      median = median(efficacySet),
      CILow = CILow,
      CIHigh = CIHigh
    )
  )
}

#' @title PoDBAY efficacy summary at three confidence levels
#'
#' @description
#' Function summarizes PoDBAY efficacy statistics (mean, median, confidence intervals) at 80\%, 90\% and user-defined confidence levels, based on the set of estimated efficacies. (Set of efficacies is a vector obtained by number of replications specified by repeatCount. These replications are performed for calculation of a confidence interval. For more details, see the supplementary material of the article).
#'
#' @param efficacySet numeric vector: estimated PoDBAY efficacies from \code{PoDBAYEfficacy} function.
#' @param ci numeric: value from (0, 1) interval, confidence level of interest
#'
#' @return named list: mean, median, CILow, CIHigh
#'
#' @usage
#' EfficacyCICoverage(efficacySet, ci = 0.95)
#'
#' @examples
#' ## Data preparation
#' data(efficacySet)
#'
#' ## Example 1
#' EfficacyCICoverage(efficacySet, ci = 0.95)
#'
#' @details
#' Confidence intervals are calculated using quantiles of estimated efficacies.
#'
#' @export
EfficacyCICoverage <- function(efficacySet, ci = 0.95) {
  efficacySet <- efficacySet[!is.na(efficacySet)]
  sortedEfficacy <- sort(efficacySet)
  CILow95 <- quantile(sortedEfficacy, (1 - ci) / 2, names = F)
  CIHigh95 <- quantile(sortedEfficacy, ci + (1 - ci) / 2, names = F)
  
  CILow90 <- quantile(sortedEfficacy, (1 - 0.9) / 2, names = F)
  CIHigh90 <- quantile(sortedEfficacy, 0.9 + (1 - 0.9) / 2, names = F)
  
  CILow80 <- quantile(sortedEfficacy, (1 - 0.8) / 2, names = F)
  CIHigh80 <- quantile(sortedEfficacy, 0.8 + (1 - 0.8) / 2, names = F)
  
  return(
    list(
      mean = mean(sortedEfficacy),
      median = median(sortedEfficacy),
      CILow95 = CILow95,
      CIHigh95 = CIHigh95,
      CILow90 = CILow90,
      CIHigh90 = CIHigh90,
      CILow80 = CILow80,
      CIHigh80 = CIHigh80
    )
  )
}

#' @title Population mean jittering
#'
#' @description
#' Function jitters the mean of the population.
#'
#' Jittering is adding noise to the mean. The jittered mean is sampled from the distribution with the population mean and population standard deviation divided by the number of subjects in the population. The input population is provided in the form of population class objects (see the \code{Population-class} function for more details).
#'
#' \deqn{Mean_{jitter} \sim N(mean, \frac{sd}{N} )}{ MeanJitter ~ N(mean, sd/N)}
#'
#' @param blindPopulation \code{Population-class} object with N, mean, stdDev attributes
#'
#' @return
#' Jittered mean, numeric value
#'
#' @usage
#' JitterMean(blindPopulation)
#'
#' @examples
#' ## Data preparation
#' data(vaccinated)
#'
#' ## Example 1
#' vaccinated$mean
#' JitterMean(vaccinated)
#'
#' @export
JitterMean <- function(blindPopulation) {
  if (!is(blindPopulation, "Population")) {incorrectPopulationInput("blindPopulation")}
  
  blindSample <- rnorm(
    1,
    blindPopulation$mean,
    blindPopulation$stdDev / sqrt(round(blindPopulation$N))
  )
  
  return(blindSample)
}

#' @title PoDBAY efficacy equation
#'
#' @description
#' Function calculates the PoDBAY efficacy based on the PoD curve parameters and titer distribution parameters (mean, sd) for vaccinated and control groups.
#'
#' @param PoDParameters named data frame ("pmax", "slope", "et50"): PoD curve parameters
#' @param means named list ("vaccinated", "control"): mean values of vaccinated and control subjects titers
#' @param standardDeviations named list ("vaccinated", "control"): standard deviations of vaccinated and control subjects titers
#' @param adjustTiters boolean: set to TRUE if titer values should be adjusted, for details see \code{PoD} function
#' @param adjustFrom numeric: value specifying the detection limit, all values below the detection limit will be adjusted to adjustTo value
#' @param adjustTo numeric: value to which titers below the detection limit will be adjusted
#'
#' @return
#' efficacy: numeric value
#'
#' @usage
#' efficacyComputation(PoDParameters, 
#'                     means = NA, 
#'                     standardDeviations = NA,
#'                     adjustTiters = FALSE,
#'                     adjustFrom = NA,
#'                     adjustTo = NA)
#'
#' @examples
#' ## Data preparation
#' data(vaccinated)
#' data(control)
#' data(PoDParams)
#'
#' ## Example 1
#' means <- list(vaccinated = vaccinated$mean, control = control$mean)
#'
#' standardDeviations <- list(vaccinated = vaccinated$stdDev, control = control$stdDev)
#'
#' efficacyComputation(PoDParams, means, standardDeviations)
#'
#'
#' @details
#' \deqn{Efficacy = 1 - \frac{E[PoD_{vaccinated}]}{{E[PoD_{control}]} } }{ Efficacy = 1 - E[PoD|vaccinated]/ E[PoD|control]}.
#'
#' E[PoD] for each group is calculated as integral from -Inf to Inf of (titer density function) * (PoD Curve); for further details see Example2 and\code{ExpectedPoD} function.
#'
#' @export
efficacyComputation <- function(PoDParameters,
                                means = NA,
                                standardDeviations = NA,
                                adjustTiters = FALSE,
                                adjustFrom = NA,
                                adjustTo = NA) {
  
  if(any(is.na(match(names(PoDParameters), c("pmax", "et50", "slope") )))){
    stop(paste("The input value for PoDParameters is incorrect. 'PoDParameters' parameter has wrong names."))
  }
  
  if(any(is.na(match(names(means), c("vaccinated", "control") )))){
    stop(paste("The input value for means is incorrect. 'means' parameter has wrong names."))
  }
  
  if(any(is.na(match(names(standardDeviations), c("vaccinated", "control") )))){
    stop(paste("The input value for standardDeviations is incorrect. 'standardDeviations' parameter has wrong names."))
  }
  
  # PoD curve function
  funPoD <- function(x) PoD(x,
                            pmax = PoDParameters$pmax,
                            et50 = PoDParameters$et50,
                            slope = PoDParameters$slope,
                            adjustTiters = adjustTiters,
                            adjustFrom = adjustFrom,
                            adjustTo = adjustTo)
  
  if (standardDeviations$vaccinated == 0) {
    aucVaccinated <- funPoD(means$vaccinated)
  } else {
    funVaccinated <- function(x) dnorm(x, mean = means$vaccinated, sd = standardDeviations$vaccinated)
    aucVaccinated <- ExpectedPoD(funPoD, funVaccinated)
  }

  if (standardDeviations$control == 0) {
    aucControl <- funPoD(means$control)
  } else {
    funControl <- function(x) dnorm(x, mean = means$control, sd = standardDeviations$control)
    aucControl <- ExpectedPoD(funPoD, funControl)
  }
  
  efficacy <- 1 - aucVaccinated / aucControl
  
  return(efficacy)
}

#' @title Expected probability of disease
#'
#' @description Function calculates the integral of multiplication of two functions: PoD curve and titer probability density function.
#'
#' @param f.pod function(x):  PoD curve, estimated sigmoid function relating titers to a probability of disease
#' @param f.titer function(x): titer probability density function, distribution of titer values in a group.
#'
#' @return Value of the integral of the multiplication of the two functions
#'
#' @usage
#' ExpectedPoD(f.pod, f.titer)
#'
#' @examples
#' # Example 1 
#' data(vaccinated)
#' data(control)
#' data(PoDParams)
#'
#' # Defining the PoD curve
#' funPoD <- function(x) PoD(x, pmax = PoDParams$pmax, et50 = PoDParams$et50, slope = PoDParams$slope)
#'
#' # Defining the titer distribution for vaccinated and control groups
#' funVaccinated <- function(x) dnorm(x, mean = vaccinated$mean, sd = vaccinated$stdDev)
#' funControl <- function(x) dnorm(x, mean = control$mean, sd = control$stdDev)
#'
#' # Calculating the expected probability of disease 
#' aucVaccinated <- ExpectedPoD(funPoD, funVaccinated)
#' aucControl <- ExpectedPoD(funPoD, funControl)
#'
#' # PoDBAY efficacy estimation
#' efficacy <- 1 - aucVaccinated/aucControl
#'
#' @details Function calculates integral from -Inf to +Inf of titer probability density function multiplied by the PoD curve.
#'
#' It is used mainly in the PoDBAY efficacy calculation \code{efficacyComputation}.
#'
#' @export
ExpectedPoD <- function(f.pod, f.titer) {
  AUC <- integrate(function(x) f.pod(x) * f.titer(x), 
                   -Inf, 
                   Inf, 
                   abs.tol = 1e-10)$value
  return(AUC)
}