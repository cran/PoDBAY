#' @title PoD curve parameters estimation
#'
#' @description
#' Function estimates the PoD curve parameters (pmax, slope, et50) using \code{PoDMLE} function. Number of PoD curves estimated equals to the repeatCount input parameter.
#'
#' The estimation is performed using provided diseased and non-diseased subject level data.
#'
#' @param diseasedTiters numeric vector: all diseased titers, subject level data
#' @param nondiseasedTiters numeric vector: non-diseased titers from immunogenicity subset, subject level data
#' @param nondiseasedGenerationCount numeric: total number of non-diseased subjects in the clinical trial
#' @param repeatCount numeric: how many times is the dataset bootstrapped and the PoD curve parameter estimation performed
#' @param adjustTiters boolean: set to TRUE if titer values should be adjusted, for details see \code{PoD} function
#' @param adjustFrom numeric: value specifying the detection limit, all values below the detection limit will be adjusted to adjustTo value
#' @param adjustTo numeric: value to which titers below the detection limit will be adjusted
#'
#' @return
#' results: PoD curve parameters after resetting the disease status, named data.frame of estimated PoD curve parameters (pmax, slope, et50); see details for more information
#'
#' resultsPriorReset: PoD curve parameters prior to resetting the status, named data.frame of estimated PoD curve parameters (pmax, slope, et50); see details for more information
#'
#' failcount: number of iterations in which MLE failed to converge; see details for more information
#'
#' @usage
#' PoDParamEstimation(diseasedTiters,
#'                    nondiseasedTiters,
#'                    nondiseasedGenerationCount,
#'                    repeatCount = 500,
#'                    adjustTiters = FALSE,
#'                    adjustFrom = log2(10),
#'                    adjustTo = log2(5))
#'
#' @examples
#' ## Data preparation
#' data(diseased)
#' data(nondiseased)
#'
#' ## Example 1
#' # Creating imunogenicity subset, method = "Full"
#' NondiseasedImmunogenicitySubset <- 
#'     ImmunogenicitySubset(diseased, 
#'                          nondiseased, 
#'                          method = list(name = "Full", 
#'                                        value = "NA"))
#'
#' # Number of all non-diseased subjects in the clinical trial
#' nondiseasedGenerationCount <- nondiseased$N
#'
#' PoDParamEstimation(diseased$titers,
#'                    NondiseasedImmunogenicitySubset$titers,
#'                    nondiseasedGenerationCount,
#'                    repeatCount = 10)
#'
#' ## Example 2
#' # Creating imunogenicity subset, method = "Ratio", value = 4
#' NondiseasedImmunogenicitySubset <- 
#'     ImmunogenicitySubset(diseased, 
#'                          nondiseased, 
#'                          method = list(name = "Ratio", 
#'                                        value = 4))
#'                                        
#' # Number of all non-diseased subjects in the clinical trial
#' nondiseasedGenerationCount <- nondiseased$N
#'
#' PoDParamEstimation(diseased$titers,
#'                    NondiseasedImmunogenicitySubset$titers,
#'                    nondiseasedGenerationCount,
#'                    repeatCount = 10)
#'
#' @details
#'
#' diseasedTiters: subject level titers of all diseased in the clinical trial
#'
#' nondiseasedTiters: subject level titers of non-diseased subjects in the immunogenicity subset
#'
#' There are two possible scenarios
#'
#' \itemize{
#'   \item Full: Full information about non-diseased titers is available, i.e subject level data for all non-diseased subjects from the clinical trial (nondiseasedGenerationCount = number of all non-diseased subjects in the clinical trial).
#'
#'   \item Ratio or Fixed: Information about non-diseased titers is available only for the immunogenicity subset. In order to compensate for these missing titers we upsampling of this subset to the total number of non-diseased (nondiseasedGenerationCount) in the trial is needed. 
#'
#' }
#'
#' nondiseasedGenerationCount: number of all non-diseased subjects in the clinical trial
#'
#' NOTE: Number of estimated parameters can be lower than repeatCount as MLE does not necessary converge in all estimations; failcount (number of iterations in which MLE failed to converge) is also returned; for details see \code{MLE} function.
#'
#'
#' Function steps
#'
#' \itemize{
#'   \item Upsample non-diseased if needed (needed for methods Ratio and Fixed) - from immunogenicity subset size (N = NondiseasedImmunogenicitySubset$N) to the whole trial size (N = nondiseasedGenerationCount). For details see \code{GenerateNondiseased} function.
#'
#'   \item Estimate PoD curve: resultsPriorReset
#'
#'   \item Reset disease status: the purpose is to estimate the confidence intervals of the PoD curve and its parameters
#'
#'   Part of the reset disease status procedure is the non-parametric bootstrap: titers of diseased and non-diseased subjects are pooled, and associated PoDs are calculated using their titer values and estimated PoD curve. Based on the subject level probabilities (PoDs), the disease status is reestimated.
#'
#'   \item Re-estimate PoD curve: new diseased and non-diseased titers are used to reestimate the PoD curve
#' }
#'
#' @export
PoDParamEstimation <- function(diseasedTiters,
                               nondiseasedTiters,
                               nondiseasedGenerationCount,
                               repeatCount = 500,
                               adjustTiters = FALSE,
                               adjustFrom = log2(10),
                               adjustTo = log2(5)) {
  
  if (floor(nondiseasedGenerationCount) < length(nondiseasedTiters)) {
    warning(paste("The input value for \"nondiseasedGenerationCount\" is lower than number of nondiseasedTiters"))
  }
  
  results <- data.frame("pmax" = 0, "slope" = 0, "et50" = 0)
  resultsPriorReset <- data.frame("pmax" = 0, "slope" = 0, "et50" = 0)
  
  failCount <- 0
  for (i in 1:repeatCount) {
    print(paste("iteration", i))
    
    # upsample non-diseased immunogenicity subset to the size of the trial
    newNondiseased <- GenerateNondiseased(nondiseasedTiters,
                                         nondiseasedGenerationCount)
    
    # estimate PoD curve parameters
    estimatedParameters <- PoDMLE(
      nondiseasedTiters = newNondiseased,
      diseasedTiters = diseasedTiters,
      adjustTiters = adjustTiters,
      adjustFrom = adjustFrom,
      adjustTo = adjustTo
    )
    
    #if the estimation fails, note the failure and skip to the next iteration
    if (is.null(estimatedParameters)) {
      failCount <- failCount + 1
      next
    }
    
    ## reset status ##
    
    # pool non-diseased and diseased titers
    completePopulation <- generatePopulation(0)
    completePopulationTiters <- c(newNondiseased, diseasedTiters)
    
    # bootstrap all titers
    completePopulation$titers <- sample(sample(completePopulationTiters),
                                        length(completePopulationTiters),
                                        replace = T)
    
    # assign to each subject (titer) its probability of disease
    completePopulation$assignPoD(
      PoD(
        completePopulation$titers,
        estimatedParameters$pmax,
        estimatedParameters$et50,
        estimatedParameters$slope,
        adjustTiters = adjustTiters,
        adjustFrom = adjustFrom,
        adjustTo = adjustTo
      ))
    
    # re-assign the disease status based on the probabilities
    dStatus <- mapply(rbinom, n = 1, size = 1, prob = completePopulation$PoDs)
    completePopulation$diseaseStatus <- sapply(dStatus, numToBool)
    
    # new diseased and non-diseased populations
    completePopulationDiseased <- generatePopulation(0)
    completePopulationDiseased$N <- completePopulation$getDiseasedCount()
    completePopulationDiseased$titers <- completePopulation$getDiseasedTiters()
    completePopulationDiseased$diseaseStatus <- rep(TRUE, completePopulationDiseased$N)
    
    completePopulationNondiseased <- generatePopulation(0)
    completePopulationNondiseased$N <- completePopulation$getNondiseasedCount()
    completePopulationNondiseased$titers <- completePopulation$getNondiseasedTiters()
    completePopulationNondiseased$diseaseStatus <- rep(FALSE, completePopulationNondiseased$N)
    
    # create non-diseased immunogenicity sample of the new population
    # choose appropriate method for the immunogenicty sample creation = Full or Ratio based on the input to the function
    # if input FULL SAMPLE -> method = "FULL"
    # otherwise method = "Ratio", method value = calculated based on the data inputs #Nondiseased/#Diseased
    
    if (round(nondiseasedGenerationCount) == length(nondiseasedTiters)) {
      method <- list(name = "Full",
                     value = NA)
    } else {
      method <- list(name = "Ratio",
                     value = length(nondiseasedTiters) / length(diseasedTiters))
    }

    immunogenicitySample <- ImmunogenicitySubset(completePopulationDiseased,
                                                 completePopulationNondiseased,
                                                 method = method)
    
    immunogenicitySample$diseaseStatus <- grepl("TRUE", names(immunogenicitySample$titers))
    
    nondiseasedSampleTiters <- immunogenicitySample$getNondiseasedTiters()
    
    # upsample the new non-diseased immunogenicity subset to the full non-diseased population
    completePopulationNewNondiseased <- GenerateNondiseased(nondiseasedSampleTiters, completePopulationNondiseased$N )
    
    # estimate new PoD curve parameters
    estimatedParameters2 <- PoDMLE(
      nondiseasedTiters = completePopulationNewNondiseased,
      diseasedTiters = completePopulationDiseased$titers,
      adjustTiters = adjustTiters,
      adjustFrom = adjustFrom,
      adjustTo = adjustTo
    )
    
    if (is.null(estimatedParameters2)) {
      failCount <- failCount + 1
      next
    }
    
    resultsPriorReset <- rbind(resultsPriorReset,
                               c(estimatedParameters$pmax,
                                 estimatedParameters$slope,
                                 estimatedParameters$et50)
    )
    
    results <- rbind(results,
                     c(estimatedParameters2$pmax,
                       estimatedParameters2$slope,
                       estimatedParameters2$et50)
    )
  }
  return(list(results = results[-1, ],
              resultsPriorReset = resultsPriorReset[-1, ],
              failCount = failCount)
  )
}

#' @title Generation of upsampled non-diseased subjects titers
#'
#' @description Function upsamples (by random sampling with replacement) titers from the immunogenicity subset to the required size.
#'
#' If the size of the immunogenicity subset matches the required size, nothing happens and the original titers from the immunogenicity subset are returned.
#'
#' @param blindNondiseasedTiters numeric vector: vector of non-diseased subjects titer values
#' @param nondiseasedCount numeric: total number of non-diseased subjects, required size of the non-diseased population
#'
#' @return
#' nondiseasedTiters: numeric vector of all non-diseased subjects titers 
#'
#' @usage
#' GenerateNondiseased(blindNondiseasedTiters, nondiseasedCount)
#'
#' @examples
#' ## Data preparation
#' data(nondiseased)
#'
#' ## Example 1
#' # Creating imunogenicity subset, method = "Full"
#' NondiseasedImmunogenicitySubset <- 
#'     ImmunogenicitySubset(diseased, 
#'                          nondiseased, 
#'                          method = list(name = "Full", 
#'                                        value = "NA"))
#'
#' # Number of all non-diseased subjects in the clinical trial
#' nondiseasedGenerationCount <- nondiseased$N
#'
#' # Upsampling of non-diseased titers
#' GenerateNondiseased(NondiseasedImmunogenicitySubset$titers, nondiseasedGenerationCount)
#'
#' @details
#' The inputs should come from immunogenicity subset. "nondiseasedCount" represents number of all non-diseased patients in the clinical trial.
#'
#' Immunogenicity subset populations are obtained from function \code{BlindSampling}. Immunogenicity subset represents a sample from the non-diseased population.
#'
#' In this function, sampling with replacement to the required "nondiseasedCount" of the immunogenecitry subset is performed. The function is used inside \code{PoDParamEstimation} function.
#'
#' @export
GenerateNondiseased <- function(blindNondiseasedTiters, nondiseasedCount) {
  blindNondiseasedTiters <- unlist(blindNondiseasedTiters)
  nondiseasedTiters <-  if (round(length(blindNondiseasedTiters)) == floor(nondiseasedCount))
  {blindNondiseasedTiters}
  else {sample(blindNondiseasedTiters, nondiseasedCount, replace = T)}
  
  return(nondiseasedTiters = nondiseasedTiters)
}

#' @title Setup for the maximum likelihood estimation (MLE)
#'
#' @description Function estimates the optimal PoD curve parameters (pmax, et50, slope) using diseased and non-diseased titers. Initial guess of the slope parameter needs to be provided as an input to the optimization, as well as the lowTiterPercent parameter, which is needed for initial guess of the pmax parameter calculation.
#'
#' @param nondiseasedTiters numeric vector: non-diseased subjects titers
#' @param diseasedTiters numeric vector: diseased subjects titers
#' @param initialSlope numeric: initial guess of the slope parameter for the optimization function
#' @param lowTiterPercent numeric: value in the interval (0,1) - it represents a fraction of bottom titer values of the whole clinical trial used for calculation of inital guess of the pmax parameter.
#' @param adjustTiters boolean: set to TRUE if titer values should be adjusted, for details see \code{PoD} function
#' @param adjustFrom numeric: value specifying the detection limit, all values below the detection limit will be adjusted to adjustTo value
#' @param adjustTo numeric: value to which titers below the detection limit will be adjusted
#'
#' @return
#' list("et50", "slope", "pmax"), named list of PoD paraters: if MLE converges.
#'
#' Null: if MLE does not converge.
#'
#' @usage PoDMLE(nondiseasedTiters,
#'               diseasedTiters,
#'               adjustTiters = FALSE,
#'               adjustFrom = log2(10),
#'               adjustTo = log2(5),
#'               initialSlope = 6,
#'               lowTiterPercent = 0.2)
#'
#' @details
#'
#' Initial guess of pmax = (number of diseased in the bottom titers + 0.5) / (number of non-diseased and diseased in the bottom titers + 0.5),
#' Initial et50 = intersection point of distributions of non-diseased and diseased groups. If L-BFGS-B optimization fails to converge, a new et50 initial guess is set to median value of all titers.
#'
#' PoDMLE function estimates the PoD curve parameters by maximizing the likelihood value (see \code{MLE} function for details) based on the provided titers for diseased and non-diseased groups.
#'
#' The \code{optim} function is used for optimization with method = "L-BFGS-B", 500 maximum iterations, (0.1,Inf) boundaries for et50, (1e-6,1) boundaries for pmax and (-slopeBoundary, slopeBoundary) boundaries for slope.
#'
#' NOTE: The reason for slope boundary settings is because from certain value of slope parameter the shape of the PoD curve and the corresponding PoD values for given titers are almost identical. This parameter is expected to limit the resulting slope value and help MLE to converge to optimal parameters.
#' The value of "slopeBoundaries" is calculated as described by Dunning, 2015 (https://doi.org/10.1186/s12874-015-0096-9).
#'
#' @examples
#' ## EXAMPLE 1:
#' # Data preparation
#' data(diseased)
#' data(nondiseased)
#'
#' # PoD curve parameter estimation
#' PoDMLE(nondiseased$titers,
#'        diseased$titers)
#'
#' ## EXAMPLE 2:
#' ## initialSlope and lowTiterPercent variables are adjusted.
#' PoDMLE(nondiseased$titers,
#'        diseased$titers,
#'        initialSlope = 5,
#'        lowTiterPercent = 0.3)
#'
#' @export
PoDMLE <- function(nondiseasedTiters,
                   diseasedTiters,
                   adjustTiters = FALSE,
                   adjustFrom = log2(10),
                   adjustTo = log2(5),
                   initialSlope = 6,
                   lowTiterPercent = 0.2) {
  
  # initial guess of pmax 
  titers <- c(nondiseasedTiters, diseasedTiters)
  
  diseaseStatus <- c(rep(0, length(nondiseasedTiters)),
                     rep(1, length(diseasedTiters)))
  newOrder <- order(titers)
  titers <- titers[newOrder]
  diseaseStatus <- diseaseStatus[newOrder]
  
  lowTitersStatus <- diseaseStatus[1:round(length(diseaseStatus) * lowTiterPercent)]
  initialPmax <- (sum(lowTitersStatus) + 0.5) / (length(lowTitersStatus) + 0.5)
  
  # slope boundaries 
  coeff <- 9.1902
  range <- max(titers) - min(titers)
  slopeUpperBound <- coeff * 50 / range
  
  tryError <- F
  tryCatch(
    {  # initial guess of et50 
      RootFunction <- function(x) dnorm(x, 
                                        mean = mean(diseasedTiters), 
                                        sd = sd(diseasedTiters)) - dnorm(x, 
                                                                         mean = mean(nondiseasedTiters), 
                                                                         sd = sd(nondiseasedTiters))
      RootSearch <- uniroot(RootFunction, interval = c(mean(diseasedTiters), mean(nondiseasedTiters)))
      initialEt50 <- RootSearch$root
      
      initialParams <- c(initialEt50, initialSlope, initialPmax)
      names(initialParams) <- c("et50", "slope", "pmax")
      
      estimatedParameters <- optim(
        par = initialParams,
        fn = cppMLE,
        nondiseasedTiters = nondiseasedTiters,
        diseasedTiters = diseasedTiters,
        adjustTiters = adjustTiters,
        adjustFrom = adjustFrom,
        adjustTo = adjustTo,
        lower = c(0.1, -slopeUpperBound, 1e-6),
        upper = c(Inf, slopeUpperBound, 1),
        method = "L-BFGS-B",
        control = list(fnscale = -1, maxit = 500))
    },
    error = function(e) {
      tryError <<- T
    }
  )
  
  # if L-BFGS-B fails to converge or other error (e.g. no intersection of inf/uninf is found)
  if (tryError) {
    initialEt50 <- median(titers[titers > 0])
    initialParams <- c(initialEt50, initialSlope, initialPmax)
    names(initialParams) <- c("et50", "slope", "pmax")
    tryError <- F
    tryCatch(
      estimatedParameters <- optim(
        par = initialParams,
        fn = cppMLE,
        nondiseasedTiters = nondiseasedTiters,
        diseasedTiters = diseasedTiters,
        adjustTiters = adjustTiters,
        adjustFrom = adjustFrom,
        adjustTo = adjustTo,
        lower = c(0.1, -slopeUpperBound, 1e-6),
        upper = c(Inf, slopeUpperBound, 1),
        method = "L-BFGS-B",
        control = list(fnscale = -1, maxit = 500)),
      error = function(e) {
        tryError <<- T
      }
    )
  }
  
  if (tryError) {return(NULL)}
  
  return(list("et50"  = estimatedParameters$par[1],
              "slope" = estimatedParameters$par[2],
              "pmax"  = estimatedParameters$par[3]))
}

#' @title Maximum likelihood estimation: cpp
#'
#' @description Function calculates the log likelihood value which is used after the initial guesses of the parameters are set in the \code{PoDMLE} function.
#'
#' @param params named numeric vector: PoD curve parameters ("et50", "slope", "pmax")
#' @param nondiseasedTiters numeric vector: non-diseased subjects titers
#' @param diseasedTiters numeric vector: diseased subjects titers
#' @param adjustTiters boolean: set to TRUE if titer values should be adjusted, for details see \code{PoD} function
#' @param adjustFrom numeric: value specifying the detection limit, all values below the detection limit will be adjusted to adjustTo value
#' @param adjustTo numeric: value to which titers below the detection limit will be adjusted
#'
#' @return log likelihood, numeric value
#'
#' @usage
#' cppMLE(params,
#'        nondiseasedTiters,
#'        diseasedTiters,
#'        adjustTiters = FALSE,
#'        adjustFrom = log2(10),
#'        adjustTo = log2(5))
#'
#' @examples
#' # Data preparation
#' data(diseased)
#' data(nondiseased)
#' data(PoDParams)
#'
#' # MLE calculation
#' cppMLE(PoDParams, nondiseased$titers, diseased$titers)
#'
#' @details cppMLE function is used inside of PoDMLE function and estimates the PoD curve paramers.
#'
#' Based on the provided titers for diseased and non-diseased groups the PoD curve parameters which maximize the log likelihood are chosen as optimal.
#'
#' Difference between MLE and cppMLE is only that cppMLE use cppPoD function instead of PoD. This step significantly improves the computation speed and provides the same results.
#'
#' @export
cppMLE <- function(params,
                   nondiseasedTiters,
                   diseasedTiters,
                   adjustTiters = FALSE,
                   adjustFrom = log2(10),
                   adjustTo = log2(5)) {
  
  if (!adjustTiters) {
    adjustFrom <- NA
    adjustTo <- NA
  }
  
  if(any(is.na(match(names(params), c("pmax", "et50", "slope") )))){
    stop(paste("The input value for params is incorrect. 'params' parameter has wrong names."))
  }
  
  pmax <- as.double(params["pmax"])
  et50 <- as.double(params["et50"])
  slope <- as.double(params["slope"])
  
  probDiseased <- cppPoD(titer = diseasedTiters,
                         pmax, et50, slope,
                         adjustTiters = adjustTiters,
                         adjustFrom = adjustFrom,
                         adjustTo = adjustTo)
  
  probNondiseased <- cppPoD(titer = nondiseasedTiters,
                           pmax, et50, slope,
                           adjustTiters = adjustTiters,
                           adjustFrom = adjustFrom,
                           adjustTo = adjustTo)
  
  logLikelihood <- sum(log(probDiseased)) + sum(log(1 - probNondiseased))
  return(logLikelihood)
}

#' @title Maximum Likelihood estimation
#'
#' @description Function calculates the log likelihood value which is used after the initial guesses of the parameters are set in the \code{PoDMLE} function.
#'
#' @param params named numeric vector: PoD curve parameters (et50, slope, pmax)
#' @param nondiseasedTiters numeric vector: non-diseased subjects titers
#' @param diseasedTiters numeric vector: diseased subjects titers
#' @param adjustTiters boolean: set to TRUE if titer values should be adjusted, for details see \code{PoD} function
#' @param adjustFrom numeric: value specifying the detection limit, all values below the detection limit will be adjusted to adjustTo value
#' @param adjustTo numeric: value to which titers below the detection limit will be adjusted
#'
#' @return log likelihood, numeric value
#'
#' @usage
#' MLE(params,
#'     nondiseasedTiters,
#'     diseasedTiters,
#'     adjustTiters = FALSE,
#'     adjustFrom = log2(10),
#'     adjustTo = log2(5))
#'
#' @examples
#' # Data preparation
#' data(diseased)
#' data(nondiseased)
#' data(PoDParams)
#'
#' # MLE calculation
#' MLE(PoDParams, nondiseased$titers, diseased$titers)
#'
#' @details MLE function is used inside of PoDMLE function and esimates the PoD curve parameters.
#'
#' Based on the provided titers for diseased and non-diseased subjects the PoD curve parameters which maximize the log likelihood are chosen as optimal estimates of parameters.
#'
#' @export
MLE <- function(params,
                nondiseasedTiters,
                diseasedTiters,
                adjustTiters = FALSE,
                adjustFrom = log2(10),
                adjustTo = log2(5)) {
  
  if(any(is.na(match(names(params), c("pmax", "et50", "slope") )))){
    stop(paste("The input value for params is incorrect. 'params' parameter has wrong names."))
  }
  
  pmax <- as.double(params["pmax"])
  et50 <- as.double(params["et50"])
  slope <- as.double(params["slope"])
  
  probDiseased <- PoD(titer = diseasedTiters,
                      pmax, et50, slope,
                      adjustTiters = adjustTiters,
                      adjustFrom = adjustFrom,
                      adjustTo = adjustTo)
  
  probNondiseased <- PoD(titer = nondiseasedTiters,
                        pmax, et50, slope,
                        adjustTiters = adjustTiters,
                        adjustFrom = adjustFrom,
                        adjustTo = adjustTo)
  
  logLikelihood <- sum(log(probDiseased)) + sum(log(1 - probNondiseased))
  return(logLikelihood)
}

#' @title Confidence intervals of PoD curve parameters
#'
#' @description Function calculates confidence intervals of the PoD curve parameters (pmax, et50, slope) at user-defined confidence level.
#'
#' @param estimatedParameters output of \code{PoDParamEstimation} function
#' @param ci numeric: value from (0, 1) interval, confidence level of interest
#'
#' @return CI of all PoD curve parameters
#'
#' @usage
#' PoDParamsCI(estimatedParameters, ci = 0.95)
#'
#' @export
PoDParamsCI <- function(estimatedParameters, ci = 0.95) {
  
  if(any(is.na(match(names(estimatedParameters), c("pmax", "et50", "slope") )))){
    stop(paste("The input value for estimatedParameters is incorrect. 'estimatedParameters' parameter has wrong names."))
  }
  
  sortedPmax <- sort(estimatedParameters[["pmax"]])
  sortedEt50 <- sort(estimatedParameters[["et50"]])
  sortedSlope <- sort(estimatedParameters[["slope"]])
  
  PmaxCILow <- quantile(sortedPmax, (1 - ci) / 2, names = F)
  PmaxCIHigh <- quantile(sortedPmax, ci + (1 - ci) / 2, names = F)
  
  Et50CILow <- quantile(sortedEt50, (1 - ci) / 2, names = F)
  Et50CIHigh <- quantile(sortedEt50, ci + (1 - ci) / 2, names = F)
  
  SlopeCILow <- quantile(sortedSlope, (1 - ci) / 2, names = F)
  SlopeCIHigh <- quantile(sortedSlope, ci + (1 - ci) / 2, names = F)
  
  return(
    list(
      PmaxCILow = PmaxCILow,
      PmaxCIHigh = PmaxCIHigh,
      Et50CILow = Et50CILow,
      Et50CIHigh = Et50CIHigh,
      SlopeCILow = SlopeCILow,
      SlopeCIHigh = SlopeCIHigh
    )
  )
}

#' @title Confidence intervals of PoD curve parameters at three confidence levels
#'
#' @description Function calculates confidence intervals (80\%, 90\% and user-defined) of the PoD curve parameters (pmax, et50, slope).
#'
#' @param estimatedParameters output of \code{PoDParamEstimation} function
#' @param ci numeric: value from (0, 1) interval, confidence level of interest
#'
#' @return CI of all PoD curve parameters
#'
#' @usage
#' PoDParamsCICoverage(estimatedParameters, ci = 0.95)
#'
#' @export
PoDParamsCICoverage <- function(estimatedParameters, ci = 0.95) {
  
  if(any(is.na(match(names(estimatedParameters), c("pmax", "et50", "slope") )))){
    stop(paste("The input value for estimatedParameters is incorrect. 'estimatedParameters' parameter has wrong names."))
  }
  
  sortedPmax <- sort(estimatedParameters[["pmax"]])
  sortedEt50 <- sort(estimatedParameters[["et50"]])
  sortedSlope <- sort(estimatedParameters[["slope"]])
  
  PmaxCILow95 <- quantile(sortedPmax, (1 - ci) / 2, names = F)
  PmaxCIHigh95 <- quantile(sortedPmax, ci + (1 - ci) / 2, names = F)
  PmaxCILow90 <- quantile(sortedPmax, (1 - 0.9) / 2, names = F)
  PmaxCIHigh90 <- quantile(sortedPmax, 0.9 + (1 - 0.9) / 2, names = F)
  PmaxCILow80 <- quantile(sortedPmax, (1 - 0.8) / 2, names = F)
  PmaxCIHigh80 <- quantile(sortedPmax, 0.8 + (1 - 0.8) / 2, names = F)
  
  Et50CILow95 <- quantile(sortedEt50, (1 - ci) / 2, names = F)
  Et50CIHigh95 <- quantile(sortedEt50, ci + (1 - ci) / 2, names = F)
  Et50CILow90 <- quantile(sortedEt50, (1 - 0.9) / 2, names = F)
  Et50CIHigh90 <- quantile(sortedEt50, 0.9 + (1 - 0.9) / 2, names = F)
  Et50CILow80 <- quantile(sortedEt50, (1 - 0.8) / 2, names = F)
  Et50CIHigh80 <- quantile(sortedEt50, 0.8 + (1 - 0.8) / 2, names = F)
  
  SlopeCILow95 <- quantile(sortedSlope, (1 - ci) / 2, names = F)
  SlopeCIHigh95 <- quantile(sortedSlope, ci + (1 - ci) / 2, names = F)
  SlopeCILow90 <- quantile(sortedSlope, (1 - 0.9) / 2, names = F)
  SlopeCIHigh90 <- quantile(sortedSlope, 0.9 + (1 - 0.9) / 2, names = F)
  SlopeCILow80 <- quantile(sortedSlope, (1 - 0.8) / 2, names = F)
  SlopeCIHigh80 <- quantile(sortedSlope, 0.8 + (1 - 0.8) / 2, names = F)
  
  return(
    list(
      PmaxCILow95 = PmaxCILow95,
      PmaxCIHigh95 = PmaxCIHigh95,
      PmaxCILow90 = PmaxCILow90,
      PmaxCIHigh90 = PmaxCIHigh90,
      PmaxCILow80 = PmaxCILow80,
      PmaxCIHigh80 = PmaxCIHigh80,
      Et50CILow95 = Et50CILow95,
      Et50CIHigh95 = Et50CIHigh95,
      Et50CILow90 = Et50CILow90,
      Et50CIHigh90 = Et50CIHigh90,
      Et50CILow80 = Et50CILow80,
      Et50CIHigh80 = Et50CIHigh80,
      SlopeCILow95 = SlopeCILow95,
      SlopeCIHigh95 = SlopeCIHigh95,
      SlopeCILow90 = SlopeCILow90,
      SlopeCIHigh90 = SlopeCIHigh90,
      SlopeCILow80 = SlopeCILow80,
      SlopeCIHigh80 = SlopeCIHigh80
    )
  )
}

#' @title PoD curve point estimate
#'
#' @description Function returns PoD curve parameters corresponding to the point estimate of PoD curve.
#'
#' @param resultsPriorReset named data frame ("pmax", "slope", "et50"): set of estimated PoD curve parameters before resetting the disease status; for further details see \code{PoDParamEstimation} function.
#' @param titers numeric vector: a grid of titers for PoD curve point estimate calculation
#' @param optim_titers logical: TRUE for a predefined sequence of titers 
#'
#' @return
#' paramsPointEstimate: named data frame of PoD curve parameters corresponding to the PoD curve point estimate
#'
#' @usage
#' PoDParamPointEstimation(resultsPriorReset, 
#'                         titers = seq(from = 0, to = 20, by = 0.01), 
#'                         optim_titers = FALSE)
#'
#' @details
#' For each of estimated PoD curves in resultsPriorReset, the function values (probabilities of disease, PoD) for provided grid of titers are calculated.
#'
#' Median of function values (PoDs) at each provided titer is calculated.
#'
#' Subsequently, the PoD curve model is fitted to the median datapoins using \code{fitPoD} function, in order to get PoD curve parameters close to this median curve.
#'
#' 
#'
#' @examples
#' ## Data preparation
#' data(estimatedParameters)
#'
#' ## Example 1
#' # titers for which we want to optimize the functional values
#' titers <- seq(from = 0, to = 20, by = 0.01)
#'
#' # Point estimate of PoD curve
#' PoDParamPointEstimation(estimatedParameters$resultsPriorReset, titers)
#'
#' @export
PoDParamPointEstimation <- function(resultsPriorReset,
                                    titers = seq(from = 0, to = 20, by = 0.01),
                                    optim_titers = FALSE){
  
  if(any(is.na(match(names(resultsPriorReset), c("pmax", "et50", "slope") )))){
    stop(paste("The input value for resultsPriorReset is incorrect. 'resultsPriorReset' parameter has wrong names."))
  }
  
  pmax <- resultsPriorReset["pmax"]
  et50 <- resultsPriorReset["et50"]
  slope <- resultsPriorReset["slope"]
  
  if (optim_titers) {
    et50_med <- round(median(unlist(et50)),1)
    titers <- c(seq(0, et50_med - 1, by = 0.1),
                seq(et50_med - 0.99, et50_med + 0.99, by = 0.01),
                seq(et50_med + 1, 2 * et50_med , by = 0.1))
  }
  
  # number of estimated PoD curves
  repeatCount <- nrow(resultsPriorReset) 
  
  # for each estimated PoD curve, calculate functional values
  functionValues <- matrix(NA, nrow = repeatCount, ncol = length(titers))
  for (i in 1:repeatCount) { 
    functionValues[i, ] <- PoD(titers, 
                               pmax = resultsPriorReset[i, 1], 
                               et50 = resultsPriorReset[i, 3], 
                               slope = resultsPriorReset[i, 2], 
                               adjustTiters = FALSE) # we are fitting the curve so there is no need to adjust titers
  }
  
  # median of functional values
  curveMedian <- apply(functionValues, 2, median)
  
  # median PoD curve parameters
  paramsMedian <- apply(resultsPriorReset, 2, median)
  
  # set initial guess of PoD curve oarameters for optimization
  paramsInitial <- (1 - (paramsMedian["slope"] / 100)) * paramsMedian

  # find optimal PoD curve
  paramsPointEstimate <- optim(
    par = paramsInitial,
    fn = fitPoD,
    TitersInput = titers,
    CurveTitersMedian = curveMedian,
    lower = c(1e-6, -100, 0.1),
    upper = c( 1, 100, Inf),
    method = "L-BFGS-B",
    control = list(fnscale = -1, maxit = 1000, factr = 1e4)
  )
  
  return(paramsPointEstimate = data.frame("pmax" = paramsPointEstimate$par[1],
                                          "slope" = paramsPointEstimate$par[2],
                                          "et50" = paramsPointEstimate$par[3]))
  
}

#' @title PoD curve: fitting function
#'
#' @description Function calculates the root mean squared error (RMSE) between provided PoD values and calculated PoD values. The latter are calculated using for provided titers and provided PoD curve parameters.
#'
#' By using the input titers \code{PoDParamPointEstimation} function and median of the estimated set of PoD curve parameters (output of \code{PoDParamEstimation} function), the point estimate of PoD curve can be obtained (for details see \code{PoDParamPointEstimation} function).
#'
#' @param params named data frame ("pmax", "slope", "et50"): provided PoD curve parameters
#' @param TitersInput numeric vector: provided titers 
#' @param CurveTitersMedian numeric vector: provided PoD values 
#'
#' @return
#' negative RMSE
#'
#' @usage
#' fitPoD(params, TitersInput, CurveTitersMedian)
#'
#' @details
#'
#' \deqn{RMSE = \sqrt{\frac{\sum_{i}^{N} (PoD_{median}(titers) - PoD_{optimized}(titers))^2}{N}}}{ RMSE = sqrt( mean( (PoDmedian(titers) - PoDoptimized(titers))2 ) )}
#'
#' @examples
#' ## Data preparation
#' data(estimatedParameters)
#' data(PoDParams)
#'
#' ## Example 1
#'
#' # grid of titers
#' TitersInput <- seq(from = 0, to = 20, by = 0.01)
#'
#' # for each estimated PoD curve calculate functional values
#' functionValues <- 
#'   matrix(NA, 
#'          nrow = nrow(estimatedParameters$resultsPriorReset), 
#'          ncol = length(TitersInput))
#' 
#' for (i in 1:nrow(estimatedParameters$resultsPriorReset)) { 
#'   functionValues[i,] <- PoD(TitersInput,
#'   pmax = estimatedParameters$resultsPriorReset[i,1], 
#'   et50 = estimatedParameters$resultsPriorReset[i,3], 
#'   slope = estimatedParameters$resultsPriorReset[i,2], adjustTiters = FALSE)
#' }
#'
#' # functional values corresponding to the median of the estimated PoD curve parameters
#' CurveTitersMedian <- apply(functionValues, 2, median)
#'
#' # squared error of CurveTitersMedian and functional values of "params" curve
#' fitPoD(PoDParams, TitersInput, CurveTitersMedian)
#'
#' @export
fitPoD <- function(params, TitersInput, CurveTitersMedian) {
  
  if(any(is.na(match(names(params), c("pmax", "et50", "slope") )))){
    stop(paste("The input value for params is incorrect. 'params' parameter has wrong names."))
  }
  
  pmax  <- as.double(params["pmax"])
  slope <- as.double(params["slope"])
  et50  <- as.double(params["et50"])
  
  # functional values
  CurveTitersCalculated <- PoD(titer = TitersInput,
                               pmax = pmax,
                               slope = slope,
                               et50 = et50)
  
  # RMSE
  error <- - sqrt(mean( (CurveTitersMedian - CurveTitersCalculated) ^ 2) )
  
  return(error)
}

#' @title PoD curve confidence ribbon
#'
#' @description Supplementary function for \code{PoDCurvePlot} function. Function calculates the confidence ribbon around the PoD curve.
#' 
#' @param data numeric vector for which we the confidence intervals should be calculated
#' @param ci numeric: required confidence level
#'
#' @return
#' lower bound of CI
#' median value
#' upper bound of CI
#'
#' @usage
#' PoDCI(data, ci = 0.95)
#'
#' @examples
#' ## Data preparation
#' data <- 0:100
#'
#' ## Example 1
#' PoDCI(data,
#'       ci = 0.95)
#'
#' @export
PoDCI <- function(data, ci = 0.95){
  lower <- quantile(data, (1-ci)/2, na.rm = TRUE)
  upper <- quantile(data, ci + (1-ci)/2, na.rm  = TRUE )
  median <- median(data, na.rm  = TRUE)
  
  return(list(lower = lower,
              median = median,
              upper = upper)
  )
}

#' @title PoD curve: plot
#'
#' @description Supplementary function for plotting the PoD curve with the confidence ribbon (of a required level). Input values are related to PoDBAY package structure. 
#' See \code{vignette("EfficacyEstimation", package = "PoDBAY")} for an example of application of this function.
#' 
#' @param titers numeric vector: grid of titers at which the confidence ribbon should be calculated
#' @param estimatedParameters estimatedParameters named data frame (pmax, slope, et50): set of estimated PoD curve parameters, output of \code{PoDParamEstimation} function.
#' @param ci numeric, required confidence level
#'
#' @return
#' PoD curve plot 
#'
#' @usage
#' PoDCurvePlot(titers,
#'              estimatedParameters,
#'              ci = 0.95)
#'
#' @examples
#' ## Data preparation
#' library(ggplot2)
#' data(PoDParams)
#' data(estimatedParameters)
#'
#' ## Example 1
#' # titers for which we want calculate the confidence intervals
#' titers <- seq(from = 0, to = 15, by = 0.01)
#'
#' # squared error of CurveTitersMedian and functional values of "params" curve
#' PoDCurvePlot(titers,
#'              estimatedParameters,
#'              ci = 0.95)
#'
#' @export
PoDCurvePlot <- function(titers,
                         estimatedParameters,
                         ci = 0.95){
  
  lowerCI <- NULL
  upperCI <- NULL
  
  # plot format setting
  size <- 12
  themePlot <-   
    theme(
      plot.title = element_text(size=size,hjust = 0.5),
      axis.text.y   = element_text(size=size),
      axis.text.x   = element_text(size=size),
      axis.title.y  = element_text(size=size),
      axis.title.x  = element_text(size=size),
      legend.text = element_text(size = size),
      panel.background = element_blank(),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
    ) 
  
  # calculate PoDs for the point estimate
  PoDParamsPointEst <- PoDParamPointEstimation(estimatedParameters$resultsPriorReset, titers)
  
  PoDParamsPointEstPoD <- PoD(titers, pmax = PoDParamsPointEst$pmax, et50 = PoDParamsPointEst$et50, slope = PoDParamsPointEst$slope)
  
  # calculate PoDs for all PoD curves in estimatedParameters$results
  PoDMatrix <- matrix(0, ncol = length(titers), nrow = nrow(estimatedParameters$results))
  for (i in 1:nrow(estimatedParameters$results)) {
    PoDMatrix[i,] <- PoD(titers, pmax = estimatedParameters$results$pmax[i], et50 = estimatedParameters$results$et50[i], slope = estimatedParameters$results$slope[i])
  }
  
  # calculate CI of functional values 
  PoDCIMatrix <- matrix(0, ncol = 3, nrow = ncol(PoDMatrix))
  for (i in 1:ncol(PoDMatrix)) {
    PoDCIMatrix[i,] <- unlist(PoDCI(PoDMatrix[,i], ci = ci))
  }
  
  # prepare data for ggplot
  PoDCIdf <- as.data.frame(PoDCIMatrix)
  colnames(PoDCIdf) <- c("lowerCI", "median", "upperCI")
  PoDCIdf$titers <- titers
  
  # ggplot
  PoDCurvePlot <- ggplot(PoDCIdf) + 
    geom_line(aes(x = titers, y = lowerCI), linetype = 5) +
    geom_line(aes(x = titers, y = PoDParamsPointEstPoD )) +
    geom_line(aes(x = titers, y = upperCI), linetype = 5) +
    ylab("PoD") + 
    ggtitle(paste0("PoD curve with ",ci*100,"% confidence intervals")) +
    themePlot
  
  return(PoDCurvePlot)
}

#' @title Optimization objective function: efficacy squared error
#'
#' @description
#' Function calculates squared difference between input (reference value, or for example true in the simulation setup) efficacy and
#' efficacy calculated based on input parameters of PoD curve and input titer 
#' distributions of vaccinated and control groups.
#' 
#' @param params numeric vector: vector of et50 and slope; efficacy calculation is independent of Pmax and thus Pmax is excluded
#' @param TrueEfficacy numeric value: input efficacy value
#' @param titerFun list: list of probability density functions for vaccinated and control groups 
#' @param adjustTiters boolean: set to TRUE if titer values should be adjusted, for details see \code{PoD} function
#' @param adjustFrom numeric: value specifying the detection limit, all values below the detection limit will be adjusted to adjustTo value
#' @param adjustTo numeric: value to which titers below the detection limit will be adjusted
#'
#' @return Squared difference between calculated and reference efficacy
#'
#' @usage
#' efficacySquaredError(params, 
#'                      TrueEfficacy, 
#'                      titerFun,
#'                      adjustTiters = FALSE,
#'                      adjustFrom = 0,
#'                      adjustTo = 0)
#'
#' @details
#' Function is used inside the \code{PoDEfficacySquaredError} function for calculation of the PoD parameters.
#'
#' @examples
#'
#' ## Example 1
#' data(vaccinated)
#' data(control)
#' data(PoDParams)
#'
#' # Choosing et50 and slope as the inputs
#' params <- list("et50" = 4, "slope" = 6)
#'
#' # Using probability density function from the populations
#' 
#' titerFun <- 
#'   list(
#'       function(x) {dnorm(x, mean = vaccinated$mean, sd = vaccinated$stdDev)},
#'       function(x) {dnorm(x, mean = control$mean, sd = control$stdDev)}
#'       )
#' 
#' # Assigning true efficacy
#' TrueEfficacy <- 0.53
#'
#' # Sqaured difference between true and calcuated efficacy
#' efficacySquaredError(params, TrueEfficacy, titerFun)
#'
#' @export
efficacySquaredError <- function(params, 
                                 TrueEfficacy, 
                                 titerFun, 
                                 adjustTiters = FALSE, 
                                 adjustFrom = 0, 
                                 adjustTo = 0) {
  
  et50 <- as.double(params["et50"])
  slope <- as.double(params["slope"])
  
  funVaccinated <- titerFun[[1]]
  funControl <- titerFun[[2]]
  
  funPoD <- function(x) PoD(x, 
                            pmax = 1, 
                            et50 = et50, 
                            slope = slope, 
                            adjustTiters = adjustTiters,
                            adjustFrom = adjustFrom,
                            adjustTo = adjustTo)
  
  effEst <- 1 - ExpectedPoD(f.pod = funPoD, f.titer = funVaccinated) /
    ExpectedPoD(f.pod = funPoD, f.titer = funControl)
  
  effDiff <- (effEst - TrueEfficacy) ^ 2
  
  return(effDiff)
}

#' @title Optimization function: finds PoD curve paramaters (et50, slope) 
#'
#' @description
#' Function finds PoD curve parameters (et50, slope) using population summary statistics (mean, sd) and input (reference value, or for example true in the simulation setup) efficacy. 
#' Efficacy is independent of pmax parameter thus pmax is estimated separately using \code{PmaxEstimation} function.
#'
#' @param TrueEfficacy numeric: input reference efficacy
#' @param vaccinated \code{Population-class} object: vaccinated group (mean, sd)
#' @param control \code{Population-class} object: control group (mean, sd)
#' @param initialSlope numeric: initial slope parameter for the optimization function
#' @param adjustTiters boolean: set to TRUE if titer values should be adjusted, for details see \code{PoD} function
#' @param adjustFrom numeric: value specifying the detection limit, all values below the detection limit will be adjusted to adjustTo value
#' @param adjustTo numeric: value to which titers below the detection limit will be adjusted
#'
#' @return PoD curve parameters (et50, slope)
#'
#' @usage
#' PoDEfficacySquaredError(TrueEfficacy, 
#'                         vaccinated, 
#'                         control,
#'                         initialSlope = 6,
#'                         adjustTiters = FALSE, 
#'                         adjustFrom = NA, 
#'                         adjustTo = NA)
#'
#' @details
#' Function returns et50 and slope PoD curve parameters obtained using \code{efficacySquaredError} 
#' i.e. the opimal (output) parameters et50 and slope correspond to the minimal squared difference between input reference efficacy and calculated efficacy. 
#' 
#' Pmax parameter is not obtained as efficacy is independent on pmax. 
#'
#' The \code{optim} function is used for optimization with method = "L-BFGS-B", 1000 maximum itiretations, (0.1,Inf) boundaries for et50 and (-slopeBoundary, slopeBoundary) boundaries for slope.
#'
#' NOTE: The reason for slope boundary settings is because from certain value of slope parameter the shape of the PoD curve and the corresponding PoD values for given titers are almost identical. 
#' This parameter is supposed to limit the resulting slope value and help MLE to converge to optimal parameters.
#' The value of "slopeBoundaries" is calculated from data according to Dunning, 2015 (https://doi.org/10.1186/s12874-015-0096-9).
#' 
#' @examples
#'
#' ## Example 1
#' data(vaccinated)
#' data(control)
#'
#' # Assigning reference efficacy 
#' TrueEfficacy <- 0.53
#'
#' # PoD curve parameter estimation
#' PoDEfficacySquaredError(TrueEfficacy, vaccinated, control)
#'
#' @export
PoDEfficacySquaredError <- function(TrueEfficacy, 
                                    vaccinated, 
                                    control,
                                    initialSlope = 6,
                                    adjustTiters = FALSE, 
                                    adjustFrom = NA, 
                                    adjustTo = NA) {
  
  if (!is(vaccinated, "Population")) {
    incorrectPopulationInput("vaccinated")
  }
  if (!is(control, "Population")) {
    incorrectPopulationInput("control")
  }
  
  if(vaccinated$mean < control$mean) {
    warning(paste("Mean of vaccinated population is lower than mean of control population"))
  }
  
  coeff <- 9.1902
  range <- (vaccinated$mean + 1.96 * vaccinated$stdDev) - (control$mean - 1.96 * control$stdDev)
  slopeUpperBound <- coeff * 50 / range
  
  initialEt50 <- (vaccinated$mean + control$mean) / 2
  initialParams <- c(initialEt50, initialSlope)
  names(initialParams) <- c("et50", "slope")
  
  tryError <- F
  
  tryCatch(
    params <- optim(par = initialParams,
                    fn = efficacySquaredError,
                    TrueEfficacy = TrueEfficacy,
                    titerFun = list(vaccinated$popFun(), control$popFun()),
                    adjustTiters = adjustTiters,
                    adjustFrom = adjustFrom,
                    adjustTo = adjustTo,
                    upper = c(Inf, slopeUpperBound),
                    lower = c(0.1, -slopeUpperBound),
                    method = "L-BFGS-B",
                    control = list(trace = F, maxit = 1000)),
    error = function(e) {
      tryError <<- T
    }
  )
  
  return(params$par)
}

#' @title PoD curve paramater, pmax,  estimation
#'
#' @description
#' Function finds the pmax parameter of the PoD curve using control subjects summary statistics (mean, sd), observed incidence rate and previsouly estimated et50 and slope by \code{PoDEfficacySquaredError} function. 
#'
#' @param IncidenceRate numeric: observed incidence rate in overall (control) subjects
#' @param params numeric vector: et50 and slope
#' @param control \code{Population-class} object: control subjects (mean, sd)
#' @param adjustTiters boolean: set to TRUE if titer values should be adjusted, for details see \code{PoD} function
#' @param adjustFrom numeric: value specifying the detection limit, all values below the detection limit will be adjusted to adjustTo value
#' @param adjustTo numeric: value to which titers below the detection limit will be adjusted
#'
#' @return PoD curve parameter pmax
#'
#' @usage
#' PmaxEstimation(IncidenceRate,
#'                params, 
#'                control,
#'                adjustTiters = FALSE, 
#'                adjustFrom = NA, 
#'                adjustTo = NA)
#'
#' @examples
#'
#' ## Example 1
#' data(vaccinated)
#' data(control)
#' 
#' # Assigning true efficacy 
#' TrueEfficacy <- 0.53
#'
#' # PoD curve parameters (et50, slope) estimation
#' params <- PoDEfficacySquaredError(TrueEfficacy, vaccinated, control)
#'
#' # Assigning incidence rate (observed incidence rate)
#' IncidenceRate <- 0.2
#' 
#' # pmax estimation
#' pmax <- PmaxEstimation(IncidenceRate, params, control)
#' 
#' # combining PoD curve parameters
#' PoDParams <- unlist(c(params, pmax))
#' 
#' @export
PmaxEstimation <- function(IncidenceRate,
                           params,
                           control, 
                           adjustTiters = FALSE, 
                           adjustFrom = NA, 
                           adjustTo = NA){
  
  if (!is(control, "Population")) {
    incorrectPopulationInput("control")
  }
  
  if(any(is.na(match(names(params), c("et50", "slope") )))){
    stop(paste("The input value for params is incorrect. 'params' parameter has wrong names."))
  }
  
  et50 <- as.double(params["et50"])
  slope <- as.double(params["slope"])
  
  funPoD <-  function(x) PoD(x,
                             pmax = 1,
                             et50 = et50,
                             slope = slope, 
                             adjustTiters = adjustTiters,
                             adjustFrom = adjustFrom,
                             adjustTo = adjustTo)
  
  funControl <- function(x) dnorm(x, mean = control$mean, sd = control$stdDev)
  
  aucControl = ExpectedPoD(f.pod = funPoD,
                                f.titer = funControl)
  
  pmax = IncidenceRate / aucControl
  
  return(list(pmax = pmax))
}