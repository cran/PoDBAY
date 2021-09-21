## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(dplyr)
library(ggplot2)
set.seed(1234)

## ---- echo = FALSE------------------------------------------------------------
library(PoDBAY)

## -----------------------------------------------------------------------------
vaccinated <- list()
vaccinated$N = 2000
vaccinated$mean = 8
vaccinated$stdDev = 2

control <- list()
control$N = 1000
control$mean = 5
control$stdDev = 2

PoDParams$pmax = 0.03
PoDParams$et50 = 7
PoDParams$slope = 7
  
#methods: "Full", "Ratio", "Fixed"
  method <- list(name = "Full",
                 value = NA)
  
  # method <- list(name = "Ratio",
  #                value = 4)
  
  # method <- list(name = "Fixed",
  #                value = 300)
  
parameters <- list(vaccinated = vaccinated,
                   control = control,
                   PoDParams = PoDParams,
                   method = method,
                   repeatCount = 50,
                   adjustTiters = FALSE,
                   adjustFrom = NA,
                   adjustTo = NA)

## -----------------------------------------------------------------------------
PoDParams <- parameters$PoDParams

means <- list(vaccinated = vaccinated$mean, 
              control = control$mean)

standardDeviations <- list(vaccinated = vaccinated$stdDev, 
                           control = control$stdDev) 

TrueEfficacy <- efficacyComputation(PoDParams, 
                                    means, 
                                    standardDeviations)
TrueEfficacy

## -----------------------------------------------------------------------------
vaccinated <- generatePopulation(parameters$vaccinated$N,
                                 parameters$vaccinated$mean,
                                 parameters$vaccinated$stdDev)
  
control <- generatePopulation(parameters$control$N,
                                   parameters$control$mean,
                                   parameters$control$stdDev)

str(vaccinated)
str(control)

## -----------------------------------------------------------------------------
vaccinated$assignPoD(
  PoD(titer = vaccinated$popX(),
      pmax = PoDParams$pmax,
      et50 = PoDParams$et50,
      slope = PoDParams$slope,
      adjustTiters = parameters$adjustTiters,
      adjustFrom = parameters$adjustFrom,
      adjustTo = parameters$adjustTo)
)

control$assignPoD(
  PoD(titer = control$popX(),
      pmax = PoDParams$pmax,
      et50 = PoDParams$et50,
      slope = PoDParams$slope,
      adjustTiters = parameters$adjustTiters,
      adjustFrom = parameters$adjustFrom,
      adjustTo = parameters$adjustTo)
)

str(vaccinated)
str(control)

## -----------------------------------------------------------------------------
CaseCountEfficacy <- ClinicalTrialCoverage(vaccinated,
                                           control)

list(CaseCountEfficacy    = CaseCountEfficacy$efficacy,
     confidenceInterval95 = unlist(CaseCountEfficacy$confidenceInterval95),
     confidenceInterval90 = unlist(CaseCountEfficacy$confidenceInterval90),
     confidenceInterval80 = unlist(CaseCountEfficacy$confidenceInterval80))

str(vaccinated)
str(control)

## -----------------------------------------------------------------------------
diseasedAll <- ExtractDiseased(CaseCountEfficacy$vaccinated, 
                               CaseCountEfficacy$control)
  
nondiseasedAll <- ExtractNondiseased(CaseCountEfficacy$vaccinated, 
                                   CaseCountEfficacy$control)

str(diseasedAll)
str(nondiseasedAll)

## -----------------------------------------------------------------------------
ImmunogenicitySample <- BlindSampling(diseasedAll,
                                      nondiseasedAll,
                                      parameters$method)

## -----------------------------------------------------------------------------
str(ImmunogenicitySample$ImmunogenicityNondiseased)  

## -----------------------------------------------------------------------------
str(ImmunogenicitySample$ImmunogenicityVaccinated)

## -----------------------------------------------------------------------------
str(ImmunogenicitySample$ImmunogenicityControl)

## ---- results = "hide", cache = TRUE------------------------------------------
estimatedParameters <- PoDParamEstimation(diseasedAll$titers,
                                            ImmunogenicitySample$ImmunogenicityNondiseased$titers,
                                            nondiseasedAll$N,
                                            parameters$repeatCount,
                                            adjustTiters = parameters$adjustTiters,
                                            adjustFrom = parameters$adjustFrom,
                                            adjustTo = parameters$adjustTo)

## -----------------------------------------------------------------------------
# Confidence intervals
PoDParamsCI <- PoDParamsCICoverage(estimatedParameters$results)
unlist(PoDParamsCI)

# Point estimate
PoDParamsPointEst <- PoDParamPointEstimation(estimatedParameters$resultsPriorReset)
unlist(PoDParamsPointEst)

## -----------------------------------------------------------------------------
# Point estimate
meansBlind <- list("vaccinated" = ImmunogenicitySample$ImmunogenicityVaccinated$mean,
                   "control" = ImmunogenicitySample$ImmunogenicityControl$mean)

standardDeviationsBlind  <- list("vaccinated" = ImmunogenicitySample$ImmunogenicityVaccinated$stdDev,
                                 "control" = ImmunogenicitySample$ImmunogenicityControl$stdDev)
  
EfficacyPointEst <- efficacyComputation(PoDParamsPointEst, meansBlind, standardDeviationsBlind)

# PoDBAY efficacy set
efficacySet <- PoDBAYEfficacy(estimatedParameters$results, 
                             ImmunogenicitySample$ImmunogenicityVaccinated, 
                             ImmunogenicitySample$ImmunogenicityControl,
                             adjustTiters = parameters$adjustTiters,
                             adjustFrom = parameters$adjustFrom,
                             adjustTo = parameters$adjustTo)
  
# PoDBAY efficacy confidence intervals 
CI <- EfficacyCICoverage(efficacySet)

EfficacyPointEst
unlist(CI)

## ---- fig.width = 5-----------------------------------------------------------
result <- list(
  TrueEfficacy = TrueEfficacy,
  CaseCountEfficacy    = CaseCountEfficacy$efficacy,
  confidenceInterval95 = unlist(CaseCountEfficacy$confidenceInterval95),
  confidenceInterval90 = unlist(CaseCountEfficacy$confidenceInterval90),
  confidenceInterval80 = unlist(CaseCountEfficacy$confidenceInterval80),
  EfficacyPointEst = EfficacyPointEst,
  efficacyCI = unlist(CI),
  PoDParamsPointEst = unlist(PoDParamsPointEst),
  PoDParamsCI = unlist(PoDParamsCI))
  
result

## -----------------------------------------------------------------------------
vaccinated <- list()
vaccinated$N = 2000
vaccinated$mean = 8
vaccinated$stdDev = 2

control <- list()
control$N = 1000
control$mean = 5
control$stdDev = 2

PoDParams$pmax = 0.03
PoDParams$et50 = 7
PoDParams$slope = 7
  
#methods: "Full", "Ratio", "Fixed"
  # method <- list(name = "Full",
  #                value = NA)
  
  # method <- list(name = "Ratio",
  #                value = 4)
  
  method <- list(name = "Fixed",
                 value = 300)
  
parameters <- list(vaccinated = vaccinated,
                   control = control,
                   PoDParams = PoDParams,
                   method = method,
                   repeatCount = 50,
                   adjustTiters = FALSE,
                   adjustFrom = NA,
                   adjustTo = NA)

