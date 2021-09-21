## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(dplyr)
library(ggplot2)
set.seed(1234)

## ----TrialAData---------------------------------------------------------------
library(PoDBAY)
data(diseased)
data(nondiseased)
str(diseased)
str(nondiseased)

## ----TrialAPoDEst, results = "hide", cache = TRUE-----------------------------
estimatedParameters <- PoDParamEstimation(diseasedTiters = diseased$titers,
                                          nondiseasedTiters = nondiseased$titers, 
                                          nondiseasedGenerationCount = nondiseased$N,
                                          repeatCount = 50)

## ---- echo = FALSE------------------------------------------------------------
as_tibble(estimatedParameters$resultsPriorReset)

## ---- echo = FALSE------------------------------------------------------------
as_tibble(estimatedParameters$results)

## -----------------------------------------------------------------------------
PoDParamsPointEst <- PoDParamPointEstimation(estimatedParameters$resultsPriorReset)
PoDParamsPointEst

## -----------------------------------------------------------------------------
PoDParametersCI <- PoDParamsCI(estimatedParameters$results, ci = 0.95)
unlist(PoDParametersCI)

## ---- fig.width = 5-----------------------------------------------------------
titers <- seq(from = 0, to = 15, by = 0.01)

PoDCurve <- PoDCurvePlot(titers,
                         estimatedParameters,
                         ci = 0.95)

PoDCurve

## ----TrialBData---------------------------------------------------------------
data(vaccinated)
data(control)
str(vaccinated)
str(control)

## ----EfficacyPointEstimate----------------------------------------------------
means <- list("vaccinated" = vaccinated$mean,
              "control" = control$mean)
  
standardDeviations  <- list("vaccinated" = vaccinated$stdDev,
                            "control" = control$stdDev)
  
EfficacyPointEst <- efficacyComputation(PoDParamsPointEst, 
                                        means, 
                                        standardDeviations)
EfficacyPointEst

## ----EfficacySet--------------------------------------------------------------
efficacySet <- PoDBAYEfficacy(estimatedParameters$results, 
                              vaccinated, 
                              control)

## ----EfficacyCI---------------------------------------------------------------
CI <- EfficacyCI(efficacySet, ci = 0.95)
unlist(CI)

## ---- fig.width = 5-----------------------------------------------------------
result <- list(
  EfficacyPointEst = EfficacyPointEst,
  efficacyCI = unlist(CI),
  PoDParamsPointEst = PoDParamsPointEst,
  PoDParametersCI = unlist(PoDParametersCI),
  PoDCurve = PoDCurve
)
  
result

## -----------------------------------------------------------------------------
data(diseased)
data(nondiseased)

# Immunogenicity sample created
ImmunogenicitySample <- BlindSampling(diseased, nondiseased, method = list(name = "Fixed", value = 200))
nondiseasedImmunogenicitySample <- ImmunogenicitySample$ImmunogenicityNondiseased

str(diseased)
str(nondiseasedImmunogenicitySample)

## ----ApxPoDEst, results = "hide", cache = TRUE--------------------------------
estimatedParametersAP <- PoDParamEstimation(diseasedTiters = diseased$titers,
                                            nondiseasedTiters = nondiseasedImmunogenicitySample$titers, 
                                            nondiseasedGenerationCount = nondiseased$N,
                                            repeatCount = 50)

## ---- echo = FALSE------------------------------------------------------------
as_tibble(estimatedParametersAP$resultsPriorReset)

## ---- echo = FALSE------------------------------------------------------------
as_tibble(estimatedParametersAP$results)

## -----------------------------------------------------------------------------
PoDParamsPointEst <- PoDParamPointEstimation(estimatedParametersAP$resultsPriorReset)
PoDParamsPointEst

## -----------------------------------------------------------------------------
PoDParametersCI <- PoDParamsCI(estimatedParametersAP$results)
unlist(PoDParametersCI)

## ---- fig.width = 5-----------------------------------------------------------
titers <- seq(from = 0, to = 15, by = 0.01)

PoDCurve <- PoDCurvePlot(titers,
                         estimatedParametersAP,
                         ci = 0.95)

PoDCurve

## ----AppendixTrialBData-------------------------------------------------------
# Immunogenicity sample - vaccinated
str(ImmunogenicitySample$ImmunogenicityVaccinated)

# Immunogenicity sample - control
str(ImmunogenicitySample$ImmunogenicityControl)

## ----AppendixEfficacyPointEstimate--------------------------------------------
means <- list("vaccinated" = ImmunogenicitySample$ImmunogenicityVaccinated$mean,
                   "control" = ImmunogenicitySample$ImmunogenicityControl$mean)
  
standardDeviations  <- list("vaccinated" = ImmunogenicitySample$ImmunogenicityVaccinated$stdDev,
                                 "control" = ImmunogenicitySample$ImmunogenicityControl$stdDev)
  
EfficacyPointEst <- efficacyComputation(PoDParamsPointEst, 
                                        means, 
                                        standardDeviations)
EfficacyPointEst

## ----AppendixEfficacySet------------------------------------------------------
efficacySet <- PoDBAYEfficacy(estimatedParametersAP$results, 
                              ImmunogenicitySample$ImmunogenicityVaccinated, 
                              ImmunogenicitySample$ImmunogenicityControl)

## ----AppendixEfficacyCI-------------------------------------------------------
CI <- EfficacyCICoverage(efficacySet)
unlist(CI)

## ---- fig.width = 5-----------------------------------------------------------
result <- list(
  EfficacyPointEst = EfficacyPointEst,
  efficacyCI = unlist(CI),
  PoDParamsPointEst = PoDParamsPointEst,
  PoDParametersCI = unlist(PoDParametersCI),
  PoDCurve = PoDCurve
)
  
result

