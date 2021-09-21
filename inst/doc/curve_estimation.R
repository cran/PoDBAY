## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(dplyr)
library(ggplot2)
library(PoDBAY)
set.seed(1234)

## -----------------------------------------------------------------------------
# Mockup vaccinated and control population class objects
data(vaccinated)
data(control)

# Observed vaccine efficacy
TrueEfficacy <- 0.53

# PoD curve parameter estimation
params_et50_slope <- PoDEfficacySquaredError(TrueEfficacy, 
                                             vaccinated, 
                                             control,
                                             initialSlope = 6)
params_et50_slope

## -----------------------------------------------------------------------------
# Incidence rate for low titer population 
IncidenceRate <- 0.02

# pmax estimation
pmax <- PmaxEstimation(IncidenceRate, params_et50_slope, control)
 
# combining PoD curve parameters
PoDParams <- unlist(c(params_et50_slope, pmax))

PoDParams

