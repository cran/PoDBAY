## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(dplyr)
library(PoDBAY)
set.seed(1)

## -----------------------------------------------------------------------------
populationEmpty <- generatePopulation()
populationEmpty

## -----------------------------------------------------------------------------
vaccinated <- generatePopulation()
vaccinated$N <- 1000
vaccinated$mean <- 7
vaccinated$stdDev <- 2

str(vaccinated)

control <- generatePopulation()
control$N <- 1000
control$mean <- 5
control$stdDev <- 2

str(control)

## -----------------------------------------------------------------------------
vaccinated <- generatePopulation(N = 1000,
                                 mean = 7,
                                 stdDev = 2)
str(vaccinated)

control <- generatePopulation(N = 1000,
                                   mean = 5,
                                   stdDev = 2)

str(control)

## -----------------------------------------------------------------------------
vaccinated <- generatePopulation()
vaccinated$N <- 1000
vaccinated$mean <- 7
vaccinated$stdDev <- 2

str(vaccinated)

control <- generatePopulation()
control$N <- 1000
control$mean <- 5
control$stdDev <- 2

str(control)

## ---- results = "hide"--------------------------------------------------------
vaccinated$getTiters()
control$getTiters()

## -----------------------------------------------------------------------------
str(vaccinated)
str(control)

## -----------------------------------------------------------------------------
# Define PoD curve parameters
PoDParams <- data.frame("pmax" = 0.05, "et50" = 5, "slope" = 7)

# Assign PoD
vaccinated$assignPoD(
  PoD(titer = vaccinated$titers,
      pmax = PoDParams$pmax,
      et50 = PoDParams$et50,
      slope = PoDParams$slope,
      adjustTiters = FALSE
  ))

control$assignPoD(
  PoD(titer = control$titers,
      pmax = PoDParams$pmax,
      et50 = PoDParams$et50,
      slope = PoDParams$slope,
      adjustTiters = FALSE
  ))

str(vaccinated)
str(control)

## ---- eval = FALSE------------------------------------------------------------
#  vaccinated$assignPoD(
#    rep(0.05, vaccinated$N))
#  
#  control$assignPoD(
#    rep(0.1, control$N))

## -----------------------------------------------------------------------------
CaseCount <- ClinicalTrial(vaccinated, control, CI = 0.95)

str(vaccinated)
str(control)

## ---- warning = FALSE---------------------------------------------------------
vaccinated$getDiseasedCount()
vaccinated$getNondiseasedCount()
vaccinated$getDiseasedTiters()
as_tibble(vaccinated$getNondiseasedTiters())

## -----------------------------------------------------------------------------
dataTrial <- data.frame("patno" = 1:20,
                        "treatment" = c(rep(FALSE, 10), rep(TRUE, 10)),
                        "titers" = as.numeric(c(rnorm(10, 5, 2),
                                                rnorm(10, 7, 2)))
                        )

dataTrial

## -----------------------------------------------------------------------------
# vaccinated
vacc <- dataTrial %>% filter(treatment)

vaccinated <- generatePopulation()
vaccinated$N <- nrow(vacc)
vaccinated$mean <- mean(vacc$titers)
vaccinated$stdDev <- sd(vacc$titers)
vaccinated$titers <- vacc$titers

vaccinated

# control
ctrl <- dataTrial %>% filter(!treatment)

control <- generatePopulation()
control$N <- nrow(ctrl)
control$mean <- mean(ctrl$titers)
control$stdDev <- sd(ctrl$titers)
control$titers <- ctrl$titers


control


