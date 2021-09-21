#' @title Population class
#'
#' @description Population reference class which provides summary and subject level information about the population
#'
#' @field N numeric: number of subjects in the population
#' @field mean numeric: mean value of titers
#' @field stdDev numeric: standard deviation of titers 
#' @field unknownDistribution logical: TRUE if titer distribution is not normally /log-normally distributed; titer disrtibution function needs to be defined by user
#' @field UDFunction function: user-defined titer distribution
#' @field titers numeric: subject level titers, generated with \code{getTiters} method
#' @field PoDs numeric: subject level probability of disease (PoD), generated with \code{assginPoD} method
#' @field diseaseStatus logical: subject level disease status (TRUE if diseased), generated with \code{ClinicilaTrial} function
#'
#' @exportClass Population
population <- setRefClass(
  "Population",
  fields = list(
    N = "numeric",
    mean = "numeric",
    stdDev = "numeric",
    unknownDistribution = "logical",
    UDFunction = "function", # function with one parameter that defines the length of generated sequence
    titers = "numeric", # should be generated automatically
    PoDs = "numeric", # probability of disease for each subject
    diseaseStatus = "logical"
  ),
  methods = list(
    initialize = function() {
      mean <<- NA_real_
      stdDev <<- NA_real_
      unknownDistribution <<- FALSE
    },
    getTiters = function() {
      if (length(titers) == 0) {
        titers <<- rnorm(N, mean, stdDev)
      }
      return(titers)
    },
    popFun = function() {
      if (unknownDistribution) {
        toRet <- approxfun(
          density(rnorm(1e6, mean = mean, sd = stdDev) + getUnknown(1e6)),
          yleft = 0,
          yright = 0)
        return(toRet)
      }
      else {
        toRet <- function(x) {
          dnorm(x, mean = mean, sd = stdDev)
        }
        return(toRet)
      }
    },
    popX = function() {
      if (unknownDistribution) {
        return(titers + getUnknown(N))
      } else {
        return(titers)
      }
    },
    getUnknown = function(n) {
      return(UDFunction(n))
    },
    assignPoD = function(x){
      PoDs <<- x
    },
    getDiseasedCount = function() {
      return(length(which(diseaseStatus)))
    },
    getNondiseasedCount = function() {
      return(length(which(!diseaseStatus)))
    },
    getDiseasedTiters = function() {
      return(titers[which(diseaseStatus)])
    },
    getNondiseasedTiters = function() {
      return(titers[which(!diseaseStatus)])
    }
  )
)

#' @title Subject level titers
#'
#' @name getTiters
#'
#' @description Returns subject level titers. If titers are not yet generated, the function generates them based on \code{Population-class} object attributes: N, mean, stdDev.
#'
#' @details Inputs into the function (N, mean, stdDev) are taken from the \code{Population-class} object attributes.
#'
#' @return Subject level titers
NULL
population$methods(
  getTiters = function() {
    if (length(titers) == 0) {
      titers <<- rnorm(N, mean, stdDev)
    }
    return(titers)
  }
)

#' @title Population function
#'
#' @name popFun
#'
#' @description Function describing the titer distribution of the population: mean, standard deviation and an additional unknown factor affecting the shape of the distribution (e.g. mixture of two normals or other shapes defined by user).
#'
#' @details Inputs into the function (mean, stdDev, Unknowndistribution) and getUnknown method are taken from the \code{Population-class} object.
#'
#' @return Titer distribution function
NULL
population$methods(
  popFun = function() {
    if (unknownDistribution) {
      toRet <- approxfun(
        density(rnorm(1e6, mean = mean, sd = stdDev) + getUnknown(1e6)),
        yleft = 0,
        yright = 0)
      return(toRet)
    }
    else {
      toRet <- function(x) {
        dnorm(x, mean = mean, sd = stdDev)
      }
      return(toRet)
    }
  }
)

#' @title Add noise to population titers
#'
#' @name popX
#'
#' @description Function adds noise to population titers accounting for an unknown factor affecting the titer distibution.
#'
#' @details Inputs into the function: N, unknownDistribution and getUnknown() method are taken from the \code{Population-class} object.
#'
#' @return subject level titers
NULL
population$methods(
  popX = function() {
    if (unknownDistribution) {
      return(titers + getUnknown(N))
    } else {
      return(titers)
    }
  }
)

#' @title Generate unknown
#'
#' @name getUnknown
#'
#' @description Function generates unknown part of the titers which is eventually added to the original titers in \code{popX} and to the original titer distribution in \code{popFun}.
#'
#' @param n numeric: number of subjects in the population
#'
#' @details Input into the function: UDFunction is taken from the \code{Population-class} object. UDFunction is used for generating the unknown part of the titer distribution.
#'
#' @return unknown part of the titers
NULL
population$methods(
  getUnknown = function(n) {
    return(UDFunction(n))
  }
)

#' @title Assign probability of disease (PoD)
#'
#' @name assignPoD
#'
#' @description Function assigns subject-level probability of disease based on PoD curve and subject level titer.
#'
#' @param x numeric vector - vector of estimated PoD values
#'
#' @details The input into the function is either calculated using \code{PoD} function or if the PoD curve is unknown the same arbitrary PoD can be assigned to the whole population.
#'
#' @return Subject level probability of disease for the population
NULL
population$methods(
  assignPoD = function(x){
    PoDs <<- x
  }
)

#' @title Diseased count
#'
#' @name getDiseasedCount
#'
#' @description Function calculates the number of diseased subjects (disease status = TRUE) in the \code{Population-class} object.
#'
#' @details Input into the function, "diseaseStatus", is taken from the \code{Population-class} object attribute. Information about disease status is written into the \code{Population-class} object by the \code{ClinicalTrial()} function.
#'
#' @return numeric: number of the diseased subjects in the \code{Population-class} object
NULL
getDiseasedCount <- function() {
  return(length(which(diseaseStatus)))
}

#' @title Non-diseased count
#'
#' @name getNondiseasedCount
#'
#' @description Function calculates the number of non-diseased subjects (disease status = FALSE) in the \code{Population-class} object.
#'
#' @details Input into the function, "diseaseStatus", is taken from the \code{Population-class} object attribute. Information about disease status is written into the \code{Population-class} object by the \code{ClinicalTrial()} function.
#'
#' @return numeric: number of the non-diseased subjects in the \code{Population-class} object
NULL
getNondiseasedCount <- function() {
  return(length(which(!diseaseStatus)))
}

#' @title Diseased titers
#'
#' @name getDiseasedTiters
#'
#' @description Function returns titers of diseased subjects (disease status = TRUE) in the \code{Population-class} object.
#'
#' @details Input into the function, "diseaseStatus", is taken from the \code{Population-class} object attribute. Information about disease status is written into the \code{Population-class} object by the \code{ClinicalTrial()} function.
#'
#' @return numeric vector: titers of diseased subjects in the \code{Population-class} object
NULL
getDiseasedTiters <- function() {
  return(titers[which(diseaseStatus)])
}

#' @title Non-diseased titers
#'
#' @name getNondiseasedTiters
#'
#' @description Function returns titers of non-diseased subjects (disease status = FALSE) in the \code{Population-class} object.
#'
#' @details Input into the function, "diseaseStatus", is taken from the \code{Population-class} object attribute. Information about disease status is written into the \code{Population-class} object by the \code{ClinicalTrial()} function.
#'
#' @return numeric vector: titers of non-diseased subjects in the \code{Population-class} object
NULL
getNondiseasedTiters <- function() {
  return(titers[which(!diseaseStatus)])
}