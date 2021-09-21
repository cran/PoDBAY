#' @title Immunogenicity subset: vaccinated, control, non-diseased
#'
#' @description
#' Function creates non-diseased immunogenicity subset, and vaccinated and control immunogenicity subsets based on chosen method. The immunogenicity subsets are provided in the form of population class objects (see the \code{Population-class} function for more details).
#'
#' @param diseased \code{Population-class} object: diseased subjects, created using \code{ExtractDiseased} function
#' @param nondiseased \code{Population-class} object: non-diseased subjects, created using \code{ExtractNondiseased} function
#' @param method named list: "name" possible inputs "Full", "Ratio", "Fixed";
#'
#' "value" = numeric value
#'
#' @return
#' \itemize{
#'   \item ImmunogenicityVaccinated: vaccinated subjects in the immunogenicity subset, \code{Population-class} object (N, mean, stdDev, titers)
#'
#'   \item ImmunogenicityControl: control subjects in the immunogenicity subset, \code{Population-class} object (N, mean, stdDev, titers)
#'
#'   \item ImmunogenicityNondiseased: non-diseased subjects in the immunogenicity subset, \code{Population-class} object (N, mean, stdDev, titers)
#' }
#'
#' @usage
#' BlindSampling(diseased, 
#'               nondiseased,  
#'               method = list(name = "Full", value = NA))
#'
#' @details
#'
#' For details about the method parameter see \code{ImmunogenicitySubset} function.
#'
#' @examples
#' # Data preparation
#' data(diseased)
#' data(nondiseased)
#'
#' ## Example 1
#' # Creating immunogenicity subset, method = "Full"
#' ImmunogenicitySubsetFull <- 
#'     BlindSampling(diseased, 
#'                   nondiseased, 
#'                   method = list(name = "Full", 
#'                                 value = NA))
#'
#' ## Example 2
#' # Creating of immunogenicity subset, method = "Ratio"
#' ImmunogenicitySubsetRatio <- 
#'     BlindSampling(diseased, 
#'                   nondiseased, 
#'                   method = list(name = "Ratio", 
#'                                 value = 4))
#'
#' ## Example 3
#' # Creating of immunogenicity subset, method = "Fixed"
#' ImmunogenicitySubsetFixed <- 
#'     BlindSampling(diseased, 
#'                   nondiseased, 
#'                   method = list(name = "Fixed", 
#'                                 value = 100))
#'
#' @export
BlindSampling <- function(diseased, nondiseased, method = list(name = "Full",
                                                              value = NA)) {

  if (!is(diseased, "Population")) {
    incorrectPopulationInput("diseased")
  }
  if (!is(nondiseased, "Population")) {
    incorrectPopulationInput("nondiseased")
  }
  
  # create non-diseased immunogenicity sample
  ImmunogenicitySample <- ImmunogenicitySubset(diseased = diseased,
                                               nondiseased = nondiseased,
                                               method = method)

  # vaccinated and control groups in the non-diseased immunogenicty sample
  ImmunogenicityVaccinatedTiters <- ImmunogenicitySample$titers[!grepl("control", names(ImmunogenicitySample$titers))]
  ImmunogenicityControlTiters <- ImmunogenicitySample$titers[grepl("control", names(ImmunogenicitySample$titers))]
 
  ImmunogenicityVaccinated          <- generatePopulation(0)
  ImmunogenicityVaccinated$N        <- length(ImmunogenicityVaccinatedTiters)
  ImmunogenicityVaccinated$stdDev   <- sd(ImmunogenicityVaccinatedTiters)
  ImmunogenicityVaccinated$mean     <- mean(ImmunogenicityVaccinatedTiters)
  ImmunogenicityVaccinated$titers   <- ImmunogenicityVaccinatedTiters
  ImmunogenicityVaccinated$diseaseStatus <- as.logical(substr(names(ImmunogenicityVaccinatedTiters), 
                                                              6, 
                                                              length(names(ImmunogenicityVaccinatedTiters))))

  ImmunogenicityControl        <- generatePopulation(0)
  ImmunogenicityControl$mean   <- mean(ImmunogenicityControlTiters)
  ImmunogenicityControl$stdDev <- sd(ImmunogenicityControlTiters)
  ImmunogenicityControl$N      <- length(ImmunogenicityControlTiters)
  ImmunogenicityControl$titers <- ImmunogenicityControlTiters
  ImmunogenicityControl$diseaseStatus <- as.logical(substr(names(ImmunogenicityControlTiters), 
                                                                9, 
                                                                length(names(ImmunogenicityControlTiters))))

  ImmunogenicityNondiseased <- ExtractNondiseased(ImmunogenicityVaccinated, ImmunogenicityControl)
  
  return(list(
    ImmunogenicityNondiseased   = ImmunogenicityNondiseased,
    ImmunogenicityVaccinated   = ImmunogenicityVaccinated,
    ImmunogenicityControl = ImmunogenicityControl
  ))
}

#' @title Immunogenicity subset
#'
#' @description
#' Function creates the immunogenicity subset based on the chosen method.
#'
#' @param diseased \code{Population-class} object: diseased subjects with assigned vaccination status
#' @param nondiseased \code{Population-class} object: non-diseased subjects with assigned vacination status
#' @param method named list: a selected method for creating the immunogenicity subset
#'
#' method$name
#'
#' \itemize{
#'   \item Full: subject level titer information is available for all diseased and all non-diseased subjects, i.e. immunogenicity subset is the full clinical trial
#'
#'   \item Ratio: subject level titer information is available for all diseased and some non-diseased subjects.
#'
#'   \item Fixed: subject level titer information is available for all diseased and some non-diseased subjects.
#' }
#'
#'
#' method$value
#'
#' \itemize{
#'   \item Full: value = NA; immunogenicity sample is the full clinical trial (non-diseased subset contains all non-diseased in the trial; diseased subset contains all disease cases in the trial)
#'
#'   \item Ratio: value = number of non-diseased divided by number of diseased subjects; ratio of diseased vs. non-diseased subjects in the immunogenicity subset (non-diseased subset contains only non-diseased subjects, as the selection is done in the end of the study, when the disease status is known; diseased subset contains all disease cases in the trial)
#'
#'   \item Fixed: value = size of the immunogenicity subset, pre-defined number of subjects assayed for titers independently of their future disease status (non-diseased subset could rarely contain some diseased subjects, as the selection is done at the enrollment and prior the knowledge of future disease status; diseased subset contains all disease cases in the trial)
#' }
#'
#' @return
#' Immunogenicity subset with subject level information about vaccination status and disease status, provided in the form of \code{Population-class} object
#'
#' @usage
#' ImmunogenicitySubset(diseased, 
#'                      nondiseased, 
#'                      method = list(name = "Full", value = NA))
#'
#' @details
#' The total immunogenicity subset consists of the diseased immunogenicity subset and non-diseased immunogenicity subset. 
#' For all three methods implemented, we assume that the diseased immunogenicity subset contains all disease cases in the trial.
#' Based on the chosen method, the the size of the non-diseaded immunogenicity subset can be derived as follows:
#'
#' Size = number of subjects in the non-diseased immunogenicity subset
#'
#' Titers = values of titers from which we want to sample in order to simulate the non-diseased immunogenicity subset
#'
#' #Diseased = total number of diseased in the clinical trial
#'
#' #Nondiseased = total number of non-diseased in the clinical trial
#'
#'  \itemize{
#'   \item method$name = "Full"
#'   
#'   Size = #Nondiseased
#'   
#'   Titers = Nondiseased Titers 
#'
#'   \item method$name = "Ratio"
#'   
#'   Size =  method$value * #Diseased
#'   
#'   Titers = Nondiseased Titers
#'   
#'   \item method$name = "Fixed"
#'   
#'   Size = method$value
#'   
#'   Titers = Nondiseased Titers + Diseased Titers
#'   
#' }
#'
#' @examples
#' ## Example 1
#' # Data preparation
#' data(diseased)
#' data(nondiseased)
#'
#' ImmunogenicitySubset(diseased,
#'                      nondiseased,
#'                      method = list(name = "Ratio",
#'                                    value = 4))
#'
#' @export
ImmunogenicitySubset <- function(diseased,
                                 nondiseased,
                                 method = list(name = "Full",
                                               value = NA)) {
  if (!is(diseased, "Population")) {
    incorrectPopulationInput("diseased")
  }
  if (!is(nondiseased, "Population")) {
    incorrectPopulationInput("nondiseased")
  }
  
  if(any(is.na(match(names(method), c("name", "value") )))){
    stop(paste("The input value for method is incorrect. 'method' parameter has wrong names."))
  }
    
  if (is.na(
    match(method$name, c("Full", "Ratio", "Fixed")))
  ) { stop(paste("The input value for method$name is incorrect. 'method$name' needs to be either \"Full\", \"Ratio\" or \"Fixed\".")) }

  if (method$name == "Full") {
    ForImmunogenicitySample <- generatePopulation(0)
    ForImmunogenicitySample$N <- nondiseased$N + diseased$N
    ForImmunogenicitySample$titers <- c(nondiseased$titers, diseased$titers)
    ForImmunogenicitySample$diseaseStatus <- c(nondiseased$diseaseStatus, diseased$diseaseStatus)

  }  else if (method$name == "Ratio") {

    ImmunogenicityRatio <- method$value
    if (!is.numeric(ImmunogenicityRatio)| ImmunogenicityRatio <= 0 | is.na(ImmunogenicityRatio)) { 
      stop(paste("The input value for method$value is incorrect. 'method$value' needs to be positive numeric value."))
    }
    if(nondiseased$N < (ImmunogenicityRatio * diseased$N)) {
      warning(paste("nondiseased$N is less than method$value * diseased$N"))
    }

    ImmunogenicitySubsetSize <- ifelse(nondiseased$N >= (ImmunogenicityRatio * diseased$N),
                                       ImmunogenicityRatio * diseased$N,
                                       nondiseased$N)
    ForImmunogenicitySample <- generatePopulation(0)
    ForImmunogenicitySample$N <- ImmunogenicitySubsetSize
    ForImmunogenicitySample$titers <- nondiseased$titers
    ForImmunogenicitySample$diseaseStatus <- nondiseased$diseaseStatus

  } else if (method$name == "Fixed") {

    ImmunogenicitySize <- method$value
    if (!is.numeric(ImmunogenicitySize) | ImmunogenicitySize <= 0 | is.na(ImmunogenicitySize)) { 
      stop(paste("The input value for method$value is incorrect. 'method$value' needs to be positive numeric value."))
    }

    ForImmunogenicitySample <- generatePopulation(0)
    ForImmunogenicitySample$N <- ImmunogenicitySize
    ForImmunogenicitySample$titers <- c(nondiseased$titers, diseased$titers)
    ForImmunogenicitySample$diseaseStatus <- c(nondiseased$diseaseStatus, diseased$diseaseStatus)
  }

  names(ForImmunogenicitySample$titers) <- paste(names(ForImmunogenicitySample$titers), ForImmunogenicitySample$diseaseStatus, sep = "_")
  
  ImmunogenicitySampleTiters  <- sample(x = sample(ForImmunogenicitySample$titers),
                                        round(ForImmunogenicitySample$N),
                                        replace = FALSE)

  ImmunogenicitySample        <- generatePopulation(0)
  ImmunogenicitySample$N      <- length(ImmunogenicitySampleTiters)
  ImmunogenicitySample$mean   <- mean(ImmunogenicitySampleTiters)
  ImmunogenicitySample$stdDev <- sd(ImmunogenicitySampleTiters)
  ImmunogenicitySample$titers <- ImmunogenicitySampleTiters

  return(ImmunogenicitySample)
}

#' @title Diseased subjects extraction
#'
#' @description
#' Function extracts diseased subjects from vaccinated and control groups if the data have assigned disease status (for example using \code{ClinicalTrial} function). The vaccinated and control data are provided in the form of population class objects (see the \code{Population-class} function for more details).
#'
#' @param vaccinated \code{Population-class} object: vaccinated subjects with assigned disease status
#' @param control \code{Population-class} object: control subjects with assigned disease status
#'
#' @return
#' diseased subjects, \code{Population-class} object: a subset of control and vaccinated subjects with disease status = TRUE.
#'
#' @usage
#' ExtractDiseased(vaccinated, control)
#'
#' @examples
#' ## Example 1
#' # Data preparation
#' data(vaccinated)
#' data(control)
#'
#' # Estimating the disease status and case-count efficacy with CI
#' ClinicalTrial(vaccinated, control, CI = 0.95)
#'
#' # Extracting the disease cases
#' ExtractDiseased(vaccinated, control)
#'
#' @export
ExtractDiseased <- function(vaccinated, control) {
  if (!is(vaccinated, "Population")) {
    incorrectPopulationInput("vaccinated")
  }
  if (!is(control, "Population")) {
    incorrectPopulationInput("control")
  }

  diseased <- generatePopulation(0)

  titers <- c(vaccinated$getDiseasedTiters(), control$getDiseasedTiters())
  diseased$titers <- titers

  names(diseased$titers) <- c(rep("vacc", vaccinated$getDiseasedCount()),
                              rep("control", control$getDiseasedCount()))

  diseased$N <- length(titers)
  diseased$mean <- mean(titers)
  diseased$stdDev <- sd(titers)
  diseased$diseaseStatus <- rep(TRUE, length(titers))
  return(diseased)
}

#' @title Non-diseased subjects extraction
#'
#' @description
#' Function extracts non-diseased subjects from vaccinated and control groups if the data have assigned disease status (for example using \code{ClinicalTrial} function). The vaccinated and control data are provided in the form of population class objects (see the \code{Population-class} function for more details).
#'
#' @param vaccinated \code{Population-class} object: vaccinated subjects with assigned disease status
#' @param control \code{Population-class} object: control subjects with assigned disease status
#'
#' @return
#' non-diseased subjects, \code{Population-class} object: a subset of control and vaccinated subjects with disease status = FALSE.
#'
#' @usage
#' ExtractNondiseased(vaccinated, control)
#'
#' @examples
#' ## Example 1
#' # Data preparation
#' data(vaccinated)
#' data(control)
#'
#' # Estimating the disease status and case-count efficacy with CI
#' ClinicalTrial(vaccinated, control, CI = 0.95)
#'
#' # Extracting the non-diseased subjects
#' ExtractNondiseased(vaccinated, control)
#'
#' @export
ExtractNondiseased <- function(vaccinated, control) {
  if (!is(vaccinated, "Population")) {
    incorrectPopulationInput("vaccinated")
  }
  if (!is(control, "Population")) {
    incorrectPopulationInput("control")
  }
  
  nondiseased <- generatePopulation(0)

  titers <- c(vaccinated$getNondiseasedTiters(), control$getNondiseasedTiters())
  nondiseased$titers <- titers

  names(nondiseased$titers) <- c(rep("vacc", vaccinated$getNondiseasedCount()),
                                rep("control", control$getNondiseasedCount()))

  nondiseased$N <- length(titers)
  nondiseased$mean <- mean(titers)
  nondiseased$stdDev <- sd(titers)
  nondiseased$diseaseStatus <- rep(FALSE, length(titers))

  return(nondiseased)
}