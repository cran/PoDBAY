#' @title Population class object generation
#'
#' @description Function generates the population class object using provided summary statistics.
#'
#' @param N numeric: number of subjects in the population
#' @param mean numeric: mean of titers
#' @param stdDev numeric: standard deviation of titers 
#' @param unknownDistribution logical: TRUE if there is an unknown factor affacting the shape of titer distribution
#' @param UDFunction function: function defining the unknown factor affecting the shape of titer distribution
#'
#' @return
#' generated population class object with all its characteristics defined in the input parameters
#'
#' @usage
#' generatePopulation(N, mean, stdDev, unknownDistribution = FALSE, UDFunction = NULL)
#'
#' @examples
#'
#' # Example 1: empty population
#' population0 <- generatePopulation()
#'
#' # Example 2
#' population1 <- generatePopulation(N = 100,
#'                                   mean = 5,
#'                                   stdDev = 2)
#'
#' @export
generatePopulation <- function(N = 0, 
                               mean = NA_real_, 
                               stdDev = NA_real_, 
                               unknownDistribution = FALSE, 
                               UDFunction = NULL) {
  if (is.na(mean))          {mean <- NA_real_ }
  if (is.na(stdDev))        {stdDev <- NA_real_ }

  generated <- population$new()
  generated$N <- N
  generated$mean <- mean
  generated$stdDev <- stdDev
  generated$getTiters()
  generated$unknownDistribution <- unknownDistribution
  if (unknownDistribution) { generated$UDFunction <- UDFunction}

  return(generated)
}

#' @title Probability of disease calculation
#'
#' @description Function calculates probability of disease (PoD) corresponding to given titers according to a sigmoid PoD curve.
#'
#' @param titer numeric vector: subject level titers
#' @param pmax numeric: maximum PoD
#' @param et50 numeric: titer values corresponding to pmax/2 value, PoD(et50) = pmax/2
#' @param slope numeric: slope of the PoD curve
#' @param adjustTiters boolean: set to TRUE if titer values should be adjusted, for details see \code{PoD} function
#' @param adjustFrom numeric: value specifying the detection limit, all values below the detection limit will be adjusted to adjustTo value
#' @param adjustTo numeric: value to which titers below the detection limit will be adjusted
#'
#' @return vector of PoDs
#'
#' @usage
#' PoD(titer, pmax, et50, slope, adjustTiters = FALSE, adjustFrom = 0, adjustTo = 0)
#'
#' @examples
#' data(vaccinated)
#' data(PoDParams)
#'
#' PoD(vaccinated$titers, pmax = PoDParams$pmax, et50 = PoDParams$et50, slope = PoDParams$slope)
#'
#' @details PoD is calculated as: \deqn{ PoD = p_{max} \frac{ (\frac{et50}{titer})^{\gamma} }{ 1 + (\frac{et50}{titer})^{\gamma}}, \ for \ titers \ > 0}{ PoD = pmax * (et50/titer)^(slope) / (1+ (et50/titer)^(slope), for titers > 0} and 
#'\deqn{ PoD = pmax, \ for \ titers \ <= 0}{PoD = pmax for titers <= 0}.
#'
#' @export
PoD <- function(titer, pmax, et50, slope, adjustTiters = FALSE, adjustFrom = 0, adjustTo = 0) {
  
  if (adjustTiters & (adjustFrom < adjustTo) ) {warning(paste("The input value for \"adjustFrom\" is lower than \"adjustTo\" "))}
  if (adjustTiters) {titer[titer < adjustFrom] <- adjustTo}
  
  probDisease <- ifelse( titer > 0, pmax - pmax / ( 1 + ( et50 / titer ) ^ slope), pmax)
  return(probDisease)
}

# AssignPoD - see RefClassPopulation
