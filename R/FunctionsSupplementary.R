#' @title Numeric to boolean
#'
#' @description  Converts numeric format to boolean format.
#'
#' @details If the function is supposed to be used on a vector, the form \code{sapply("vector", numToBool)} needs to be applied.
#'
#' @param x numeric value (0, 1)
#'
#' @return boolean value (T, F)
#'
#' @usage
#' numToBool(x)
#'
#' @examples
#' dStatus <- c(0,0,1,1,0,1)
#' sapply(dStatus, numToBool)
#'
#' @export
numToBool <- function(x) {if (x != 0) { return(TRUE) } else { return(FALSE) }}

#' @title Error message
#'
#' @description  Error meassage: the input value for "name" is incorrent
#'
#' @param name name of the input value
#'
#' @return error message: "the input value for "name" is incorrect"
#'
#' @export
incorrectInput <- function(name) {
  stop(paste("The input value for", name, "is incorrect."))
}

#' @title Population class error message
#'
#' @description  Error meassage: the input value for "name" is incorrect.
#'
#' @param name name of the input value
#'
#' @return error message: "The input value for "name" is incorrect. Input needs to be a population class object."
#'
#' @export
incorrectPopulationInput <- function(name) {
  stop(paste("Input value for", name," is incorrect. Input needs to be a population class object. For more details see 'population' vignette - vignette(\"population\", package = \"PoDBA\")."))
}

#' @importFrom methods is new
#' @importFrom stats dnorm integrate median optim quantile rbinom rnorm sd uniroot 
#' @importFrom ggplot2 aes element_blank element_line element_rect element_text geom_line ggplot ggtitle theme ylab
NULL
