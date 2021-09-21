// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;
//' 
//' @title Probability of disease calculation
//'
//' @description Function calculates probability of disease (PoD) for given titers according to a PoD curve.
//'
//' @param titer numeric vector: vector of subject level titers
//' @param pmax numeric: maximum PoD
//' @param et50 numeric: titer value corresponding to pmax/2 value, PoD(et50) = pmax/2
//' @param slope numeric: slope of the PoD curve
//' @param adjustTiters boolean: set to TRUE if titer values should be adjusted, for details see \code{PoD} function
//' @param adjustFrom numeric: value specifying the detection limit, all values below the detection limit will be adjusted to adjustTo value
//' @param adjustTo numeric: value to which titers below the detection limit will be adjusted
//'
//' @return vector of PoDs
//'
//' @usage
//' cppPoD(titer, pmax, et50, slope, adjustTiters = FALSE, adjustFrom = 0, adjustTo = 0)
//' 
//' @details 
//' See \code{PoD} function for more details. These two functions are equivalent. Usage of cppPoD significantly improves the computation speed over the \code{PoD} function.
//' 
//' @export
// [[Rcpp::export]]
NumericVector cppPoD(NumericVector titer, double pmax, double et50, double slope, bool adjustTiters = false, double adjustFrom = 0, double adjustTo = 0) {
  int titerSize = titer.size();
  double titerVal;
  double probDisease;
  NumericVector toReturn = NumericVector(titerSize);
  
  for (int i = 0; i != titerSize; i++) {
      titerVal = titer[i];
      if (adjustTiters && titerVal < adjustFrom) {
        titerVal = adjustTo;
      }
      if (titerVal <= 0) {
        probDisease = pmax;
      } else {
        probDisease = pmax - pmax / (1 + pow((et50/titerVal), slope));
      }
      toReturn[i] = probDisease;
    }
    
    return  toReturn;
    
  }
