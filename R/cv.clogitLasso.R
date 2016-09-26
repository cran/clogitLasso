#' Plot coefficients from a \code{clogitLasso} object
#' 
#' Plot the parameter profile associated \code{clogitLasso} object
#' 
#' @param objclogitLasso An objet of type \code{clogitLasso}
#' @param K The number of folds used in cross validation
#' @param gpe A list of group defined by the user.
#' @return An object of type \code{cv.clogitLasso} with the following components: 
#' 
#' \item{lambda}{Vector of regularisation parameter}
#' 
#' \item{mean_cv}{vector of mean deviances for each value of the regularisation parameter}
#' 
#' \item{beta}{Vector of estimated coefficients with optimal regularisation parameter}
#' 
#' \item{lambdaopt}{Optimal regularisation parameter}
#' @author Marius Kwemou and Marta Avalos
#' @references Avalos, M., Pouyes, H., Grandvalet, Y., Orriols, L., & Lagarde, E. (2015). \emph{Sparse conditional logistic 
#'  regression for analyzing large-scale matched data from epidemiological studies: a simple algorithm.} BMC bioinformatics, 16(6), 1.
#'  @examples
#'  \dontrun{
#'  #generate data
#'  y = rep(c(1,0), 100)
#'  X = matrix (rnorm(20000, 0, 1), ncol = 100) # pure noise
#'  strata = sort(rep(1:100, 2))
#'  
#'  fitLasso = clogitLasso(X,y,strata,log=TRUE)
#'  # Cross validation 
#'  cv.fit = cv.clogitLasso(fitLasso)
#'  }
#' @export

cv.clogitLasso = function(objclogitLasso,K=10,gpe=NULL){
  
  if(table(objclogitLasso$arg$strata)[1]==2){
    res = crossvalidation(objclogitLasso,K=10,gpe=NULL)
  }else{
    res = crossvalidation1M(objclogitLasso,K=10,gpe=NULL)
  }
  class(res) = "cv.clogitLasso"
  return(res)
}