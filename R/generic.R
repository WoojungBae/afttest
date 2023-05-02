#' ##############################################################################
#' ## Summary
#' ##############################################################################
#' 
#' #' summary.afttest
#' #'
#' #' It gives plot for checking the aft model assumptions.
#' #' @param object is a \code{afttest} fit
#' #' @return \code{summary.afttest} returns a summary of a \code{afttest} fit:
#' #' \describe{
#' #'    \item{Class}{"afttest"}
#' #'    \item{Call}{Object Call}
#' #'    \item{Test Type}{Test Type}
#' #'    \item{Path}{The number of sample path generated}
#' #'    \item{Coeffcient}{beta coeffcient based on aftsrr function}
#' #'    \item{p-value based on standardized test}{p-value based on standardized test}
#' #'    \item{p-value based on unstandardized test}{p-value based on unstandardized test}
#' #' }
#' #' 
#' #' @export
#' summary.afttest = function(object, ...) {
#'   if (!inherits(object,"afttest")) stop("Must be afttest class")
#'   
#'   cat("Class:\n")
#'   print(class(object))
#'   cat("Call:\n")
#'   print(object$call)
#'   cat("test type:\n")
#'   print(object$TestType)
#'   cat("the number of sample path generated:\n")
#'   print(object$path)
#'   cat("beta coeffcient based on aftsrr function:\n")
#'   print(object$beta)  
#'   cat("\n p-value based on unstandardized test:\n")
#'   print(object$p_value)
#'   cat("\n p-value based on standardized test:\n")
#'   print(object$p_std_value)
#' }
#' 
#' #' print.summary.afttest
#' #' 
#' #' @param object is a \code{afttest} fit
#' #' @return \code{print.summary.afttest} returns a summary of a \code{afttest} fit:
#' #'    This function links to \code{summary.afttest}. 
#' #'    See \code{\link[afttest]{summary.afttest}}.
#' #'    
#' #' @export
#' print.summary.afttest = function(object, ...) {
#'   summary.afttest(object)
#' }
#' 
#' #' print.afttest
#' #' 
#' #' @param object is a \code{afttest} fit
#' #' @return \code{print.afttest} returns a summary of a \code{afttest} fit:
#' #'    This function links to \code{summary.afttest}. 
#' #'    See \code{\link[afttest]{summary.afttest}}.
#' #' 
#' #' @export
#' print.afttest = function(object, ...) {
#'   summary.afttest(object)
#' }
#' 
#' ##############################################################################
#' ## Plot
#' ##############################################################################
#' 
#' #' plot.afttest
#' #' 
#' #' @param object is a \code{afttest} fit
#' #' @return \code{plot.afttest} returns a plot of a \code{afttest} fit:
#' #'    This function links to \code{afttestplot}. 
#' #'    See \code{\link[afttest]{afttestplot}}.
#' #' 
#' #' @export
#' plot.afttest = function(object, ...) {
#'   afttestplot(object)
#' }
#' 
#' ##############################################################################
#' ## etc
#' ##############################################################################
#' 
#' 
#' 
#' 
