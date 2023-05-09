##############################################################################
## Summary
##############################################################################

#' summary.afttest
#'
#' @param object is a \code{afttest} fit.
#' @param ... other options.
#' @return \code{summary.afttest} returns a summary of a \code{afttest} fit:
#' 
#' @example inst/examples/ex_afttest.R
#' @export
summary.afttest = function(object, ...) {
  if (!inherits(object,"afttest")) stop("Must be afttest class")
  
  out = list(call = object$call,
             path = object$path,
             eqType = object$eqType,
             testType = object$testType,
             optimType = object$optimType,
             coefficients = object$beta,
             p_value = object$p_value,
             p_std_value = object$p_std_value,
             missingmessage = object$missingmessage)
  
  class(out) = "summary.afttest"
  out
}

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

#' ##############################################################################
#' ## etc
#' ##############################################################################
#' 
#' #' print.afttest
#' #'
#' #' @param object is a \code{afttest} fit
#' #' @param ... other options.
#' #' @return \code{print.afttest} returns a summary of a \code{afttest} fit:
#' #'    This function links to \code{summary.afttest}.
#' #'    See \code{\link[afttest]{summary.afttest}}.
#' #' 
#' #' @example inst/examples/ex_afttest.R
#' #' @export
#' print.afttest = function(object, ...) {
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

#' #' print.summary.afttest
#' #'
#' #' @param object is a \code{afttest} fit.
#' #' @param ... other options.
#' #' @return \code{print.summary.afttest} returns a summary of a \code{afttest} fit:
#' #'    This function links to \code{summary.afttest}.
#' #'    See \code{\link[afttest]{summary.afttest}}.
#' #' 
#' #' @example inst/examples/ex_afttest.R
#' #' @export
#' print.summary.afttest = function(object, ...) {
#'   summary.afttest(object)
#' }