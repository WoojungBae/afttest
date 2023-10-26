##############################################################################
## print
##############################################################################

#' print.afttest
#'
#' @param x is a \code{afttest} fit.
#' @param ... other options.
#' @return \code{print.afttest} returns a summary of a \code{afttest} fit:
#'
#' @example inst/examples/ex_afttest.R
#' @export
print.afttest <- function(x, ...) {
  if (!inherits(x,"afttest")) stop("Must be afttest class")
  
  print.default(x)
  invisible(x)
}

##############################################################################
## summary
##############################################################################

#' summary.afttest
#'
#' @param object is a \code{afttest} fit.
#' @param ... other options.
#' @return \code{summary.afttest} returns a summary of a \code{afttest} fit:
#'
#' @example inst/examples/ex_afttest.R
#' @export
summary.afttest <- function(object, ...) {
  if (!inherits(object,"afttest")) stop("Must be afttest class")

  out <- list(call = object$call,
             path = object$path,
             eqType = object$eqType,
             testType = object$testType,
             optimType = object$optimType,
             coefficients = object$beta,
             p_value = object$p_value,
             p_std_value = object$p_std_value,
             missingmessage = object$missingmessage)

  class(out) <- "summary.afttest"
  out
}

#' ##############################################################################
#' ## plot
#' ##############################################################################
#' 
#' #' plot.afttest
#' #'
#' #' @param x is a \code{afttest} fit
#' #' @return \code{plot.afttest} returns a plot of a \code{afttest} fit:
#' #'    This function links to \code{afttestplot}.
#' #'    See \code{\link[afttest]{afttestplot}}.
#' #'
#' #' @export
#' plot.afttest <- function(x, ...) {
#'   afttestplot(x)
#' }
