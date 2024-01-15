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
  
  cat("\n Call: \n")
  print(x$call)
  
  cat("\n p.value: \n")
  p.valueTAB <- data.frame(t(c(x$p_value, x$p_std_value)))
  rownames(p.valueTAB) <- ""
  colnames(p.valueTAB) <- c("p.value", "std.p.value")
  print(p.valueTAB)
}
# print.afttest <- function(x, ...) {
#   if (!inherits(x,"afttest")) stop("Must be afttest class")
#   
#   print.default(x)
#   invisible(x)
# }

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
  
  cat("\n Call: \n")
  print(object$call)
  
  cat("\n p.value: \n")
  p.valueTAB <- data.frame(t(c(object$p_value, object$p_std_value)))
  rownames(p.valueTAB) <- ""
  colnames(p.valueTAB) <- c("p.value", "std.p.value")
  print(p.valueTAB)
  
  cat("\n Coefficients (estimated by aftgee::aftsrr): \n")
  coefTAB <- data.frame(t(object$beta))
  rownames(coefTAB) <- ""
  colnames(coefTAB) <- object$names[-c(1:2)]
  print(coefTAB)
}
# summary.afttest <- function(object, ...) {
#   if (!inherits(object,"afttest")) stop("Must be afttest class")
#   
#   out <- list(call = object$call,
#              path = object$path,
#              eqType = object$eqType,
#              testType = object$testType,
#              optimType = object$optimType,
#              coefficients = object$beta,
#              p_value = object$p_value,
#              p_std_value = object$p_std_value,
#              missingmessage = object$missingmessage)
# 
#   class(out) <- "summary.afttest"
#   out
# }

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
