#' afttest: Diagnostics for Semiparametric Accelerated Failure Time Models
#'
#' @description
#' The \pkg{afttest} package provides diagnostic tools for semiparametric
#' accelerated failure time (AFT) models. It supports formal model checking
#' and graphical assessment for fitted AFT models, providing a unified
#' framework for evaluating model adequacy and interpreting results.
#'
#' @details
#' The main interface is the generic function \code{\link{afttest}}, which
#' accepts either a formula or fitted model objects from the \pkg{aftgee}
#' package, including \code{aftsrr} and \code{aftgee}.
#'
#' Implemented diagnostic procedures include:
#' \itemize{
#'   \item \strong{Omnibus test} (\code{testType = "omnibus"}), for assessing overall model adequacy;
#'   \item \strong{Link function test} (\code{testType = "link"}), for evaluating the specified link function;
#'   \item \strong{Functional form test} (\code{testType = "covForm"}), for checking whether a continuous covariate is correctly specified.
#' }
#'
#' The package supports diagnostic inference for:
#' \itemize{
#'   \item rank-based estimators with non-smooth (\code{eqType = "ns"}) and induced-smoothed (\code{eqType = "is"}) estimating equations;
#'   \item least-squares estimators (\code{estMethod = "ls"}) fit through generalized estimating equations.
#' }
#'
#' For computational efficiency, \pkg{afttest} implements a fast multiplier
#' bootstrap procedure. When \code{linApprox = TRUE}, an asymptotic linear
#' approximation is used to avoid repeated iterative model fitting, substantially
#' reducing computation time relative to standard resampling with
#' \code{linApprox = FALSE}.
#'
#' @references
#' Bae, W., Choi, D., Yan, J., & Kang, S. (2025).
#' \emph{afttest: Model diagnostics for semiparametric accelerated failure time models in R}.
#' arXiv preprint arXiv:2511.09823.
#'
#' Choi, D., Bae, W., Yan, J., & Kang, S. (2024).
#' A general model-checking procedure for semiparametric accelerated failure time models.
#' \emph{Statistics and Computing}, 34(3), 117.
#'
#' @name afttest-package
#' @keywords internal
"_PACKAGE"