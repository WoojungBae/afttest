##############################################################################
## User's Main Function
##############################################################################

.afttest_design_formula <- function(covnames) {
  terms <- paste0("`", gsub("`", "\\\\`", covnames, fixed = TRUE), "`")
  stats::as.formula(paste("survival::Surv(time, delta) ~",
                          paste(terms, collapse = " + ")))
}

.afttest_validate_args <- function(npath, testType, covTested, covnames,
                                   npathsave, linApprox, seed) {
  if (!is.numeric(npath) || length(npath) != 1 || is.na(npath) ||
      !is.finite(npath) || npath < 50 || npath != as.integer(npath)) {
    stop("npath must be a single integer of at least 50.")
  }
  npath <- as.integer(npath)

  if (!is.character(testType) || length(testType) != 1 || is.na(testType) ||
      !testType %in% c("omnibus", "link", "covForm")) {
    stop("testType must be one of 'omnibus', 'link', or 'covForm'.")
  }

  if (!is.numeric(npathsave) || length(npathsave) != 1 || is.na(npathsave) ||
      !is.finite(npathsave) || npathsave < 0 ||
      npathsave != as.integer(npathsave)) {
    stop("npathsave must be a single nonnegative integer.")
  }
  npathsave <- as.integer(npathsave)
  if (npathsave > npath) {
    warning("npathsave exceeds npath; only npath paths will be saved.")
    npathsave <- npath
  }

  if (!is.logical(linApprox) || length(linApprox) != 1 || is.na(linApprox)) {
    stop("linApprox must be a single logical value.")
  }

  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1 ||
                         is.na(seed) || !is.finite(seed) ||
                         seed != as.integer(seed))) {
    stop("seed must be a single finite integer.")
  }

  covTested.num <- 1L
  if (testType == "covForm") {
    if (is.numeric(covTested) && length(covTested) == 1 &&
        !is.na(covTested) && is.finite(covTested) &&
        covTested == as.integer(covTested) && covTested >= 1 &&
        covTested <= length(covnames)) {
      covTested.num <- as.integer(covTested)
    } else if (is.character(covTested) && length(covTested) == 1 &&
               !is.na(covTested) && covTested %in% covnames) {
      covTested.num <- match(covTested, covnames)
    } else {
      stop("covTested must identify one design-matrix column for covForm.")
    }
  }

  list(npath = npath, testType = testType, npathsave = npathsave,
       linApprox = linApprox, covTested.num = covTested.num)
}

#' Model Diagnostics for Semiparametric AFT Models
#'
#' @description
#' Performs model-checking procedures for a semiparametric AFT model.
#' This is a generic function with methods for formulas and fitted objects
#' from the \pkg{aftgee} package.
#'
#' @param object A formula or a fitted model object (e.g., from \code{aftsrr} or 
#'   \code{aftgee}).
#' @param ... Other arguments passed to methods. See the documentation for
#'   \code{afttest.formula}, \code{afttest.aftsrr}, and \code{afttest.aftgee} 
#'   for details.
#'
#' @return An object of class \code{afttest} or \code{htest}.
#'   An object is a list containing at least the following components:
#' \describe{
#'   \item{beta}{a vector of beta estimates based on \code{estMethod}}
#'   \item{hypothesis}{null hypothesis for each \code{testType}}
#'   \item{SE_process}{estimated standard error of the observed process}
#'   \item{obs_process}{observed process}
#'   \item{apprx_process}{approximated process}
#'   \item{obs_std_process}{standardized observed process}
#'   \item{apprx_std_process}{standardized approximated processes}
#'   \item{p_value}{obtained by the unstandardized test}
#'   \item{p_std_value}{obtained by the standardized test}
#'   \item{DF}{a data frame of observed failure time, right censoring indicator, 
#'   covariates (scaled), time-transformed residual based on beta estimates}
#'   \item{npath}{the number of sample paths}
#'   \item{testType}{testType}
#'   \item{eqType}{eqType}
#'   \item{estMethod}{estMethod}
#'   \item{npathsave}{npathsave}
#' }
#' 
#'   For an omnibus test, the observed process and the realizations are composed 
#'   of the n by n matrix where rows represent the t and columns represent the x 
#'   in the time-transformed residual order. The observed process and the 
#'   simulated processesfor checking a functional form and a link function are 
#'   given by the n by 1 vectorwhich is a function of x in the time-transformed 
#'   residual order. 
#' 
#' @importFrom stats optim get_all_vars as.formula model.matrix model.frame
#' @importFrom aftgee aftsrr aftgee
#' @importFrom survival Surv
#'
#' @example inst/examples/ex_afttest.R
#' @export
afttest <- function(object, ...) {
  UseMethod("afttest")
}

#' Model Diagnostics for AFT Models using Formulas
#'
#' @param object A formula expression, of the form \code{response ~ predictors}.
#'   The \code{response} is a \code{Surv} object with right censoring.
#'   See the documentation of \code{lm}, \code{coxph} and \code{formula} for details.
#' @param data An optional data frame in which to interpret the variables occurring 
#'   in the formula.
#' @param npath An integer value specifying the number of approximated processes.
#'   The default is given by 200.
#' @param testType A character string specifying the type of the test.
#'   The following are permitted:
#'   \describe{
#'     \item{\code{omnibus}}{an omnibus test}
#'     \item{\code{link}}{a link function test}
#'     \item{\code{covForm}}{a functional form of a covariate}
#' }
#' @param estMethod A character string specifying the type of the estimator used.
#'   The readers are referred to the \pkg{aftgee} package for details.
#'   The following are permitted:
#'   \describe{
#'     \item{\code{ls}}{Least-Squares Approach for Accelerated Failure Time 
#'     with Generalized Estimating Equation}
#'     \item{\code{rr}}{Accelerated Failure Time with Smooth Rank Regression}
#' }
#' @param eqType A character string specifying the type of the 
#'   estimating equation used to obtain the regression parameters.
#'   The readers are referred to the \pkg{aftgee} package for details.
#'   The following are permitted:
#'   \describe{
#'     \item{\code{ns}}{Regression parameters are estimated by directly solving 
#'     the nonsmooth estimating equations.}
#'     \item{\code{is}}{Regression parameters are estimated by directly solving 
#'     the induced-smoothing estimating equations.}
#' }
#' @param covTested A character string specifying the covariate which will be tested.
#'   The argument \code{covTested} is necessary only if \code{testType} is 
#'   \code{covForm}. The default option for \code{covTested} is given by "1", which 
#'   represents the first covariate in the formula argument.
#' @param npathsave An integer value specifying the number of paths saved among all the paths.
#'   The default is given by 50. Note that it requires a lot of memory if saving all
#'   sampled paths (N by N matrix for each npath, and so npath*N*N elements).
#' @param linApprox A logical value. If \code{TRUE}, the multiplier bootstrap is 
#'   computed using the asymptotic linear approximation, which is significantly 
#'   faster. If \code{FALSE}, the estimating equations are solved numerically for 
#'   each bootstrap replication. Defaults to \code{TRUE}.
#' @param seed An optional integer specifying the random seed for reproducibility.
#' @param ... Other arguments passed to methods.
#'  
#' @export
afttest.formula <- function(object, data, npath = 200, testType = "omnibus", 
                            estMethod = "rr", eqType = "ns", 
                            covTested = 1, npathsave = 50, linApprox = TRUE,
                            seed = NULL, ...) {
  eqType_supplied <- !missing(eqType)
  
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1 || is.na(seed) ||
        !is.finite(seed) || seed != as.integer(seed)) {
      stop("seed must be a single finite integer.")
    }
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
    } else {
      on.exit(rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    }
    set.seed(seed)
  }
  
  scall <- match.call()
  mf <- model.frame(object, data)
  Y <- mf[[1]]
  mt <- attr(mf, "terms")
  X <- model.matrix(mt, mf)
  if (sum(colnames(X) == "(Intercept)") > 0) {
    X <- X[, -which(colnames(X) == "(Intercept)"), drop = FALSE]
  }
  if (ncol(X) == 0) {
    stop("No covariates found in the formula. Intercept-only models are not supported.")
  }
  covnames <- colnames(X)
  cov.length <- length(covnames)
  args <- .afttest_validate_args(npath, testType, covTested, covnames,
                                 npathsave, linApprox, seed)
  npath <- args$npath
  testType <- args$testType
  npathsave <- args$npathsave
  linApprox <- args$linApprox
  covTested.num <- args$covTested.num
  DF <- data.frame(unclass(Y), X)
  colnames(DF) <- c("time", "delta", covnames)
  
  # check&delete NA, -Inf, Inf, ...
  missingmessage <- NA
  DF[is.infinite(as.matrix(DF))] <- NA
  whichNA_DF <- which(apply(is.na(DF), 1, sum) > 0)
  nNA_DF <- length(whichNA_DF)
  if (nNA_DF > 0) {
    missingmessage <- paste0("(", nNA_DF, " observations deleted due to missingness out of ", nrow(DF), ")")
    DF <- DF[-whichNA_DF, ]
  } else {
    missingmessage <- paste0("(No missingness observed)")
  }
  
  if (any(DF$time <= 0)) {
    return(warning("time must be positive number"))
  }
  
  # beta coefficients from aftsrr function (aftgee package) - with original covariates
  formula <- .afttest_design_formula(covnames)
  if (length(estMethod) != 1 || !estMethod %in% c("ls", "rr")) {
    stop("estMethod must be either 'ls' or 'rr'.")
  }
  if (estMethod == "ls") {
    if (eqType_supplied && !identical(eqType, "ls")) {
      warning("eqType is ignored when estMethod = 'ls'; eqType = 'ls' is used.")
    }
    eqType <- "ls"
    beta <- - aftgee::aftgee(formula, data = DF)$coef.res[-1]
  } else {
    if (length(eqType) != 1 || !eqType %in% c("ns", "is")) {
      stop("eqType must be either 'ns' or 'is' when estMethod = 'rr'.")
    }
    beta <- - aftgee::aftsrr(formula, data = DF, eqType = eqType, rankWeights = "gehan")$beta
  }
  
  # Covariate Scaling
  time <- DF$time
  delta <- DF$delta
  covariates <- scale(as.matrix(DF[, -(1:2), drop = FALSE]))
  DF <- data.frame(time = time, delta = delta, covariates)
  
  # npath
  if (length(npath) > 1){
    return(warning("npath needs to be an integer."))
  } else {
    if (!is.numeric(npath)) {
      npath <- 200
    } else {
      npath <- max(npath,50)
    }
  }
  
  # testType
  if (length(testType) > 1){
    return(warning("testType needs to be one of 'omnibus', 'link', or 'covForm'"))
  } else {
    if (!testType %in% c("omnibus","link","covForm")) {
      testType <- "omnibus"
    }
  }
  
  # npathsave
  if (!is.numeric(npathsave) || length(npathsave) != 1) {
    warning("'npathsave' must be a single numeric integer. Defaulting to npathsave = 50.")
    npathsave <- 50L
  } else {
    npathsave <- as.integer(npathsave)
  }
  
  # linApprox
  if (length(linApprox) > 1) {
    warning("linApprox needs to be a single logical value (TRUE or FALSE). Using default (TRUE).")
    linApprox <- TRUE
  } else {
    if (!is.logical(linApprox)) {
      warning("linApprox needs to be logical (TRUE or FALSE). Using default (TRUE).")
      linApprox <- TRUE
    }
  }
  
  # covTested
  # beta coefficients from aftsrr function (aftgee package) - with scaled covariates
  formula <- .afttest_design_formula(covnames)
  if (estMethod == "ls") {
    b <- - aftgee::aftgee(formula, data = DF)$coef.res[-1]
  } else if (estMethod == "rr") {
    b <- - aftgee::aftsrr(formula, data = DF, eqType = eqType, rankWeights = "gehan")$beta
  } else {
    return(warning("estMethod needs to be one of 'ls' and 'rr'"))
  }
  
  # This function contains the core logic (the C++ calls)
  out <- .afttest_worker(b, time, delta, covariates, npath, testType,
                         eqType, covTested.num, npathsave, linApprox)
  out$call <- scall
  out$beta <- beta
  # out$DF <- data
  out$DF <- DF
  out$seed <- seed
  out$estMethod <- estMethod # It's an aftsrr object
  out$missingmessage <- missingmessage
  if (testType == "covForm") {
    out$covTested <- covTested
  }
  
  return(out)
}

#' Model Diagnostics for Smooth Rank Regression (aftsrr) Objects
#'
#' @param object A fitted model object of class \code{aftsrr} from the \pkg{aftgee} package.
#' @param data An optional data frame in which to interpret the variables occurring 
#'   in the formula.
#' @param npath An integer value specifying the number of approximated processes.
#'   The default is given by 200.
#' @param testType A character string specifying the type of the test.
#'   The following are permitted:
#'   \describe{
#'     \item{\code{omnibus}}{an omnibus test}
#'     \item{\code{link}}{a link function test}
#'     \item{\code{covForm}}{a functional form of a covariate}
#' }
#' @param eqType An optional character string specifying the type of the
#'   estimating equation. For a fitted \code{aftsrr} object, the estimating
#'   equation used in the fitted object is retained. If it cannot be identified
#'   from the fitted object, \code{"ns"} is used. The permitted values are
#'   \code{"ns"} and \code{"is"}. If a different value from that used in the
#'   fitted object is supplied, it is ignored with a warning.
#' @param covTested A character string specifying the covariate which will be tested.
#'   The argument \code{covTested} is necessary only if \code{testType} is 
#'   \code{covForm}. The default option for \code{covTested} is given by "1", which 
#'   represents the first covariate in the formula argument.
#' @param npathsave An integer value specifying the number of paths saved among all the paths.
#'   The default is given by 50. Note that it requires a lot of memory if saving all
#'   sampled paths (N by N matrix for each npath, and so npath*N*N elements).
#' @param linApprox A logical value. If \code{TRUE}, the multiplier bootstrap is 
#'   computed using the asymptotic linear approximation, which is significantly 
#'   faster. If \code{FALSE}, the estimating equations are solved numerically for 
#'   each bootstrap replication. Defaults to \code{TRUE}.
#' @param seed An optional integer specifying the random seed for reproducibility.
#' @param ... Other arguments passed to methods. 
#' 
#' @export
afttest.aftsrr <- function(object, data, npath = 200, testType = "omnibus", eqType = NULL, 
                           covTested = 1, npathsave = 50, linApprox = TRUE,
                           seed = NULL, ...) {
  dots <- list(...)
  
  if ("estMethod" %in% names(dots)) {
    warning("estMethod = '", dots$estMethod, "' is ignored for an aftsrr object; ", "estMethod = 'rr' is used.")
  }
  
  if (is.null(object$call$eqType)) {
    fitted_eqType <- "ns"
  } else {
    fitted_eqType <- as.character(object$call$eqType)
  }
  
  if (length(fitted_eqType) != 1 || !fitted_eqType %in% c("ns", "is")) {
    stop("The fitted aftsrr object must use eqType = 'ns' or 'is'.")
  }
  
  if (!is.null(eqType)) {
    if (length(eqType) != 1 || !eqType %in% c("ns", "is")) {
      stop("eqType must be either 'ns' or 'is'.")
    }
    if (eqType != fitted_eqType) {
      warning("eqType = '", eqType, "' is ignored for this aftsrr object; eqType = '", fitted_eqType, "' used in the fitted object is used.")
    }
  }
  
  eqType <- fitted_eqType
  
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1 || is.na(seed) ||
        !is.finite(seed) || seed != as.integer(seed)) {
      stop("seed must be a single finite integer.")
    }
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
    } else {
      on.exit(rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    }
    set.seed(seed)
  }
  
  scall <- match.call()
  mf <- model.frame(object, data)
  Y <- mf[[1]]
  mt <- attr(mf, "terms")
  X <- model.matrix(mt, mf)
  if (sum(colnames(X) == "(Intercept)") > 0) {
    X <- X[, -which(colnames(X) == "(Intercept)"), drop = FALSE]
  }
  if (ncol(X) == 0) {
    stop("No covariates found in the formula. Intercept-only models are not supported.")
  }
  covnames <- colnames(X)
  cov.length <- length(covnames)
  args <- .afttest_validate_args(npath, testType, covTested, covnames,
                                 npathsave, linApprox, seed)
  npath <- args$npath
  testType <- args$testType
  npathsave <- args$npathsave
  linApprox <- args$linApprox
  covTested.num <- args$covTested.num
  DF <- data.frame(unclass(Y), X)
  colnames(DF) <- c("time", "delta", covnames)
  
  # check&delete NA, -Inf, Inf, ...
  missingmessage <- NA
  DF[is.infinite(as.matrix(DF))] <- NA
  whichNA_DF <- which(apply(is.na(DF), 1, sum) > 0)
  nNA_DF <- length(whichNA_DF)
  if (nNA_DF > 0) {
    missingmessage <- paste0("(", nNA_DF, " observations deleted due to missingness out of ", nrow(DF), ")")
    DF <- DF[-whichNA_DF, ]
  } else {
    missingmessage <- paste0("(No missingness observed)")
  }
  
  if (any(DF$time <= 0)) {
    return(warning("time must be positive number"))
  }
  
  # Covariate Scaling
  time <- DF$time
  delta <- DF$delta
  covariates <- scale(as.matrix(DF[, -(1:2), drop = FALSE]))
  DF <- data.frame(time = time, delta = delta, covariates)
  
  # npath
  if (length(npath) > 1){
    return(warning("npath needs to be an integer."))
  } else {
    if (!is.numeric(npath)) {
      npath <- 200
    } else {
      npath <- max(npath,50)
    }
  }
  
  # testType
  if (length(testType) > 1){
    return(warning("testType needs to be one of 'omnibus', 'link', or 'covForm'"))
  } else {
    if (!testType %in% c("omnibus","link","covForm")) {
      testType <- "omnibus"
    }
  }
  
  # npathsave
  if (!is.numeric(npathsave) || length(npathsave) != 1) {
    warning("'npathsave' must be a single numeric integer. Defaulting to npathsave = 50.")
    npathsave <- 50L
  } else {
    npathsave <- as.integer(npathsave)
  }
  
  # linApprox
  if (length(linApprox) > 1) {
    warning("linApprox needs to be a single logical value (TRUE or FALSE). Using default (TRUE).")
    linApprox <- TRUE
  } else {
    if (!is.logical(linApprox)) {
      warning("linApprox needs to be logical (TRUE or FALSE). Using default (TRUE).")
      linApprox <- TRUE
    }
  }
  
  # covTested
  # beta coefficients from aftsrr function (aftgee package)
  formula <- .afttest_design_formula(covnames)
  b <- - aftgee::aftsrr(formula, data = DF, eqType = eqType, rankWeights = "gehan")$beta
  
  # This function contains the core logic (the C++ calls)
  out <- .afttest_worker(b, time, delta, covariates, npath, testType,
                         eqType, covTested.num, npathsave, linApprox)
  out$beta <- - object$beta
  out$call <- scall
  # out$DF <- data
  out$DF <- DF
  out$seed <- seed
  out$estMethod <- "rr"
  out$missingmessage <- missingmessage
  if (testType == "covForm") {
    out$covTested <- covTested
  }
  
  return(out)
}

#' Model Diagnostics for Generalized Estimating Equation (aftgee) Objects
#'
#' @param object A fitted model object of class \code{aftgee} from the \pkg{aftgee} package.
#' @param data An optional data frame in which to interpret the variables occurring 
#'   in the formula.
#' @param npath An integer value specifying the number of approximated processes.
#'   The default is given by 200.
#' @param testType A character string specifying the type of the test.
#'   The following are permitted:
#'   \describe{
#'     \item{\code{omnibus}}{an omnibus test}
#'     \item{\code{link}}{a link function test}
#'     \item{\code{covForm}}{a functional form of a covariate}
#' }
#' @param eqType The estimating-equation type used for the diagnostic
#'   procedure. For a fitted \code{aftgee} object, this is fixed to
#'   \code{"ls"}. Any other supplied value is ignored with a warning.
#' @param covTested A character string specifying the covariate which will be tested.
#'   The argument \code{covTested} is necessary only if \code{testType} is 
#'   \code{covForm}. The default option for \code{covTested} is given by "1", which 
#'   represents the first covariate in the formula argument.
#' @param npathsave An integer value specifying the number of paths saved among all the paths.
#'   The default is given by 50. Note that it requires a lot of memory if saving all
#'   sampled paths (N by N matrix for each npath, and so npath*N*N elements).
#' @param linApprox A logical value. If \code{TRUE}, the multiplier bootstrap is 
#'   computed using the asymptotic linear approximation, which is significantly 
#'   faster. If \code{FALSE}, the estimating equations are solved numerically for 
#'   each bootstrap replication. Defaults to \code{TRUE}.
#' @param seed An optional integer specifying the random seed for reproducibility.
#' @param ... Other arguments passed to methods. 
#' 
#' @export
afttest.aftgee <- function(object, data, npath = 200, testType = "omnibus", eqType = "ls", 
                           covTested = 1, npathsave = 50, linApprox = TRUE,
                           seed = NULL, ...) {
  dots <- list(...)
  eqType_supplied <- !missing(eqType)
  
  if ("estMethod" %in% names(dots)) {
    warning("estMethod = '", dots$estMethod, "' is ignored for an aftgee object; ", "estMethod = 'ls' is used.")
  }
  
  if (eqType_supplied && !identical(eqType, "ls")) {
    warning("eqType is ignored for an aftgee object; eqType = 'ls' is used.")
  }
  
  eqType <- "ls"
  
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1 || is.na(seed) ||
        !is.finite(seed) || seed != as.integer(seed)) {
      stop("seed must be a single finite integer.")
    }
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
    } else {
      on.exit(rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    }
    set.seed(seed)
  }
  
  scall <- match.call()
  mf <- model.frame(object, data)
  Y <- mf[[1]]
  mt <- attr(mf, "terms")
  X <- model.matrix(mt, mf)
  if (sum(colnames(X) == "(Intercept)") > 0) {
    X <- X[, -which(colnames(X) == "(Intercept)"), drop = FALSE]
  }
  if (ncol(X) == 0) {
    stop("No covariates found in the formula. Intercept-only models are not supported.")
  }
  covnames <- colnames(X)
  cov.length <- length(covnames)
  args <- .afttest_validate_args(npath, testType, covTested, covnames,
                                 npathsave, linApprox, seed)
  npath <- args$npath
  testType <- args$testType
  npathsave <- args$npathsave
  linApprox <- args$linApprox
  covTested.num <- args$covTested.num
  DF <- data.frame(unclass(Y), X)
  colnames(DF) <- c("time", "delta", covnames)
  
  # check&delete NA, -Inf, Inf, ...
  missingmessage <- NA
  DF[is.infinite(as.matrix(DF))] <- NA
  whichNA_DF <- which(apply(is.na(DF), 1, sum) > 0)
  nNA_DF <- length(whichNA_DF)
  if (nNA_DF > 0) {
    missingmessage <- paste0("(", nNA_DF, " observations deleted due to missingness out of ", nrow(DF), ")")
    DF <- DF[-whichNA_DF, ]
  } else {
    missingmessage <- paste0("(No missingness observed)")
  }
  
  if (any(DF$time <= 0)) {
    return(warning("time must be positive number"))
  }
  
  # Covariate Scaling
  time <- DF$time
  delta <- DF$delta
  covariates <- scale(as.matrix(DF[, -(1:2), drop = FALSE]))
  DF <- data.frame(time = time, delta = delta, covariates)
  
  # estMethod
  estMethod = "ls"
  
  # npath
  if (length(npath) > 1){
    return(warning("npath needs to be an integer."))
  } else {
    if (!is.numeric(npath)) {
      npath <- 200
    } else {
      npath <- max(npath,50)
    }
  }
  
  # testType
  if (length(testType) > 1){
    return(warning("testType needs to be one of 'omnibus', 'link', or 'covForm'"))
  } else {
    if (!testType %in% c("omnibus","link","covForm")) {
      testType <- "omnibus"
    }
  }
  
  # npathsave
  if (!is.numeric(npathsave) || length(npathsave) != 1) {
    warning("'npathsave' must be a single numeric integer. Defaulting to npathsave = 50.")
    npathsave <- 50L
  } else {
    npathsave <- as.integer(npathsave)
  }
  
  # linApprox
  if (length(linApprox) > 1) {
    warning("linApprox needs to be a single logical value (TRUE or FALSE). Using default (TRUE).")
    linApprox <- TRUE
  } else {
    if (!is.logical(linApprox)) {
      warning("linApprox needs to be logical (TRUE or FALSE). Using default (TRUE).")
      linApprox <- TRUE
    }
  }
  
  # covTested
  # beta coefficients from aftsrr function (aftgee package)
  formula <- .afttest_design_formula(covnames)
  b <- - aftgee::aftgee(formula, data = DF)$coef.res[-1]
  
  # This function contains the core logic (the C++ calls)
  out <- .afttest_worker(b, time, delta, covariates, npath, testType,
                         eqType, covTested.num, npathsave, linApprox)
  out$beta <- - object$coef.res[-1]
  out$call <- scall
  # out$DF <- data
  out$DF <- DF
  out$seed <- seed
  out$estMethod <- "ls"
  out$missingmessage <- missingmessage
  if (testType == "covForm") {
    out$covTested <- covTested
  }
  
  return(out)
}

#' Internal worker function for afttest
#' @noRd
.afttest_worker <- function(b, time, delta, covariates, npath, testType,
                            eqType, covTested, npathsave, linApprox) {
  
  if (linApprox) {
    sigma_est <- diag(ncol(covariates)) 
    omega_res <- getOmega(beta = b, 
                          Y = time, 
                          X = covariates, 
                          delta = delta, 
                          weights = NULL,
                          gw = NULL,
                          eqType = eqType, 
                          sigma = sigma_est, 
                          B = 500)
    Omega <- omega_res$Omega
    invOmega <- omega_res$invOmega
  } else {
    Omega <-  matrix(NA)
    invOmega <-  matrix(NA)
  }
  
  # C++ functions
  if (testType == "omnibus") {
    out <- .Call("_afttest_omni_cpp", npath, b, time, delta, covariates, 
                 npathsave, eqType, linApprox, invOmega)
  } else if (testType == "link") {
    out <- .Call("_afttest_link_cpp", npath, b, time, delta, covariates, 
                 npathsave, eqType, linApprox, invOmega)
  } else if (testType == "covForm") {
    out <- .Call("_afttest_form_cpp", npath, b, time, delta, covariates, covTested, 
                 npathsave, eqType, linApprox, invOmega)
  }
  
  class(out) <- c("afttest", "htest")
  out$betascaled <- b
  out$npath <- npath
  out$eqType <- eqType
  out$testType <- testType
  out$npathsave <- npathsave
  out$linApprox <- linApprox
  out$Omega <- Omega
  out$invOmega <- invOmega
  
  return(out)
}

#' Internal worker function for linApprox = TRUE
#' @noRd
#' @importFrom stats rexp rnorm var
getOmega <- function(beta, Y, X, delta, weights = NULL, gw = NULL,
                     eqType = "is", sigma = diag(ncol(X)), B = 1e3) {
  
  X <- as.matrix(X); p <- ncol(X); n <- nrow(X)
  if (is.null(weights)) weights <- rep(1, n)
  if (is.null(gw)) gw <- rep(1, n)
  
  viEmp <- function(beta, Y, delta, X, id, weights = rep(1, nrow(X)), B = 1e3,
                    mb = TRUE, zbeta = FALSE, smooth = TRUE,
                    rankWeights = "gehan", gw = NULL,
                    sigma = diag(ncol(X))) {
    
    X <- as.matrix(X); p <- ncol(X); n <- nrow(X)
    if (is.null(gw)) gw <- rep(1, n)
    
    UnV <- matrix(0, ncol = B, nrow = p)
    zmat <- matrix(0, ncol = B, nrow = p)
    
    for (i in 1:B) {
      if (mb) Z <- rexp(n) else Z <- rep(1, n)
      
      if (zbeta) {
        zb <- rnorm(p)
        # Perturbation scale is n^-0.5
        newbeta <- beta + (n^(-0.5)) * zb
        zmat[, i] <- zb
      } else {
        newbeta <- beta
      }
      
      total_weights <- gw * Z
      
      if (smooth) {
        score_vec <- .Call("_afttest_score_gehan_is_cpp", newbeta, Y, X, delta, sigma, total_weights)
      } else {
        score_vec <- .Call("_afttest_score_gehan_ns_cpp", newbeta, Y, X, delta, total_weights)
      }
      UnV[, i] <- score_vec
    }
    vi <- var(t(UnV))
    return(list(vi = vi, zmat = zmat, UnV = UnV))
  }
  
  if (eqType == "is") {
    # Method 1: Induced Smoothing
    Omega <- .Call("_afttest_abar_gehan_cpp", beta, Y, X, delta, sigma, weights, gw) * n^{2}
  } else if (eqType == "ns") {
    # Method 2: Non-Smooth Resampling
    resamp <- viEmp(beta, Y, delta, X, id = 1:n, weights = weights, 
                    B = B, mb = FALSE, zbeta = TRUE, smooth = FALSE, 
                    rankWeights = "gehan", gw = gw)
    UnV <- resamp$UnV
    zmat <- resamp$zmat
    ZZt <- tcrossprod(zmat)        
    UZt <- tcrossprod(UnV, zmat)
    
    Omega <- UZt %*% .Call("_afttest_inv_cpp", ZZt) * n^{2}
  } else if (eqType == "ls") {
    X_weighted <- X * sqrt(weights)
    Omega <- - crossprod(X_weighted) * n^{2}
  } else {
    stop("Invalid eqType")
  }
  
  invOmega <- .Call("_afttest_inv_cpp", Omega)
  
  return(list(Omega = Omega, invOmega = invOmega))
}
