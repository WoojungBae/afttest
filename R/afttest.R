##############################################################################
## User's Main Function
##############################################################################

#' Model Diagnostics for Semiparametric AFT Models
#'
#' @description
#' Performs model-checking procedures for a semiparametric AFT model.
#' This is a generic function with methods for formulas and fitted objects
#' from the `aftgee` package.
#'
#' @param object A formula or a fitted model object (e.g., from `aftsrr` or `aftgee`).
#' @param ... Other arguments passed to methods. See the documentation for
#'   `afttest.formula` and `afttest.aftsrr` for details.
#'
#' @return An object of class \code{afttest} or \code{htest}.
#'    An object is a list containing at least the following components:
#' \describe{
#'    \item{beta}{a vector of beta estimates based on \code{estMethod}}
#'    \item{hypothesis}{null hypothesis for each \code{testType}}
#'    \item{SE_process}{estimated standard error of the observed process}
#'    \item{obs_process}{observed process}
#'    \item{apprx_process}{approximated process}
#'    \item{obs_std_process}{standardized observed process}
#'    \item{apprx_std_process}{standardized approximated processes}
#'    \item{p_value}{obtained by the unstandardized test}
#'    \item{p_std_value}{obtained by the standardized test}
#'    \item{DF}{a data frame of observed failure time, right censoring indicator, covariates (scaled), 
#'    time-transformed residual based on beta estimates}
#'    \item{npath}{the number of sample paths}
#'    \item{testType}{testType}
#'    \item{eqType}{eqType}
#'    \item{estMethod}{estMethod}
#'    \item{npathsave}{npathsave}
#' }
#' 
#'    For an omnibus test, the observed process and the realizations are composed of the 
#'    n by n matrix that rows represent the t and columns represent the x in the 
#'    time-transformed residual order.The observed process and the simulated processes
#'    for checking a functional form and a link function are given by the n by 1 vector
#'    which is a function of x in the time-transformed residual order. 
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

#' @param object A formula expression, of the form \code{response ~ predictors}.
#'    The \code{response} is a \code{Surv} object with right censoring.
#'    See the documentation of \code{lm}, \code{coxph} and \code{formula} for details.
#' @param data An optional data frame in which to interpret the variables occurring 
#'    in the formula.
#' @param npath An integer value specifies the number of approximated processes.
#'    The default is given by 200.
#' @param testType A character string specifying the type of the test.
#'    The following are permitted:
#'    \describe{
#'      \item{\code{omnibus}}{an omnibus test}
#'      \item{\code{link}}{a link function test}
#'      \item{\code{covForm}}{a functional form of a covariate}
#' }
#' @param estMethod A character string specifying the type of the estimator used.
#'    The readers are refered to the \pkg{aftgee} package for details.
#'    The following are permitted:
#'    \describe{
#'      \item{\code{ls}}{Least-Squares Approach for Accelerated Failure Time 
#'      with Generalized Estimating Equation}
#'      \item{\code{rr}}{Accelerated Failure Time with Smooth Rank Regression}
#' }
#' @param eqType A character string specifying the type of the 
#'    estimating equation used to obtain the regression parameters.
#'    The readers are refered to the \pkg{aftgee} package for details.
#'    The following are permitted:
#'    \describe{
#'      \item{\code{ns}}{Regression parameters are estimated by directly solving 
#'      the nonsmooth estimating equations.}
#'      \item{\code{is}}{Regression parameters are estimated by directly solving 
#'      the induced-smoothing estimating equations.}
#' }
#' @param covTested A character string specifying the covariate which will be tested.
#'    The argument \code{covTested} is necessary only if \code{testType} is 
#'    \code{covForm}.The default option for \code{covTested} is given by "1", which 
#'    represents the first covariate in the formula argument.
#' @param npathsave An integer value specifies he number of paths saved among all the paths.
#'    The default is given by 50. Note that it requires a lot of memory if save all
#'    sampled paths (N by N matrix for each npath andso npath*N*N elements)
#' @param linApprox A logical value. If \code{TRUE}, the multiplier bootstrap is 
#'    computed using the asymptotic linear approximation, which is significantly 
#'    faster. If \code{FALSE}, the estimating equations are solved numerically for 
#'    each bootstrap replication. Defaults to \code{TRUE}.
#' @param seed An optional integer specifying the random seed for reproducibility.
#' @param ... Other arguments passed to methods.
#'  
#' @rdname afttest
#' @export
afttest.formula <- function(object, data, npath = 200, testType = "omnibus", 
                            estMethod = "rr", eqType = "ns", 
                            covTested = 1, npathsave = 50, linApprox = TRUE,
                            seed = NULL, ...) {
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv))
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
  covnames <- all.vars(object[[3]])
  cov.length <- length(covnames)
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
  if (cov.length==1 && length(unique(covariates))==1) {
    return(warning("Intercept-only model detected; The semiparametric AFT model is unable to handle an intercept-only model"))
  }
  
  # beta coefficients from aftsrr function (aftgee package) - with original covariates
  formula <- stats::as.formula(paste0("Surv(time,delta)~",paste(covnames, collapse="+")))
  if (estMethod == "ls") {
    beta <- - aftgee::aftgee(formula, data = DF)$coef.res[-1]
  } else if (estMethod == "rr") {
    beta <- - aftgee::aftsrr(formula, data = DF, eqType = eqType, rankWeights = "gehan")$beta
  } else {
    return(warning("estMethod needs to be one of 'ls' and 'rr'"))
  }
  
  # Covariate Scaling
  time <- DF$time
  delta <- DF$delta
  covariates <- scale(as.matrix(DF[, -(1:2)]))
  DF <- data.frame(time = time, delta = delta, covariates)
  
  # eqType
  if (length(eqType) > 1){
    return(warning("testType needs to be one of 'ns' and 'is'"))
  }
  
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
  if (length(npathsave) > 1){
    return(warning("npathsave needs to be an integer."))
  } else {
    if (!is.numeric(npathsave)) {
      npathsave <- 50
    }
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
  if (testType == "covForm") {
    if (length(covTested) > 1){
      return(warning("the length if covTested needs to be exactly 1."))
    } else {
      if (is.numeric(covTested)) {
        if (covTested %%1 != 0 || covTested > cov.length) {
          return(warning("covTested needs to be postivie integer and less than the lenght of covariates."))
        }
      } else if (is.character(covTested)) {
        if (!covTested %in% covnames) {
          return(warning("covTested needs to specified the one of the covariates in the formula."))
        } else {
          covTested.num <- which(covTested == covnames)
        }
      } else {
        return(warning("covTested needs to be specified correctly."))
      }
    }
  }
  
  # beta coefficients from aftsrr function (aftgee package) - with scaled covariates
  formula <- stats::as.formula(paste0("Surv(time,delta)~",paste(covnames, collapse="+")))
  if (estMethod == "ls") {
    b <- - aftgee::aftgee(formula, data = DF)$coef.res[-1]
    if (!eqType == "ls") {
      warning("eqType must be 'ls' when estMethod is 'ls'")
    }
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

#' @param object A formula expression, of the covTested \code{response ~ predictors}.
#'    The \code{response} is a \code{Surv} object with right censoring.
#'    See the documentation of \code{lm}, \code{coxph} and \code{formula} for details.
#' @param data An optional data frame in which to interpret the variables occurring 
#'    in the formula.
#' @param npath An integer value specifies the number of approximated processes.
#'    The default is given by 200.
#' @param testType A character string specifying the type of the test.
#'    The following are permitted:
#'    \describe{
#'      \item{\code{omnibus}}{an omnibus test}
#'      \item{\code{link}}{a link function test}
#'      \item{\code{covForm}}{a functional covTested of a covariate}
#' }
#' @param eqType A character string specifying the type of the 
#'    estimating equation used to obtain the regression parameters.
#'    The readers are refered to the \pkg{aftgee} package for details.
#'    The following are permitted:
#'    \describe{
#'      \item{\code{ns}}{Regression parameters are estimated by directly solving 
#'      the nonsmooth estimating equations.}
#'      \item{\code{is}}{Regression parameters are estimated by directly solving 
#'      the induced-smoothing estimating equations.}
#' }
#' @param covTested A character string specifying the covariate which will be tested.
#'    The argument \code{covTested} is necessary only if \code{testType} is 
#'    \code{covForm}.The default option for \code{covTested} is given by "1", which 
#'    represents the first covariate in the formula argument.
#' @param npathsave An integer value specifies he number of paths saved among all the paths.
#'    The default is given by 50. Note that it requires a lot of memory if save all
#'    sampled paths (N by N matrix for each npath andso npath*N*N elements)
#' @param linApprox A logical value. If \code{TRUE}, the multiplier bootstrap is 
#'    computed using the asymptotic linear approximation, which is significantly 
#'    faster. If \code{FALSE}, the estimating equations are solved numerically for 
#'    each bootstrap replication. Defaults to \code{TRUE}.
#' @param seed An optional integer specifying the random seed for reproducibility.
#' @param ... Other arguments passed to methods. 
#' 
#' @rdname afttest
#' @export
afttest.aftsrr <- function(object, data, npath = 200, testType = "omnibus", eqType = "ns", 
                           covTested = 1, npathsave = 50, linApprox = TRUE,
                           seed = NULL, ...) {
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv))
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
  if (cov.length==1 && length(unique(covariates))==1) {
    return(warning("Intercept-only model detected; The semiparametric AFT model is unable to handle an intercept-only model"))
  }
  
  # Covariate Scaling
  time <- DF$time
  delta <- DF$delta
  covariates <- scale(as.matrix(DF[, -(1:2)]))
  DF <- data.frame(time = time, delta = delta, covariates)
  
  # eqType
  if (length(eqType) > 1){
    return(warning("testType needs to be one of 'ns' and 'is'"))
  }
  
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
  if (length(npathsave) > 1){
    return(warning("npathsave needs to be an integer."))
  } else {
    if (!is.numeric(npathsave)) {
      npathsave <- 50
    }
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
  if (testType == "covForm") {
    if (length(covTested) > 1){
      return(warning("the length if covTested needs to be exactly 1."))
    } else {
      if (is.numeric(covTested)) {
        if (covTested %%1 != 0 || covTested > cov.length) {
          return(warning("covTested needs to be postivie integer and less than the lenght of covariates."))
        }
      } else if (is.character(covTested)) {
        if (!covTested %in% covnames) {
          return(warning("covTested needs to specified the one of the covariates in the formula."))
        } else {
          covTested.num <- which(covTested == covnames)
        }
      } else {
        return(warning("covTested needs to be specified correctly."))
      }
    }
  }
  
  # beta coefficients from aftsrr function (aftgee package)
  formula <- stats::as.formula(paste0("Surv(time,delta)~",paste(covnames, collapse="+")))
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

#' @param object A formula expression, of the form \code{response ~ predictors}.
#'    The \code{response} is a \code{Surv} object with right censoring.
#'    See the documentation of \code{lm}, \code{coxph} and \code{formula} for details.
#' @param data An optional data frame in which to interpret the variables occurring 
#'    in the formula.
#' @param npath An integer value specifies the number of approximated processes.
#'    The default is given by 200.
#' @param testType A character string specifying the type of the test.
#'    The following are permitted:
#'    \describe{
#'      \item{\code{omnibus}}{an omnibus test}
#'      \item{\code{link}}{a link function test}
#'      \item{\code{covForm}}{a functional form of a covariate}
#' }
#' @param eqType A character string specifying the type of the 
#'    estimating equation used to obtain the regression parameters.
#'    The readers are refered to the \pkg{aftgee} package for details.
#'    The following are permitted:
#'    \describe{
#'      \item{\code{ns}}{Regression parameters are estimated by directly solving 
#'      the nonsmooth estimating equations.}
#'      \item{\code{is}}{Regression parameters are estimated by directly solving 
#'      the induced-smoothing estimating equations.}
#' }
#' @param covTested A character string specifying the covariate which will be tested.
#'    The argument \code{covTested} is necessary only if \code{testType} is 
#'    \code{covForm}.The default option for \code{covTested} is given by "1", which 
#'    represents the first covariate in the formula argument.
#' @param npathsave An integer value specifies he number of paths saved among all the paths.
#'    The default is given by 50. Note that it requires a lot of memory if save all
#'    sampled paths (N by N matrix for each npath andso npath*N*N elements)
#' @param linApprox A logical value. If \code{TRUE}, the multiplier bootstrap is 
#'    computed using the asymptotic linear approximation, which is significantly 
#'    faster. If \code{FALSE}, the estimating equations are solved numerically for 
#'    each bootstrap replication. Defaults to \code{TRUE}.
#' @param seed An optional integer specifying the random seed for reproducibility.
#' @param ... Other arguments passed to methods. 
#' 
#' @rdname afttest
#' @export
afttest.aftgee <- function(object, data, npath = 200, testType = "omnibus", eqType = "ls", 
                           covTested = 1, npathsave = 50, linApprox = TRUE,
                           seed = NULL, ...) {
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv))
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
  if (cov.length==1 && length(unique(covariates))==1) {
    return(warning("Intercept-only model detected; The semiparametric AFT model is unable to handle an intercept-only model"))
  }
  
  # Covariate Scaling
  time <- DF$time
  delta <- DF$delta
  covariates <- scale(as.matrix(DF[, -(1:2)]))
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
  if (length(npathsave) > 1){
    return(warning("npathsave needs to be an integer."))
  } else {
    if (!is.numeric(npathsave)) {
      npathsave <- 50
    }
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
  if (testType == "covForm") {
    if (length(covTested) > 1){
      return(warning("the length if covTested needs to be exactly 1."))
    } else {
      if (is.numeric(covTested)) {
        if (covTested %%1 != 0 || covTested > cov.length) {
          return(warning("covTested needs to be postivie integer and less than the lenght of covariates."))
        }
      } else if (is.character(covTested)) {
        if (!covTested %in% covnames) {
          return(warning("covTested needs to specified the one of the covariates in the formula."))
        } else {
          covTested.num <- which(covTested == covnames)
        }
      } else {
        return(warning("covTested needs to be specified correctly."))
      }
    }
  }
  
  # beta coefficients from aftsrr function (aftgee package)
  formula <- stats::as.formula(paste0("Surv(time,delta)~",paste(covnames, collapse="+")))
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
