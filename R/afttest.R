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
#'      \item{\code{covform}}{a functional form of a covariate}
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
#' @param cov.tested A character string specifying the covariate which will be tested.
#'    The argument \code{cov.tested} is necessary only if \code{testType} is 
#'    \code{covform}.The default option for \code{cov.tested} is given by "1", which 
#'    represents the first covariate in the formula argument.
#' @param npathsave An integer value specifies he number of paths saved among all the paths.
#'    The default is given by 50. Note that it requires a lot of memory if save all
#'    sampled paths (N by N matrix for each npath andso npath*N*N elements)
#' @param ... Other arguments passed to methods.
#'  
#' @rdname afttest
#' @export
afttest.formula <- function(object, data, npath = 200, testType = "omnibus", 
                            estMethod = "rr", eqType = "ns", 
                            cov.tested = 1, npathsave = 50, ...) {
  
  # Extract variable names and dimensions
  # varnames <- noquote(all.vars(object))
  # var.length <- length(varnames)
  # covnames <- varnames[3:var.length]
  # cov.length <- length(covnames)
  
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
  DF[DF == "-Inf" | DF == "Inf"] <- NA
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
  
  # unique_Delta <- unique(delta)
  # if (length(unique_Delta)==2){
  #   if (any(c(0,1) == sort(unique_Delta))){
  #     delta <- ifelse(delta == unique_Delta[1], 0, 1)
  #     warning(paste0(unique_Delta[1], "=0 is assumed to be observed and ",
  #                    unique_Delta[2], "=1 is assumed to be censored"))
  #   } else {
  #     return(warning("delta must have 2 statuses (0=observed and 1=censored)"))
  #   }
  # }
  
  # eqType
  if (length(eqType) > 1){
    return(warning("testType needs to be one of 'ns' and 'is'"))
  } else {
    if (!eqType %in% c("ns","is")) {
      print(warning("'ns' is used by default"))
      eqType <- "ns"
    }
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
    return(warning("testType needs to be one of 'omnibus', 'link', or 'covform'"))
  } else {
    if (!testType %in% c("omnibus","link","covform")) {
      testType <- "omnibus"
    }
  }
  
  optimMethod <- "DFSANE"
  # # optimMethod
  # if (length(optimMethod) > 1){
  #   return(warning("optimMethod needs to be one of 'DFSANE', 'Nelder-Mead', 'BFGS', 'CG', 'L-BFGS-B', 'SANN', 'Brent'"))
  # } else {
  #   if (!optimMethod %in% c("DFSANE","Nelder-Mead","BFGS","CG","L-BFGS-B","SANN","Brent")) {
  #     optimMethod <- "DFSANE"
  #   }
  # }
  
  # npathsave
  if (length(npathsave) > 1){
    return(warning("npathsave needs to be an integer."))
  } else {
    if (!is.numeric(npathsave)) {
      npathsave <- 50
    }
  }
  
  # cov.tested
  if (testType == "covform") {
    if (length(cov.tested) > 1){
      return(warning("the length if cov.tested needs to be exactly 1."))
    } else {
      if (is.numeric(cov.tested)) {
        if (cov.tested %%1 != 0 || cov.tested > cov.length) {
          return(warning("cov.tested needs to be postivie integer and less than the lenght of covariates."))
        }
      } else if (is.character(cov.tested)) {
        if (!cov.tested %in% covnames) {
          return(warning("cov.tested needs to specified the one of the covariates in the formula."))
        } else {
          cov.tested <- which(cov.tested == covnames)
        }
      } else {
        return(warning("cov.tested needs to be specified correctly."))
      }
    }
  }
  
  # beta coefficients from aftsrr function (aftgee package) - with scaled covariates
  formula <- stats::as.formula(paste0("Surv(time,delta)~",paste(covnames, collapse="+")))
  if (estMethod == "ls") {
    b <- - aftgee::aftgee(formula, data = DF)$coef.res[-1]
  } else if (estMethod == "rr") {
    b <- - aftgee::aftsrr(formula, data = DF, eqType = eqType, rankWeights = "gehan")$beta
  } else {
    return(warning("estMethod needs to be one of 'ls' and 'rr'"))
  }
  
  # This function contains the core logic (the C++ calls)
  out <- .afttest_worker(b, time, delta, covariates, npath, testType,
                         eqType, optimMethod, cov.tested, npathsave)
  
  out$call <- scall
  out$beta <- beta
  # out$DF <- data
  out$DF <- DF
  out$estMethod <- estMethod # It's an aftsrr object
  out$missingmessage <- missingmessage
  
  return(out)
}

#' @param object A formula expression, of the cov.tested \code{response ~ predictors}.
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
#'      \item{\code{covform}}{a functional cov.tested of a covariate}
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
#' @param cov.tested A character string specifying the covariate which will be tested.
#'    The argument \code{cov.tested} is necessary only if \code{testType} is 
#'    \code{covform}.The default option for \code{cov.tested} is given by "1", which 
#'    represents the first covariate in the formula argument.
#' @param npathsave An integer value specifies he number of paths saved among all the paths.
#'    The default is given by 50. Note that it requires a lot of memory if save all
#'    sampled paths (N by N matrix for each npath andso npath*N*N elements)
#' @param ... Other arguments passed to methods. 
#' 
#' @rdname afttest
#' @export
afttest.aftsrr <- function(object, data, npath = 200, testType = "omnibus", eqType = "ns", 
                           cov.tested = 1, npathsave = 50, ...) {
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
  DF[DF == "-Inf" | DF == "Inf"] <- NA
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
  
  # unique_Delta <- unique(delta)
  # if (length(unique_Delta)==2){
  #   if (any(c(0,1) == sort(unique_Delta))){
  #     delta <- ifelse(delta == unique_Delta[1], 0, 1)
  #     warning(paste0(unique_Delta[1], "=0 is assumed to be observed and ",
  #                    unique_Delta[2], "=1 is assumed to be censored"))
  #   } else {
  #     return(warning("delta must have 2 statuses (0=observed and 1=censored)"))
  #   }
  # }
  
  # eqType
  if (length(eqType) > 1){
    return(warning("testType needs to be one of 'ns' and 'is'"))
  } else {
    if (!eqType %in% c("ns","is")) {
      print(warning("'ns' is used by default"))
      eqType <- "ns"
    }
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
    return(warning("testType needs to be one of 'omnibus', 'link', or 'covform'"))
  } else {
    if (!testType %in% c("omnibus","link","covform")) {
      testType <- "omnibus"
    }
  }
  
  optimMethod <- "DFSANE"
  # # optimMethod
  # if (length(optimMethod) > 1){
  #   return(warning("optimMethod needs to be one of 'DFSANE', 'Nelder-Mead', 'BFGS', 'CG', 'L-BFGS-B', 'SANN', 'Brent'"))
  # } else {
  #   if (!optimMethod %in% c("DFSANE","Nelder-Mead","BFGS","CG","L-BFGS-B","SANN","Brent")) {
  #     optimMethod <- "DFSANE"
  #   }
  # }
  
  # npathsave
  if (length(npathsave) > 1){
    return(warning("npathsave needs to be an integer."))
  } else {
    if (!is.numeric(npathsave)) {
      npathsave <- 50
    }
  }
  
  # cov.tested
  if (testType == "covform") {
    if (length(cov.tested) > 1){
      return(warning("the length if cov.tested needs to be exactly 1."))
    } else {
      if (is.numeric(cov.tested)) {
        if (cov.tested %%1 != 0 || cov.tested > cov.length) {
          return(warning("cov.tested needs to be postivie integer and less than the lenght of covariates."))
        }
      } else if (is.character(cov.tested)) {
        if (!cov.tested %in% covnames) {
          return(warning("cov.tested needs to specified the one of the covariates in the formula."))
        } else {
          cov.tested <- which(cov.tested == covnames)
        }
      } else {
        return(warning("cov.tested needs to be specified correctly."))
      }
    }
  }
  
  # beta coefficients from aftsrr function (aftgee package)
  formula <- stats::as.formula(paste0("Surv(time,delta)~",paste(covnames, collapse="+")))
  b <- - aftgee::aftsrr(formula, data = DF, eqType = eqType, rankWeights = "gehan")$beta
  
  # --- CALL THE INTERNAL WORKER FUNCTION ---
  # This function contains the core logic (the C++ calls)
  out <- .afttest_worker(b, time, delta, covariates, npath, testType,
                         eqType, optimMethod, cov.tested, npathsave)
  out$beta <- -object$beta
  out$call <- scall
  # out$DF <- data
  out$DF <- DF
  out$estMethod <- "rr"
  out$missingmessage <- missingmessage
  
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
#'      \item{\code{covform}}{a functional form of a covariate}
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
#' @param cov.tested A character string specifying the covariate which will be tested.
#'    The argument \code{cov.tested} is necessary only if \code{testType} is 
#'    \code{covform}.The default option for \code{cov.tested} is given by "1", which 
#'    represents the first covariate in the formula argument.
#' @param npathsave An integer value specifies he number of paths saved among all the paths.
#'    The default is given by 50. Note that it requires a lot of memory if save all
#'    sampled paths (N by N matrix for each npath andso npath*N*N elements)
#' @param ... Other arguments passed to methods. 
#' 
#' @rdname afttest
#' @export
afttest.aftgee <- function(object, data, npath = 200, testType = "omnibus", eqType = "ls", 
                           cov.tested = 1, npathsave = 50, ...) {
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
  DF[DF == "-Inf" | DF == "Inf"] <- NA
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
  
  # unique_Delta <- unique(delta)
  # if (length(unique_Delta)==2){
  #   if (any(c(0,1) == sort(unique_Delta))){
  #     delta <- ifelse(delta == unique_Delta[1], 0, 1)
  #     warning(paste0(unique_Delta[1], "=0 is assumed to be observed and ",
  #                    unique_Delta[2], "=1 is assumed to be censored"))
  #   } else {
  #     return(warning("delta must have 2 statuses (0=observed and 1=censored)"))
  #   }
  # }
  
  # estMethod
  estMethod = "ls"
  
  # eqType
  if (estMethod == "ls") {
    eqType <- "ns"
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
    return(warning("testType needs to be one of 'omnibus', 'link', or 'covform'"))
  } else {
    if (!testType %in% c("omnibus","link","covform")) {
      testType <- "omnibus"
    }
  }
  
  optimMethod <- "DFSANE"
  # # optimMethod
  # if (length(optimMethod) > 1){
  #   return(warning("optimMethod needs to be one of 'DFSANE', 'Nelder-Mead', 'BFGS', 'CG', 'L-BFGS-B', 'SANN', 'Brent'"))
  # } else {
  #   if (!optimMethod %in% c("DFSANE","Nelder-Mead","BFGS","CG","L-BFGS-B","SANN","Brent")) {
  #     optimMethod <- "DFSANE"
  #   }
  # }
  
  # npathsave
  if (length(npathsave) > 1){
    return(warning("npathsave needs to be an integer."))
  } else {
    if (!is.numeric(npathsave)) {
      npathsave <- 50
    }
  }
  
  # cov.tested
  if (testType == "covform") {
    if (length(cov.tested) > 1){
      return(warning("the length if cov.tested needs to be exactly 1."))
    } else {
      if (is.numeric(cov.tested)) {
        if (cov.tested %%1 != 0 || cov.tested > cov.length) {
          return(warning("cov.tested needs to be postivie integer and less than the lenght of covariates."))
        }
      } else if (is.character(cov.tested)) {
        if (!cov.tested %in% covnames) {
          return(warning("cov.tested needs to specified the one of the covariates in the formula."))
        } else {
          cov.tested <- which(cov.tested == covnames)
        }
      } else {
        return(warning("cov.tested needs to be specified correctly."))
      }
    }
  }
  
  # beta coefficients from aftsrr function (aftgee package)
  formula <- stats::as.formula(paste0("Surv(time,delta)~",paste(covnames, collapse="+")))
  b <- - aftgee::aftgee(formula, data = DF)$coef.res[-1]
  
  # This function contains the core logic (the C++ calls)
  out <- .afttest_worker(b, time, delta, covariates, npath, testType,
                         eqType, optimMethod, cov.tested, npathsave)
  out$beta <- -object$coef.res[-1]
  out$call <- scall
  # out$DF <- data
  out$DF <- DF
  out$estMethod <- "ls"
  out$missingmessage <- missingmessage
  
  return(out)
}

#' Internal worker function for afttest
#' @noRd
.afttest_worker <- function(b, time, delta, covariates, npath, testType,
                            eqType, optimMethod, cov.tested, npathsave) {
  # C++ functions
  if (optimMethod != "DFSANE"){
    if (eqType=="ns"){
      if (testType == "omnibus") {
        out <- .Call(`_afttest_omni_mns_optim`, npath, b, time, delta, covariates, optimMethod, npathsave)
      } else if (testType == "link") {
        out <- .Call(`_afttest_link_mns_optim`, npath, b, time, delta, covariates, optimMethod, npathsave)
      } else if (testType == "covform") {
        out <- .Call(`_afttest_form_mns_optim`, npath, b, time, delta, covariates, optimMethod, cov.tested, npathsave)
      }
    } else if (eqType=="is"){
      if (testType == "omnibus") {
        out <- .Call(`_afttest_omni_mis_optim`, npath, b, time, delta, covariates, optimMethod, npathsave)
      } else if (testType == "link") {
        out <- .Call(`_afttest_link_mis_optim`, npath, b, time, delta, covariates, optimMethod, npathsave)
      } else if (testType == "covform") {
        out <- .Call(`_afttest_form_mis_optim`, npath, b, time, delta, covariates, optimMethod, cov.tested, npathsave)
      }
    }
  } else if (optimMethod == "DFSANE"){
    if (eqType=="ns"){
      if (testType == "omnibus") {
        out <- .Call(`_afttest_omni_mns_DFSANE`, npath, b, time, delta, covariates, npathsave)
      } else if (testType == "link") {
        out <- .Call(`_afttest_link_mns_DFSANE`, npath, b, time, delta, covariates, npathsave)
      } else if (testType == "covform") {
        out <- .Call(`_afttest_form_mns_DFSANE`, npath, b, time, delta, covariates, cov.tested, npathsave)
      }
    } else if (eqType=="is"){
      if (testType == "omnibus") {
        out <- .Call(`_afttest_omni_mis_DFSANE`, npath, b, time, delta, covariates, npathsave)
      } else if (testType == "link") {
        out <- .Call(`_afttest_link_mis_DFSANE`, npath, b, time, delta, covariates, npathsave)
      } else if (testType == "covform") {
        out <- .Call(`_afttest_form_mis_DFSANE`, npath, b, time, delta, covariates, cov.tested, npathsave)
      }
    }
  } else {
    return(warning("Check your code"))
  }
  
  class(out) <- c("afttest", "htest")
  out$betascaled <- b
  out$npath <- npath
  out$eqType <- eqType
  out$testType <- testType
  out$optimMethod <- optimMethod
  out$npathsave <- npathsave
  if (testType == "covform") {
    out$cov.tested <- cov.tested
  }
  
  return(out)
}
