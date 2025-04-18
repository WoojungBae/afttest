##############################################################################
## User's Main Function
##############################################################################

#' afttest
#' 
#' @param formula A formula expression, of the form \code{response ~ predictors}.
#'    The \code{response} is a \code{Surv} object object with right censoring.
#'    See the documentation of \code{lm}, \code{coxph} and \code{formula} for details.
#' @param data An optional data frame in which to interpret the variables occurring 
#'    in the formula.
#' @param npath An integer value specifies the number of approximated processes.
#'    The default is given by 200.
#' @param testType A character string specifying the type of the test.
#'    The following are permitted:
#'    \describe{
#'      \item{\code{omni}}{an omnibus test}
#'      \item{\code{link}}{a link function test}
#'      \item{\code{form}}{a functional form}
#' }
#' @param eqType A character string specifying the type of the 
#'    estimating equation used to obtain the regression parameters.
#'    The readers are refered to the \pkg{aftgee} package for details.
#'    The following are permitted:
#'    \describe{
#'      \item{\code{mns}}{Regression parameters are estimated by iterating 
#'      the monotonic non-smoothed Gehan-based estimating equations.}
#'      \item{\code{mis}}{Regression parameters are estimated by iterating 
#'      the monotonic smoothed Gehan-based estimating equations.}
#' }
#' @param estType A character string specifying the type of the estimator used.
#'    The readers are refered to the \pkg{aftgee} package for details.
#'    The following are permitted:
#'    \describe{
#'      \item{\code{aftgee}}{Least-Squares Approach for Accelerated Failure Time 
#'      with Generalized Estimating Equation}
#'      \item{\code{aftsrr}}{Accelerated Failure Time with Smooth Rank Regression}
#' }
#' @param optimType A character string specifying the type of the optimization method.
#'    The following are permitted:
#'    \describe{
#'      \item{\code{DFSANE}}{See the documentation of \pkg{BB} packages for details.}
#'      \item{\code{Nelder-Mead}}{See the documentation of \code{optim} for details.}
#'      \item{\code{BFGS}}{See the documentation of \code{optim} for details.}
#'      \item{\code{CG}}{See the documentation of \code{optim} for details.}
#'      \item{\code{L-BFGS-B}}{See the documentation of \code{optim} for details.}
#'      \item{\code{SANN}}{See the documentation of \code{optim} for details.}
#'      \item{\code{Brent}}{See the documentation of \code{optim} for details.}
#' }
#' @param form A character string specifying the covariate which will be tested.
#'    The argument form is necessary only if \code{testType} is \code{form}.
#'    The default option for \code{form} is given by "1", which represents the 
#'    first covariate in the formula argument.
#' @param npathsave An integer value specifies he number of paths saved among all the paths.
#'    The default is given by 50. Note that it requires a lot of memory if save all
#'    sampled paths (N by N matrix for each npath andso npath*N*N elements)
#' @return \code{afttest} returns an object of class \code{afttest}.
#'    An object of class \code{afttest} is a list containing at least the following components:
#' \describe{
#'    \item{beta}{a vector of beta estimates based on \code{aftsrr}}
#'    \item{SE_process}{estimated standard error of the observed process}
#'    \item{obs_process}{observed process}
#'    \item{app_process}{approximated process}
#'    \item{obs_std_process}{standardized observed process}
#'    \item{app_std_process}{standardized approximated processes}
#'    \item{p_value}{obtained by the unstandardized test}
#'    \item{p_std_value}{obtained by the standardized test}
#'    \item{DF}{a data frame of observed failure time, right censoring indicator, covariates (scaled), 
#'    time-transformed residual based on beta estimates}
#'    \item{npath}{the number of sample paths}
#'    \item{testType}{testType}
#'    \item{eqType}{eqType}
#'    \item{estType}{estType}
#'    \item{optimType}{optimType}
#'    \item{npathsave}{npathsave}
#' }
#'    For an omnibus test, the observed process and the realizations are composed of the 
#'    n by n matrix that rows represent the t and columns represent the x in the 
#'    time-transformed residual order.The observed process and the simulated processes
#'    for checking a functional form and a link function are given by the n by 1 vector
#'    which is a function of x in the time-transformed residual order. 
#' 
#' @importFrom stats optim get_all_vars as.formula model.matrix
#' @importFrom aftgee aftsrr aftgee
#' @importFrom survival Surv
#' 
#' @example inst/examples/ex_afttest.R
#' @export
afttest <- function(formula, data, npath = 200, testType = "omni", 
                    eqType = "mns", estType = "aftsrr", optimType = "DFSANE", 
                    form = 1, npathsave = 50) {
  # , "aftgee"
  # Data Frame
  # DF <- stats::get_all_vars(formula)
  # varnames <- noquote(all.vars(formula))
  # var.length <- length(DF)
  # cov.length <- var.length - 2
  # covnames <- varnames[3:var.length]
  # colnames(DF) <- c("Time", "Delta", paste0("Covari", 1:cov.length))
  
  varnames <- noquote(all.vars(formula))
  var.length <- length(varnames)
  covnames <- varnames[3:var.length]
  cov.length <- length(covnames)
  
  scall <- match.call()
  mnames <- c("", "formula", "data")
  cnames <- names(scall)
  cnames <- cnames[match(mnames, cnames, 0)]
  mcall <- scall[cnames]
  mcall[[1]] <- as.name("model.frame")
  m <- eval(mcall, parent.frame())
  mterms <- attr(m, "terms")
  obj <- unclass(m[, 1])
  formula[[2]] <- NULL
  if (formula == ~1) {
    DF <- cbind(obj, zero = 0)
  } else {
    DF <- cbind(obj, model.matrix(mterms, m))
    if (sum(colnames(DF) == "(Intercept)") > 0) {
      DF <- DF[, -which(colnames(DF) == "(Intercept)")]
    }
  }
  DF <- as.data.frame(DF)
  colnames(DF) <- c("Time", "Delta", paste0("Covari", 1:cov.length))
  
  # check&delete NA, -Inf, Inf, ...
  missingmessage <- NA
  DF[DF=="-Inf" | DF=="Inf"] <- NA
  whichNA_DF <- which(apply(is.na(DF), 1, sum)>0)
  nNA_DF <- length(whichNA_DF)
  if (nNA_DF > 0){
    missingmessage <- paste0("(",nNA_DF, " observations deleted due to missingness out of ", nrow(DF), ")")
    DF <- DF[-whichNA_DF,]
  } else {
    missingmessage <- paste0("(No missing observed)")
  }
  
  if (any(DF$Time <= 0)) {
    return(warning("Time must be positive number"))
  }
  if (cov.length==1 && length(unique(Covari))==1) {
    return(warning("Intercept-only model detected; The semiparametric AFT model is unable to handle an intercept-only model"))
  }
  
  # Covariate Scaling
  Time <- DF$Time
  Delta <- DF$Delta
  Covari <- scale(as.matrix(DF[, 3:var.length]))
  DF[, 3:var.length] <- Covari
  
  # unique_Delta <- unique(Delta)
  # if (length(unique_Delta)==2){
  #   if (any(c(0,1) == sort(unique_Delta))){
  #     Delta <- ifelse(Delta == unique_Delta[1], 0, 1)
  #     warning(paste0(unique_Delta[1], "=0 is assumed to be observed and ", unique_Delta[2], "=1 is assumed to be censred"))
  #   }
  #  else {
  #   return(warning("Delta must have 2 statuses (0=observed and 1=censored)"))
  # }
  
  # beta coefficients from aftsrr function (aftgee package)
  formula <- stats::as.formula(paste0("Surv(Time,Delta)~",paste(paste0("Covari", 1:cov.length), collapse="+")))
  if (estType == "aftgee") {
    b <- - aftgee::aftgee(formula, data = DF)$coef.res[-1]
  } else if (estType == "aftsrr") {
    b <- - aftgee::aftsrr(formula, data = DF, eqType = eqType, rankWeights = "gehan")$beta
  } else {
    return(warning("estType needs to be one of 'aftgee' and 'aftsrr'"))
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
    return(warning("testType needs to be one of 'omni', 'link', or 'form'"))
  } else {
    if (!testType %in% c("omni","link","form")) {
      testType <- "omni"
    }
  }
  
  # eqType
  if (length(eqType) > 1){
    return(warning("testType needs to be one of 'mns' and 'mis'"))
  } else {
    if (!eqType %in% c("mns","mis")) {
      eqType <- "mns"
    }
  }
  
  # optimType
  if (length(optimType) > 1){
    return(warning("optimType needs to be one of 'DFSANE', 'Nelder-Mead', 'BFGS', 'CG', 'L-BFGS-B', 'SANN', 'Brent'"))
  } else {
    if (!optimType %in% c("DFSANE","Nelder-Mead","BFGS","CG","L-BFGS-B","SANN","Brent")) {
      optimType <- "DFSANE"
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
  
  # form
  if (testType == "form") {
    if (length(form) > 1){
      return(warning("the length if form needs to be exactly 1."))
    } else {
      if (is.numeric(form)) {
        if (form%%1 != 0 || form > cov.length) {
          return(warning("form needs to be postivie integer and less than the lenght of covariates."))
        }
      } else if (is.character(form)) {
        if (!form %in% covnames) {
          return(warning("form needs to specified the one of the covariates in the formula."))
        } else {
          form <- which(form == covnames)
        }
      } else {
        return(warning("form needs to be specified correctly."))
      }
    }
  }
  
  if (optimType != "DFSANE"){
    if (eqType=="mns"){
      if (testType == "omni") {
        out <- .Call(`_afttest_omni_mns_optim`, npath, b, Time, Delta, Covari, optimType, npathsave)
      } else if (testType == "link") {
        out <- .Call(`_afttest_link_mns_optim`, npath, b, Time, Delta, Covari, optimType, npathsave)
      } else if (testType == "form") {
        out <- .Call(`_afttest_form_mns_optim`, npath, b, Time, Delta, Covari, optimType, form, npathsave)
      }
    } else if (eqType=="mis"){
      if (testType == "omni") {
        out <- .Call(`_afttest_omni_mis_optim`, npath, b, Time, Delta, Covari, optimType, npathsave)
      } else if (testType == "link") {
        out <- .Call(`_afttest_link_mis_optim`, npath, b, Time, Delta, Covari, optimType, npathsave)
      } else if (testType == "form") {
        out <- .Call(`_afttest_form_mis_optim`, npath, b, Time, Delta, Covari, optimType, form, npathsave)
      }
    }
  } else if (optimType == "DFSANE"){
    if (eqType=="mns"){
      if (testType == "omni") {
        out <- .Call(`_afttest_omni_mns_DFSANE`, npath, b, Time, Delta, Covari, npathsave)
      } else if (testType == "link") {
        out <- .Call(`_afttest_link_mns_DFSANE`, npath, b, Time, Delta, Covari, npathsave)
      } else if (testType == "form") {
        out <- .Call(`_afttest_form_mns_DFSANE`, npath, b, Time, Delta, Covari, form, npathsave)
      }
    } else if (eqType=="mis"){
      if (testType == "omni") {
        out <- .Call(`_afttest_omni_mis_DFSANE`, npath, b, Time, Delta, Covari, npathsave)
      } else if (testType == "link") {
        out <- .Call(`_afttest_link_mis_DFSANE`, npath, b, Time, Delta, Covari, npathsave)
      } else if (testType == "form") {
        out <- .Call(`_afttest_form_mis_DFSANE`, npath, b, Time, Delta, Covari, form, npathsave)
      }
    }
  } else {
    return(warning("Check your code"))
  }
  
  class(out) <- "afttest"
  out$names <- varnames
  out$call <- match.call()
  out$missingmessage <- missingmessage
  
  out$DF <- DF
  out$beta <- b
  out$npath <- npath
  out$eqType <- eqType
  out$testType <- testType
  out$estType <- estType
  out$optimType <- optimType
  out$npathsave <- npathsave
  if (testType == "form") {out$form <- form}
  
  return(out)
}