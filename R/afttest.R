##############################################################################
## User's Main Function
##############################################################################

#' afttest
#' 
#' @param formula A formula expression, of the form \code{response ~ predictors}.
#'    The \code{response} is a \code{Surv} object object with right censoring.
#'    See the documentation of \code{lm}, \code{coxph} and \code{formula} for details.
#' @param path A numeric value specifies the approximated processes number.
#'    The default is given by 200.
#' @param testtype A character string specifying the type of the test.
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
#'      \item{\code{mis}}{Regression parameters are estimated by iterating 
#'      the monotonic smoothed Gehan-based estimating equations.}
#'      \item{\code{mns}}{Regression parameters are estimated by iterating 
#'      the monotonic non-smoothed Gehan-based estimating equations.}
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
#'    The argument form is necessary only if testtype \code{form}.
#'    The default option for testtype is given by "1", which represents the 
#'    first covariate in the formula argument.
#' @param pathsave A numeric value specifies he number of paths saved among all the paths.
#'    The default is given by 100. Note that it requires a lot of memory if save all
#'    sampled paths (N by N matrix for each path andso path*N*N elements)
#' @return \code{afttest} returns an object of class \code{afttest}.
#'    An object of class \code{afttest} is a list containing at least the following components:
#' \describe{
#'    \item{beta}{a vector of beta estimates based on \code{aftsrr}}
#'    \item{path}{the number of sample paths}
#'    \item{Time}{observed failure time}
#'    \item{Delta}{right censoring indicator}
#'    \item{Covari}{covariates}
#'    \item{Resid}{time-transformed residual based on beta estimates}
#'    \item{SE_process}{estimated standard error of the observed process}
#'    \item{obs_process}{observed process}
#'    \item{app_process}{approximated process}
#'    \item{obs_std_process}{standardized observed process}
#'    \item{app_std_process}{standardized approximated processes}
#'    \item{p_value}{obtained by the unstandardized test}
#'    \item{p_std_value}{obtained by the standardized test}
#' }
#'    For an omnibus test, the observed process and the realizations are composed of the 
#'    n by n matrix that rows represent the t and columns represent the x in the 
#'    time-transformed residual order.The observed process and the simulated processes
#'    for checking a functional form and a link function are given by the n by 1 vector
#'    which is a function of x in the time-transformed residual order. 
#'    
#' @example inst/examples/ex_afttest.R
#' @export
afttest = function(formula, path = 200, testtype = "omni", eqType = "mis", 
                   optimType = "DFSANE", form = 1, pathsave = 100) {
  
  DF = get_all_vars(formula)
  varnames = noquote(all.vars(formula))
  var.length = ncol(DF)
  cov.length = var.length - 2
  
  colnames(DF) = c("Time", "Delta", paste0("Covari", 1:cov.length))
  
  Time = DF$Time
  Delta = DF$Delta
  Covari = scale(as.matrix(DF[, 3:var.length]))
  
  # path
  if (!is.numeric(path)) {
    path = 200
  } else {
    path = max(path,10)
  }
  
  # testtype
  if (!testtype %in% c("omni","link","form")) {
    testtype = "omni"
  } else {
    testtype = match.arg(testtype, c("omni","link","form"))
  }
  
  # eqType
  if (!eqType %in% c("mis","mns")) {
    eqType = "mis"
  } else {
    eqType = match.arg(eqType, c("mis","mns"))
  }
  
  # optimType
  if (!optimType %in% c("DFSANE","Nelder-Mead","BFGS","CG","L-BFGS-B","SANN","Brent")) {
    optimType = "DFSANE"
  } else {
    optimType = match.arg(optimType, c("DFSANE","Nelder-Mead","BFGS","CG","L-BFGS-B","SANN","Brent"))
  }
  
  # pathsave
  if (!is.numeric(pathsave)) {
    pathsave = 200
  } else {
    pathsave = max(path,10)
  }
  
  # form
  if (testtype == "form") {
    if (length(form) > 1){
      stop("the length if form needs to be exactly 1." )
    } else {
      if (is.character(form)) {
        if (!form %in% varnames) {
          form = "mis"
        } else {
          form = match.arg(form, varnames) - 2
        }
      } else if (is.numeric(form)) {
        if (form > cov.length){
          stop("form is greater than the lenght of covariates. form needs to be specified correctly." )
        } 
      } else {
        stop("form needs to be specified correctly." )
      }
    }
  }
  
  # beta coefficients from aftsrr function (aftgee package)
  b = - aftsrr(formula, eqType = eqType)$beta
  
  if (optimType != "DFSANE"){
    if (eqType=="mns"){
      if (testtype == "omni") {
        out = .Call(`_afttest_omni_mns_optim`, path, b, Time, Delta, Covari, optimType, pathsave)
      } else if (testtype == "link") {
        out = .Call(`_afttest_link_mns_optim`, path, b, Time, Delta, Covari, optimType, pathsave)
      } else if (testtype == "form") {
        out = .Call(`_afttest_form_mns_optim`, path, b, Time, Delta, Covari, optimType, form, pathsave)
      }
    } else if (eqType=="mis"){
      if (testtype == "omni") {
        out = .Call(`_afttest_omni_mis_optim`, path, b, Time, Delta, Covari, optimType, pathsave)
      } else if (testtype == "link") {
        out = .Call(`_afttest_link_mis_optim`, path, b, Time, Delta, Covari, optimType, pathsave)
      } else if (testtype == "form") {
        out = .Call(`_afttest_form_mis_optim`, path, b, Time, Delta, Covari, optimType, form, pathsave)
      }
    }
  } else if (optimType == "DFSANE"){
    if (eqType=="mns"){
      if (testtype == "omni") {
        out = .Call(`_afttest_omni_mns_DFSANE`, path, b, Time, Delta, Covari, pathsave)
      } else if (testtype == "link") {
        out = .Call(`_afttest_link_mns_DFSANE`, path, b, Time, Delta, Covari, pathsave)
      } else if (testtype == "form") {
        out = .Call(`_afttest_form_mns_DFSANE`, path, b, Time, Delta, Covari, form, pathsave)
      }
    } else if (eqType=="mis"){
      if (testtype == "omni") {
        out = .Call(`_afttest_omni_mis_DFSANE`, path, b, Time, Delta, Covari, pathsave)
      } else if (testtype == "link") {
        out = .Call(`_afttest_link_mis_DFSANE`, path, b, Time, Delta, Covari, pathsave)
      } else if (testtype == "form") {
        out = .Call(`_afttest_form_mis_DFSANE`, path, b, Time, Delta, Covari, form, pathsave)
      }
    }
  } else {
    stop("Check your code")
  }
  
  class(out) = "afttest"
  out$names = varnames
  out$call = match.call()
  
  return(out)
}