#' afttest
#'
#' It gives several test statistics for cheking the aft model assumptions.
#' 
#' @param formula The argument formula specifies the model to be fitted with
#' the variables coming with data. The expression of the formula argument
#' is equivalent to the Surv in the survival package. The object Surv
#' consists of two columns. The first one is the observed failure time and
#'  the second one is the indicator variable, specifying right censoring.
#' @param path The argument path determines the number of simulations of
#' the approximated process. The default is given by 200.
#' @param testtype The argument testtype includes the aforementioned an
#' omnibus test ("omni"), a functional form ("form") and a link function
#' ("linkftn"). The rank weight in the package is the Gehan"s weight and
#' each weight of the test statistics is determined by the testtype
#' arguments. The default option for testtype is given by "omni".
#' @param eqType The argument eqType determines the equation type to estimate
#' the regression parameter while generating approximated process. The following
#' are permitted. Regression parameters are estimated by directly solving
#' the monotonic nonsmooth estimating equations ("mns"). Regression parameters
#' are estimated by directly solving the monotonic induced-smoothing
#' estimating equations.
#' @param optimType The argument optimType determines the algorithm to the
#' objective function be minimized. User can choose one of the following algorithms: 
#' "DFSANE", "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", and "Brent". The
#' default option is "DFSANE".
#' @param form The argument form is necessary only if testtype
#' is given as "form" and it determines a covariate which will be tested.
#' It needs to be specified the name of covariates in the formula argument
#' and the default option is "1, which represents the first covariate
#' in the formula argument.
#' 
#' @return The function afttest gives the list as a result. The result 
#' consists of the number of paths ($path), the estimated beta ($beta), 
#' the observed failure time ($Time), the right censoring indicator ($Delta), 
#' the covariates ($Covari), the time-transformed residual ($Resid), the 
#' estimated standard error of the observed process ($SE_process), the 
#' observed process ($obs_process), a number of the simulated processes 
#' ($app_process), the standardized observed process ($obs_std_process), 
#' the standardized processes of realizations ($app_std_process) and two 
#' kinds of the p-value obtained by the unstandardized test and the 
#' standardized test ($p_value and $p_std_value). Now, we offer two types of 
#' p-values for all tests even though the p-value for the standardized test 
#' is only used for an omnibus test. For an omnibus test, the observed process 
#' and the realizations are composed of the n by n matrix that rows represent 
#' the t and columns represent the x in the time-transformed residual order. 
#' The observed process and the simulated processes for checking a functional 
#' form and a link function are given by the n by 1 vector which is a function 
#' of x in the time-transformed residual order. 
#' 
#' @examples 
#' library(afttest)
#' library(survival)
#' 
#' set.seed(1)
#' 
#' cgd_data = subset(cgd,enum==1,c(tstop,status,treat,age,steroids))
#' 
#' trt = ifelse(cgd_data$treat=="placebo",0,1)
#' age = cgd_data$age ; log_age = log(age)
#' str = cgd_data$steroids
#' 
#' X_cgd = cgd_data$tstop
#' D_cgd = cgd_data$status
#' table(D_cgd)
#' 
#' #---------------------------------------------------------
#' 
#' path=30
#' 
#' result_omni = afttest(Surv(X_cgd,D_cgd)~trt+str+age,path=path,testtype="omni")
#' 
#' # graph_omni_unstd = afttestplot(result_omni,path=10,std="unstd")
#' # graph_omni_std = afttestplot(result_omni,path=10,std="std")
#' 
#' # result_link = afttest(Surv(X_cgd,D_cgd)~trt+str+age,path=path,testtype="link",eqType="DFSANE")
#' #   
#' # graph_link_unstd = afttestplot(result_link,path=10,std="unstd")
#' # graph_link_std = afttestplot(result_link,path=10,std="std")
#' # 
#' # result_form = afttest(Surv(X_cgd,D_cgd)~trt+str+age,path=path,testtype="form",eqType="DFSANE")
#' # 
#' # graph_form_unstd = afttestplot(result_form,path=10,std="unstd")
#' # graph_form_std = afttestplot(result_form,path=10,std="std")
#'  
#' @importFrom stats optim get_all_vars
#' @importFrom aftgee aftsrr
#' @importFrom survival Surv
#' 
#' @export
afttest = function(formula, path = 200, testtype = c("omni","link","form"), eqType = c("mis","mns"), 
                   optimType = c("DFSANE","Nelder-Mead","BFGS","CG","L-BFGS-B","SANN","Brent"), form = 1) {
  
  if(length(testtype) != 1){testtype = "omni"}
  if(length(eqType) != 1){eqType = "mis"}
  if(length(optimType) != 1){optimType = "DFSANE"}
  
  dataset = get_all_vars(formula)
  varnames = noquote(all.vars(formula))
  var.length = ncol(dataset)
  cov.length = var.length - 2
  
  colnames(dataset) = c("Time", "Delta", paste0("Covari", 1:cov.length))
  
  Time = dataset$Time
  Delta = dataset$Delta
  Covari = as.matrix(dataset[, 3:var.length])
  
  if (length(which(varnames == form[1])) != 0) {
    form = which(varnames == form[1]) - 2
  }
  
  b = - aftgee::aftsrr(formula, eqType = eqType)$beta
  
  if (optimType != "DFSANE"){
    if (eqType=="mns"){
      if (testtype == "omni") {
        return(.Call(`_afttest_omni_mns_optim`, path, b, Time, Delta, Covari, optimType))
      } else if (testtype == "link") {
        return(.Call(`_afttest_link_mns_optim`, path, b, Time, Delta, Covari, optimType))
      } else if (testtype == "form") {
        return(.Call(`_afttest_form_mns_optim`, path, b, Time, Delta, Covari, optimType, form))
      }
    } else if (eqType=="mis"){
      if (testtype == "omni") {
        return(.Call(`_afttest_omni_mis_optim`, path, b, Time, Delta, Covari, optimType))
      } else if (testtype == "link") {
        return(.Call(`_afttest_link_mis_optim`, path, b, Time, Delta, Covari, optimType))
      } else if (testtype == "form") {
        return(.Call(`_afttest_form_mis_optim`, path, b, Time, Delta, Covari, optimType, form))
      }
    }
  } else if (optimType == "DFSANE"){
    if (eqType=="mns"){
      if (testtype == "omni") {
        return(.Call(`_afttest_omni_mns_DFSANE`, path, b, Time, Delta, Covari))
      } else if (testtype == "link") {
        return(.Call(`_afttest_link_mns_DFSANE`, path, b, Time, Delta, Covari))
      } else if (testtype == "form") {
        return(.Call(`_afttest_form_mns_DFSANE`, path, b, Time, Delta, Covari, form))
      }
    } else if (eqType=="mis"){
      if (testtype == "omni") {
        return(.Call(`_afttest_omni_mis_DFSANE`, path, b, Time, Delta, Covari))
      } else if (testtype == "link") {
        return(.Call(`_afttest_link_mis_DFSANE`, path, b, Time, Delta, Covari))
      } else if (testtype == "form") {
        return(.Call(`_afttest_form_mis_DFSANE`, path, b, Time, Delta, Covari, form))
      }
    }
  }
  
  return(print("Check your code"))
}

