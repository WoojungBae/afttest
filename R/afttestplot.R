#' afttestplot
#'
#' It gives plot for cheking the aft model assumptions.
#' @param result For function afttestplot, the only required argument is
#' afttestresult based on the result from the afttest function. Whatever
#' the testtype of the result is, it automatically gives the corresponding
#' graph.
#' @param path The argument path is for the number of the simulation paths
#' that is plotted in the graph. Therefore it needs to be equal or less than
#' the number of paths used in by afttest function, otherwise it is given as
#' the number of paths used in by afttest function. The default is set to be 100.
#' @param std The option for the argument std is "unstd" and "std". In this 
#' argument, "std" is the default.
#' @return Basically, a graph from the afttestplot is based on the packages
#' ggplot2 (Wickham, 2009) and gridExtra (Auguie, 2017). It offers a graph that
#' y-axis is the test statistics and x-axis represents the rank of the subjects
#' ordered by time transformed residual. Since the result of the omnibus test
#' is the form of n by n matrix, some quantiles of x, which are used in weight,
#' are plotted for graphs, i.e. 0%, 10%, 25%, 40%, 50%, 60%, 75%, 90%, and 100%
#' are used.
#' 
#' @include source_r.R
#' 
#' @examples 
#' library(afttest)
#' library(survival)
#' 
#' set.seed(1)
#' path = 5
#' 
#' cgd_data = subset(cgd,enum==1)
#' D_cgd = cgd_data$status
#' X_cgd = cgd_data$tstop - cgd_data$tstart
#' X_cgd = X_cgd + runif(length(X_cgd))/1e4
#' 
#' trt = ifelse(cgd_data$treat=="placebo",0,1)
#' str = cgd_data$steroids
#' age = cgd_data$age; log_age = log(age)
#' wei = cgd_data$weight; log_wei = log(wei)
#' 
#' trt = trt - 0.5
#' str = str - 0.5
#' age = scale(age);log_age = scale(log_age)
#' wei = scale(wei);log_wei = scale(log_wei)
#' 
#' result01_afttest_omni_mns=afttest(Surv(X_cgd,D_cgd)~
#' trt+str+age+wei,path=path,testtype="omni",eqType="mns")
#' afttestplot(result01_afttest_omni_mns,std="std")
#' 
#' @importFrom ggplot2 ggplot geom_step theme theme_minimal ggtitle ylab xlab aes element_text
#' @importFrom gridExtra grid.arrange
#' @importFrom stats quantile
#' 
#' @export
afttestplot=function(result,path=100,std="std"){
  
  testtype=result$TestType
  
  if(testtype=="Omni"){
    return(plotting_omni(result,path,std))
  }
  if(testtype=="Form"){
    return(plotting_form(result,path,std))
  }
  if(testtype=="Link"){
    return(plotting_link(result,path,std))
  }
}

