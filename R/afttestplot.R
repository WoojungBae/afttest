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
#' the number of paths used in by afttest function. The default option for
#' argument path is set to be 100.
#' @param std The option for the argument std is "unstd" and "std". In this 
#' argument, "unstd" is the default.
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
#' # result_link = afttest(Surv(X_cgd,D_cgd)~trt+str+age,path=path,testtype="link")
#' #   
#' # graph_link_unstd = afttestplot(result_link,path=10,std="unstd")
#' # graph_link_std = afttestplot(result_link,path=10,std="std")
#' # 
#' # result_form = afttest(Surv(X_cgd,D_cgd)~trt+str+age,path=path,testtype="form")
#' # 
#' # graph_form_unstd = afttestplot(result_form,path=10,std="unstd")
#' # graph_form_std = afttestplot(result_form,path=10,std="std")
#' 
#' @importFrom ggplot2 ggplot geom_step theme theme_minimal ggtitle ylab xlab aes element_text
#' @importFrom gridExtra grid.arrange
#' @importFrom stats quantile
#' 
#' @export
afttestplot=function(result,path=100,std="unstd"){
  
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

