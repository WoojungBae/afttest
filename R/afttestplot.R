##############################################################################
## User's Main Function
##############################################################################

#' afttestplot
#'
#' @param object is a \code{afttest} fit
#' @param path A numeric value specifies the number of approximated processes plotted
#'    The default is set to be 100.
#' @param stdType A character string specifying if the graph is based on the 
#'    unstandardized test statistics or standardized test statistics
#'    The default is set to be "std".
#' @return \code{afttestplot} returns a plot based on the \code{testType}:
#' \describe{
#'    \item{omni}{an object of the omnibus test is the form of n by n matrix, 
#'    some quantiles of x, which are used in weight, are plotted for graphs, 
#'    i.e. 10\%, 25\%, 50\%, 75\%, and 90\% are used.}
#'    \item{link}{an object of the link function test is the form of n by 1 matrix}
#'    \item{form}{an object of the functional form test is the form of n by 1 matrix}
#' }
#'    See the documentation of \pkg{ggplot2} and \pkg{gridExtra} for details.
#'    
#' @example inst/examples/ex_afttestplot.R
#' @export
##############################################################################
## User's Main Function
##############################################################################

#' afttestplot
#'
#' @param object is a \code{afttest} fit
#' @param path A numeric value specifies the number of approximated processes plotted
#'    The default is set to be 100.
#' @param stdType A character string specifying if the graph is based on the 
#'    unstandardized test statistics or standardized test statistics
#'    The default is set to be "std".
#' @param quantile A numeric vector specifies 5 of five quantiles within the range [0,1]. 
#'    The default is set to be c(0.1,0.25,0.5,0.75,0.9).
#' @return \code{afttestplot} returns a plot based on the \code{testType}:
#' \describe{
#'    \item{omni}{an object of the omnibus test is the form of n by n matrix, 
#'    some quantiles of x, which are used in weight, are plotted for graphs, 
#'    i.e. 0\%, 10\%, 25\%, 40\%, 50\%, 60\%, 75\%, 90\%, and 100\% are used.}
#'    \item{link}{an object of the link function test is the form of n by 1 matrix}
#'    \item{form}{an object of the functional form test is the form of n by 1 matrix}
#' }
#'    See the documentation of \pkg{ggplot2} and \pkg{gridExtra} for details.
#'    
#' @example inst/examples/ex_afttestplot.R
#' @export
afttestplot = function(object, path = 50, stdType = "std", quantile = NULL){
  
  # class
  if (!inherits(object,"afttest")) return(warning("Must be afttest class"))
  # eqType
  eqType = object$eqType
  # testType
  testType = object$testType
  # stdType
  if (!stdType %in% c("std","unstd")) {
    stdType = "std"
  }
  # path
  if (!is.numeric(path)) {
    path = 50
  } else {
    path = max(min(path,length(object$app_std_path)),10)
  }
  
  stdTypeQuote = ifelse(stdType=="std","standardized","unstandardized")
  testTypeQuote = ifelse(eqType=="mns","non-smooth","induced-smoothed")
  
  x_axis = 1:nrow(object$DF)
  
  defaultQ = c(0.1,0.25,0.5,0.75,0.9)
  lengthdefaultQ = length(defaultQ)
  
  quantile = sort(quantile)
  lengthquantile = length(quantile)
  if (is.null(quantile)){
    Q = defaultQ
  } else if (!is.numeric(quantile) || !lengthquantile == 5 || 
             min(quantile)<0 || max(quantile)>1) {
    return(warning("quantile needs to be numeric vector of 5 quantiles in [0,1]."))
  } else {
    Q = quantile
  }
  Q = round(quantile(x_axis,Q))
  K = length(Q)
  
  if(testType=="omni"){
    resid = c(NA)
    app = matrix(NA)
    obs = matrix(NA)
    
    Figure = list(NA)
    if (stdType == "std") {
      for(k in 1:K){
        Q_k = Q[k]
        
        # DF_app
        DF_app=data.frame()
        for (group in 1:path){
          temp = object$app_std_path[[group]][,Q_k]
          temp = data.frame(group,resid=x_axis,app=temp)
          DF_app = rbind(DF_app,temp)
        }
        
        # DF_obs
        DF_obs = data.frame(group,resid=x_axis,obs=object$obs_std_path[,Q_k])
        
        # Figure
        if (k==((K+1)/2)){
          Figure_k =
            ggplot() +
            geom_step(data=DF_app,aes(x=resid,y=app,group=group),colour="grey",alpha=0.5) +
            geom_step(data=DF_obs,aes(x=resid,y=obs),colour="tomato",lwd=0.25) +
            # ylab("Test Statistic") + xlab("Residuals") +
            labs(y = NULL) + labs(x = NULL) +
            ggtitle(paste0("Omnibus test ","(",stdTypeQuote,")"," with ",names(Q)[k], " percentile for z")) +
            # ggtitle(paste0("Omnibus test ","(",stdTypeQuote,", ",testTypeQuote,")"," with ",names(Q)[k], " percentile for z")) + 
            scale_y_continuous(breaks = round(seq(min(c(DF_app$app,DF_obs$obs)), max(c(DF_app$app,DF_obs$obs)), length.out = 5),1)) +
            theme(plot.title=element_text(hjust=0.5))
        } else {
          Figure_k =
            ggplot() +
            geom_step(data=DF_app,aes(x=resid,y=app,group=group),colour="grey",alpha=0.5) +
            geom_step(data=DF_obs,aes(x=resid,y=obs),colour="tomato",lwd=0.25) +
            # ylab("Test Statistic") + xlab("Residuals") +
            labs(y = NULL) + labs(x = NULL) +
            ggtitle(paste0(names(Q)[k])) + 
            scale_y_continuous(breaks = round(seq(min(c(DF_app$app,DF_obs$obs)), max(c(DF_app$app,DF_obs$obs)), length.out = 3),1)) +
            theme(plot.title=element_text(hjust=0.5))
        } 
        Figure[[k]] = Figure_k
      }
    } else {
      for(k in 1:K){
        Q_k = Q[k]
        
        #DF_app
        DF_app = data.frame()
        for (group in 1:path){
          temp = object$app_path[[group]][,Q_k]
          temp = data.frame(group,resid=x_axis,app=temp)
          DF_app = rbind(DF_app,temp)
        }
        
        #DF_obs
        DF_obs = data.frame(group,resid=x_axis,obs=object$obs_path[,Q_k])
        
        # Figure
        if (k==((K+1)/2)){
          Figure_k =
            ggplot() +
            geom_step(data=DF_app,aes(x=resid,y=app,group=group),colour="grey",alpha=0.5) +
            geom_step(data=DF_obs,aes(x=resid,y=obs),colour="tomato",lwd=0.25) +
            # ylab("Test Statistic") + xlab("Residuals") +
            labs(y = NULL) + labs(x = NULL) +
            ggtitle(paste0("Omnibus test ","(",stdTypeQuote,")"," with ",names(Q)[k], " percentile for z")) + 
            # ggtitle(paste0("Omnibus test ","(",stdTypeQuote,", ",testTypeQuote,")"," with ",names(Q)[k], " percentile for z")) + 
            scale_y_continuous(breaks = round(seq(min(c(DF_app$app,DF_obs$obs)), max(c(DF_app$app,DF_obs$obs)), length.out = 5),1)) +
            theme(plot.title=element_text(hjust=0.5))
        } else {
          Figure_k =
            ggplot() +
            geom_step(data=DF_app,aes(x=resid,y=app,group=group),colour="grey",alpha=0.5) +
            geom_step(data=DF_obs,aes(x=resid,y=obs),colour="tomato",lwd=0.25) +
            # ylab("Test Statistic") + xlab("Residuals") +
            labs(y = NULL) + labs(x = NULL) +
            ggtitle(paste0(names(Q)[k])) + 
            scale_y_continuous(breaks = round(seq(min(c(DF_app$app,DF_obs$obs)), max(c(DF_app$app,DF_obs$obs)), length.out = 3),1)) +
            theme(plot.title=element_text(hjust=0.5))
        }
        Figure[[k]] = Figure_k
      }
    }
    
    lay = rbind(c(1,1),c(1,1),c(2,3),c(4,5))
    return(grid.arrange(Figure[[3]],
                        Figure[[1]],Figure[[2]],
                        Figure[[4]],Figure[[5]],
                        left = "Test Statistic",
                        bottom = "Residuals",
                        layout_matrix=lay))
    
    # lay = rbind(c(1,1,1,1),c(1,1,1,1),c(2,3,4,5))
    # return(grid.arrange(Figure[[3]],
    #                     Figure[[1]],Figure[[2]],
    #                     Figure[[4]],Figure[[5]],
    #                     layout_matrix=lay))
    
  } else if(testType=="link"){
    resid = c(NA)
    app = c(NA)
    obs = c(NA)
    if (stdType == "std"){
      # DF_app
      DF_app = data.frame()
      for (group in 1:path){
        temp = object$app_std_path[[group]]
        temp = data.frame(group,resid=x_axis,app=temp)
        DF_app = rbind(DF_app,temp)
      }
      
      # DF_obs
      DF_obs = data.frame(group,resid=x_axis,obs=object$obs_std_path)
      
      # Figure
      Figure =
        ggplot() +
        geom_step(data=DF_app,aes(x=resid,y=app,group=group),colour="grey",alpha=0.5) +
        geom_step(data=DF_obs,aes(x=resid,y=obs),colour="tomato",lwd=0.25) +
        ylab("Test Statistic")+xlab("Residuals")+
        ggtitle(paste0("Link Function ","(",stdTypeQuote,")","")) + 
        # ggtitle(paste0("Link Function ","(",stdTypeQuote,", ",testTypeQuote,")","")) + 
        scale_y_continuous(breaks = round(seq(min(c(DF_app$app,DF_obs$obs)), max(c(DF_app$app,DF_obs$obs)), length.out = 5),1)) +
        theme(plot.title=element_text(hjust=0.5))
    } else {
      # DF_app
      DF_app = data.frame()
      for (group in 1:path){
        temp = object$app_path[[group]]
        temp = data.frame(group,resid=x_axis,app=temp)
        DF_app = rbind(DF_app,temp)
      }
      
      # DF_obs
      DF_obs = data.frame(group,resid=x_axis,obs=object$obs_path)
      
      # Figure
      Figure =
        ggplot() +
        geom_step(data=DF_app,aes(x=resid,y=app,group=group),colour="grey",alpha=0.5) +
        geom_step(data=DF_obs,aes(x=resid,y=obs),colour="tomato",lwd=0.25) +
        ylab("Test Statistic")+xlab("Residuals")+
        ggtitle(paste0("Link Function ","(",stdTypeQuote,")","")) + 
        # ggtitle(paste0("Link Function ","(",stdTypeQuote,", ",testTypeQuote,")","")) + 
        scale_y_continuous(breaks = round(seq(min(c(DF_app$app,DF_obs$obs)), max(c(DF_app$app,DF_obs$obs)), length.out = 5),1)) +
        theme(plot.title=element_text(hjust=0.5))
    }
    
    return(Figure)
    
  } else if(testType=="form"){
    resid = c(NA)
    app = c(NA)
    obs = c(NA)
    if (stdType == "std"){
      # DF_app
      DF_app = data.frame()
      for (group in 1:path){
        temp = object$app_std_path[[group]]
        temp = data.frame(group,resid=x_axis,app=temp)
        DF_app = rbind(DF_app,temp)
      }
      
      # DF_obs
      DF_obs = data.frame(group,resid=x_axis,obs=object$obs_std_path)
      
      # Figure
      Figure =
        ggplot() +
        geom_step(data=DF_app,aes(x=resid,y=app,group=group),colour="grey",alpha=0.5) +
        geom_step(data=DF_obs,aes(x=resid,y=obs),colour="tomato",lwd=0.25) +
        ylab("Test Statistic")+xlab("Residuals") +
        ggtitle(paste0("Functional Form ","(",stdTypeQuote,")","")) + 
        # ggtitle(paste0("Functional Form ","(",stdTypeQuote,", ",testTypeQuote,")","")) + 
        scale_y_continuous(breaks = round(seq(min(c(DF_app$app,DF_obs$obs)), max(c(DF_app$app,DF_obs$obs)), length.out = 5),1)) +
        theme(plot.title=element_text(hjust=0.5))
    } else {
      # DF_app
      DF_app = data.frame()
      for (group in 1:path){
        temp = object$app_path[[group]]
        temp = data.frame(group,resid=x_axis,app=temp)
        DF_app = rbind(DF_app,temp)
      }
      
      # DF_obs
      DF_obs = data.frame(group,resid=x_axis,obs=object$obs_path)
      
      # Figure
      Figure =
        ggplot() +
        geom_step(data=DF_app,aes(x=resid,y=app,group=group),colour="grey",alpha=0.5) +
        geom_step(data=DF_obs,aes(x=resid,y=obs),colour="tomato",lwd=0.25) +
        ylab("Test Statistic")+xlab("Residuals") +
        ggtitle(paste0("Functional Form ","(",stdTypeQuote,")","")) + 
        # ggtitle(paste0("Functional Form ","(",stdTypeQuote,", ",testTypeQuote,")",")","")) + 
        scale_y_continuous(breaks = round(seq(min(c(DF_app$app,DF_obs$obs)), max(c(DF_app$app,DF_obs$obs)), length.out = 5))) +
        theme(plot.title=element_text(hjust=0.5))
    }
    
    return(Figure)
    
  } else {
    return(warning("Check your code"))
  }
}