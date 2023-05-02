##############################################################################
## User's Main Function
##############################################################################

#' afttestplot
#'
#' @param object is a \code{afttest} fit
#' @param path A numeric value specifies the number of approximated processes plotted
#'    The default is set to be 100.
#' @param std A character string specifying if the graph is based on the 
#'    unstandardized test statistics or standardized test statistics
#'    The default is set to be "std".
#' @return \code{afttestplot} returns a plot based on the \code{testtype}:
#' \describe{
#'    \item{omni}{an object of the omnibus test is the form of n by n matrix, 
#'    some quantiles of x, which are used in weight, are plotted for graphs, 
#'    i.e. 0\%, 10\%, 25\%, 40\%, 50\%, 60\%, 75\%, 90\%, and 100\% are used.}
#'    \item{link}{an object of the link function test is the form of n by 1 matrix}
#'    \item{form}{an object of the functional form test is the form of n by 1 matrix}
#' }
#'    See the documentation of \pkg{ggplot2} and \pkg{gridExtra} for details.
#'    
#' @export
#' @example inst/examples/ex_afttestplot.R
afttestplot = function(object, path = 50, std = "std"){
  
  # class
  if (!inherits(object,"afttest")) stop("Must be afttest class")
  # testtype
  testtype = object$TestType
  # std
  if (!std %in% c("std","unstd")) {
    std = "std"
  }
  # path
  if (!is.numeric(path)) {
    path = 50
  } else {
    path = max(min(path,length(object$app_std_path)),10)
  }
  
  x_axis = 1:length(object$Resid)
  Q = c(0,0.1,0.25,0.4,0.5,0.6,0.75,0.9,1)
  Q = round(quantile(x_axis,Q))
  K = length(Q)
  
  if(testtype=="Omni"){
    resid = c(NA)
    app = matrix(NA)
    obs = matrix(NA)
    
    Figure = list(NA)
    if (std == "std") {
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
        if (k==5){
          Figure_k =
            ggplot() +
            geom_step(data=DF_app,aes(x=resid,y=app,group=group),colour="grey",alpha=0.5) +
            geom_step(data=DF_obs,aes(x=resid,y=obs),colour="tomato",lwd=0.25) +
            ylab("Test Statistic")+xlab("Residuals") +
            ggtitle(paste0("Standardized Omnibus: quantile(z): ",names(Q)[k])) + 
            scale_y_continuous(breaks = round(seq(min(c(DF_app$app,DF_obs$obs)), max(c(DF_app$app,DF_obs$obs)), length.out = 5),1)) +
            theme(plot.title=element_text(hjust=0.5))
        } else {
          Figure_k =
            ggplot() +
            geom_step(data=DF_app,aes(x=resid,y=app,group=group),colour="grey",alpha=0.5) +
            geom_step(data=DF_obs,aes(x=resid,y=obs),colour="tomato",lwd=0.25) +
            ylab("Test Statistic")+xlab("Residuals") +
            ggtitle(paste0("quantile(z): ",names(Q)[k])) + 
            scale_y_continuous(breaks = round(seq(min(c(DF_app$app,DF_obs$obs)), max(c(DF_app$app,DF_obs$obs)), length.out = 5),1)) +
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
        if (k==5){
          Figure_k =
            ggplot() +
            geom_step(data=DF_app,aes(x=resid,y=app,group=group),colour="grey",alpha=0.5) +
            geom_step(data=DF_obs,aes(x=resid,y=obs),colour="tomato",lwd=0.25) +
            ylab("Test Statistic")+xlab("Residuals") +
            ggtitle(paste0("Unstandardized Omnibus: quantile(z): ",names(Q)[k])) + 
            scale_y_continuous(breaks = round(seq(min(c(DF_app$app,DF_obs$obs)), max(c(DF_app$app,DF_obs$obs)), length.out = 5),1)) +
            theme(plot.title=element_text(hjust=0.5))
        } else {
          Figure_k =
            ggplot() +
            geom_step(data=DF_app,aes(x=resid,y=app,group=group),colour="grey",alpha=0.5) +
            geom_step(data=DF_obs,aes(x=resid,y=obs),colour="tomato",lwd=0.25) +
            ylab("Test Statistic")+xlab("Residuals") +
            ggtitle(paste0("quantile(z): ",names(Q)[k])) + 
            scale_y_continuous(breaks = round(seq(min(c(DF_app$app,DF_obs$obs)), max(c(DF_app$app,DF_obs$obs)), length.out = 5),1)) +
            theme(plot.title=element_text(hjust=0.5))
        }
        Figure[[k]] = Figure_k
      }
    }
    
    lay = rbind(c(1,1,1,1),c(1,1,1,1),c(2,3,4,5),c(6,7,8,9))
    return(grid.arrange(Figure[[5]],
                        Figure[[1]],Figure[[2]],Figure[[3]],Figure[[4]],
                        Figure[[6]],Figure[[7]],Figure[[8]],Figure[[9]],
                        layout_matrix=lay))
    
  } else if(testtype=="Link"){
    resid = c(NA)
    app = c(NA)
    obs = c(NA)
    if (std=="std"){
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
        ylab("Test Statistic")+xlab("Residuals")+ggtitle("Standardized Link Function") + 
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
        ylab("Test Statistic")+xlab("Residuals")+ggtitle("Untandardized Link Function") + 
        scale_y_continuous(breaks = round(seq(min(c(DF_app$app,DF_obs$obs)), max(c(DF_app$app,DF_obs$obs)), length.out = 5),1)) +
        theme(plot.title=element_text(hjust=0.5))
    }
    
    return(Figure)
    
  } else if(testtype=="Form"){
    resid = c(NA)
    app = c(NA)
    obs = c(NA)
    if (std=="std"){
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
        ylab("Test Statistic")+xlab("Residuals")+ggtitle("Standardized Functional Form") + 
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
        ylab("Test Statistic")+xlab("Residuals")+ggtitle("Untandardized Functional Form") + 
        scale_y_continuous(breaks = round(seq(min(c(DF_app$app,DF_obs$obs)), max(c(DF_app$app,DF_obs$obs)), length.out = 5))) +
        theme(plot.title=element_text(hjust=0.5))
    }
    
    return(Figure)
    
  } else {
    stop("Check your code")
  }
  
}