##############################################################################
## User's Main Function
##############################################################################
#' afttestplot
#'
#' It gives plot for checking the aft model assumptions.
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
afttestplot=function(object,path=100,std="std"){
  
  testtype=object$TestType
  
  if(testtype=="Omni"){
    return(afttestplot_omni(object,path,std))
  }
  if(testtype=="Form"){
    return(afttestplot_form(object,path,std))
  }
  if(testtype=="Link"){
    return(afttestplot_link(object,path,std))
  }
}

##############################################################################
## Background functions
##############################################################################
#' afttestplot_omni
#'
#' It gives plot for cheking the aft model assumptions.
#' @param object is a \code{afttest} fit
#' @param path A numeric value specifies the number of approximated processes plotted
#'    The default is set to be 100.
#' @param std A character string specifying if the graph is based on the 
#'    unstandardized test statistics or standardized test statistics
#'    The default is set to be "std".
#' @return \code{afttestplot_omni} returns a plot based on the \code{omni}:
#'    An object of the omnibus test is the form of n by n matrix, 
#'    some quantiles of x, which are used in weight, are plotted for graphs, 
#'    i.e. 0\%, 10\%, 25\%, 40\%, 50\%, 60\%, 75\%, 90\%, and 100\% are used.
#'    See the documentation of \pkg{ggplot2} and \pkg{gridExtra} for details.
#' 
#' @keywords internal
afttestplot_omni=function(object,path,std){
  
  e_i = c(NA)
  std_What = matrix(NA)
  std_W = matrix(NA)
  What = matrix(NA)
  W = matrix(NA)
  
  object_resid=object$Resid
  n=length(object_resid)
  xaxix=(1:n)
  
  pathsave = length(object$app_std_path)
  if (pathsave<path){
    path = pathsave
  }
  
  # if (xaxix=="rank"){xaxix=(1:n)[order(object_Time)]}
  # else if (xaxix=="real"){xaxix=object_Time}
  # else (xaxix=object_Time)
  if (std=="std"){
    Figure=list(NA)
    
    for(k in 1:9){
      
      quant=round(quantile(1:n,c(0,0.1,0.25,0.4,0.5,0.6,0.75,0.9,1)))
      
      Q=quant[k]
      
      dataset_std_What=data.frame()
      
      for (i in 1:path){
        group=i
        A=object$app_std_path[[i]][,Q]
        AA=data.frame(group,e_i=xaxix,std_What=A)
        dataset_std_What=rbind(dataset_std_What,AA)
      }
      #dataset_std_What
      
      dataset_std_W=data.frame(group,e_i=xaxix,std_W=object$obs_std_path[,Q])
      #dataset_std_W
      
      Figure_std_W=
        ggplot()+
        geom_step(data=dataset_std_What,aes(x=e_i,y=std_What,group=group),colour="grey",alpha=0.5)+
        geom_step(data=dataset_std_W,aes(x=e_i,y=std_W),colour="tomato",lwd=0.25)+
        ylab("Test Statistic")+xlab("Residuals")+
        ggtitle(paste("Quantile of z",names(quant)[k]))+theme(plot.title=element_text(hjust=0.5))
      #Figure_W
      Figure[[k]]=Figure_std_W
    }
    
    # Figure[[1]]
    
    lay=rbind(c(1,1,1,1),c(1,1,1,1),c(2,3,4,5),c(6,7,8,9))
    
    return(grid.arrange(Figure[[5]],Figure[[1]],Figure[[2]],Figure[[3]],Figure[[4]],
                        Figure[[6]],Figure[[7]],Figure[[8]],Figure[[9]],layout_matrix=lay))
  } else {
    Figure=list(NA)
    
    for(k in 1:9){
      
      quant=round(quantile(1:n,c(0,0.1,0.25,0.4,0.5,0.6,0.75,0.9,1)))
      
      Q=quant[k]
      
      dataset_What=data.frame()
      
      for (i in 1:path){
        group=i
        A=object$app_path[[i]][,Q]
        AA=data.frame(group,e_i=xaxix,What=A)
        dataset_What=rbind(dataset_What,AA)
      }
      #dataset_What
      
      dataset_W=data.frame(group,e_i=xaxix,W=object$obs_path[,Q])
      #dataset_W
      
      Figure_W=
        ggplot()+
        geom_step(data=dataset_What,aes(x=e_i,y=What,group=group),colour="grey",alpha=0.5)+
        geom_step(data=dataset_W,aes(x=e_i,y=W),colour="tomato",lwd=0.25)+
        ylab("Test Statistic")+xlab("Residuals")+
        ggtitle(paste("Quantile of z",names(quant)[k]))+theme(plot.title=element_text(hjust=0.5))
      #Figure_W
      Figure[[k]]=Figure_W
    }
    
    # Figure[[1]]
    
    lay=rbind(c(1,1,1,1),c(1,1,1,1),c(2,3,4,5),c(6,7,8,9))
    
    return(grid.arrange(Figure[[5]],Figure[[1]],Figure[[2]],Figure[[3]],Figure[[4]],
                        Figure[[6]],Figure[[7]],Figure[[8]],Figure[[9]],layout_matrix=lay))
  }
}

#' afttestplot_link
#'
#' It gives plot for cheking the aft model assumptions.
#' @param object is a \code{afttest} fit
#' @param path A numeric value specifies the number of approximated processes plotted
#'    The default is set to be 100.
#' @param std A character string specifying if the graph is based on the 
#'    unstandardized test statistics or standardized test statistics
#'    The default is set to be "std".
#' @return \code{afttestplot_link} returns a plot based on the \code{link}:
#'    An object of the functional form test is the form of n by 1 matrix.
#'    See the documentation of \pkg{ggplot2} and \pkg{gridExtra} for details.
#' 
#' @keywords internal
afttestplot_link=function(object,path,std){
  
  e_i = c(NA)
  std_What = c(NA)
  std_W = c(NA)
  What = c(NA)
  W = c(NA)
  
  n=length(object$Time)
  xaxix=(1:n)
  
  pathsave = length(object$app_std_path)
  if (pathsave<path){
    path = pathsave
  }
  # object_Covari=object$Covari
  # if (xaxix=="rank"){xaxix=(1:n)[order(object_Covari)]}
  # else {xaxix=object_Covari}
  if (std=="std"){
    
    dataset_std_What=data.frame()
    
    for (i in 1:path){
      group=i
      A=object$app_std_path[[i]]
      AA=data.frame(group,e_i=xaxix,std_What=A)
      dataset_std_What=rbind(dataset_std_What,AA)
    }
    #dataset_std_What
    
    dataset_std_W=data.frame(group,e_i=xaxix,std_W=object$obs_std_path)
    #dataset_std_W
    
    Figure_std_W=
      ggplot()+
      geom_step(data=dataset_std_What,aes(x=e_i,y=std_What,group=group),colour="grey",alpha=0.5)+
      geom_step(data=dataset_std_W,aes(x=e_i,y=std_W),colour="tomato",lwd=0.25)+
      ylab("Test Statistic")+xlab("Residuals")+
      theme_minimal()
    #Figure_std_W
    
    return(Figure_std_W)
    
  } else {
    
    dataset_What=data.frame()
    
    for (i in 1:path){
      group=i
      A=object$app_path[[i]]
      AA=data.frame(group,e_i=xaxix,What=A)
      dataset_What=rbind(dataset_What,AA)
    }
    #dataset_What
    
    dataset_W=data.frame(group,e_i=xaxix,W=object$obs_path)
    #dataset_W
    
    Figure_W=
      ggplot()+
      geom_step(data=dataset_What,aes(x=e_i,y=What,group=group),colour="grey",alpha=0.5)+
      geom_step(data=dataset_W,aes(x=e_i,y=W),colour="tomato",lwd=0.25)+
      ylab("Test Statistic")+xlab("Residuals")+
      theme_minimal()
    #Figure_W
    
    return(Figure_W)
    
  }
}

#' afttestplot_form
#'
#' It gives plot for cheking the aft model assumptions.
#' @param object is a \code{afttest} fit
#' @param path A numeric value specifies the number of approximated processes plotted
#'    The default is set to be 100.
#' @param std A character string specifying if the graph is based on the 
#'    unstandardized test statistics or standardized test statistics
#'    The default is set to be "std".
#' @return \code{afttestplot_form} returns a plot based on the \code{form}:
#'    An object of the functional form test is the form of n by 1 matrix.
#'    See the documentation of \pkg{ggplot2} and \pkg{gridExtra} for details.
#' 
#' @keywords internal
afttestplot_form=function(object,path,std){
  
  e_i = c(NA)
  std_What = c(NA)
  std_W = c(NA)
  What = c(NA)
  W = c(NA)
  
  n=length(object$Time)
  xaxix=(1:n)
  
  pathsave = length(object$app_std_path)
  if (pathsave<path){
    path = pathsave
  }
  # object_Covari=object$Covari
  # if (xaxix=="rank"){xaxix=(1:n)[order(object_Covari)]}
  # else {xaxix=object_Covari}
  if (std=="std"){
    
    dataset_std_What=data.frame()
    
    for (i in 1:path){
      group=i
      A=object$app_std_path[[i]]
      AA=data.frame(group,e_i=xaxix,std_What=A)
      dataset_std_What=rbind(dataset_std_What,AA)
    }
    #dataset_std_What
    
    dataset_std_W=data.frame(group,e_i=xaxix,std_W=object$obs_std_path)
    #dataset_std_W
    
    Figure_std_W=
      ggplot()+
      geom_step(data=dataset_std_What,aes(x=e_i,y=std_What,group=group),colour="grey",alpha=0.5)+
      geom_step(data=dataset_std_W,aes(x=e_i,y=std_W),colour="tomato",lwd=0.25)+
      ylab("Test Statistic")+xlab("Residuals")+
      theme_minimal()
    #Figure_std_W
    
    return(Figure_std_W)
    
  } else {
    
    dataset_What=data.frame()
    
    for (i in 1:path){
      group=i
      A=object$app_path[[i]]
      AA=data.frame(group,e_i=xaxix,What=A)
      dataset_What=rbind(dataset_What,AA)
    }
    #dataset_What
    
    dataset_W=data.frame(group,e_i=xaxix,W=object$obs_path)
    #dataset_W
    
    Figure_W=
      ggplot()+
      geom_step(data=dataset_What,aes(x=e_i,y=What,group=group),colour="grey",alpha=0.5)+
      geom_step(data=dataset_W,aes(x=e_i,y=W),colour="tomato",lwd=0.25)+
      ylab("Test Statistic")+xlab("Residuals")+
      theme_minimal()
    #Figure_W
    
    return(Figure_W)
    
  }
}