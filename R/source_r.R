#' plotting_omni
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
#' @importFrom ggplot2 ggplot geom_step theme theme_minimal ggtitle ylab xlab aes element_text
#' @importFrom gridExtra grid.arrange
#' @importFrom stats quantile
#' 
plotting_omni=function(result,path,std){
  
  e_i = c(NA)
  std_What = matrix(NA)
  std_W = matrix(NA)
  What = matrix(NA)
  W = matrix(NA)
  
  result_resid=result$Resid
  n=length(result_resid)
  xaxix=(1:n)
  
  # if (xaxix=="rank"){xaxix=(1:n)[order(result_Time)]}
  # else if (xaxix=="real"){xaxix=result_Time}
  # else (xaxix=result_Time)
  if (std=="std"){
    Figure=list(NA)
    
    for(k in 1:9){
      
      quant=round(quantile(1:n,c(0,0.1,0.25,0.4,0.5,0.6,0.75,0.9,1)))
      
      Q=quant[k]
      
      dataset_std_What=data.frame()
      
      for (i in 1:path){
        group=i
        A=result$app_std_path[[i]][,Q]
        AA=data.frame(group,e_i=xaxix,std_What=A)
        dataset_std_What=rbind(dataset_std_What,AA)
      }
      #dataset_std_What
      
      dataset_std_W=data.frame(group,e_i=xaxix,std_W=result$obs_std_path[,Q])
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
        A=result$app_path[[i]][,Q]
        AA=data.frame(group,e_i=xaxix,What=A)
        dataset_What=rbind(dataset_What,AA)
      }
      #dataset_What
      
      dataset_W=data.frame(group,e_i=xaxix,W=result$obs_path[,Q])
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

#' plotting_link
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
#' ordered by time transformed residual.
#' 
#' @importFrom ggplot2 ggplot geom_step theme theme_minimal ggtitle ylab xlab aes
#' @importFrom gridExtra grid.arrange
#' 
plotting_link=function(result,path,std){
  
  e_i = c(NA)
  std_What = c(NA)
  std_W = c(NA)
  What = c(NA)
  W = c(NA)
  
  n=length(result$Time)
  xaxix=(1:n)
  
  # result_Covari=result$Covari
  # if (xaxix=="rank"){xaxix=(1:n)[order(result_Covari)]}
  # else {xaxix=result_Covari}
  
  if (std=="std"){
    
    dataset_std_What=data.frame()
    
    for (i in 1:path){
      group=i
      A=result$app_std_path[[i]]
      AA=data.frame(group,e_i=xaxix,std_What=A)
      dataset_std_What=rbind(dataset_std_What,AA)
    }
    #dataset_std_What
    
    dataset_std_W=data.frame(group,e_i=xaxix,std_W=result$obs_std_path)
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
      A=result$app_path[[i]]
      AA=data.frame(group,e_i=xaxix,What=A)
      dataset_What=rbind(dataset_What,AA)
    }
    #dataset_What
    
    dataset_W=data.frame(group,e_i=xaxix,W=result$obs_path)
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

#' plotting_form
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
#' ordered by time transformed residual.
#' 
#' @include source_r.R
#' 
#' @importFrom ggplot2 ggplot geom_step theme theme_minimal ggtitle ylab xlab aes
#' @importFrom gridExtra grid.arrange
#' 
plotting_form=function(result,path,std){
  
  e_i = c(NA)
  std_What = c(NA)
  std_W = c(NA)
  What = c(NA)
  W = c(NA)
  
  n=length(result$Time)
  xaxix=(1:n)
  
  # result_Covari=result$Covari
  # if (xaxix=="rank"){xaxix=(1:n)[order(result_Covari)]}
  # else {xaxix=result_Covari}
  
  if (std=="std"){
    
    dataset_std_What=data.frame()
    
    for (i in 1:path){
      group=i
      A=result$app_std_path[[i]]
      AA=data.frame(group,e_i=xaxix,std_What=A)
      dataset_std_What=rbind(dataset_std_What,AA)
    }
    #dataset_std_What
    
    dataset_std_W=data.frame(group,e_i=xaxix,std_W=result$obs_std_path)
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
      A=result$app_path[[i]]
      AA=data.frame(group,e_i=xaxix,What=A)
      dataset_What=rbind(dataset_What,AA)
    }
    #dataset_What
    
    dataset_W=data.frame(group,e_i=xaxix,W=result$obs_path)
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
