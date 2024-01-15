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
#'    See the documentation of \pkg{ggplot2} and \pkg{gridExtra} for details.\
#' 
#' @importFrom ggplot2 ggplot geom_step geom_label theme theme_minimal ggtitle aes unit labs ylab xlab scale_y_continuous element_text
#' @importFrom gridExtra grid.arrange
#' @importFrom stats quantile
#' 
#' @example inst/examples/ex_afttestplot.R
#' @export
afttestplot <- function(object, path = 50, stdType = "std", quantile = NULL){
  
  # class
  if (!inherits(object,"afttest")) return(warning("Must be afttest class"))
  # pathsave
  if ((object$pathsave<1)) return(warning("afttest is conduced with pathsave=0"))
  # eqType
  eqType <- object$eqType
  # testType
  testType <- object$testType
  # stdType
  if (!stdType %in% c("std","unstd")) {
    stdType <- "std"
  }
  # path
  if (length(path) > 1){
    return(warning("path needs to be an integer."))
  } else {
    if (!is.numeric(path)) {
      path <- 50
    } else {
      path <- min(path,object$pathsave)
    }
  }
  
  stdTypeQuote <- ifelse(stdType=="std","standardized","unstandardized")
  testTypeQuote <- ifelse(eqType=="mns","non-smooth","induced-smoothed")
  
  x_axis <- 1:nrow(object$DF)
  
  defaultQ <- c(0.1,0.25,0.5,0.75,0.9)
  lengthdefaultQ <- length(defaultQ)
  
  quantile <- sort(quantile)
  lengthquantile <- length(quantile)
  if (is.null(quantile)){
    Q <- defaultQ
  } else if (!is.numeric(quantile) || !lengthquantile == 5 || 
             min(quantile)<0 || max(quantile)>1) {
    return(warning("quantile needs to be numeric vector of 5 quantiles in [0,1]."))
  } else {
    Q <- quantile
  }
  Q <- round(stats::quantile(x_axis,Q))
  # names(Q) <- paste0(Q*100,"%")
  K <- length(Q)
  
  if(testType=="omni"){
    resid <- c(NA)
    app <- matrix(NA)
    obs <- matrix(NA)
    
    Figure <- list(NA)
    for(k in 1:K){
      if (stdType == "std") {
        # DF_app
        DF_app=data.frame()
        for (group in 1:path){
          temp <- object$app_std_path[[group]][,Q[k]]
          temp <- data.frame(group,resid=x_axis,app=temp)
          DF_app <- rbind(DF_app,temp)
        }
        # DF_obs
        DF_obs <- data.frame(group,resid=x_axis,obs=object$obs_std_path[,Q[k]])
        
      } else {
        #DF_app
        DF_app <- data.frame()
        for (group in 1:path){
          temp <- object$app_path[[group]][,Q[k]]
          temp <- data.frame(group,resid=x_axis,app=temp)
          DF_app <- rbind(DF_app,temp)
        }
        #DF_obs
        DF_obs <- data.frame(group,resid=x_axis,obs=object$obs_path[,Q[k]])
      }
      breaks <- c(DF_app$app,DF_obs$obs)
      breaks <- breaks[which(is.finite(breaks))]
      y_breaksMIN <- min(breaks, na.rm = TRUE)
      y_breaksMAX <- max(breaks, na.rm = TRUE)
      
      # Figure
      if (k==((K+1)/2)){
        y_breaks <- round(seq(y_breaksMIN, y_breaksMAX, length.out = 5),1)
        y_labels <- format(y_breaks, nsmall = 1)
        Figure_k <-
          ggplot() +
          geom_step(data=DF_app,aes(x=resid,y=app,group=group),colour="grey",alpha=0.5) +
          geom_step(data=DF_obs,aes(x=resid,y=obs),colour="tomato",lwd=0.25) +
          labs(y = NULL) + labs(x = NULL) +
          ggtitle(paste0("Omnibus test ","(",stdTypeQuote,")"," with ",names(Q[k]), " percentile for z")) +
          scale_y_continuous(breaks = y_breaks, labels = y_labels) +
          theme(plot.title=element_text(hjust=0.5),
                plot.margin = rep(unit(0,"null"),4),
                panel.spacing = unit(0,"null"))
      } else {
        y_breaks <- round(seq(y_breaksMIN, y_breaksMAX, length.out = 3),1)
        y_labels <- format(y_breaks, nsmall = 1)
        Figure_k <-
          ggplot() +
          geom_step(data=DF_app,aes(x=resid,y=app,group=group),colour="grey",alpha=0.5) +
          geom_step(data=DF_obs,aes(x=resid,y=obs),colour="tomato",lwd=0.25) +
          labs(y = NULL) + labs(x = NULL) + 
          scale_y_continuous(breaks = y_breaks, labels = y_labels) +
          geom_label(aes(x=-Inf,y=Inf),label=paste0(names(Q[k])),fill="darkgrey",label.size=NA,size=3,hjust=-0.0,vjust=1.0) +
          theme(plot.margin = rep(unit(0,"null"),4),
                panel.spacing = unit(0,"null")) 
      }
      Figure[[k]] <- Figure_k
    }
    
    lay <- rbind(c(1,1),c(1,1),c(2,3),c(4,5))
    return(gridExtra::grid.arrange(Figure[[3]],
                                   Figure[[1]],Figure[[2]],
                                   Figure[[4]],Figure[[5]],
                                   left = "Test Statistic",
                                   bottom = "Residuals",
                                   layout_matrix=lay))
    
  } else if(testType=="link"){
    resid <- c(NA)
    app <- c(NA)
    obs <- c(NA)
    if (stdType == "std"){
      # DF_app
      DF_app <- data.frame()
      for (group in 1:path){
        temp <- object$app_std_path[[group]]
        temp <- data.frame(group,resid=x_axis,app=temp)
        DF_app <- rbind(DF_app,temp)
      }
      # DF_obs
      DF_obs <- data.frame(group,resid=x_axis,obs=object$obs_std_path)
    } else {
      # DF_app
      DF_app <- data.frame()
      for (group in 1:path){
        temp <- object$app_path[[group]]
        temp <- data.frame(group,resid=x_axis,app=temp)
        DF_app <- rbind(DF_app,temp)
      }
      # DF_obs
      DF_obs <- data.frame(group,resid=x_axis,obs=object$obs_path)
    }
    breaks <- c(DF_app$app,DF_obs$obs)
    breaks <- breaks[which(is.finite(breaks))]
    y_breaksMIN <- min(breaks, na.rm = TRUE)
    y_breaksMAX <- max(breaks, na.rm = TRUE)
    
    # Figure
    Figure <-
      ggplot() +
      geom_step(data=DF_app,aes(x=resid,y=app,group=group),colour="grey",alpha=0.5) +
      geom_step(data=DF_obs,aes(x=resid,y=obs),colour="tomato",lwd=0.25) +
      ylab("Test Statistic")+xlab("Residuals")+
      ggtitle(paste0("Link Function ","(",stdTypeQuote,")","")) + 
      scale_y_continuous(breaks = round(seq(y_breaksMIN, y_breaksMAX, length.out = 5),1)) +
      theme(plot.title=element_text(hjust=0.5))
    
    return(Figure)
    
  } else if(testType=="form"){
    resid <- c(NA)
    app <- c(NA)
    obs <- c(NA)
    if (stdType == "std"){
      # DF_app
      DF_app <- data.frame()
      for (group in 1:path){
        temp <- object$app_std_path[[group]]
        temp <- data.frame(group,resid=x_axis,app=temp)
        DF_app <- rbind(DF_app,temp)
      }
      # DF_obs
      DF_obs <- data.frame(group,resid=x_axis,obs=object$obs_std_path)
      
    } else {
      # DF_app
      DF_app <- data.frame()
      for (group in 1:path){
        temp <- object$app_path[[group]]
        temp <- data.frame(group,resid=x_axis,app=temp)
        DF_app <- rbind(DF_app,temp)
      }
      # DF_obs
      DF_obs <- data.frame(group,resid=x_axis,obs=object$obs_path)
      
    }
    breaks <- c(DF_app$app,DF_obs$obs)
    breaks <- breaks[which(is.finite(breaks))]
    y_breaksMIN <- min(breaks, na.rm = TRUE)
    y_breaksMAX <- max(breaks, na.rm = TRUE)
    
    # Figure
    Figure <-
      ggplot() +
      geom_step(data=DF_app,aes(x=resid,y=app,group=group),colour="grey",alpha=0.5) +
      geom_step(data=DF_obs,aes(x=resid,y=obs),colour="tomato",lwd=0.25) +
      ylab("Test Statistic")+xlab("Residuals") +
      ggtitle(paste0("Functional Form ","(",stdTypeQuote,")","")) + 
      scale_y_continuous(breaks = round(seq(y_breaksMIN, y_breaksMAX, length.out = 5),1)) +
      theme(plot.title=element_text(hjust=0.5))
    
    return(Figure)
    
  } else {
    return(warning("Check your code"))
  }
}