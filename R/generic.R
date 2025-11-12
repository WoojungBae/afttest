##############################################################################
## print
##############################################################################

#' print.afttest
#'
#' @param x is a \code{afttest} fit.
#' @param ... other options.
#' @return \code{print.afttest} returns a summary of a \code{afttest} fit:
#'
#' @example inst/examples/ex_afttest.R
#' @export
print.afttest <- function(x, ...) {
  if (!inherits(x,"afttest")) stop("Must be afttest class")
  
  cat("Call: \n")
  print(x$call)
  
  cat("\n p-values: \n")
  p_valueTAB <- data.frame(t(c(x$p_value, x$p_std_value)))
  rownames(p_valueTAB) <- ""
  colnames(p_valueTAB) <- c("unstandardized", "standardized")
  print(p_valueTAB)
  
  if (x$testType == "omnibus") {
    cat(paste0("\n H0: The assumed semiparametric AFT model fits the data adequately. \n"))
  } else if (x$testType == "link") {
    cat(paste0("\n H0: The relationship between covariates and the log survival time "))
    cat(paste0("\n     is correctly specified. \n"))
  } else if (x$testType == "covform") {
    cat(paste0("\n H0: The relationship between the log survival time and the specific "))
    cat(paste0("\n     covariate chosen by the cov.tested argument is correctly specified. \n"))
  }
  
  # if (x$estMethod=="ls"){
  #   cat(paste0("\n Coefficients (estimated by aftgee::aftgee):"))
  # } else if (x$estMethod=="rr"){
  #   cat(paste0("\n Coefficients (estimated by aftgee::aftsrr):"))
  # }
  # coefTAB <- data.frame(t(-x$beta))
  # rownames(coefTAB) <- ""
  # colnames(coefTAB) <- x$names[-c(1:2)]
  # print(coefTAB)
  
  invisible(x)
}
# print.afttest <- function(x, ...) {
#   if (!inherits(x,"afttest")) stop("Must be afttest class")
#   
#   print.default(x)
#   invisible(x)
# }

##############################################################################
## summary
##############################################################################

#' summary.afttest
#'
#' @param object is a \code{afttest} fit.
#' @param ... other options.
#' @return \code{summary.afttest} returns a summary of a \code{afttest} fit:
#'
#' @example inst/examples/ex_afttest.R
#' @export
summary.afttest <- function(object, ...) {
  if (!inherits(object,"afttest")) stop("Must be afttest class")
  
  cat("Call: \n")
  print(object$call)
  
  cat("\n p-values: \n")
  p.valueTAB <- data.frame(t(c(object$p_value, object$p_std_value)))
  rownames(p.valueTAB) <- ""
  colnames(p.valueTAB) <- c("unstandardized", "standardized")
  print(p.valueTAB)
  
  if (object$testType == "omnibus") {
    cat(paste0("\n H0: The assumed semiparametric AFT model fits the data adequately. \n"))
  } else if (object$testType == "link") {
    cat(paste0("\n H0: The relationship between covariates and the log survival time "))
    cat(paste0("\n     is correctly specified. \n"))
  } else if (object$testType == "covform") {
    cat(paste0("\n H0: The relationship between the log survival time and the specific "))
    cat(paste0("\n     covariate chosen by the cov.tested argument is correctly specified. \n"))
  }
  
  # if (object$estMethod=="ls"){
  #   cat(paste0("\n Coefficients (estimated by aftgee::aftgee):"))
  # } else if (object$estMethod=="rr"){
  #   cat(paste0("\n Coefficients (estimated by aftgee::aftsrr):"))
  # }
  # coefTAB <- data.frame(t(-object$beta))
  # rownames(coefTAB) <- ""
  # colnames(coefTAB) <- object$names[-c(1:2)]
  # print(coefTAB)
  
  invisible(object)
}

##############################################################################
## plot
##############################################################################

#' plot.afttest
#'
#' @param x is a \code{afttest} fit
#' @param npath A numeric value specifies the number of approximated processes plotted.
#'    The default is set to be 100.
#' @param std A character string specifying if the graph is based on 
#'    the unstandardized test statistics or standardized test statistics
#'    The default is set to be "std".
#' @param quantile A numeric vector specifies 5 of five quantiles within the range [0,1]. 
#'    The default is set to be c(0.1,0.25,0.5,0.75,0.9).
#' @param ... for future extension
#' @return \code{plot.afttest} returns a plot based on the \code{testType}:
#' \describe{
#'    \item{omnibus}{an x of the omnibus test is the form of n by n matrix, 
#'    some quantiles of x, which are used in weight, are plotted for graphs, 
#'    i.e. 0\%, 10\%, 25\%, 40\%, 50\%, 60\%, 75\%, 90\%, and 100\% are used.}
#'    \item{link}{an x of the link function test is the form of n by 1 matrix}
#'    \item{covform}{an x of the functional form test is the form of n by 1 matrix}
#' }
#'    See the documentation of \pkg{ggplot2} and \pkg{gridExtra} for details.
#' 
#' @importFrom ggplot2 ggplot geom_step geom_label theme theme_minimal ggtitle aes unit labs ylab xlab scale_y_continuous element_text
#' @importFrom gridExtra grid.arrange
#' @importFrom stats quantile
#' 
#' @example inst/examples/ex_plot.afttest.R
#' @export
plot.afttest <- function(x, npath = 50, std = TRUE, quantile = NULL, ...){
  
  # class
  if (!inherits(x,"afttest")) return(warning("Must be afttest class"))
  # npathsave
  if ((x$npathsave<1)) return(warning("afttest is conduced with npathsave=0"))
  # eqType
  eqType <- x$eqType
  # testType
  testType <- x$testType
  # std
  if (isTRUE(std)) {
    std <- "std"
  } else if (!std %in% c(T, F)) {
    return(warning("std needs to be logical."))
  }
  # npath
  if (!is.numeric(npath) || !(length(npath)==1)){
    return(warning("npath needs to be a positive integer."))
  } else {
    if (!is.numeric(npath)) {
      npath <- 50
    } else {
      npath <- min(npath,x$npathsave)
    }
  }
  
  stdTypeQuote <- ifelse(std=="std","standardized","unstandardized")
  testTypeQuote <- ifelse(eqType=="ns","non-smooth","induced-smoothed")
  
  x_axis <- 1:nrow(x$DF)
  
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
  
  if(testType=="omnibus"){
    resid <- c(NA)
    app <- matrix(NA)
    obs <- matrix(NA)
    
    Figure <- list(NA)
    for(k in 1:K){
      if (std == "std") {
        # DF_app
        DF_app=data.frame()
        for (group in 1:npath){
          temp <- x$apprx_std_npath[[group]][,Q[k]]
          temp <- data.frame(group,resid=x_axis,app=temp)
          DF_app <- rbind(DF_app,temp)
        }
        # DF_obs
        DF_obs <- data.frame(group,resid=x_axis,obs=x$obs_std_npath[,Q[k]])
        
      } else {
        #DF_app
        DF_app <- data.frame()
        for (group in 1:npath){
          temp <- x$apprx_npath[[group]][,Q[k]]
          temp <- data.frame(group,resid=x_axis,app=temp)
          DF_app <- rbind(DF_app,temp)
        }
        #DF_obs
        DF_obs <- data.frame(group,resid=x_axis,obs=x$obs_npath[,Q[k]])
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
          ggtitle(paste0("Omnibus test ","(",stdTypeQuote,")")) +
          # ggtitle(paste0("Omnibus test ","(",stdTypeQuote,")"," with ",names(Q[k]), " percentile for z")) +
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
    if (std == "std"){
      # DF_app
      DF_app <- data.frame()
      for (group in 1:npath){
        temp <- x$apprx_std_npath[[group]]
        temp <- data.frame(group,resid=x_axis,app=temp)
        DF_app <- rbind(DF_app,temp)
      }
      # DF_obs
      DF_obs <- data.frame(group,resid=x_axis,obs=x$obs_std_npath)
    } else {
      # DF_app
      DF_app <- data.frame()
      for (group in 1:npath){
        temp <- x$apprx_npath[[group]]
        temp <- data.frame(group,resid=x_axis,app=temp)
        DF_app <- rbind(DF_app,temp)
      }
      # DF_obs
      DF_obs <- data.frame(group,resid=x_axis,obs=x$obs_npath)
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
    
  } else if(testType=="covform"){
    resid <- c(NA)
    app <- c(NA)
    obs <- c(NA)
    if (std == "std"){
      # DF_app
      DF_app <- data.frame()
      for (group in 1:npath){
        temp <- x$apprx_std_npath[[group]]
        temp <- data.frame(group,resid=x_axis,app=temp)
        DF_app <- rbind(DF_app,temp)
      }
      # DF_obs
      DF_obs <- data.frame(group,resid=x_axis,obs=x$obs_std_npath)
      
    } else {
      # DF_app
      DF_app <- data.frame()
      for (group in 1:npath){
        temp <- x$apprx_npath[[group]]
        temp <- data.frame(group,resid=x_axis,app=temp)
        DF_app <- rbind(DF_app,temp)
      }
      # DF_obs
      DF_obs <- data.frame(group,resid=x_axis,obs=x$obs_npath)
      
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

##############################################################################
## Surv
##############################################################################
#' \code{Surv} function imported from \code{survival}
#'
#' This function is imported from the \code{survival} package. See
#' \code{\link[survival]{Surv}}.
#'
#' @importFrom survival Surv
#' @name export_Surv
#' @aliases Surv 
#' @export Surv
NULL
