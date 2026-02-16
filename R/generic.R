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
  cat("\n\tAFT Goodness-of-Fit Test\n\n")
  cat("Call: \n")
  print(x$call)
  
  cat("\n--- Null Hypothesis (H0) ---\n")
  if (x$testType %in% c("omnibus", "omni")) {
    cat(paste0("The assumed semiparametric AFT model fits the data adequately. \n"))
  } else if (x$testType == "link") {
    cat(paste0("The relationship between covariates and the log survival time is \n"))
    cat(paste0("correctly specified. \n"))
  } else if (x$testType %in% c("covForm", "form")) {
    cov_name <- ifelse(!is.null(x$covTested), x$covTested, "selected covariate")
    cat(paste0("The functional form for covariate '", cov_name, "' is correctly specified.\n"))
    # cat(paste0("The relationship between the log survival time and the specific \n"))
    # cat(paste0("covariate chosen by the covTested argument is correctly specified. \n"))
  }
  
  cat("\n--- p-values ---\n")
  p_vals <- c(x$p_value, x$p_std_value)
  p_vals_fmt <- format.pval(p_vals, digits = 3, eps = 0.001)
  p_valueTAB <- data.frame(t(p_vals_fmt))
  rownames(p_valueTAB) <- ""
  colnames(p_valueTAB) <- c("unstandardized", "standardized")
  print(p_valueTAB)
  
  invisible(x)
}

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
  cat("\n\tAFT Goodness-of-Fit Test Summary\n\n")
  cat("Call:\n")
  print(object$call)
  
  cat("\n--- Null Hypothesis (H0) ---\n")
  if (object$testType %in% c("omnibus", "omni")) {
    cat(paste0("The assumed semiparametric AFT model fits the data adequately. \n"))
  } else if (object$testType == "link") {
    cat(paste0("The relationship between covariates and the log survival time is \n"))
    cat(paste0("correctly specified. \n"))
  } else if (object$testType %in% c("covForm", "form")) {
    cov_name <- ifelse(!is.null(object$covTested), object$covTested, "selected covariate")
    cat(paste0("The functional form for covariate '", cov_name, "' is correctly specified.\n"))
    # cat(paste0("The relationship between the log survival time and the specific \n"))
    # cat(paste0("covariate chosen by the covTested argument is correctly specified. \n"))
  }
  
  cat("\n--- Simulation Parameters ---\n")
  cat(paste("Test Type:       ", object$testType, "\n"))
  cat(paste("Bootstrap:       ", ifelse(object$linApprox, "Linear Approximation", "Standard Resampling"), "\n"))
  cat(paste("Resampling Paths:", object$npath, "\n"))
  cat(paste("Random Seed:     ", object$seed, "\n"))
  
  cat("\n--- p-values ---\n")
  p_vals <- c(object$p_value, object$p_std_value)
  p_vals_fmt <- format.pval(p_vals, digits = 3, eps = 0.001)
  
  p.valueTAB <- data.frame(t(p_vals_fmt))
  rownames(p.valueTAB) <- ""
  colnames(p.valueTAB) <- c("unstandardized", "standardized")
  print(p.valueTAB)
  
  method_name <- ifelse(object$estMethod == "ls", "Least Squares", "Rank-Based")
  cat(paste0("\n--- Coefficients (Method: ", method_name, ") ---"))
  if (object$estMethod=="ls"){
    cat(paste0("\n Coefficients (estimated by aftgee::aftgee):"))
  } else if (object$estMethod=="rr"){
    cat(paste0("\n Coefficients (estimated by aftgee::aftsrr):"))
  }
  coefTAB <- data.frame(t(object$beta))
  rownames(coefTAB) <- ""
  colnames(coefTAB) <- object$names[-c(1:2)]
  print(coefTAB)
  
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
#'    \item{covForm}{an x of the functional form test is the form of n by 1 matrix}
#' }
#'    See the documentation of \pkg{ggplot2} and \pkg{gridExtra} for details.
#' 
#' @importFrom ggplot2 ggplot geom_step geom_label theme theme_minimal ggtitle aes unit labs ylab xlab scale_y_continuous element_text
#' @importFrom gridExtra grid.arrange
#' @importFrom stats quantile
#' 
#' @example inst/examples/ex_afttest.R
#' @export
plot.afttest <- function(x, npath = 50, std = TRUE, quantile = NULL, ...){
  # eqType
  eqType <- x$eqType
  
  # testType
  testType <- x$testType
  
  # std
  if (isFALSE(std)) {
    std <- "unstd"
  } else if (isTRUE(std)) {
    std <- "std"
  } else {
    warning("Argument 'std' must be a logical value (TRUE or FALSE). Defaulting to TRUE.")
    std <- "std"
  }
  
  # npathsave
  if (is.null(x$npathsave) || x$npathsave == 0) {
    stop("Cannot create plot: 'afttest' was conducted with npathsave = 0. Please rerun the model with npathsave > 0 to generate paths.")
  }
  
  # npath
  if (!is.numeric(npath) || length(npath) != 1 || npath <= 0) {
    stop("Argument 'npath' must be a single positive integer.")
  }
  npath <- as.integer(npath)
  if (npath > x$npathsave) {
    warning(sprintf("Requested npath (%d) exceeds stored paths (%d). Plotting all %d available paths.", 
                    npath, x$npathsave, x$npathsave))
    npath <- x$npathsave
  }
  
  # eqType
  if (eqType=="ns") {
    testTypeQuote <- "non-smooth"
  } else if (eqType=="ns") {
    testTypeQuote <- "induced-smoothed"
  } else if (eqType=="ns") {
    testTypeQuote <- "least-squares"
  }
  
  stdTypeQuote <- ifelse(std=="std","standardized","unstandardized")
  
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
          ggplot2::ggplot() +
          ggplot2::geom_step(data=DF_app,ggplot2::aes(x=resid,y=app,group=group),colour="grey",alpha=0.5) +
          ggplot2::geom_step(data=DF_obs,ggplot2::aes(x=resid,y=obs),colour="tomato",lwd=0.25) +
          ggplot2::labs(y = NULL) + ggplot2::labs(x = NULL) +
          ggplot2::ggtitle(paste0("Omnibus test ","(",stdTypeQuote,")")) +
          # ggtitle(paste0("Omnibus test ","(",stdTypeQuote,")"," with ",names(Q[k]), " percentile for z")) +
          ggplot2::scale_y_continuous(breaks = y_breaks, labels = y_labels) +
          ggplot2::theme(plot.title=element_text(hjust=0.5),
                plot.margin = rep(unit(0,"null"),4),
                panel.spacing = unit(0,"null"))
      } else {
        y_breaks <- round(seq(y_breaksMIN, y_breaksMAX, length.out = 3),1)
        y_labels <- format(y_breaks, nsmall = 1)
        Figure_k <-
          ggplot2::ggplot() +
          ggplot2::geom_step(data=DF_app,ggplot2::aes(x=resid,y=app,group=group),colour="grey",alpha=0.5) +
          ggplot2::geom_step(data=DF_obs,ggplot2::aes(x=resid,y=obs),colour="tomato",lwd=0.25) +
          ggplot2::labs(y = NULL) + ggplot2::labs(x = NULL) + 
          ggplot2::scale_y_continuous(breaks = y_breaks, labels = y_labels) +
          ggplot2::geom_label(aes(x=-Inf,y=Inf),label=paste0(names(Q[k])),fill="darkgrey",label.size=NA,size=3,hjust=-0.0,vjust=1.0) +
          ggplot2::theme(plot.margin = rep(unit(0,"null"),4),
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
    
  } else if(testType=="covForm"){
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
