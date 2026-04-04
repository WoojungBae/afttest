# Define functions -----------------------------------------------------------------
# Define function to generate datasets
generate_data = function(n, gamma0, Scenario) {
  if (Scenario==11){
    # --------------------------------- Scenario 1 ---------------------------------
    # True Coefficient: beta_0 = -4, beta_1 = 1
    beta_0 = -4
    beta_1 = 1
    if(gamma0 == 0){
      tau = 78.900
    } else if(gamma0 == 0.1){
      tau = 61.030
    } else if(gamma0 == 0.2){
      tau = 48.900
    } else if(gamma0 == 0.3){
      tau = 39.950
    } else if(gamma0 == 0.4){
      tau = 33.090
    } else if(gamma0 == 0.5){
      tau = 28.000
    }
    
    # ------------------------------------------------------------------------------
    # Generate covariate Z
    Z = rnorm(n,2,1)
    
    # ------------------------------------------------------------------------------
    # censoring rate 20% 
    # T: rnorm(n,0,1)
    # C: runif(n,0,tau)
    
    # T_data: true event time
    # C_data: true censoring time
    # X_data: observed time
    # D_data: observed indicator
    # Z_data: covariates fitted
    T_data = exp(-beta_0-beta_1*Z-gamma0*(Z^{2})+rnorm(n,0,1))
    C_data = runif(n,0,tau)
    X_data = C_data*(T_data>C_data)+T_data*(T_data<=C_data)
    D_data = 0*(T_data>C_data)+1*(T_data<=C_data)
    Z_data = Z
    
    return(data.frame(X=X_data, D=D_data, Z=Z))
  } else if (Scenario==12){
    # --------------------------------- Scenario 1 ---------------------------------
    # True Coefficient: beta_0 = -4, beta_1 = 1
    beta_0 = -4
    beta_1 = 1
    if(gamma0 == 0){
      tau = 28.200
    } else if(gamma0 == 0.1){
      tau = 19.530
    } else if(gamma0 == 0.2){
      tau = 13.900
    } else if(gamma0 == 0.3){
      tau = 10.040
    } else if(gamma0 == 0.4){
      tau = 7.310
    } else if(gamma0 == 0.5){
      tau = 5.350
    }
    
    # ------------------------------------------------------------------------------
    # Generate covariate Z
    Z = rnorm(n,2,1)
    
    # ------------------------------------------------------------------------------
    # censoring rate 40% 
    # T: rnorm(n,0,1)
    # C: runif(n,0,tau)
    
    # T_data: true event time
    # C_data: true censoring time
    # X_data: observed time
    # D_data: observed indicator
    # Z_data: covariates fitted
    T_data = exp(-beta_0-beta_1*Z-gamma0*(Z^{2})+rnorm(n,0,1))
    C_data = runif(n,0,tau)
    X_data = C_data*(T_data>C_data)+T_data*(T_data<=C_data)
    D_data = 0*(T_data>C_data)+1*(T_data<=C_data)
    Z_data = Z
    
    return(data.frame(X=X_data, D=D_data, Z=Z))
  } else if (Scenario==21){
    # --------------------------------- Scenario 2 ---------------------------------
    # True Coefficient: beta_0 = -4, beta_1 = 1, beta_2 = 1
    beta_0 = -4
    beta_1 = 1
    beta_2 = 1
    if(gamma0 == 0){
      tau = 51.350
    } else if(gamma0 == 0.1){
      tau = 39.280
    } else if(gamma0 == 0.2){
      tau = 31.095
    } else if(gamma0 == 0.3){
      tau = 25.270
    } else if(gamma0 == 0.4){
      tau = 20.840
    } else if(gamma0 == 0.5){
      tau = 17.415
    }
    
    # ------------------------------------------------------------------------------
    # Generate covariate Z
    Z1 = rbinom(n,1,0.5)
    Z2 = rnorm(n,2,1)
    
    # ------------------------------------------------------------------------------
    # censoring rate 20% 
    # T: rnorm(n,0,1)
    # C: runif(n,0,tau)
    
    # T_data: true event time
    # C_data: true censoring time
    # X_data: observed time
    # D_data: observed indicator
    # Z_data: covariates fitted
    T_data = exp(- beta_0 - beta_1 * Z1 - beta_2 * Z2 - gamma0 * (Z2^{2}) + rnorm(n,0,1))
    C_data = runif(n,0,tau)
    X_data = C_data * (T_data > C_data) + T_data * (T_data <= C_data)
    D_data = 0 * (T_data > C_data) + 1 * (T_data <= C_data)
    Z_data = cbind(Z1,Z2)
    
    return(data.frame(X=X_data, D=D_data, Z1=Z1, Z2=Z2))
  } else if (Scenario==22){
    # --------------------------------- Scenario 2 ---------------------------------
    # True Coefficient: beta_0 = -4, beta_1 = 1, beta_2 = 1
    beta_0 = -4
    beta_1 = 1
    beta_2 = 1
    if(gamma0 == 0){
      tau = 17.550
    } else if(gamma0 == 0.1){
      tau = 12.063
    } else if(gamma0 == 0.2){
      tau = 8.515
    } else if(gamma0 == 0.3){
      tau = 6.130
    } else if(gamma0 == 0.4){
      tau = 4.445
    } else if(gamma0 == 0.5){
      tau = 3.255
    }
    
    # ------------------------------------------------------------------------------
    # Generate covariate Z
    Z1 = rbinom(n,1,0.5)
    Z2 = rnorm(n,2,1)
    
    # ------------------------------------------------------------------------------
    # censoring rate 40% 
    # T: rnorm(n,0,1)
    # C: runif(n,0,tau)
    
    # T_data: true event time
    # C_data: true censoring time
    # X_data: observed time
    # D_data: observed indicator
    # Z_data: covariates fitted
    T_data = exp(- beta_0 - beta_1 * Z1 - beta_2 * Z2 - gamma0 * (Z2^{2}) + rnorm(n,0,1))
    C_data = runif(n,0,tau)
    X_data = C_data * (T_data > C_data) + T_data * (T_data <= C_data)
    D_data = 0 * (T_data > C_data) + 1 * (T_data <= C_data)
    Z_data = cbind(Z1,Z2)
    
    return(data.frame(X=X_data, D=D_data, Z1=Z1, Z2=Z2))
  } else if (Scenario==31){
    # --------------------------------- Scenario 3 ---------------------------------
    # True Coefficient: beta_0 = -4, beta_1 = 1, beta_2 = 1
    beta_0 = -4
    beta_1 = 1
    beta_2 = 1
    if(gamma0 == 0){
      tau = 51.335
    } else if(gamma0 == 0.1){
      tau = 50.610
    } else if(gamma0 == 0.2){
      tau = 49.960
    } else if(gamma0 == 0.3){
      tau = 49.500
    } else if(gamma0 == 0.4){
      tau = 49.100
    } else if(gamma0 == 0.5){
      tau = 48.600
    }
    
    # ------------------------------------------------------------------------------
    # Generate covariate Z
    Z1 = rbinom(n,1,0.5)
    Z2 = rnorm(n,2,1)
    
    # ------------------------------------------------------------------------------
    # censoring rate 20% 
    # T: rnorm(n,0,1)
    # C: runif(n,0,tau)
    
    # T_data: true event time
    # C_data: true censoring time
    # X_data: observed time
    # D_data: observed indicator
    # Z_data: covariates fitted
    T_data = exp(- beta_0 - beta_1 * Z1 - beta_2 * Z2 - gamma0 * lgamma(1+0.3*Z2^{2}) + rnorm(n,0,1))
    C_data = runif(n,0,tau)
    X_data = C_data * (T_data > C_data) + T_data * (T_data <= C_data)
    D_data = 0 * (T_data > C_data) + 1 * (T_data <= C_data)
    Z_data = cbind(Z1,Z2)
    
    return(data.frame(X=X_data, D=D_data, Z1=Z1, Z2=Z2))
  } else if (Scenario==32){
    # --------------------------------- Scenario 3 ---------------------------------
    # True Coefficient: beta_0 = -4, beta_1 = 1, beta_2 = 1
    beta_0 = -4
    beta_1 = 1
    beta_2 = 1
    if(gamma0 == 0){
      tau = 17.550
    } else if(gamma0 == 0.1){
      tau = 17.035
    } else if(gamma0 == 0.2){
      tau = 16.585
    } else if(gamma0 == 0.3){
      tau = 16.180
    } else if(gamma0 == 0.4){
      tau = 15.855
    } else if(gamma0 == 0.5){
      tau = 15.535
    }
    
    # ------------------------------------------------------------------------------
    # Generate covariate Z
    Z1 = rbinom(n,1,0.5)
    Z2 = rnorm(n,2,1)
    
    # ------------------------------------------------------------------------------
    # censoring rate 40% 
    # T: rnorm(n,0,1)
    # C: runif(n,0,tau)
    
    # T_data: true event time
    # C_data: true censoring time
    # X_data: observed time
    # D_data: observed indicator
    # Z_data: covariates fitted
    T_data = exp(- beta_0 - beta_1 * Z1 - beta_2 * Z2 - gamma0 * lgamma(1+0.3*Z2^{2}) + rnorm(n,0,1))
    C_data = runif(n,0,tau)
    X_data = C_data * (T_data > C_data) + T_data * (T_data <= C_data)
    D_data = 0 * (T_data > C_data) + 1 * (T_data <= C_data)
    Z_data = cbind(Z1,Z2)
    
    return(data.frame(X=X_data, D=D_data, Z1=Z1, Z2=Z2))
  }
  # return(list(X=X_data, D=D_data, Z=Z_data))
}

rejectionratio = function(Scenario) {
  # Type 1 error control
  alpha = 0.05
  
  # ------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------
  S1 = as.numeric(paste0(Scenario,1))
  {
    N = 100
    {
      # ------------------------------------------------------------------------------
      # ------------------------------------------------------------------------------
      gamma0 = 0
      txt.title = paste0("afttest","Scn",S1,"N",N,"gamma",gamma0*10,"_result.txt")
      
      gamma0_result = read.table(txt.title, header = TRUE, sep = "", dec = ".", fill = T)
      gamma0_sim = nrow(gamma0_result)
      gamma0_rejectionratio = matrix(apply(gamma0_result[,-1], 2, function(l) sum(l<=alpha, na.rm = T)),ncol=3)/gamma0_sim
      gamma0_rejectionratio = round(gamma0_rejectionratio,4)
      
      # ------------------------------------------------------------------------------
      gamma0 = 0.1
      txt.title = paste0("afttest","Scn",S1,"N",N,"gamma",gamma0*10,"_result.txt")
      
      gamma01_result = read.table(txt.title, header = TRUE, sep = "", dec = ".", fill = T)
      gamma01_sim = nrow(gamma01_result)
      gamma01_rejectionratio = matrix(apply(gamma01_result[,-1], 2, function(l) sum(l<=alpha, na.rm = T)),ncol=3)/gamma01_sim
      gamma01_rejectionratio = round(gamma01_rejectionratio,4)
      
      # ------------------------------------------------------------------------------
      gamma0 = 0.2
      txt.title = paste0("afttest","Scn",S1,"N",N,"gamma",gamma0*10,"_result.txt")
      
      gamma02_result = read.table(txt.title, header = TRUE, sep = "", dec = ".", fill = T)
      gamma02_sim = nrow(gamma02_result)
      gamma02_rejectionratio = matrix(apply(gamma02_result[,-1], 2, function(l) sum(l<=alpha, na.rm = T)),ncol=3)/gamma02_sim
      gamma02_rejectionratio = round(gamma02_rejectionratio,4)
      
      # ------------------------------------------------------------------------------
      gamma0 = 0.3
      txt.title = paste0("afttest","Scn",S1,"N",N,"gamma",gamma0*10,"_result.txt")
      
      gamma03_result = read.table(txt.title, header = TRUE, sep = "", dec = ".", fill = T)
      gamma03_sim = nrow(gamma03_result)
      gamma03_rejectionratio = matrix(apply(gamma03_result[,-1], 2, function(l) sum(l<=alpha, na.rm = T)),ncol=3)/gamma03_sim
      gamma03_rejectionratio = round(gamma03_rejectionratio,4)
      
      # ------------------------------------------------------------------------------
      gamma0 = 0.4
      txt.title = paste0("afttest","Scn",S1,"N",N,"gamma",gamma0*10,"_result.txt")
      
      gamma04_result = read.table(txt.title, header = TRUE, sep = "", dec = ".", fill = T)
      gamma04_sim = nrow(gamma04_result)
      gamma04_rejectionratio = matrix(apply(gamma04_result[,-1], 2, function(l) sum(l<=alpha, na.rm = T)),ncol=3)/gamma04_sim
      gamma04_rejectionratio = round(gamma04_rejectionratio,4)
      
      # ------------------------------------------------------------------------------
      gamma0 = 0.5
      txt.title = paste0("afttest","Scn",S1,"N",N,"gamma",gamma0*10,"_result.txt")
      
      gamma05_result = read.table(txt.title, header = TRUE, sep = "", dec = ".", fill = T)
      gamma05_sim = nrow(gamma05_result)
      gamma05_rejectionratio = matrix(apply(gamma05_result[,-1], 2, function(l) sum(l<=alpha, na.rm = T)),ncol=3)/gamma05_sim
      gamma05_rejectionratio = round(gamma05_rejectionratio,4)
      
      # ------------------------------------------------------------------------------
      # ------------------------------------------------------------------------------
      Scn1n100 = cbind(gamma0_rejectionratio,gamma01_rejectionratio,gamma02_rejectionratio,
                       gamma03_rejectionratio,gamma04_rejectionratio,gamma05_rejectionratio)
    }
    
    N = 300
    {
      # ------------------------------------------------------------------------------
      # ------------------------------------------------------------------------------
      gamma0 = 0
      txt.title = paste0("afttest","Scn",S1,"N",N,"gamma",gamma0*10,"_result.txt")
      
      gamma0_result = read.table(txt.title, header = TRUE, sep = "", dec = ".", fill = T)
      gamma0_sim = nrow(gamma0_result)
      gamma0_rejectionratio = matrix(apply(gamma0_result[,-1], 2, function(l) sum(l<=alpha, na.rm = T)),ncol=3)/gamma0_sim
      gamma0_rejectionratio = round(gamma0_rejectionratio,4)
      
      # ------------------------------------------------------------------------------
      gamma0 = 0.1
      txt.title = paste0("afttest","Scn",S1,"N",N,"gamma",gamma0*10,"_result.txt")
      
      gamma01_result = read.table(txt.title, header = TRUE, sep = "", dec = ".", fill = T)
      gamma01_sim = nrow(gamma01_result)
      gamma01_rejectionratio = matrix(apply(gamma01_result[,-1], 2, function(l) sum(l<=alpha, na.rm = T)),ncol=3)/gamma01_sim
      gamma01_rejectionratio = round(gamma01_rejectionratio,4)
      
      # ------------------------------------------------------------------------------
      gamma0 = 0.2
      txt.title = paste0("afttest","Scn",S1,"N",N,"gamma",gamma0*10,"_result.txt")
      
      gamma02_result = read.table(txt.title, header = TRUE, sep = "", dec = ".", fill = T)
      gamma02_sim = nrow(gamma02_result)
      gamma02_rejectionratio = matrix(apply(gamma02_result[,-1], 2, function(l) sum(l<=alpha, na.rm = T)),ncol=3)/gamma02_sim
      gamma02_rejectionratio = round(gamma02_rejectionratio,4)
      
      # ------------------------------------------------------------------------------
      gamma0 = 0.3
      txt.title = paste0("afttest","Scn",S1,"N",N,"gamma",gamma0*10,"_result.txt")
      
      gamma03_result = read.table(txt.title, header = TRUE, sep = "", dec = ".", fill = T)
      gamma03_sim = nrow(gamma03_result)
      gamma03_rejectionratio = matrix(apply(gamma03_result[,-1], 2, function(l) sum(l<=alpha, na.rm = T)),ncol=3)/gamma03_sim
      gamma03_rejectionratio = round(gamma03_rejectionratio,4)
      
      # ------------------------------------------------------------------------------
      gamma0 = 0.4
      txt.title = paste0("afttest","Scn",S1,"N",N,"gamma",gamma0*10,"_result.txt")
      
      gamma04_result = read.table(txt.title, header = TRUE, sep = "", dec = ".", fill = T)
      gamma04_sim = nrow(gamma04_result)
      gamma04_rejectionratio = matrix(apply(gamma04_result[,-1], 2, function(l) sum(l<=alpha, na.rm = T)),ncol=3)/gamma04_sim
      gamma04_rejectionratio = round(gamma04_rejectionratio,4)
      
      # ------------------------------------------------------------------------------
      gamma0 = 0.5
      txt.title = paste0("afttest","Scn",S1,"N",N,"gamma",gamma0*10,"_result.txt")
      
      gamma05_result = read.table(txt.title, header = TRUE, sep = "", dec = ".", fill = T)
      gamma05_sim = nrow(gamma05_result)
      gamma05_rejectionratio = matrix(apply(gamma05_result[,-1], 2, function(l) sum(l<=alpha, na.rm = T)),ncol=3)/gamma05_sim
      gamma05_rejectionratio = round(gamma05_rejectionratio,4)
      
      # ------------------------------------------------------------------------------
      # ------------------------------------------------------------------------------
      Scn1n300 = cbind(gamma0_rejectionratio,gamma01_rejectionratio,gamma02_rejectionratio,
                       gamma03_rejectionratio,gamma04_rejectionratio,gamma05_rejectionratio)
    }
    
    N = 500
    {
      # ------------------------------------------------------------------------------
      # ------------------------------------------------------------------------------
      gamma0 = 0
      txt.title = paste0("afttest","Scn",S1,"N",N,"gamma",gamma0*10,"_result.txt")
      
      gamma0_result = read.table(txt.title, header = TRUE, sep = "", dec = ".", fill = T)
      gamma0_sim = nrow(gamma0_result)
      gamma0_rejectionratio = matrix(apply(gamma0_result[,-1], 2, function(l) sum(l<=alpha, na.rm = T)),ncol=3)/gamma0_sim
      gamma0_rejectionratio = round(gamma0_rejectionratio,4)
      
      # ------------------------------------------------------------------------------
      gamma0 = 0.1
      txt.title = paste0("afttest","Scn",S1,"N",N,"gamma",gamma0*10,"_result.txt")
      
      gamma01_result = read.table(txt.title, header = TRUE, sep = "", dec = ".", fill = T)
      gamma01_sim = nrow(gamma01_result)
      gamma01_rejectionratio = matrix(apply(gamma01_result[,-1], 2, function(l) sum(l<=alpha, na.rm = T)),ncol=3)/gamma01_sim
      gamma01_rejectionratio = round(gamma01_rejectionratio,4)
      
      # ------------------------------------------------------------------------------
      gamma0 = 0.2
      txt.title = paste0("afttest","Scn",S1,"N",N,"gamma",gamma0*10,"_result.txt")
      
      gamma02_result = read.table(txt.title, header = TRUE, sep = "", dec = ".", fill = T)
      gamma02_sim = nrow(gamma02_result)
      gamma02_rejectionratio = matrix(apply(gamma02_result[,-1], 2, function(l) sum(l<=alpha, na.rm = T)),ncol=3)/gamma02_sim
      gamma02_rejectionratio = round(gamma02_rejectionratio,4)
      
      # ------------------------------------------------------------------------------
      gamma0 = 0.3
      txt.title = paste0("afttest","Scn",S1,"N",N,"gamma",gamma0*10,"_result.txt")
      
      gamma03_result = read.table(txt.title, header = TRUE, sep = "", dec = ".", fill = T)
      gamma03_sim = nrow(gamma03_result)
      gamma03_rejectionratio = matrix(apply(gamma03_result[,-1], 2, function(l) sum(l<=alpha, na.rm = T)),ncol=3)/gamma03_sim
      gamma03_rejectionratio = round(gamma03_rejectionratio,4)
      
      # ------------------------------------------------------------------------------
      gamma0 = 0.4
      txt.title = paste0("afttest","Scn",S1,"N",N,"gamma",gamma0*10,"_result.txt")
      
      gamma04_result = read.table(txt.title, header = TRUE, sep = "", dec = ".", fill = T)
      gamma04_sim = nrow(gamma04_result)
      gamma04_rejectionratio = matrix(apply(gamma04_result[,-1], 2, function(l) sum(l<=alpha, na.rm = T)),ncol=3)/gamma04_sim
      gamma04_rejectionratio = round(gamma04_rejectionratio,4)
      
      # ------------------------------------------------------------------------------
      gamma0 = 0.5
      txt.title = paste0("afttest","Scn",S1,"N",N,"gamma",gamma0*10,"_result.txt")
      
      gamma05_result = read.table(txt.title, header = TRUE, sep = "", dec = ".", fill = T)
      gamma05_sim = nrow(gamma05_result)
      gamma05_rejectionratio = matrix(apply(gamma05_result[,-1], 2, function(l) sum(l<=alpha, na.rm = T)),ncol=3)/gamma05_sim
      gamma05_rejectionratio = round(gamma05_rejectionratio,4)
      
      # ------------------------------------------------------------------------------
      # ------------------------------------------------------------------------------
      Scn1n500 = cbind(gamma0_rejectionratio,gamma01_rejectionratio,gamma02_rejectionratio,
                       gamma03_rejectionratio,gamma04_rejectionratio,gamma05_rejectionratio)
    }
  }
  
  # ------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------
  result = t(rbind(Scn1n100,Scn1n300,Scn1n500))
  
  row2 = 1:18*2; row1 = row2 - 1
  resultdataframe = matrix(nrow=36, ncol=9)
  for (cc in 1:9) {
    ccc = cc*2
    cccc = ccc - 1
    resultdataframe[row1,cc] = result[,ccc]  # odd: standardized test 
    resultdataframe[row2,cc] = result[,cccc] # even: unstandardized test
  }
  colnames(resultdataframe) = rep(c("ns", "is", "ls"), 3)
  rownames(resultdataframe) = rep(c("omni", "link", "form"), 6, each = 2)
  
  resultdataframe = data.frame(gamma0 = rep(paste0("gamma",0:5), each = 6),
                               rep(c("omni", "link", "form"),6, each = 2),
                               resultdataframe)
  colnames(resultdataframe) = NULL
  rownames(resultdataframe) = NULL
  
  return(resultdataframe)
}
# End function definitions ---------------------------------------------------------

# Load R packages
# library(Rcpp)
# library(RcppArmadillo)
# library(survival)
# library(ggplot2)
# library(gridExtra)
library(aftgee)
library(afttest)

# Load R code
# source("https://raw.githubusercontent.com/WoojungBae/afttest_analysis/main/source/afttest_source_r.R")
# setwd("/Users/woojung/Documents/Rproject/afttest_package_manuscript/analysis/simulation/Sim2")
# source("afttest_source.R")

Scenario = 21

# Type 2 error check
# gamma_0 = 0
# gamma_0 = 0.1
# gamma_0 = 0.2
# gamma_0 = 0.3
# gamma_0 = 0.4
# gamma_0 = 0.5

# The number of the approximated paths for each simulation is 200
npath = 200

gamma_0s = c(0, 0.5)
Ns = c(100, 500)
# Ns = c(300)

# Extract ID for simulated dataset (specific to LSF computing cluster)
# Note: The LSB_JOBINDEX is specified in the bsub command using the -J
# option
# run_ID = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# if (!dir.exists("Results")) {
#   dir.create("Results")
# }

# Extract ID for simulated dataset (specific to LSF computing cluster)
# Note: The LSB_JOBINDEX is specified in the bsub command using the -J
# option
run_ID = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

for (gamma_0 in gamma_0s) {
  for (N in Ns) {
    # ------ Define of constants ------
    # Define number of observations for each dataset
    txt.title = paste0("Results/afttest_TIME_","Scn",Scenario,"N",N,"gamma",gamma_0*10,"_result.txt")
    
    # 1. SETUP OUTPUT FILE HEADERS
    if (run_ID == 1) {
      df = data.frame(matrix(ncol = 19, nrow = 0)) # 1 run_ID + 18 time columns
      df_col_names = c("run_ID", 
                       # Omnibus
                       "omni_ns_linT", "omni_ns_linF",
                       "omni_is_linT", "omni_is_linF",
                       "omni_ls_linT", "omni_ls_linF",
                       # Link
                       "link_ns_linT", "link_ns_linF",
                       "link_is_linT", "link_is_linF",
                       "link_ls_linT", "link_ls_linF",
                       # Functional Form
                       "form_ns_linT", "form_ns_linF",
                       "form_is_linT", "form_is_linF",
                       "form_ls_linT", "form_ls_linF")
      colnames(df) = df_col_names
      write.table(df, file = txt.title, sep = "\t", row.names = FALSE, col.names = TRUE)
    }
    
    # A. Generate Data
    temp_data = generate_data(N, gamma_0, Scenario)
    
    # B. OPTIMIZATION: Fit Models ONCE (Not timed, or timed separately if desired)
    # 1. Rank-Based Non-Smoothed (RR-NS)
    fit_ns <- aftsrr(Surv(X, D) ~ Z1 + Z2, data = temp_data, eqType = "ns", rankWeights = "gehan")
    # 2. Rank-Based Induced-Smoothed (RR-IS)
    fit_is <- aftsrr(Surv(X, D) ~ Z1 + Z2, data = temp_data, eqType = "is", rankWeights = "gehan")
    # 3. Least Squares (LS)
    fit_ls <- aftgee(Surv(X, D) ~ Z1 + Z2, data = temp_data)
    
    # C. Run Tests and Record Time (using "elapsed" time)
    
    # --- 1. OMNIBUS TESTS ---
    # NS
    omni_ns_linT <- system.time(afttest(fit_ns, data = temp_data, npath = npath, testType = 'omnibus', npathsave = 0, linApprox = TRUE, seed = 1))[["elapsed"]]
    omni_ns_linF <- system.time(afttest(fit_ns, data = temp_data, npath = npath, testType = 'omnibus', npathsave = 0, linApprox = FALSE, seed = 1))[["elapsed"]]
    # IS
    omni_is_linT <- system.time(afttest(fit_is, data = temp_data, npath = npath, testType = 'omnibus', npathsave = 0, linApprox = TRUE, seed = 1))[["elapsed"]]
    omni_is_linF <- system.time(afttest(fit_is, data = temp_data, npath = npath, testType = 'omnibus', npathsave = 0, linApprox = FALSE, seed = 1))[["elapsed"]]
    # LS
    omni_ls_linT <- system.time(afttest(fit_ls, data = temp_data, npath = npath, testType = 'omnibus', npathsave = 0, linApprox = TRUE, seed = 1))[["elapsed"]]
    omni_ls_linF <- system.time(afttest(fit_ls, data = temp_data, npath = npath, testType = 'omnibus', npathsave = 0, linApprox = FALSE, seed = 1))[["elapsed"]]
    
    # --- 2. LINK FUNCTION TESTS ---
    # NS
    link_ns_linT <- system.time(afttest(fit_ns, data = temp_data, npath = npath, testType = 'link', npathsave = 0, linApprox = TRUE, seed = 1))[["elapsed"]]
    link_ns_linF <- system.time(afttest(fit_ns, data = temp_data, npath = npath, testType = 'link', npathsave = 0, linApprox = FALSE, seed = 1))[["elapsed"]]
    # IS
    link_is_linT <- system.time(afttest(fit_is, data = temp_data, npath = npath, testType = 'link', npathsave = 0, linApprox = TRUE, seed = 1))[["elapsed"]]
    link_is_linF <- system.time(afttest(fit_is, data = temp_data, npath = npath, testType = 'link', npathsave = 0, linApprox = FALSE, seed = 1))[["elapsed"]]
    # LS
    link_ls_linT <- system.time(afttest(fit_ls, data = temp_data, npath = npath, testType = 'link', npathsave = 0, linApprox = TRUE, seed = 1))[["elapsed"]]
    link_ls_linF <- system.time(afttest(fit_ls, data = temp_data, npath = npath, testType = 'link', npathsave = 0, linApprox = FALSE, seed = 1))[["elapsed"]]
    
    # --- 3. FUNCTIONAL FORM TESTS (Z2) ---
    # NS
    form_ns_linT <- system.time(afttest(fit_ns, data = temp_data, npath = npath, testType = 'covForm', covTested = "Z2", npathsave = 0, linApprox = TRUE, seed = 1))[["elapsed"]]
    form_ns_linF <- system.time(afttest(fit_ns, data = temp_data, npath = npath, testType = 'covForm', covTested = "Z2", npathsave = 0, linApprox = FALSE, seed = 1))[["elapsed"]]
    # IS
    form_is_linT <- system.time(afttest(fit_is, data = temp_data, npath = npath, testType = 'covForm', covTested = "Z2", npathsave = 0, linApprox = TRUE, seed = 1))[["elapsed"]]
    form_is_linF <- system.time(afttest(fit_is, data = temp_data, npath = npath, testType = 'covForm', covTested = "Z2", npathsave = 0, linApprox = FALSE, seed = 1))[["elapsed"]]
    # LS
    form_ls_linT <- system.time(afttest(fit_ls, data = temp_data, npath = npath, testType = 'covForm', covTested = "Z2", npathsave = 0, linApprox = TRUE, seed = 1))[["elapsed"]]
    form_ls_linF <- system.time(afttest(fit_ls, data = temp_data, npath = npath, testType = 'covForm', covTested = "Z2", npathsave = 0, linApprox = FALSE, seed = 1))[["elapsed"]]
    
    # D. Store Results
    allinfo = data.frame(
      "run_ID"        = run_ID,
      
      "omni_ns_linT"  = omni_ns_linT,
      "omni_ns_linF"  = omni_ns_linF,
      "omni_is_linT"  = omni_is_linT,
      "omni_is_linF"  = omni_is_linF,
      "omni_ls_linT"  = omni_ls_linT,
      "omni_ls_linF"  = omni_ls_linF,
      
      "link_ns_linT"  = link_ns_linT,
      "link_ns_linF"  = link_ns_linF,
      "link_is_linT"  = link_is_linT,
      "link_is_linF"  = link_is_linF,
      "link_ls_linT"  = link_ls_linT,
      "link_ls_linF"  = link_ls_linF,
      
      "form_ns_linT"  = form_ns_linT,
      "form_ns_linF"  = form_ns_linF,
      "form_is_linT"  = form_is_linT,
      "form_is_linF"  = form_is_linF,
      "form_ls_linT"  = form_ls_linT,
      "form_ls_linF"  = form_ls_linF
    )
    
    write.table(allinfo, file = txt.title, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
    
    # Clean up
    rm(list = c("fit_ns", "fit_is", "fit_ls"))
    gc()
  }
}

