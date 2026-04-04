# ----------------------------------------------------------------------------
# Source codes
# ----------------------------------------------------------------------------
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

generate_result_table <- function(scenario = 21, alpha = 0.05, results_dir = "Results") {
  
  # Parameters used in your simulation
  gamma_0s <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
  Ns <- c(100, 300, 500)
  tests <- c("omni", "link", "form")
  ests <- c("ns", "is", "ls")
  
  # Helper to construct column names based on your df_col_names
  get_colname <- function(test, est, std) {
    paste0(test, "_", est, ifelse(std, "_stdpvalue", "_pvalue"))
  }
  
  # Helper to format numbers to 3 decimal places with optional bolding
  fmt <- function(val, is_bold) {
    if (is.na(val)) return(ifelse(is_bold, "\\textbf{NA}", "NA"))
    str_val <- sprintf("%.3f", val)
    if (is_bold) paste0("\\textbf{", str_val, "}") else str_val
  }
  
  # 1. Read all files and store rejection rates
  res_list <- list()
  for (g in gamma_0s) {
    for (n in Ns) {
      # Reconstruct the exact file name from your simulation
      filename <- paste0("afttestScn", scenario, "N", n, "gamma", g * 10, "_result.txt")
      filepath <- file.path(results_dir, filename)
      
      if (file.exists(filepath)) {
        df <- read.table(filepath, header = TRUE, sep = "\t")
        # Calculate rejection rates (proportion < alpha)
        rej_rates <- colMeans(df[, -1] < alpha, na.rm = TRUE)
        res_list[[paste(g, n, sep = "_")]] <- rej_rates
      } else {
        warning(paste("File not found:", filepath))
        res_list[[paste(g, n, sep = "_")]] <- NA
      }
    }
  }
  # g;n
  
  # 2. Print the LaTeX Table Header
  cat("\\begin{table}[htp]\n")
  cat("    \\fontsize{10}{11}\\selectfont\n")
  cat("    \\centering\n")
  cat("    \\caption{\\label{tab:sim:result} Rejection rates with the nominal level of \n")
  cat("    0.05 under the Model \\eqref{model:sim1}. The omnibus (omni), link function \n")
  cat("    (link), and functional form (form) tests are considered using the non-smooth \n")
  cat("    estimator (ns), induced smoothed estimator (is), and the least-squares \n")
  cat("    estimator (ls) under the sample sizes of $n = 100, 300$, and $500$ with a \n")
  cat("    censoring rate of 20\\%. Each number in a cell is based on 1,000 replications \n")
  cat("    for the standardized statistic (bold) or unstandardized statistic.} \n")
  cat("    \\begin{tabular}{cccccccccccc}\n")
  cat("        \\toprule\n")
  cat("        \\multicolumn{3}{c}{ } & \\multicolumn{3}{c}{N = 100} & \\multicolumn{3}{c}{N = 300} & \\multicolumn{3}{c}{N = 500} \\\\\n")
  cat("        \\cmidrule(l{3pt}r{3pt}){4-6} \\cmidrule(l{3pt}r{3pt}){7-9} \\cmidrule(l{3pt}r{3pt}){10-12}\n")
  cat("        $\\gamma$ & Test & std & ns & is & ls & ns & is & ls & ns & is & ls\\\\\n")
  cat("        \\midrule\n")
  
  # 3. Build the table rows dynamically
  for (i in seq_along(gamma_0s)) {
    g <- gamma_0s[i]
    
    for (t_idx in seq_along(tests)) {
      test <- tests[t_idx]
      
      # Extract values for T (Standardized) and F (Unstandardized)
      vals_T <- c(); vals_F <- c()
      
      for (n in Ns) {
        rates <- res_list[[paste(g, n, sep = "_")]]
        
        for (est in ests) {
          col_T <- get_colname(test, est, std = TRUE)
          col_F <- get_colname(test, est, std = FALSE)
          
          v_T <- if (is.na(rates[1])) NA else rates[[col_T]]
          v_F <- if (is.na(rates[1])) NA else rates[[col_F]]
          
          vals_T <- c(vals_T, fmt(v_T, is_bold = TRUE))
          vals_F <- c(vals_F, fmt(v_F, is_bold = FALSE))
        }
      }
      
      # Print T row
      cat("        &  & T & ", paste(vals_T, collapse = " & "), "\\\\\n")
      cat("        \n")
      
      # Print F row (with multirows for test labels and gamma labels)
      if (test == "omni") {
        cat("        & \\multirow[t]{-2}{*}{omni} & F & ", paste(vals_F, collapse = " & "), "\\\\\n")
        cat("        \n")
      } else if (test == "link") {
        cat("        & \\multirow[t]{-2}{*}{link} & F & ", paste(vals_F, collapse = " & "), "\\\\\n")
        cat("        \n")
      } else if (test == "form") {
        # The form row closes the block for a specific gamma, so it gets the gamma multirow label
        cat("        \\multirow[t]{-6}{*}{", g, "} & \\multirow[t]{-2}{*}{form} & F & ", paste(vals_F, collapse = " & "), "\\\\\n")
        
        # Add space between gamma blocks, unless it's the very last one
        if (i < length(gamma_0s)) cat("        \n")
      }
    }
  }
  
  # 4. Print the LaTeX Table Footer
  cat("        \\bottomrule\n")
  cat("    \\end{tabular}\n")
  cat("\\end{table}\n")
}

generate_time_table <- function(scenario = 21, results_dir = "Results") {
  
  # Parameters to include in the table
  gamma_0s <- c(0, 0.5)
  Ns <- c(100, 500) # Only grabbing 100 and 500 for this specific table
  tests <- c("omni", "link", "form")
  ests <- c("ns", "is", "ls")
  lins <- c("T", "F")
  
  # 1. Read files and calculate average times
  res_list <- list()
  for (g in gamma_0s) {
    for (n in Ns) {
      # Reconstruct the exact file name
      filename <- paste0("afttest_TIME_Scn", scenario, "N", n, "gamma", g * 10, "_result.txt")
      filepath <- file.path(results_dir, filename)
      
      if (file.exists(filepath)) {
        df <- read.table(filepath, header = TRUE, sep = "\t")
        # Calculate column means (ignoring run_ID)
        avg_times <- colMeans(df[, -1], na.rm = TRUE)
        res_list[[paste(g, n, sep = "_")]] <- avg_times
      } else {
        warning(paste("File not found:", filepath))
        res_list[[paste(g, n, sep = "_")]] <- NULL
      }
    }
  }
  
  # Helper to format numbers to 1 decimal place
  fmt <- function(val) {
    if (is.null(val) || is.na(val)) return("NA")
    sprintf("%.1f", val)
  }
  
  # 2. Print LaTeX Header
  cat("{\n")
  cat("\\begin{table}[htp]\n")
  cat("    % \\fontsize{10}{9}\\selectfont\n")
  cat("    \\caption{\\label{tab:TIMEsim:result} Average running times (in seconds) for the \n")
  cat("    omnibus, link function, and functional form tests based on 1,000 simulations \n")
  cat("    (each with 500 generated resampling paths) under the Model \\eqref{model:sim1}. \n") # FIXED HERE
  cat("    The table compares the computational efficiency of the proposed asymptotic \n")
  cat("    linear approximation (\\code{linApprox = TRUE}) against the standard resampling \n")
  cat("    method (\\code{linApprox = FALSE}) for non-smooth (ns), induced-smoothed (is), and \n")
  cat("    least squares (ls) estimators across sample sizes $n \\in \\{100, 500\\}$ with \n")
  cat("    a 20\\% censoring rate. Computations were performed using parallel computing \n")
  cat("    on the University of Florida HiPerGator AI cluster (70,320 cores, 18GB RAM \n")
  cat("    per core).}\n")
  cat("    \\centering\n")
  cat("    \\begin{tabular}[t]{ccccccccc}\n")
  cat("        \\toprule\n")
  cat("        \\multicolumn{3}{c}{ } & \\multicolumn{3}{c}{$n = 100$} & \\multicolumn{3}{c}{$n = 500$} \\\\\n")
  cat("        \\cmidrule(l{3pt}r{3pt}){4-6} \\cmidrule(l{3pt}r{3pt}){7-9}\n")
  cat("        $\\gamma$ & Test & linApprox & ns & is & ls & ns & is & ls\\\\\n")
  cat("        \\midrule\n")
  
  # 3. Build Table Rows
  for (i in seq_along(gamma_0s)) {
    g <- gamma_0s[i]
    
    for (j in seq_along(tests)) {
      test <- tests[j]
      
      for (k in seq_along(lins)) {
        lin <- lins[k]
        
        # Determine if we should print the gamma label (only first row of the block)
        str_g <- if (j == 1 && k == 1) as.character(g) else " "
        
        # Determine if we should print the test label (only first row of the test pair)
        str_test <- if (k == 1) test else " "
        
        # Extract values for N=100 and N=500
        vals <- c()
        for (n in Ns) {
          rates <- res_list[[paste(g, n, sep = "_")]]
          for (est in ests) {
            col_name <- paste0(test, "_", est, "_lin", lin)
            val <- if(is.null(rates)) NA else rates[[col_name]]
            vals <- c(vals, fmt(val))
          }
        }
        
        # Print the data row
        cat("       ", str_g, "&", str_test, "&", lin, "&", paste(vals, collapse = " & "), "\\\\\n")
      }
    }
    
    # Add midrule between gamma blocks (if not the last block)
    if (i < length(gamma_0s)) {
      cat("        \\cmidrule{1-9}\n")
    }
  }
  
  # 4. Print LaTeX Footer
  cat("        \\bottomrule\n")
  cat("    \\end{tabular}\n")
  cat("\\end{table}\n")
  cat("}\n")
}

# ----------------------------------------------------------------------------
# Install packages
# ----------------------------------------------------------------------------
required_packages <- c("survival", "aftgee", "afttest", 
                       "foreach", "doParallel", "doRNG", 
                       "Rcpp", "RcppArmadillo")
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(missing_packages)) {
  install.packages(missing_packages)
}
lapply(required_packages, library, character.only = TRUE)

library(survival)
library(aftgee)
library(afttest)
library(foreach)
library(doParallel)
library(doRNG)

set.seed(0)

# directory
currentwd = getwd()
setwd(paste0(currentwd, "/Replication"))
if (!dir.exists("Results")) {
  dir.create("Results")
}

# ----------------------------------------------------------------------------
# Application: Load dataset from survival packages
# ----------------------------------------------------------------------------
pbc <- within(survival::pbc, {
  status <- ifelse(status == 2, 1, 0)
  log_bili <- log(bili)})

# ----------------------------------------------------------------------------
# Model 1: Surv(time, status) ~ bili + protime + albumin + age + edema
# ----------------------------------------------------------------------------
# estimation using aftgee::aftsrr under Model 1
pbc1_aftsrr <- aftgee::aftsrr(
  Surv(time, status) ~ bili + protime + albumin + age + edema, 
  data = pbc, eqType = "ns", rankWeights = "gehan")

# omnibus test under Model 1
system.time(
  pbc1_omni_ns <- 
    afttest(pbc1_aftsrr, data = pbc, npath = 200, testType = "omnibus", 
            estMethod = "rr", eqType = "ns", linApprox = TRUE, seed = 0))
pbc1_omni_ns
plot(pbc1_omni_ns, std = TRUE)

# link function test under Model 1
system.time(
  pbc1_link_ns <- 
    afttest(pbc1_aftsrr, data = pbc, npath = 200, testType = "link", 
            estMethod = "rr", eqType = "ns", linApprox = TRUE, seed = 0))
pbc1_link_ns
plot(pbc1_link_ns, std = TRUE)

# functional form test for "bili" under Model 1
system.time(
  pbc1_form1_ns <- afttest(
    Surv(time, status) ~ bili + protime + albumin + age + edema, 
    data = pbc, npath = 200, testType = "covForm", 
    estMethod = "rr", eqType = "ns", covTested = "bili", 
    linApprox = TRUE, seed = 0))
pbc1_form1_ns
plot(pbc1_form1_ns, std = TRUE)

# functional form test for protime" under Model 1
system.time(
  pbc1_form2_ns <- afttest(
    Surv(time, status) ~ bili + protime + albumin + age + edema, 
    data = pbc, npath = 200, testType = "covForm", 
    estMethod = "rr", eqType = "ns", covTested = "protime",
    linApprox = TRUE, seed = 0))
pbc1_form2_ns
plot(pbc1_form2_ns, std = TRUE)

# functional form test for "albumin" under Model 1
system.time(
  pbc1_form3_ns <- afttest(
    Surv(time, status) ~ bili + protime + albumin + age + edema, 
    data = pbc, npath = 200, testType = "covForm", 
    estMethod = "rr", eqType = "ns", covTested = "albumin",
    linApprox = TRUE, seed = 0))
pbc1_form3_ns
plot(pbc1_form3_ns, std = TRUE)

# functional form test for "age" under Model 1
system.time(
  pbc1_form4_ns <- afttest(
    Surv(time, status) ~ bili + protime + albumin + age + edema, 
    data = pbc, npath = 200, testType = "covForm", 
    estMethod = "rr", eqType = "ns", covTested = "age",
    linApprox = TRUE, seed = 0))
pbc1_form4_ns
plot(pbc1_form4_ns, std = TRUE)

# ----------------------------------------------------------------------------
# Model 2: Surv(time, status) ~ log_bili + protime + albumin + age + edema
# ----------------------------------------------------------------------------
# estimation using aftgee::aftsrr under Model 2
pbc2_aftsrr <- aftgee::aftsrr(
  Surv(time, status) ~ log_bili + protime + albumin + age + edema, 
  data = pbc, eqType = "ns", rankWeights = "gehan")

# omnibus test under Model 2
system.time(
  pbc2_omni_ns <- 
    afttest(pbc2_aftsrr, data = pbc, npath = 200, testType = "omnibus", 
            estMethod = "rr", linApprox = TRUE, seed = 0))
pbc2_omni_ns
plot(pbc2_omni_ns)

# link function test under Model 2
system.time(
  pbc2_link_ns <- 
    afttest(pbc2_aftsrr, data = pbc, npath = 200, testType = "link",
            estMethod = "rr", linApprox = TRUE, seed = 0))
pbc2_link_ns
plot(pbc2_link_ns)

# functional form test for "log_bili" under Model 2
system.time(
  pbc2_form1_ns <- afttest(
    Surv(time, status) ~ log_bili + protime + albumin + age + edema, 
    data = pbc, npath = 200, testType = "covForm", 
    estMethod = "rr", eqType = "ns", covTested = "log_bili",
    linApprox = TRUE, seed = 0))

pbc2_form1_ns 
plot(pbc2_form1_ns)

# functional form test for "protime" under Model 2
system.time(
  pbc2_form2_ns <- afttest(
    Surv(time, status) ~ log_bili + protime + albumin + age + edema, 
    data = pbc, npath = 200, testType = "covForm", 
    estMethod = "rr", eqType = "ns", covTested = "protime",
    linApprox = TRUE, seed = 0))
pbc2_form2_ns 
plot(pbc2_form2_ns)

# functional form test for "albumin" under Model 2
system.time(
  pbc2_form3_ns <- afttest(
    Surv(time, status) ~ log_bili + protime + albumin + age + edema, 
    data = pbc, npath = 200, testType = "covForm", 
    estMethod = "rr", eqType = "ns", covTested = "albumin",
    linApprox = TRUE, seed = 0))
pbc2_form3_ns
plot(pbc2_form3_ns)

# functional form test for "age" under Model 2
system.time(
  pbc2_form4_ns <- afttest(
    Surv(time, status) ~ log_bili + protime + albumin + age + edema, 
    data = pbc, npath = 200, testType = "covForm", 
    estMethod = "rr", eqType = "ns", covTested = "age",
    linApprox = TRUE, seed = 0))
pbc2_form4_ns
plot(pbc2_form4_ns)

# ----------------------------------------------------------------------------
# Warning: Recommend to use hpc to generate Table 1 and Table 2
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# Table 1 - Time Results
# ----------------------------------------------------------------------------
# For each sim_per_file, it takes {33 min 02 sec} on average.
# Thus, 1000 iteration would take about {550 hours}.
# The Results in the manuscript was generated by hpc (UF HiPerGator).
# See package_sim2_hpc_TIME.R and package_sim2_hpc_TIME.sbatch

# Scenario
Scenario <- 21
# Number of simulation runs per setting
sim_per_file <- 1000
# The number of the approximated paths for each simulation is 200
npath = 200

# Type 1/Type 2 errors check
gamma_0s <- c(0, 0.5)

Ns <- c(100, 300, 500)
for (N in Ns) {
  for (gamma_0 in gamma_0s) {
    
    txt.title <- paste0(
      "Results/afttest_TIME_",
      "Scn", Scenario,
      "N", N,
      "gamma", gamma_0 * 10,
      "_result.txt"
    )
    
    # Write header
    df <- data.frame(matrix(ncol = 19, nrow = 0))
    df_col_names <- c(
      "run_ID",
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
      "form_ls_linT", "form_ls_linF"
    )
    colnames(df) <- df_col_names
    write.table(df, file = txt.title, sep = "\t",
                row.names = FALSE, col.names = TRUE)
    
    cat(paste0("Running: N = ", N, ", gamma_0 = ", gamma_0, "\n"))
    
    results_df <-
      foreach(
        run_ID = 1:sim_per_file,
        .combine = rbind,
        .packages = c("survival", "aftgee", "afttest")
      ) %dorng% {
        temp_data <- generate_data(N, gamma_0, Scenario)
        fit_ns <- aftsrr(
          Surv(X, D) ~ Z1 + Z2,
          data = temp_data,
          eqType = "ns",
          rankWeights = "gehan"
        )
        fit_is <- aftsrr(
          Surv(X, D) ~ Z1 + Z2,
          data = temp_data,
          eqType = "is",
          rankWeights = "gehan"
        )
        fit_ls <- aftgee(
          Surv(X, D) ~ Z1 + Z2,
          data = temp_data
        )
        
        # 1. Omnibus
        omni_ns_linT <- system.time(
          afttest(
            fit_ns, data = temp_data, npath = npath,
            testType = "omnibus", npathsave = 0,
            linApprox = TRUE, seed = 1
          )
        )[["elapsed"]]
        omni_ns_linF <- system.time(
          afttest(
            fit_ns, data = temp_data, npath = npath,
            testType = "omnibus", npathsave = 0,
            linApprox = FALSE, seed = 1
          )
        )[["elapsed"]]
        omni_is_linT <- system.time(
          afttest(
            fit_is, data = temp_data, npath = npath,
            testType = "omnibus", npathsave = 0,
            linApprox = TRUE, seed = 1
          )
        )[["elapsed"]]
        omni_is_linF <- system.time(
          afttest(
            fit_is, data = temp_data, npath = npath,
            testType = "omnibus", npathsave = 0,
            linApprox = FALSE, seed = 1
          )
        )[["elapsed"]]
        omni_ls_linT <- system.time(
          afttest(
            fit_ls, data = temp_data, npath = npath,
            testType = "omnibus", npathsave = 0,
            linApprox = TRUE, seed = 1
          )
        )[["elapsed"]]
        omni_ls_linF <- system.time(
          afttest(
            fit_ls, data = temp_data, npath = npath,
            testType = "omnibus", npathsave = 0,
            linApprox = FALSE, seed = 1
          )
        )[["elapsed"]]
        
        # 2. Link
        link_ns_linT <- system.time(
          afttest(
            fit_ns, data = temp_data, npath = npath,
            testType = "link", npathsave = 0,
            linApprox = TRUE, seed = 1
          )
        )[["elapsed"]]
        link_ns_linF <- system.time(
          afttest(
            fit_ns, data = temp_data, npath = npath,
            testType = "link", npathsave = 0,
            linApprox = FALSE, seed = 1
          )
        )[["elapsed"]]
        link_is_linT <- system.time(
          afttest(
            fit_is, data = temp_data, npath = npath,
            testType = "link", npathsave = 0,
            linApprox = TRUE, seed = 1
          )
        )[["elapsed"]]
        link_is_linF <- system.time(
          afttest(
            fit_is, data = temp_data, npath = npath,
            testType = "link", npathsave = 0,
            linApprox = FALSE, seed = 1
          )
        )[["elapsed"]]
        link_ls_linT <- system.time(
          afttest(
            fit_ls, data = temp_data, npath = npath,
            testType = "link", npathsave = 0,
            linApprox = TRUE, seed = 1
          )
        )[["elapsed"]]
        link_ls_linF <- system.time(
          afttest(
            fit_ls, data = temp_data, npath = npath,
            testType = "link", npathsave = 0,
            linApprox = FALSE, seed = 1
          )
        )[["elapsed"]]
        
        # 3. Functional form
        form_ns_linT <- system.time(
          afttest(
            fit_ns, data = temp_data, npath = npath,
            testType = "covForm", covTested = "Z2",
            npathsave = 0, linApprox = TRUE, seed = 1
          )
        )[["elapsed"]]
        form_ns_linF <- system.time(
          afttest(
            fit_ns, data = temp_data, npath = npath,
            testType = "covForm", covTested = "Z2",
            npathsave = 0, linApprox = FALSE, seed = 1
          )
        )[["elapsed"]]
        form_is_linT <- system.time(
          afttest(
            fit_is, data = temp_data, npath = npath,
            testType = "covForm", covTested = "Z2",
            npathsave = 0, linApprox = TRUE, seed = 1
          )
        )[["elapsed"]]
        form_is_linF <- system.time(
          afttest(
            fit_is, data = temp_data, npath = npath,
            testType = "covForm", covTested = "Z2",
            npathsave = 0, linApprox = FALSE, seed = 1
          )
        )[["elapsed"]]
        form_ls_linT <- system.time(
          afttest(
            fit_ls, data = temp_data, npath = npath,
            testType = "covForm", covTested = "Z2",
            npathsave = 0, linApprox = TRUE, seed = 1
          )
        )[["elapsed"]]
        form_ls_linF <- system.time(
          afttest(
            fit_ls, data = temp_data, npath = npath,
            testType = "covForm", covTested = "Z2",
            npathsave = 0, linApprox = FALSE, seed = 1
          )
        )[["elapsed"]]
        
        data.frame(
          run_ID = run_ID,
          
          omni_ns_linT = omni_ns_linT,
          omni_ns_linF = omni_ns_linF,
          omni_is_linT = omni_is_linT,
          omni_is_linF = omni_is_linF,
          omni_ls_linT = omni_ls_linT,
          omni_ls_linF = omni_ls_linF,
          
          link_ns_linT = link_ns_linT,
          link_ns_linF = link_ns_linF,
          link_is_linT = link_is_linT,
          link_is_linF = link_is_linF,
          link_ls_linT = link_ls_linT,
          link_ls_linF = link_ls_linF,
          
          form_ns_linT = form_ns_linT,
          form_ns_linF = form_ns_linF,
          form_is_linT = form_is_linT,
          form_is_linF = form_is_linF,
          form_ls_linT = form_ls_linT,
          form_ls_linF = form_ls_linF
        )
      }
    
    write.table(
      results_df,
      file = txt.title,
      sep = "\t",
      row.names = FALSE,
      col.names = FALSE,
      append = TRUE
    )
    
    gc()
  }
}

# Read the generated .txt files to compute average times
generate_time_table()

# ----------------------------------------------------------------------------
# Table 2 - Simulation Results
# ----------------------------------------------------------------------------
# For each sim_per_file, it takes {4min 49sec} on average.
# Thus, 1000 iteration would take about {80 hours}.
# The Results in the manuscript was generated by hpc (UF HiPerGator).
# See package_sim2_hpc.R and package_sim2_hpc.sbatch

# Type 1/Type 2 errors check
gamma_0s = c(0:5)/10

# Extract ID for simulated dataset (specific to LSF computing cluster)
# Note: The LSB_JOBINDEX is specified in the bsub command using the -J
# option
# run_ID = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

Ns = c(100, 300, 500)
for (N in Ns) {
  for (gamma_0 in gamma_0s) {
    # ------ Define of constants (adjust to fit the data generating scenario) ------
    # N = 100
    
    # Define number of observations for each dataset
    txt.title = paste0("Results/afttest","Scn",Scenario,"N",N,"gamma",gamma_0*10,"_result.txt")
    df = data.frame(matrix(ncol = 19, nrow = 0))
    df_col_names = c("run_ID", 
                     "omni_ns_pvalue", "omni_ns_stdpvalue",
                     "omni_is_pvalue", "omni_is_stdpvalue",
                     "omni_ls_pvalue", "omni_ls_stdpvalue",
                     "link_ns_pvalue", "link_ns_stdpvalue",
                     "link_is_pvalue", "link_is_stdpvalue",
                     "link_ls_pvalue", "link_ls_stdpvalue",
                     "form_ns_pvalue", "form_ns_stdpvalue",
                     "form_is_pvalue", "form_is_stdpvalue",
                     "form_ls_pvalue", "form_ls_stdpvalue")
    colnames(df) = df_col_names
    write.table(df, file = txt.title, sep = "\t", row.names = FALSE, col.names = TRUE)
    
    cat(paste0("Running: N=", N, ", Gamma=", gamma_0, "\n"))
    
    # --- 4. Parallel Simulation Loop ---
    # Combine results row-by-row into a single data frame
    results_df <- 
      foreach(run_ID = 1:sim_per_file, .combine = rbind, 
              .packages = c("survival", "aftgee", "afttest")) %dorng% 
      {
        temp_data = generate_data(N, gamma_0, Scenario)
        
        # 1. Rank-Based Non-Smoothed (RR-NS)
        fit_ns <- aftsrr(Surv(X, D) ~ Z1 + Z2, data = temp_data, eqType = "ns", rankWeights = "gehan")
        # 2. Rank-Based Induced-Smoothed (RR-IS)
        fit_is <- aftsrr(Surv(X, D) ~ Z1 + Z2, data = temp_data, eqType = "is", rankWeights = "gehan")
        # 3. Least Squares (LS)
        fit_ls <- aftgee(Surv(X, D) ~ Z1 + Z2, data = temp_data)
        
        # Omnibus
        omni_ns <- afttest(fit_ns, data = temp_data, npath = npath, testType = 'omnibus', npathsave = 0, linApprox = TRUE, seed = 1)
        omni_is <- afttest(fit_is, data = temp_data, npath = npath, testType = 'omnibus', npathsave = 0, linApprox = TRUE, seed = 1)
        omni_ls <- afttest(fit_ls, data = temp_data, npath = npath, testType = 'omnibus', npathsave = 0, linApprox = TRUE, seed = 1)
        
        # Link Function
        link_ns <- afttest(fit_ns, data = temp_data, npath = npath, testType = 'link', npathsave = 0, linApprox = TRUE, seed = 1)
        link_is <- afttest(fit_is, data = temp_data, npath = npath, testType = 'link', npathsave = 0, linApprox = TRUE, seed = 1)
        link_ls <- afttest(fit_ls, data = temp_data, npath = npath, testType = 'link', npathsave = 0, linApprox = TRUE, seed = 1)
        
        # Functional Form (Z2)
        form_ns <- afttest(fit_ns, data = temp_data, npath = npath, testType = 'covForm', covTested = "Z2", npathsave = 0, linApprox = TRUE, seed = 1)
        form_is <- afttest(fit_is, data = temp_data, npath = npath, testType = 'covForm', covTested = "Z2", npathsave = 0, linApprox = TRUE, seed = 1)
        form_ls <- afttest(fit_ls, data = temp_data, npath = npath, testType = 'covForm', covTested = "Z2", npathsave = 0, linApprox = TRUE, seed = 1)
        
        # This data frame is returned by the worker and automatically combined into 'results_df'
        data.frame(
          "run_ID"              = run_ID,
          
          "omni_ns_p_value"     = omni_ns$p_value,
          "omni_ns_std_p_value" = omni_ns$p_std_value,
          "omni_is_p_value"     = omni_is$p_value,
          "omni_is_std_p_value" = omni_is$p_std_value,
          "omni_ls_p_value"     = omni_ls$p_value,
          "omni_ls_std_p_value" = omni_ls$p_std_value,
          
          "link_ns_p_value"     = link_ns$p_value,
          "link_ns_std_p_value" = link_ns$p_std_value,
          "link_is_p_value"     = link_is$p_value,
          "link_is_std_p_value" = link_is$p_std_value,
          "link_ls_p_value"     = link_ls$p_value,
          "link_ls_std_p_value" = link_ls$p_std_value,
          
          "form_ns_p_value"     = form_ns$p_value,
          "form_ns_std_p_value" = form_ns$p_std_value,
          "form_is_p_value"     = form_is$p_value,
          "form_is_std_p_value" = form_is$p_std_value,
          "form_ls_p_value"     = form_ls$p_value,
          "form_ls_std_p_value" = form_ls$p_std_value
        )
      }
  }
}

# Read the generated .txt files to compute output the LaTeX table
generate_result_table()
