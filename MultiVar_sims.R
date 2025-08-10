source("helper_functions.R")

#----------------------------------------------
#Situation A with Multiple DVs
multivar_sim <- function(k_sims, n = 20, r = 0.5, 
                  control_mu = c(0,0), exp_mu = c(0,0), sd = 1) {
  zscores_A <- numeric(k_sims)
  total_pvals <- list()
  A <- 1
  
  #Vector of all p-values generated
  pvalues_scenarioA <- numeric(k_sims)
  
  for (i in 1:k_sims) {
    st.d <- rep(sd, times = length(control_mu))
    
    control_A <- rnorm_multi(
      n, vars = length(control_mu), control_mu, st.d, r)
    exp_A <- rnorm_multi(
      n, vars = length(exp_mu), exp_mu, st.d, r)
    
    pvalue_A <- multivar_ttests(control_A, exp_A)
    min_pvalue <- min(pvalue_A)
    
    #Add p-value to pvalues_scenarioA regardless of significance
    pvalues_scenarioA[i] <- min_pvalue
    #Add every p-value
    total_pvals[[length(total_pvals) + length(pvalue_A)]] <- pvalue_A
    
    #If the smallest p-value of the three t-tests is significant at .05, convert 
    #to z-score and add to zscores_A vector
    if (min_pvalue <= 0.05) {
      zvalue_A <- pval_converter(min_pvalue)
      zscores_A[A] <- zvalue_A
      A <- A + 1
    }
  }
  
  #Trim unused cells
  zscores_A <- zscores_A[1:(A - 1)]
  
  fit_A <- zcurve(zscores_A, control = list(parallel = TRUE))
  
  total_pvals <- unlist(total_pvals)
  
  A_list <- list(fit_A = fit_A,
                 zscores_A = zscores_A,
                 pvalues_scenarioA = pvalues_scenarioA,
                 total_pvals = total_pvals)
  
  return(A_list)
}

#-------------------------------------------------------------------------------
# EDR of MultiVar sim
edr_m.sim <- function(k_sims, n = 20, r = 0.5, 
                         control_mu = c(0,0), exp_mu = c(0,0), sd = 1) {
  zscores_A <- numeric(k_sims)
  A <- 1
  
  for (i in 1:k_sims) {
    st.d <- rep(sd, times = length(control_mu))
    
    control_A <- rnorm_multi(
      n, vars = length(control_mu), control_mu, st.d, r)
    exp_A <- rnorm_multi(
      n, vars = length(exp_mu), exp_mu, st.d, r)
    
    pvalue_A <- multivar_ttests(control_A, exp_A)
    min_pvalue <- min(pvalue_A)
    
    #If the smallest p-value of the three t-tests is significant at .05, convert 
    #to z-score and add to zscores_A vector
    if (min_pvalue <= 0.05) {
      zvalue_A <- pval_converter(min_pvalue)
      zscores_A[A] <- zvalue_A
      A <- A + 1
    }
  }
  
  #Trim unused cells
  zscores_A <- zscores_A[1:(A - 1)]
  
  fit_A <- zcurve(zscores_A, control = list(parallel = TRUE))
  
  edr <- coef(fit_A)["EDR"] 
  
  return(edr)
}

#-------------------------------------------------------------------------------
#ERR of MultiVar sim
err_m.sim <- function(k_sims, n = 20, r = 0.5, 
                          control_mu = c(0,0), exp_mu = c(0,0), sd = 1) {
  zscores_A <- numeric(k_sims)
  A <- 1
  
  for (i in 1:k_sims) {
    st.d <- rep(sd, times = length(control_mu))
    
    control_A <- rnorm_multi(
      n, vars = length(control_mu), control_mu, st.d, r)
    exp_A <- rnorm_multi(
      n, vars = length(exp_mu), exp_mu, st.d, r)
    
    pvalue_A <- multivar_ttests(control_A, exp_A)
    min_pvalue <- min(pvalue_A)
    
    #If the smallest p-value of the three t-tests is significant at .05, convert 
    #to z-score and add to zscores_A vector
    if (min_pvalue <= 0.05) {
      zvalue_A <- pval_converter(min_pvalue)
      zscores_A[A] <- zvalue_A
      A <- A + 1
    }
  }
  
  #Trim unused cells
  zscores_A <- zscores_A[1:(A - 1)]
  
  fit_A <- zcurve(zscores_A, control = list(parallel = TRUE))
  
  err <- coef(fit_A)["ERR"] 
  
  return(err)
}
#-------------------------------------------------------------------------------
#Fixed number of z-scores
multivar_fixed <- function(k_sig, n = 20, r = 0.5, 
                           control_mu = c(0,0), exp_mu = c(0,0), sd = 1) {
  
  zscores_A <- numeric(k_sig)
  total_pvals <- list()
  pval_list <- list()
  A <- 1
  
  while (A <= k_sig) {
    st.d <- rep(sd, times = length(control_mu))
    
    control_A <- rnorm_multi(
      n, vars = length(control_mu), control_mu, st.d, r)
    exp_A <- rnorm_multi(
      n, vars = length(exp_mu), exp_mu, st.d, r)
    
    pvalue_A <- multivar_ttests(control_A, exp_A)
    min_pvalue <- min(pvalue_A)
    
    #Add p-value to pvalues_scenarioA regardless of significance
    pval_list[[length(pval_list) + 1]] <- min_pvalue
    #Add every p-value
    total_pvals[[length(total_pvals) + length(pvalue_A)]] <- pvalue_A
    
    #If the smallest p-value of the three t-tests is significant at .05, convert 
    #to z-score and add to zscores_A vector
    if (min_pvalue <= 0.05) {
      zvalue_A <- pval_converter(min_pvalue)
      zscores_A[A] <- zvalue_A
      A <- A + 1
    }
  }
  
  fit_A <- zcurve(zscores_A, control = list(parallel = TRUE))
  pvals <- unlist(pval_list)
  
  A_list <- list(fit_A = fit_A,
                 zscores_A = zscores_A,
                 pvals = pvals)
  
  return(A_list)
}

#-------------------------------------------------------------------------------

