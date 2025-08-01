source("helper_functions.R")

#----------------------------------------------
#Situation A with Multiple DVs
multivar_sim <- function(k_sims, n = 20, r = 0.5, 
                  control_mu = c(0,0), exp_mu = c(0,0), sd = 1) {
  zscores_A <- numeric(k_sims)
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
  
  A_list <- list(fit_A = fit_A,
                 zscores_A = zscores_A,
                 pvalues_scenarioA = pvalues_scenarioA)
  
  return(A_list)
}

#-------------------------------------------------------------------------------
#Fixed number of z-scores
multivar_fixed <- function(k_sig, n = 20, r = 0.5, 
                           control_mu = c(0,0), exp_mu = c(0,0), sd = 1) {
  
  zscores_A <- numeric(k_sig)
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
# 3 Dependent Variables under null
multi3 <- multivar_sim(5000, control_mu = c(0,0,0), exp_mu = c(0,0,0))
summary(multi3$fit_A)
multi3.plot <- plot(multi3$fit_A,
                    CI = TRUE, annotation = TRUE,
                    main = "3 DVs under Null Hypothesis")

multi3.pvals <- zcurve(p = multi3$pvalues_scenarioA,
                       control = list(parallel = TRUE))
summary(multi3.pvals)
multi3.pvals.plot <- plot(multi3.pvals, ymax = 10,
                          CI = TRUE, annotation = TRUE,
                          main = "3DVs P-vals under Null Hypothesis")
# ZCURVE 3.0
source(zcurve3)
ymax <- 1.2
TEST4HETEROGENEITY <- 0
TEST4BIAS <- TRUE
multi3_3.0 <- Zing(pval_converter(multi3$pvalues_scenarioA))

#----------------------------------------------
# 4 Dependent Variables
multi4 <- multivar_sim(5000, control_mu = c(0,0,0,0), exp_mu = c(0,0,0,0))
summary(multi4$fit_A)
multi4.plot <- plot(multi4$fit_A,
                    CI = TRUE, annotation = TRUE,
                    main = "4 DVs under Null Hypothesis")

multi4.pvals <- zcurve(p = multi4$pvalues_scenarioA,
                       control = list(parallel = TRUE))
summary(multi4.pvals)
multi4.pvals.plot <- plot(multi4.pvals, ymax = 10,
                          CI = TRUE, annotation = TRUE,
                          main = "4DVs P-vals under Null Hypothesis")
# ZCURVE 3.0
source(zcurve3)
ymax <- 1.2
TEST4HETEROGENEITY <- 0
TEST4BIAS <- TRUE
multi4_3.0 <- Zing(pval_converter(multi4$pvalues_scenarioA))












