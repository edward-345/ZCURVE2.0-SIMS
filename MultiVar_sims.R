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

#----------------------------------------------
# 5 Dependent Variables
multi5 <- multivar_sim(5000, control_mu = c(0,0,0,0,0), exp_mu = c(0,0,0,0,0))
summary(multi5$fit_A)
multi5.plot <- plot(multi5$fit_A,
                    CI = TRUE, annotation = TRUE,
                    main = "5 DVs under Null Hypothesis")

multi5.pvals <- zcurve(p = multi5$pvalues_scenarioA,
                       control = list(parallel = TRUE))
summary(multi5.pvals)
multi5.pvals.plot <- plot(multi5.pvals, ymax = 10,
                          CI = TRUE, annotation = TRUE,
                          main = "5DVs P-vals under Null Hypothesis")
# ZCURVE 3.0
source(zcurve3)
ymax <- 1.2
TEST4HETEROGENEITY <- 0
TEST4BIAS <- TRUE
multi5_3.0 <- Zing(pval_converter(multi5$pvalues_scenarioA))

#----------------------------------------------
# 5 Dependent Variables with Fixed
multi5_fixed <- multivar_fixed(5000,
                               control_mu = c(0,0,0,0,0),
                               exp_mu = c(0,0,0,0,0))
summary(multi5_fixed$fit_A)
multi5.plot <- plot(multi5_fixed$fit_A,
                    CI = TRUE, annotation = TRUE,
                    main = "5 DVs under Null Hypothesis (Fixed)")

multi5.pvals <- zcurve(p = multi5_fixed$pvals,
                       control = list(parallel = TRUE))
summary(multi5.pvals)
multi5.pvals.plot <- plot(multi5.pvals, ymax = 10,
                          CI = TRUE, annotation = TRUE,
                          main = "5DVs P-vals under Null Hypothesis")
# ZCURVE 3.0
source(zcurve3)
ymax <- 1.2
TEST4HETEROGENEITY <- 0
TEST4BIAS <- TRUE
multi5_fixed.3.0 <- Zing(pval_converter(multi5_fixed$pvals))

#----------------------------------------------
# 5 Dependent Variables True Effect
multi5_med <- multivar_sim(5000, control_mu = c(0,0,0,0,0),
                       exp_mu = c(0.5,0.5,0.5,0.5,0.5))
summary(multi5_med$fit_A)
multi5_med.plot <- plot(multi5_med$fit_A,
                    CI = TRUE, annotation = TRUE,
                    main = "5 DVs under Null Hypothesis")

multi5_med.pvals <- zcurve(p = multi5_med$pvalues_scenarioA,
                       control = list(parallel = TRUE))
summary(multi5_med.pvals)
multi5_med.pvals.plot <- plot(multi5_med.pvals, ymax = 10,
                          CI = TRUE, annotation = TRUE,
                          main = "5DVs P-vals under Null Hypothesis")
# ZCURVE 3.0
source(zcurve3)
ymax <- 1.2
TEST4HETEROGENEITY <- 0
TEST4BIAS <- TRUE
multi5_3.0 <- Zing(pval_converter(multi5_med$pvalues_scenarioA))

#----------------------------------------------
tester_5000 <- multivar_sim(1000, n = 40,
                            control_mu = rep(0, 5), exp_mu = rep(0.5, 5))
summary(tester_5000$fit_A)
tester_5000.plot <- plot(tester_5000$fit_A,
                        CI = TRUE, annotation = TRUE,
                        main = "5 DVs under Null Hypothesis")

tester_5000.pvals <- zcurve(p = tester_5000$total_pvals,
                           control = list(parallel = TRUE))
summary(tester_5000.pvals)
tester_5000.pvals.plot <- plot(tester_5000.pvals, ymax = 10,
                              CI = TRUE, annotation = TRUE,
                              main = "Total Pvals with 5DVunder Null Hypothesis")
# ZCURVE 3.0
source(zcurve3)
ymax <- 0.8
TEST4HETEROGENEITY <- 0
#TEST4BIAS <- TRUE
Title <- paste("EDR: ", round(true.edr*100), " ERR: ", round(true.err*100))
sim.z <- pval_converter(tester_5000$total_pvals)
results_tester1000 <- Zing(sim.z)

#----------------------------------------------

corr_test <- multivar_sim(3000, n = 40,
                          control_mu = rep(0, 5), exp_mu = rep(0.5, 5), r=0.7)
summary(corr_test$fit_A)
corr_test.plot <- plot(corr_test$fit_A,
                         CI = TRUE, annotation = TRUE,
                         main = "5DVs correlated by 0.9")

corr_test.pvals <- zcurve(p = corr_test$pvalues_scenarioA,
                            control = list(parallel = TRUE))
summary(corr_test.pvals)
corr_test.pvals.plot <- plot(corr_test.pvals, ymax = 10,
                               CI = TRUE, annotation = TRUE,
                               main = "5DVs correlated by 0.9 (pvals)")

corr_test.total <- zcurve(p = corr_test$total_pvals,
                          control = list(parallel = TRUE))
summary(corr_test.total)
corr_test.total.plot <- plot(corr_test.pvals, ymax = 10,
                             CI = TRUE, annotation = TRUE,
                             main = "5DVs correlated by 0.9 (total)")

# ZCURVE 3.0
source(zcurve3)
ymax <- 1.2
TEST4HETEROGENEITY <- 0
#TEST4BIAS <- TRUE
#Title <- paste("EDR: ", round(true.edr*100), " ERR: ", round(true.err*100))
sim.z <- pval_converter(corr_test$total_pvals)
#zsds = rep(0.7, 7)
results_corr_test <- Zing(sim.z)

source(zcurve3)
ymax <- 1.2
TEST4HETEROGENEITY <- 0
#TEST4BIAS <- TRUE
#Title <- paste("EDR: ", round(true.edr*100), " ERR: ", round(true.err*100))
sim.z <- pval_converter(corr_test$total_pvals)
Est.Method = "EXT"
ncz = 2
zsds = rep(5, 1)
Int.Beg = 0
new_corr <- Zing(sim.z)
