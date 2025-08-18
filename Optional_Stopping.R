source('helper_functions.R')

#----------------------------------------------
# Optional stopping
opt_stopping <- function(k_sims, n_obs = 20, extra_n = 10,
                         control_mu = 0, exp_mu = 0, st_d = 1) {
  
  #Vectors to store significant ony z-scores and all p-vals
  sig_zscores <- numeric(k_sims)
  pval_list <- list()
  
  #Counter
  B <- 1 
  
  for (i in 1:k_sims) {
    # Resample both groups until signifcant
    repeat {
      #Sample experimental and control group from normal distribution
      pre_exp <- rnorm(n = n_obs, mean = exp_mu, sd = st_d)
      pre_control <- rnorm(n = n_obs, mean = control_mu, sd = st_d)
      
      #P-value from two-sided two sample t-test
      pre_pval <- t.test(pre_exp, pre_control, var.equal = TRUE)$p.value
      #Add p-value to pval_list regardless of significance
      pval_list[[length(pval_list) + 1]] <- pre_pval
      if (pre_pval <= 0.05) break
    }
    
    # Sampling 10 additional participants per group
    extra_exp <- rnorm(n = extra_n, mean = exp_mu, sd = st_d)
    extra_control <- rnorm(n = extra_n, mean = control_mu, sd = st_d)
    
    experimental <- c(pre_exp, extra_exp)
    control <- c(pre_control, extra_control)
    
    #T-Test with 10 additional
    pval <- t.test(experimental, control, var.equal = TRUE)$p.value
    #Add p-value to pval_list regardless of significance
    pval_list[[length(pval_list) + 1]] <- pval
    
    #If increasing of sample size reduced p-value, use increased sample size
    if (pval <= pre_pval && pval <= 0.05) {
      zscore <- pval_converter(pval)
      sig_zscores[B] <- zscore
      B <- B + 1
      
      # Otherwise use original
    } else if (pval > pre_pval && pre_pval <= 0.05) {
      zscore <- pval_converter(pre_pval)
      sig_zscores[B] <- zscore
      B <- B + 1
    }
  }
  
  #Trim vector
  sig_zscores <- sig_zscores[1:(B - 1)]
  #Fit model
  model <- zcurve(sig_zscores, control = list(parallel = TRUE))
  
  all_pvals <- unlist(pval_list)
  pval_model <- zcurve(p = all_pvals, control = list(parallel = TRUE))
  
  results <- list(sig_zscores = sig_zscores,
                  model = model,
                  pval_model = pval_model,
                  all_pvals = all_pvals)
  
  return(results)
}

######################----------------------------------------------------------
# ZCURVE 2.0 RESULTS
######################----------------------------------------------------------
# Null Hypothesis true
OptSt_null <- opt_stopping(1000)

summary(OptSt_null$model)
OptSt_null.plot <- plot(OptSt_null$model,
                    CI = TRUE, annotation = TRUE,
                    main = "Optional Stopping under Null Hypothesis")

summary(OptSt_null$pval_model)
OptSt_null.plot <- plot(OptSt_null$pval_model,
                        CI = TRUE, annotation = TRUE,
                        main = "Optional Stopping under Null with Total Pvals")

#-------------------------------------------------------------------------------
#Weak ES = 0.2
OptSt_weak <- opt_stopping(1000, exp_mu = 0.2)

summary(OptSt_weak$model)
OptSt_weak.plot <- plot(OptSt_weak$model,
                        CI = TRUE, annotation = TRUE,
                        main = "Optional Stopping Weak ES = .2")

summary(OptSt_weak$pval_model)
OptSt_weak.plot <- plot(OptSt_weak$pval_model,
                        CI = TRUE, annotation = TRUE,
                        main = "Optional Stopping Weak ES = .2 with Total Pvals")

#-------------------------------------------------------------------------------
#Med ES = 0.5
OptSt_med <- opt_stopping(1000, exp_mu = 0.5)

summary(OptSt_med$model)
OptSt_med.plot <- plot(OptSt_med$model,
                        CI = TRUE, annotation = TRUE,
                        main = "Optional Stopping Med ES = .5")

summary(OptSt_med$pval_model)
OptSt_med.plot <- plot(OptSt_med$pval_model,
                        CI = TRUE, annotation = TRUE,
                        main = "Optional Stopping Med ES = .5 with Total Pvals")

#-------------------------------------------------------------------------------
#Strong ES = 0.8
OptSt_strong <- opt_stopping(1000, exp_mu = 0.8)

summary(OptSt_strong$model)
OptSt_strong.plot <- plot(OptSt_strong$model,
                       CI = TRUE, annotation = TRUE,
                       main = "Optional Stopping Strong ES = .8")

summary(OptSt_strong$pval_model)
OptSt_strong.plot <- plot(OptSt_strong$pval_model,
                       CI = TRUE, annotation = TRUE,
                       main = "Optional Stopping Strong ES = .8 with Total Pvals")


######################----------------------------------------------------------
# ZCURVE 3.0 RESULTS
zcurve3 <- "https://raw.githubusercontent.com/UlrichSchimmack/zcurve3.0/refs/heads/main/Zing.25.07.11.test.R"
source(zcurve3)

ymax <- 1.8
# TEST4HETEROGENEITY <- 0
TEST4BIAS <- TRUE
# Title <- past("title string")
Est.method <- "EXT"
# boot.iter <- 500
ncz <- c(0:6)
#W.FIXED <- TRUE
#w.fix   <- 
# Component locations (z-values at which densities are centered)
components <- length(ncz)           # Number of components
zsd <- 0.5                            # SD of standard normal z-distribution
zsds = rep(zsd,components)          # one SD for each component
# Int.Beg <- 1.96
# Int.End <- 6

######################----------------------------------------------------------
# Null Hypothesis true
ymax <- 1.9
optst.null <- Zing(OptSt_null$sig_zscores)

ymax <- 0.9
optst.null_pvals <- Zing(pval_converter(OptSt_null$all_pvals))

#-------------------------------------------------------------------------------
# Weak ES = 0.2
ymax <- 1.9
optst.weak <- Zing(OptSt_weak$sig_zscores)

ymax <- 0.9
optst.weak_pvals <- Zing(pval_converter(OptSt_weak$all_pvals))

#-------------------------------------------------------------------------------
# Med ES = 0.5
ymax <- 1.0
optst.med <- Zing(OptSt_med$sig_zscores)

ymax <- 0.5
optst.med_pvals <- Zing(pval_converter(OptSt_med$all_pvals))

#-------------------------------------------------------------------------------
# Strong ES = 0.8
ymax <- 0.8
optst.strong <- Zing(OptSt_strong$sig_zscores)

ymax <- 0.6
optst.strong_pvals <- Zing(pval_converter(OptSt_strong$all_pvals))


















