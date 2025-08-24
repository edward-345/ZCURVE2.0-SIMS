source('helper_functions.R')

#----------------------------------------------
# Optional stopping
opt_stopping <- function(k_sims, n_obs = 20,
                         control_mu = 0, exp_mu = 0, st_d = 1) {
  
  #Vectors to store significant ony z-scores and all p-vals
  sig_zscores <- numeric(k_sims)
  pval_list <- list()
  
  #Counter
  B <- 1 
  
  for (i in 1:k_sims) {
      #Sample experimental and control group from normal distribution
      
    
      exp_group <- rnorm(n = n_obs, mean = exp_mu, sd = st_d)
      control <- rnorm(n = n_obs, mean = control_mu, sd = st_d)
      
      exp <- rnorm(n = 3*n_obs, mean = exp_mu, sd = st_d) #n = 60
      control <- rnorm(n = 3*n_obs, mean = control_mu, sd = st_d)
      
      exp_1 <- exp[1:(length(exp)/3)]
      control_1 <- control[1:(length(exp)/3)]
      
      exp_2 <- exp[(length(exp)/3):(2*(length(exp)/3))]
      control_2 <- control[1:(length(exp)/3) + (length(exp)/3)]
      
      exp_3 <- exp[1:((length(exp)/3)*3)]
      control_3 <- control[1:((length(exp)/3)*3)]
      
      exp_2 <- exp[]
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


multi_sample <- function(k_sims, n_obs = 20,
                         control_mu = 0, exp_mu = 0, st_d = 1) {
  
  #Vectors to store significant ony z-scores and all p-vals
  sig_zscores <- numeric(k_sims)
  total_pvals <- list()
  
  #Counter
  K <- 1 
  
  for (i in 1:k_sims) {
    sample1 <- 1:n_obs
    sample2 <- (n_obs+1):(2*n_obs)
    sample3 <- ((2*n_obs)+1):(3*n_obs)
    samples <- list(sample1 = sample1,
                    sample2 = sample2,
                    sample3 = sample3,
                    sample12 = c(sample1, sample2),
                    sample13 = c(sample1, sample3),
                    sample23 = c(sample2, sample3),
                    sample123 = c(sample1, sample2, sample3))
    
    exp_group <- rnorm(n = 3*n_obs, mean = exp_mu, sd = st_d)
    control <- rnorm(n = 3*n_obs, mean = control_mu, sd = st_d)
    
    pval_storage <- numeric(7)
    
    for (j in 1:7) {
      sub_exp <- exp_group[unlist(samples[j])]
      sub_control <- control[unlist(samples[j])]
      
      pval <- t.test(sub_exp, sub_control, var.equal=TRUE)$p.value
      pval_storage[j] <- pval
      
      total_pvals[[length(total_pvals) + 1]] <- pval
    }
    
    if (min(pval_storage) < 0.05) {
      sig_zscores[K] <- pval_converter(min(pval_storage))
      K <- K+1
    } 
    
  }
  
  total_pvals <- unlist(total_pvals)
  sig_zscores <- sig_zscores[1:(K - 1)]
  
  model <- Zing(sig_zscores)
  pval_model <- Zing(pval_converter(total_pvals))
  
  results <- list(model = model,
                 pval_model = pval_model, 
                 sig_zscores = sig_zscores,
                 total_pvals = total_pvals)
  
  return(results)
}

#To use, store it as a variable of a function
OS_null <- multi_sample(1000)


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


















