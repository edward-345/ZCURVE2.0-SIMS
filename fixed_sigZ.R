library(zcurve)
library(faux)
library(truncnorm)
library(tidyverse)
library(ggplot2)
source("helper_functions.R")
set.seed(666)

#----------------------------------------------
#SITUATION ALPHA: Two DVs for each observation
alpha_sim <- function(k_sig, n = 20, r = 0.5, 
                      control_mu = c(0,0), exp_mu = c(0,0), sd = c(1,1)) {
  
  zscores_alpha <- numeric(k_sig)
  pval_list <- list()
  alpha <- 1
  
  while (alpha <= k_sig) {
    #Generating control and exp group from N(0,1) with 2 DVs correlated by r=0.5 
    control_alpha <- rnorm_multi(
      n, vars = 2, control_mu, sd, r,
      varnames = c("Control1","Control2"))
    exp_alpha <- rnorm_multi(
      n, vars = 2, exp_mu, sd, r,
      varnames = c("Dependent1","Dependent2"))
    
    #Helper function conducts 3 t-tests, one on each of two dependent variables 
    #and a third on the average of these two variables
    pvalue_alpha <- sitA_ttests(
      control_alpha$Control1, control_alpha$Control2,
      exp_alpha$Dependent1, exp_alpha$Dependent2)
    min_pvalue <- min(pvalue_alpha)
    
    #If the smallest p-value of the three t-tests is significant at .05, convert 
    #to z-score and add to zscores_alpha vector
    if (min_pvalue <= 0.05) {
      pval_list[[length(pval_list) + 1]] <- min_pvalue
      zvalue_alpha <- pval_converter(min_pvalue)
      zscores_alpha[alpha] <- zvalue_alpha
      alpha <- alpha + 1
    } else {
      pval_list[[length(pval_list) + 1]] <- min_pvalue
    }
  }
  
  fit_alpha <- zcurve(zscores_alpha, control = list(parallel = TRUE))
  pval_list <- unlist(pval_list)
  
  alpha_list <- list(fit_alpha = fit_alpha,
                     pval_list = pval_list)
  
  return(alpha_list)
}

# 500 studies under null hypothesis, both groups come from N(0,1)
alpha_500 <- alpha_sim(500)
summary(alpha_500$fit_alpha)
alpha_500.plot <- plot(alpha_500$fit_alpha,
                       CI = TRUE, annotation = TRUE, main = "Scenario Alpha")

alpha_500.pvalModel <- zcurve(p = alpha_500$pval_list,
                             control = list(parallel = TRUE))
alpha_500.pvalPlot <- plot(alpha_500.pvalModel,
                          CI = TRUE, annotation = TRUE, main = "Scenario Alpha")
alpha_500.pvals <- hist(alpha_500$pval_list)
#-------------------------------------------------------------------------------
#SITUATION BETA: Optional Stopping
beta_sim <- function(k_sig, n = 20, extra_n = 10,
                     control_mu = 0, exp_mu = 0) {
  
  zscores_beta <- numeric(k_sig)
  pval_list <- list()
  beta <- 1 
  
  while (beta <= k_sig) {
    #Conducting one t-test after collecting 20 observations per cell 
    control_beta <- rnorm(n, control_mu, sd = 1)
    exp_beta <- rnorm(n, exp_mu, sd = 1)
    result_beta <- t.test(control_beta, exp_beta, var.equal = TRUE)
    pvalue_beta <- result_beta$p.value
    
    if (pvalue_beta <= 0.05) {
      #If the result is significant, the researcher stops collecting data and 
      #reports the result
      pval_list[[length(pval_list) + 1]] <- pvalue_beta
      zvalue_beta <- pval_converter(pvalue_beta)
      zscores_beta[beta] <- zvalue_beta
      beta <- beta + 1
    } else {
      #If the result is non significant, the researcher collects 10 additional 
      #observations per condition
      extracontrol_beta <- rnorm(extra_n, control_mu, sd = 1)
      extraexp_beta <- rnorm(extra_n, exp_mu, sd = 1)
      extraresult_beta <- t.test(c(control_beta, extracontrol_beta),
                                 c(exp_beta, extraexp_beta), var.equal = TRUE)
      extrapvalue_beta <- extraresult_beta$p.value
      
      if (extrapvalue_beta <= 0.05) {
        #then again tests for significance
        pval_list[[length(pval_list) + 1]] <- extrapvalue_beta
        extrazvalue_beta <- pval_converter(extrapvalue_beta)
        zscores_beta[beta] <- extrazvalue_beta
        beta <- beta+1
      } else {
        pval_list[[length(pval_list) + 1]] <- extrapvalue_beta
      }
    }
  }
  
  fit_beta <- zcurve(zscores_beta, control = list(parallel = TRUE))
  pval_list <- unlist(pval_list)
  
  beta_list <- list(fit_beta = fit_beta,
                    pval_list = pval_list)
  
  return(beta_list)
}

# 500 studies under null hypothesis, both groups come from N(0,1)
beta_500 <- beta_sim(500)
summary(beta_500$fit_beta)
beta_500.plot <- plot(beta_500$fit_beta,
                      CI = TRUE, annotation = TRUE, main = "Scenario Beta")

beta_500.pvalModel <- zcurve(p = beta_500$pval_list,
                             control = list(parallel = TRUE))
beta_500.pvalPlot <- plot(beta_500.pvalModel,
                          CI = TRUE, annotation = TRUE, main = "Scenario Beta")
beta_500.pvals <- hist(beta_500$pval_list)
#-------------------------------------------------------------------------------
#Situation Gamma: Main effect or interaction term ANCOVAs
gamma_sim <- function(k_sig, n = 20, control_mu = 0, exp_mu = 0, sd = 1) {
  
  zscores_gamma <- numeric(k_sig)
  pval_list <- list()
  gamma <- 1
  
  while (gamma <= k_sig) {
    groups <- sample(
      rep(c("control", "experimental"), each = n))
    dv <- rnorm(length(groups),
                mean = ifelse(groups=="control", control_mu, exp_mu), sd)
    #Each observation was assigned a 50% probability of being female
    gender <- rbinom(n*2, size = 1, p = 0.5)
    gender <- as.factor(
      ifelse(gender == 1, "female", "male"))
    data_gamma <- data.frame(dv, gender, groups)
    
    #Results for Situation C were obtained by conducting a t-test...
    ttest_result <- t.test(dv ~ groups, data = data_gamma, var.equal = TRUE)
    ttest_pvalue <- ttest_result$p.value
    
    #...an analysis of covariance with a gender main effect..
    main_model <- lm(dv ~ groups + gender, data = data_gamma)
    main_pvalue <- summary(
      main_model)$coefficients["groupsexperimental", "Pr(>|t|)"]
    
    #...and an analysis of covariance with a gender interaction.
    int_model <- lm(dv ~ groups*gender, data = data_gamma)
    coefs <- coef(summary(int_model))
    int_term_name <- grep("group.*:gender.*", rownames(coefs), value = TRUE)
    
    #Edge case guard of sample consisting entirely of one gender for int_pvalue
    if (length(int_term_name) == 1) {
      int_pvalue <- coefs[int_term_name, "Pr(>|t|)"]
    } else {
      int_pvalue <- 1
    }
    
    pvalues_gamma <- c(ttest_pvalue, main_pvalue)
    
    minpval_gamma <- min(pvalues_gamma)
    
    #We report a significant effect if the effect of condition was significant in 
    #any of these analyses or if the Gender×Condition interaction was significant.
    if (int_pvalue <= 0.05) {
      pval_list[[length(pval_list) + 1]] <- int_pvalue
      zvalue_gamma <- pval_converter(int_pvalue)
      zscores_gamma[gamma] <- zvalue_gamma
      gamma <- gamma + 1
    } else {
      if (minpval_gamma <= 0.05) {
        pval_list[[length(pval_list) + 1]] <- minpval_gamma
        zvalue_gamma <- pval_converter(minpval_gamma)
        zscores_gamma[gamma] <- zvalue_gamma
        gamma <- gamma + 1
      } else {
        pval_list[[length(pval_list) + 1]] <- minpval_gamma
      }
    }
  }
  
  fit_gamma <- zcurve(zscores_gamma, control = list(parallel = TRUE))
  pval_list <- unlist(pval_list)
  
  gamma_list <- list(fit_gamma = fit_gamma,
                     pval_list = pval_list)
  
  return(gamma_list)
}

#Null hypothesis
gamma_500 <- gamma_sim(500)
summary(gamma_500$fit_gamma)
gamma_500.plot <- plot(gamma_500$fit_gamma,
                       CI = TRUE, annotation = TRUE, main = "Scenario Gamma")

gamma_500.pvalModel <- zcurve(p = gamma_500$pval_list,
                             control = list(parallel = TRUE))
gamma_500.pvalPlot <- plot(gamma_500.pvalModel,
                          CI = TRUE, annotation = TRUE, main = "Scenario Gamma")
gamma_500.pvals <- hist(gamma_500$pval_list)
#-------------------------------------------------------------------------------
#Situation Delta: Ordinal test conditions
delta_sim <- function(k_sig, n = 20, mu = 0, sd = 1) {
  
  zscores_delta <- numeric(k_sig)
  pval_list <- list()
  delta <- 1
  
  while (delta <= k_sig) {
    #Running three conditions (e.g., low, medium, high) 
    conditions <- sample(
      rep(c("low", "medium", "high"), length.out = n*2))
    dv <- rnorm(n*2,
                mean = ifelse(conditions == "low", -mu,
                              ifelse(conditions == "medium", 0, mu)), sd)
    data_delta <- data.frame(conditions, dv)
    
    #Conducting t tests for each of the three possible pairings of conditions 
    LowMed_pvalue <- t.test(
      dv ~ conditions,
      data = subset(data_delta, conditions %in% c("low","medium")))$p.value
    LowHigh_pvalue <-t.test(
      dv ~ conditions,
      data = subset(data_delta, conditions %in% c("low","high")))$p.value
    MedHigh_pvalue <- t.test(
      dv ~ conditions,
      data = subset(data_delta, conditions %in% c("medium","high")))$p.value
    
    #ordinary least squares regression for the linear trend of all three 
    #conditions (coding: low = –1, medium = 0, high = 1)
    lm_coding <- ifelse(data_delta$conditions == "low", -1,
                        ifelse(data_delta$conditions == "medium", 0, 1))
    model_delta <- lm(dv~lm_coding)
    model_pvalue <- coef(summary(model_delta))["lm_coding", "Pr(>|t|)"]
    
    pvalues_delta <- c(LowMed_pvalue, LowHigh_pvalue, MedHigh_pvalue,
                       model_pvalue)
    minpval_delta <- min(pvalues_delta)
    
    #Report if the lowest of all three t-tests and OLS regression is significant
    if (minpval_delta <= 0.05) {
      pval_list[[length(pval_list) + 1]] <- minpval_delta
      zvalue_delta <- pval_converter(minpval_delta)
      zscores_delta[delta] <- zvalue_delta
      delta <- delta + 1
    } else {
      pval_list[[length(pval_list) + 1]] <- minpval_delta
    }
  }
  
  fit_delta <- zcurve(zscores_delta, control = list(parallel = TRUE))
  pval_list <- unlist(pval_list)
  
  delta_list <- list(fit_delta = fit_delta,
                     pval_list = pval_list)
  
  return(delta_list)
}

delta_500 <- delta_sim(500)
summary(delta_500$fit_delta)
delta_500.plot <- plot(delta_500$fit_delta,
                       CI = TRUE, annotation = TRUE, main = "Scenario Delta")

delta_500.pvalModel <- zcurve(p = delta_500$pval_list,
                              control = list(parallel = TRUE))
delta_500.pvalPlot <- plot(delta_500.pvalModel,
                           CI = TRUE, annotation = TRUE, main = "Scenario Delta")
delta_500.pvals <- hist(delta_500$pval_list)





