library(zcurve)
library(faux)
library(truncnorm)
library(tidyverse)
library(ggplot2)
source("helper_functions.R")
set.seed(666)

#Number of significant z-scores to be produced
k_sig <- 15000
#----------------------------------------------
#SITUATION ALPHA: Two DVs for each observation
zscores_alpha <- numeric(k_sig)
alpha <- 1

while (alpha != k_sig + 1) {
  #Generating control and exp group from N(0,1) with 2 DVs correlated by r=0.5 
  control_alpha <- rnorm_multi(
    n = 20, vars = 2, mu = c(0,0),sd = c(1,1), r = 0.5,
    varnames = c("Control1","Control2"))
  exp_alpha <- rnorm_multi(
    n = 20, vars = 2, mu = c(0,0), sd = c(1,1), r = 0.5,
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
    zvalue_alpha <- pval_converter(min_pvalue)
    zscores_alpha[alpha] <- zvalue_alpha
    alpha <- alpha + 1
  }
}

zscores_alpha

fit_alpha <- zcurve(zscores_alpha, control = list(parallel = TRUE))

alpha_plot <- plot(fit_alpha, CI = TRUE, annotation = TRUE,
                   main = "Scenario Alpha")

#-------------------------------------------------------------------------------
#SITUATION BETA: Optional Stopping
zscores_beta <- numeric(k_sig)
beta <- 1 

while (beta != k_sig + 1) {
  #Conducting one t-test after collecting 20 observations per cell 
  control_beta <- rnorm(n = 20, mean = 0, sd = 1)
  exp_beta <- rnorm(n = 20, mean = 0, sd = 1)
  result_beta <- t.test(control_beta, exp_beta, var.equal = TRUE)
  pvalue_beta <- result_beta$p.value
  
  if (pvalue_beta <= 0.05) {
    #If the result is significant, the researcher stops collecting data and 
    #reports the result
    zvalue_beta <- pval_converter(pvalue_beta)
    zscores_beta[beta] <- zvalue_beta
    beta <- beta + 1
  } else {
    #If the result is non significant, the researcher collects 10 additional 
    #observations per condition
    extracontrol_beta <- rnorm(n = 10, mean = 0, sd = 1)
    extraexp_beta <- rnorm(n = 10, mean = 0, sd = 1)
    extraresult_beta <- t.test(c(control_beta, extracontrol_beta),
                            c(exp_beta, extraexp_beta), var.equal = TRUE)
    extrapvalue_beta <- extraresult_beta$p.value
    
    if (extrapvalue_beta <= 0.05) {
      #then again tests for significance
      extrazvalue_beta <- pval_converter(extrapvalue_beta)
      zscores_beta[beta] <- extrazvalue_beta
      beta <- beta+1
    }
  }
}

zscores_beta

fit_beta <- zcurve(zscores_beta, control = list(parallel = TRUE))

beta_plot <- plot(fit_beta, CI = TRUE, annotation = TRUE,
                  main = "Scenario Beta")
#-------------------------------------------------------------------------------
#Situation Gamma: Main effect or interaction term ANCOVAs
zscores_gamma <- numeric(k_sig)
gamma <- 1

while (gamma != k_sig + 1) {
  groups <- sample(
    rep(c("control", "experimental"), each = 20))
  dv <- rnorm(n = 40, mean = 0, sd = 1)
  #Each observation was assigned a 50% probability of being female
  gender <- rbinom(n = 40, size = 1, p = 0.5)
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
    zvalue_gamma <- pval_converter(int_pvalue)
    zscores_gamma[gamma] <- zvalue_gamma
    gamma <- gamma + 1
  } else {
    if (minpval_gamma <= 0.05) {
      zvalue_gamma <- pval_converter(minpval_gamma)
      zscores_gamma[gamma] <- zvalue_gamma
      gamma <- gamma + 1
    }
  }
}

zscores_gamma

fit_gamma <- zcurve(zscores_gamma, control = list(parallel = TRUE))

gamma_plot <- plot(fit_gamma, CI = TRUE, annotation = TRUE,
                   main = "Scenario Gamma")