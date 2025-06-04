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
