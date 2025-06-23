library(zcurve)
library(faux)
library(truncnorm)
library(tidyverse)
library(ggplot2)
source("helper_functions.R")
set.seed(666)

#Number of significant z-scores to be produced
k_sig <- 1500
#----------------------------------------------
#Situation A,B combined, fixed number of significant z-scores
zscores_chi <- numeric(k_sig)
chi <- 1 

while (chi <= k_sig) {
  control_chi <- rnorm_multi(
    n = 20, vars = 2, mu = c(0,0),sd = c(1,1), r = 0.5,
    varnames = c("Control1","Control2"))
  exp_chi <- rnorm_multi(
    n = 20, vars = 2, mu = c(0,0), sd = c(1,1), r = 0.5,
    varnames = c("Dependent1","Dependent2"))
  
  pvalues_chi <- sitA_ttests(control_chi$Control1, control_chi$Control2,
                           exp_chi$Dependent1, exp_chi$Dependent2)
  minpval_chi <- min(pvalues_chi)
  
  if (minpval_chi <= 0.05) {
    zvalue_chi <- pval_converter(minpval_chi)
    zscores_chi[chi] <- zvalue_chi
    chi <- chi + 1
  } else {
    extracontrol_chi <- rnorm_multi(
      n = 10, vars = 2, mu = c(0,0),sd = c(1,1), r = 0.5,
      varnames = c("Control1","Control2"))
    extraexp_chi <- rnorm_multi(
      n = 10, vars = 2, mu = c(0,0),sd = c(1,1), r = 0.5,
      varnames = c("Dependent1","Dependent2"))
    
    combined_Control <- rbind(control_chi, extracontrol_chi)
    combined_Dep <- rbind(exp_chi, extraexp_chi)
    
    ExtraPvalues_chi <- sitA_ttests(combined_Control$Control1,
                                  combined_Control$Control2,
                                  combined_Dep$Dependent1,
                                  combined_Dep$Dependent2)
    min_ExtraPvalue <- min(ExtraPvalues_chi)
    
    if (min_ExtraPvalue <= 0.05) {
      extrazvalue_chi <- pval_converter(min_ExtraPvalue)
      zscores_chi[chi] <- extrazvalue_chi
      chi <- chi+1
    }
  }
}

zscores_chi

fit_chi <- zcurve(zscores_chi, control = list(parallel = TRUE))

chi_plot <- plot(fit_chi, CI = TRUE, annotation = TRUE, main = "Scenario Chi")

#-------------------------------------------------------------------------------
#Situation A,B,C combined, fixed number of significant z-scores
zscores_psi <- numeric(k_sig)
psi <- 1

while (psi <= k_sig) {
  
}















