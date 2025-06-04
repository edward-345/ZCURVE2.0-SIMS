library(zcurve)
library(faux)
library(truncnorm)
library(tidyverse)
library(ggplot2)
source("helper_functions.R")
set.seed(666)

#Number of simulations
k_sims <- 15000
#----------------------------------------------
#General case of Z-Curve under null hypothesis
z_scores <- numeric(k_sims)
j <- 1

#Baseline is a two-condition design with 20 observations per cell
#Each of the 15,000 observations independently drawn from a normal distribution
for (i in 1:k_sims) {
  #Standard normal distribution N(0,1) since false positive occurs under 
  #the null hypothesis
  control_group <- rnorm(n = 20, mean = 0, sd = 1)
  exp_group <- rnorm(n = 20, mean = 0, sd = 1)
  
  result <- t.test(control_group, exp_group, var.equal = TRUE)
  p_value <- result$p.value
  z_value <- pval_converter(p_value)
  
  if (abs(z_value) >= qnorm(0.975)) {
    z_scores[j] <- z_value
    j <- j + 1
  }
}

z_scores <- z_scores[1:(j - 1)]

fit <- zcurve(z_scores, control = list(parallel = TRUE))

sims_plot <- plot(
  fit, CI = TRUE, annotation = TRUE, main = "Simulation under Null")

#-------------------------------------------------------------------------------
#SITUATION A: Two DVs for each observation
zscores_A <- numeric(k_sims)
A <- 1

#Vector of all 15,000 p-values generated
pvalues_scenarioA <- numeric(k_sims)

for (i in 1:k_sims) {
  #Generating control and exp group from N(0,1) with 2 DVs correlated by r=0.5 
  control_A <- rnorm_multi(
    n = 20, vars = 2, mu = c(0,0),sd = c(1,1), r = 0.5,
    varnames = c("Control1","Control2"))
  exp_A <- rnorm_multi(
    n = 20, vars = 2, mu = c(0,0), sd = c(1,1), r = 0.5,
    varnames = c("Dependent1","Dependent2"))
  
  #Helper function conducts 3 t-tests, one on each of two dependent variables 
  #and a third on the average of these two variables
  pvalue_A <- sitA_ttests(control_A$Control1, control_A$Control2,
              exp_A$Dependent1, exp_A$Dependent2)
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

A_plot <- plot(fit_A, CI = TRUE, annotation = TRUE, main = "Scenario A")

#Note that the proportion of p-values align with Simmons et al., 2011
proportions_A <- sig_pvalues(pvalues_scenarioA)

#-------------------------------------------------------------------------------
#SITUATION B: Optional Stopping
zscores_B <- numeric(k_sims)
B <- 1 

pvalues_scenarioB <- numeric(k_sims)

for (i in 1:k_sims) {
  #Conducting one t-test after collecting 20 observations per cell 
  control_B <- rnorm(n = 20, mean = 0, sd = 1)
  exp_B <- rnorm(n = 20, mean = 0, sd = 1)
  result_B <- t.test(control_B, exp_B, var.equal = TRUE)
  pvalue_B <- result_B$p.value
  
  if (pvalue_B <= 0.05) {
    #If the result is significant, the researcher stops collecting data and 
    #reports the result
    pvalues_scenarioB[i] <- pvalue_B
    zvalue_B <- pval_converter(pvalue_B)
    zscores_B[B] <- zvalue_B
    B <- B + 1
  } else {
    #If the result is non significant, the researcher collects 10 additional 
    #observations per condition
    extracontrol_B <- rnorm(n = 10, mean = 0, sd = 1)
    extraexp_B <- rnorm(n = 10, mean = 0, sd = 1)
    extraresult_B <- t.test(c(control_B, extracontrol_B),
                            c(exp_B, extraexp_B), var.equal = TRUE)
    extrapvalue_B <- extraresult_B$p.value
    
    if (extrapvalue_B <= 0.05) {
      #then again tests for significance
      pvalues_scenarioB[i] <- extrapvalue_B
      extrazvalue_B <- pval_converter(extrapvalue_B)
      zscores_B[B] <- extrazvalue_B
      B <- B+1
    } else {
      pvalues_scenarioB[i] <- extrapvalue_B
    }
  }
}

zscores_B <- zscores_B[1:(B - 1)]

fit_B <- zcurve(zscores_B, control = list(parallel = TRUE))

B_plot <- plot(fit_B, CI = TRUE, annotation = TRUE, main = "Scenario B")

#Note that the proportion of p-values align with Simmons et al., 2011
proportions_B <- sig_pvalues(pvalues_scenarioB)

#-------------------------------------------------------------------------------
#SITUATION C: Main effect or interaction term ANCOVAs
zscores_C <- numeric(k_sims)
C <- 1 

pvalues_scenarioC <- numeric(k_sims)

for (i in 1:k_sims) {
  groups <- sample(
    rep(c("control", "experimental"), each = 20))
  dv <- rnorm(n = 40, mean = 0, sd = 1)
  #Each observation was assigned a 50% probability of being female
  gender <- rbinom(n = 40, size = 1, p = 0.5)
  gender <- as.factor(
    ifelse(gender == 1, "female", "male"))
  data_C <- data.frame(dv, gender, groups)
  
  #List of p-values from ttest, main effect and interaction ANCOVAs
  sitC_pvals <- sitC_tests(data_C)
  
  min_pvalueC <- min(sitC_pvals$pvalues_C)
  int_pvalC <- sitC_pvals$int_pvalue
  
  #We report a significant effect if the effect of condition was significant in 
  #any of these analyses or if the Gender×Condition interaction was significant.
  if (int_pvalC <= 0.05) {
    pvalues_scenarioC[i] <- int_pvalC
    zvalue_C <- pval_converter(int_pvalC)
    zscores_C[C] <- zvalue_C
    C <- C + 1
  } else {
    if (min(pvalues_C) <= 0.05) {
      pvalues_scenarioC[i] <- min_pvalueC
      zvalue_C <- pval_converter(min_pvalueC)
      zscores_C[C] <- zvalue_C
      C <- C + 1
    } else {
      pvalues_scenarioC[i] <- min_pvalueC
    }
  }
  
}

zscores_C <- zscores_C[1:(C - 1)]

fit_C <- zcurve(zscores_C, control = list(parallel = TRUE))

C_plot <- plot(fit_C, CI = TRUE, annotation = TRUE, main = "Scenario C")

#Note that the proportion of p-values align with Simmons et al., 2011
proportions_C <- sig_pvalues(pvalues_scenarioC)

#-------------------------------------------------------------------------------
#SITUATION D: Ordinal test condition
zscores_D <- numeric(k_sims)
D <- 1 

pvalues_scenarioD <- numeric(k_sims)

for (i in 1:k_sims) {
  #Running three conditions (e.g., low, medium, high) 
  conditions <- sample(
    rep(c("low", "medium", "high"), length.out = 40))
  dv <- rnorm(n = 40, mean = 0, sd = 1)
  data_D <- data.frame(conditions, dv)
  
  #Conducting t tests for each of the three possible pairings of conditions 
  LowMed_pvalue <- t.test(
    dv ~ conditions,
    data = subset(data_D, conditions %in% c("low","medium")))$p.value
  LowHigh_pvalue <-t.test(
    dv ~ conditions,
    data = subset(data_D, conditions %in% c("low","high")))$p.value
  MedHigh_pvalue <- t.test(
    dv ~ conditions,
    data = subset(data_D, conditions %in% c("medium","high")))$p.value
  
  #ordinary least squares regression for the linear trend of all three 
  #conditions (coding: low = –1, medium = 0, high = 1)
  lm_coding <- ifelse(data_D$conditions == "low", -1,
                      ifelse(data_D$conditions == "medium", 0, 1))
  model_D <- lm(dv~lm_coding)
  model_pvalue <- coef(summary(model_D))["lm_coding", "Pr(>|t|)"]
  
  pvalues_D <- c(LowMed_pvalue, LowHigh_pvalue, MedHigh_pvalue,
                 model_pvalue)
  min_pvalueD <- min(pvalues_D)
  
  #Report if the lowest of all three t-tests and OLS regression is significant
  if (min_pvalueD <= 0.05) {
    pvalues_scenarioD[i] <- min_pvalueD
    zvalue_D <- pval_converter(min_pvalueD)
    zscores_D[D] <- zvalue_D
    D <- D + 1
  } else {
    pvalues_scenarioD[i] <- min_pvalueD
  }
  
}

zscores_D <- zscores_D[1:(D - 1)]

fit_D <- zcurve(zscores_D, control = list(parallel = TRUE))

D_plot <- plot(fit_D, CI = TRUE, annotation = TRUE, main = "Scenario D")

#Note that the proportion of p-values align with Simmons et al., 2011
proportions_D <- sig_pvalues(pvalues_scenarioD)

#-------------------------------------------------------------------------------











