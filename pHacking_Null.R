library(zcurve)
library(faux)
library(truncnorm)
library(tidyverse)
library(ggplot2)
source("helper_functions.R")
set.seed(666)

#ZCURVE3.0 Imports
zcurve3 <- "https://raw.githubusercontent.com/UlrichSchimmack/zcurve3.0/refs/heads/main/Zing.25.07.11.test.R"
source(zcurve3)

#----------------------------------------------
#General case of Z-Curve under null hypothesis
#Baseline is a two-condition design with 20 observations per cell (n = 20)
simulation <- function(k_sims, n) {
  z_scores <- numeric(k_sims)
  j <- 1
  for (i in 1:k_sims) {
    #Standard normal distribution N(0,1) since false positive occurs under null
    control_group <- rnorm(n, mean = 0, sd = 1)
    exp_group <- rnorm(n, mean = 0, sd = 1)
    
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
  #sim_plot <- plot(fit, CI = TRUE, annotation = TRUE, main = "Simulation under Null")
  return(summary(fit))
}


#-------------------------------------------------------------------------------
#SITUATION A: Two DVs for each observation
A_sim <- function(k_sims, n = 20, r = 0.5, 
                  control_mu = c(0,0), exp_mu = c(0,0), sd = c(1,1)) {
  zscores_A <- numeric(k_sims)
  A <- 1
  
  #Vector of all p-values generated
  pvalues_scenarioA <- numeric(k_sims)
  
  for (i in 1:k_sims) {
    #Generating control and exp group from N(0,1) with 2 DVs correlated by r=0.5 
    control_A <- rnorm_multi(
      n, vars = 2, control_mu, sd = c(1,1), r,
      varnames = c("Control1","Control2"))
    exp_A <- rnorm_multi(
      n, vars = 2, exp_mu, sd = c(1,1), r,
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
  
  A_list <- list(fit_A = fit_A,
                 zscores_A = zscores_A,
                 pvalues_scenarioA = pvalues_scenarioA)
  
  return(A_list)
}

#-------------------------------------------------------------------------------
#SITUATION B: Optional Stopping
B_sim <- function(k_sims, n = 20, extra_n = 10,
                  control_mu = 0, exp_mu = 0, sd = 1) {
  zscores_B <- numeric(k_sims)
  B <- 1 
  
  pvalues_scenarioB <- numeric(k_sims)
  
  for (i in 1:k_sims) {
    #Conducting one t-test after collecting 20 observations per cell 
    control_B <- rnorm(n, control_mu, sd)
    exp_B <- rnorm(n, exp_mu, sd)
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
      extracontrol_B <- rnorm(extra_n, control_mu, sd)
      extraexp_B <- rnorm(extra_n, exp_mu, sd)
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
  
  B_list <- list(fit_B = fit_B,
                 zscores_B = zscores_B,
                 pvalues_scenarioB = pvalues_scenarioB)
  
  return(B_list)
}

#-------------------------------------------------------------------------------
#SITUATION C: Main effect or interaction term ANCOVAs
C_sim <- function(k_sims, n = 20, control_mu = 0, exp_mu = 0, sd = 1) {
  zscores_C <- numeric(k_sims)
  C <- 1 
  
  pvalues_scenarioC <- numeric(k_sims)
  
  for (i in 1:k_sims) {
    groups <- sample(
      rep(c("control", "experimental"), each = n))
    dv <- rnorm(n*2, 
                mean = ifelse(groups=="control", control_mu, exp_mu), sd)
    #Each observation was assigned a 50% probability of being female
    gender <- rbinom(n*2, size = 1, p = 0.5)
    gender <- as.factor(
      ifelse(gender == 1, "female", "male"))
    data_C <- data.frame(dv, gender, groups)
    
    #Results for Situation C were obtained by conducting a t-test...
    ttest_result <- t.test(dv ~ groups, data = data_C, var.equal = TRUE)
    ttest_pvalue <- ttest_result$p.value
    
    #...an analysis of covariance with a gender main effect..
    main_model <- lm(dv ~ groups + gender, data = data_C)
    main_pvalue <- summary(
      main_model)$coefficients["groupsexperimental", "Pr(>|t|)"]
    
    #...and an analysis of covariance with a gender interaction.
    int_model <- lm(dv ~ groups*gender, data = data_C)
    coefs <- coef(summary(int_model))
    int_term_name <- grep("group.*:gender.*", rownames(coefs), value = TRUE)
    
    #Edge case guard of sample consisting entirely of one gender for int_pvalue
    if (length(int_term_name) == 1) {
      int_pvalue <- coefs[int_term_name, "Pr(>|t|)"]
    } else {
      int_pvalue <- 1
    }
    
    pvalues_C <- c(ttest_pvalue, main_pvalue)
    
    min_pvalueC <- min(pvalues_C)
    
    #We report a significant effect if the effect of condition was significant in 
    #any of these analyses or if the Gender×Condition interaction was significant.
    if (int_pvalue <= 0.05) {
      pvalues_scenarioC[i] <- int_pvalue
      zvalue_C <- pval_converter(int_pvalue)
      zscores_C[C] <- zvalue_C
      C <- C + 1
    } else {
      if (min_pvalueC <= 0.05) {
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
  
  C_list <- list(fit_C = fit_C,
                 zscores_C = zscores_C,
                 pvalues_scenarioC = pvalues_scenarioC)
  
  return(C_list)
}

#-------------------------------------------------------------------------------
#SITUATION D: Ordinal test condition
D_sim <- function(k_sims, n = 20, mu = 0, sd = 1) {
  zscores_D <- numeric(k_sims)
  D <- 1 
  
  pvalues_scenarioD <- numeric(k_sims)
  
  for (i in 1:k_sims) {
    #Running three conditions (e.g., low, medium, high) 
    conditions <- sample(
      rep(c("low", "medium", "high"), length.out = n*2))
    dv <- rnorm(n*2, mean = mu, sd)
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
  
  D_list <- list(fit_D = fit_D,
                 zscores_D = zscores_D, 
                 pvalues_scenarioD = pvalues_scenarioD)
  
  return(D_list)
}

#-------------------------------------------------------------------------------











