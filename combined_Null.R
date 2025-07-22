library(zcurve)
library(faux)
library(truncnorm)
library(tidyverse)
library(ggplot2)
source("helper_functions.R")
library(tidyverse)
set.seed(666)

#ZCURVE3.0 Imports
zcurve3 <- "https://raw.githubusercontent.com/UlrichSchimmack/zcurve3.0/refs/heads/main/Zing.25.07.11.test.R"
source(zcurve3)
#----------------------------------------------
#SITUATIONS A,B COMBINED
X_sim <- function(k_sims, n = 20, extra_n = 10, r = 0.5,  
                  control_mu = c(0,0), exp_mu = c(0,0), sd = c(1,1)) {
  zscores_X <- numeric(k_sims)
  X <- 1 
  
  pvalues_scenarioX <- numeric(k_sims)
  
  for (i in 1:k_sims) {
    control_X <- rnorm_multi(
      n, vars = 2, mu = control_mu, sd, r,
      varnames = c("Control1","Control2"))
    exp_X <- rnorm_multi(
      n, vars = 2, mu = exp_mu, sd, r,
      varnames = c("Dependent1","Dependent2"))
    
    pvalues_X <- sitA_ttests(control_X$Control1, control_X$Control2,
                             exp_X$Dependent1, exp_X$Dependent2)
    min_pvalueX <- min(pvalues_X)
    
    if (min_pvalueX <= 0.05) {
      pvalues_scenarioX[i] <- min_pvalueX
      zvalue_X <- pval_converter(min_pvalueX)
      zscores_X[X] <- zvalue_X
      X <- X + 1
    } else {
      extracontrol_X <- rnorm_multi(
        extra_n, vars = 2, mu = control_mu, sd, r,
        varnames = c("Control1","Control2"))
      extraexp_X <- rnorm_multi(
        extra_n, vars = 2, mu = exp_mu, sd, r,
        varnames = c("Dependent1","Dependent2"))
      
      combined_Control <- rbind(control_X, extracontrol_X)
      combined_Dep <- rbind(exp_X, extraexp_X)
      
      ExtraPvalues_X <- sitA_ttests(combined_Control$Control1,
                                    combined_Control$Control2,
                                    combined_Dep$Dependent1,
                                    combined_Dep$Dependent2)
      min_ExtraPvalue <- min(ExtraPvalues_X)
      
      if (min_ExtraPvalue <= 0.05) {
        pvalues_scenarioX[i] <- min_ExtraPvalue
        extrazvalue_X <- pval_converter(min_ExtraPvalue)
        zscores_X[X] <- extrazvalue_X
        X <- X+1
      } else {
        pvalues_scenarioX[i] <- min_ExtraPvalue
      }
    }
  }
  
  
  zscores_X <- zscores_X[1:(X - 1)]
  
  fit_X <- zcurve(zscores_X, control = list(parallel = TRUE))
  plot(fit_X, CI = TRUE, annotation = TRUE, main = "Scenario A+B")
  X_plot <- recordPlot()
  #Note that the proportion of p-values align with Simmons et al., 2011
  proportions_X <- sig_pvalues(pvalues_scenarioX)
  
  X_list <- list(fit_X = fit_X,
                 X_plot = X_plot,
                 proportions_X = proportions_X)
  
  return(X_list)
}

#----------------------------------------------
#SITUATIONS A,B,C COMBINED
Y_sim <- function(k_sims, n = 20, extra_n = 10, r = 0.5,
                  control_mu = c(0,0), exp_mu = c(0,0), sd = c(1,1)) {
  zscores_Y <- numeric(k_sims)
  Y <- 1 
  
  pvalues_scenarioY <- numeric(k_sims)
  
  for (i in 1:k_sims) {
    control_Y <- rnorm_multi(
      n, vars = 2, mu = control_mu, sd, r,
      varnames = c("DV1","DV2"))
    exp_Y <- rnorm_multi(
      n, vars = 2, mu = exp_mu, sd, r,
      varnames = c("DV1","DV2"))
    
    control_Y$group <- "Control"
    exp_Y$group <- "Experimental"
    
    data_Y <- rbind(control_Y, exp_Y)
    
    gender <- rbinom(n*2, size = 1, p = 0.5)
    gender <- as.factor(
      ifelse(gender == 1, "Female", "Male"))
    data_Y <- data_Y %>% mutate(gender = gender)
    
    pvalues_Y <- sitA_ttests(control_Y$DV1, control_Y$DV2,
                             exp_Y$DV1, exp_Y$DV2)
    
    dv1_model <- lm(DV1 ~ group + gender, data = data_Y)
    dv1_fitP <- summary(dv1_model)$coefficients[
      "groupExperimental","Pr(>|t|)"]
    pvalues_Y <- c(pvalues_Y, dv1_fitP)
    
    dv2_model <- lm(DV2 ~ group + gender, data = data_Y)
    dv2_fitP <- summary(dv2_model)$coefficients[
      "groupExperimental","Pr(>|t|)"]
    pvalues_Y <- c(pvalues_Y, dv2_fitP)
    
    dv1_intmodel <- lm(DV1 ~ group*gender, data = data_Y)
    dv1_intfitP <- summary(dv1_intmodel)$coefficients[
      "groupExperimental:genderMale","Pr(>|t|)"]
    
    dv2_intmodel <- lm(DV2 ~ group*gender, data = data_Y)
    dv2_intfitP <- summary(dv2_intmodel)$coefficients[
      "groupExperimental:genderMale","Pr(>|t|)"]
    int_pvaluesY <- c(dv1_intfitP, dv2_intfitP)
    
    min_pvalueY <- min(pvalues_Y)
    mInt_pvalueY <- min(int_pvaluesY) 
    
    if (mInt_pvalueY <= 0.05) {
      pvalues_scenarioY[i] <- mInt_pvalueY
      zvalue_Y <- pval_converter(mInt_pvalueY)
      zscores_Y[Y] <- zvalue_Y
      Y <- Y + 1
    } else if (min_pvalueY <= 0.05){
      pvalues_scenarioY[i] <- min_pvalueY
      zvalue_Y <- pval_converter(min_pvalueY)
      zscores_Y[Y] <- zvalue_Y
      Y <- Y + 1
    } else {
      extracontrol_Y <- rnorm_multi(
        extra_n, vars = 2, mu = control_mu,sd, r,
        varnames = c("DV1","DV2"))
      extraexp_Y <- rnorm_multi(
        extra_n, vars = 2, mu = exp_mu,sd, r,
        varnames = c("DV1","DV2"))
      
      extracontrol_Y$group <- "Control"
      extraexp_Y$group <- "Experimental"
      
      extra_dataY <- rbind(extracontrol_Y, extraexp_Y)
      gender <- rbinom(extra_n*2, size = 1, p = 0.5)
      gender <- as.factor(
        ifelse(gender == 1, "Female", "Male"))
      extra_dataY <- extra_dataY %>% mutate(gender = gender)
      
      combined_dataY <- rbind(data_Y, extra_dataY)
      
      extra_pvalY <- sitA_ttests(extracontrol_Y$DV1, extracontrol_Y$DV2,
                                 extraexp_Y$DV1, extraexp_Y$DV2)
      
      Exdv1_model <- lm(DV1 ~ group + gender, data = combined_dataY)
      Exdv1_fitP <- summary(Exdv1_model)$coefficients[
        "groupExperimental","Pr(>|t|)"]
      extra_pvalY <- c(extra_pvalY, Exdv1_fitP)
      Exdv2_model <- lm(DV2 ~ group + gender, data = combined_dataY)
      Exdv2_fitP <- summary(Exdv2_model)$coefficients[
        "groupExperimental","Pr(>|t|)"]
      extra_pvalY <- c(extra_pvalY, Exdv2_fitP)
      
      Exdv1_intmodel <- lm(DV1 ~ group*gender, data = combined_dataY)
      Exdv1_intfitP <- summary(Exdv1_intmodel)$coefficients[
        "groupExperimental:genderMale","Pr(>|t|)"]
      Exdv2_intmodel <- lm(DV2 ~ group*gender, data = combined_dataY)
      Exdv2_intfitP <- summary(Exdv2_intmodel)$coefficients[
        "groupExperimental:genderMale","Pr(>|t|)"]
      Extraint_pvaLY <- c(Exdv1_intfitP, Exdv2_intfitP)
      
      Exmin_pvalY <- min(extra_pvalY)
      ExmInt_pvalY <- min(Extraint_pvaLY)
      
      if (ExmInt_pvalY <= 0.05) {
        pvalues_scenarioY[i] <- ExmInt_pvalY
        extrazvalue_Y <- pval_converter(ExmInt_pvalY)
        zscores_Y[Y] <- extrazvalue_Y
        Y <- Y+1
      } else if (Exmin_pvalY <= 0.05) {
        pvalues_scenarioY[i] <- Exmin_pvalY
        extrazvalue_Y <- pval_converter(Exmin_pvalY)
        zscores_Y[Y] <- extrazvalue_Y
        Y <- Y+1
      } else {
        all_pvalY <- c(extra_pvalY, Extraint_pvaLY)
        pvalues_scenarioY[i] <- min(all_pvalY)
      }
    }
  }
  
  zscores_Y <- zscores_Y[1:(Y - 1)]
  
  fit_Y <- zcurve(zscores_Y, control = list(parallel = TRUE))
  plot(fit_Y, CI = TRUE, annotation = TRUE, main = "Scenario A+B+C")
  Y_plot <- recordPlot()
  #Note that the proportion of p-values align with Simmons et al., 2011
  proportions_Y <- sig_pvalues(pvalues_scenarioY)
  
  Y_list <- list(fit_Y = fit_Y,
                 Y_plot = Y_plot,
                 proportions_Y = proportions_Y)
  
  return(Y_list)
}

#----------------------------------------------
#SITUATIONS A,B,C,D COMBINED
Z_sim <- function(k_sims, n = 20, extra_n = 10, r = 0.5,
                  mu = c(0,0), sd = c(1,1)) {
  zscores_Z <- numeric(k_sims)
  Z <- 1 
  
  pvalues_scenarioZ <- numeric(k_sims)
  
  for (i in 1:k_sims) {
    #Generating data frame of Control group DV1, DV2 outcome values for Sit A
    data_Z <- rnorm_multi(
      n*2, vars = 2, mu, sd, r,
      varnames = c("DV1","DV2"))
    
    #Column of Sit D "conditions" variables added
    data_Z$conditions <- sample(rep(c("low", "medium", "high"), length.out = n*2))
    
    #Adding Sit C "gender" variable for each observation
    gender <- rbinom(n*2, size = 1, p = 0.5)
    gender <- as.factor(
      ifelse(gender == 1, "Female", "Male"))
    data_Z <- data_Z %>% mutate(gender = gender)
    
    #T-tests for all 3 combinations of conditions (Sit D) for DV1 and DV2 (Sit A)
    pvalues_Z <- sitZ_ttests(data_Z)
    
    #Gender main effect ANCOVA for DV1 (Sit C)
    dv1_model <- lm(DV1 ~ conditions + gender, data = data_Z)
    dv1_fitP <- summary(dv1_model)$coefficients[
      grep("^conditions", rownames(summary(dv1_model)$coefficients)),
      "Pr(>|t|)"]
    pvalues_Z <- c(pvalues_Z, dv1_fitP)
    
    #Gender main effect ANCOVA for DV2 (Sit C)
    dv2_model <- lm(DV2 ~ conditions + gender, data = data_Z)
    dv2_fitP <- summary(dv2_model)$coefficients[
      grep("^conditions", rownames(summary(dv2_model)$coefficients)),
      "Pr(>|t|)"]
    pvalues_Z <- c(pvalues_Z, dv2_fitP)
    
    #Gender*group interaction effect ANCOVA for DV1 (Sit C)
    dv1_intmodel <- lm(DV1 ~ conditions*gender, data = data_Z)
    dv1_intfitP <- summary(
      dv1_intmodel)$coefficients[grep("conditions.*:genderMale",
                                      rownames(summary(
                                        dv1_intmodel)$coefficients)), "Pr(>|t|)"]
    
    #Gender*group interaction effect ANCOVA for DV2 (Sit C)
    dv2_intmodel <- lm(DV2 ~ conditions*gender, data = data_Z)
    dv2_intfitP <- summary(
      dv2_intmodel)$coefficients[grep("conditions.*:genderMale",
                                      rownames(summary(
                                        dv2_intmodel)$coefficients)), "Pr(>|t|)"]
    
    #Interaction term ANCOVA model p-values
    int_pvaluesZ <- c(dv1_intfitP, dv2_intfitP)
    
    #OLS regression for the linear trend of conditions for DV1 (Sit D)
    lm_coding <- ifelse(data_Z$conditions == "low", -1,
                        ifelse(data_Z$conditions == "medium", 0, 1))
    
    OLS_modelDV1 <- lm(DV1~lm_coding, data = data_Z)
    DV1model_pvalue <- coef(summary(OLS_modelDV1))["lm_coding", "Pr(>|t|)"]
    
    OLS_modelDV2 <- lm(DV2~lm_coding, data = data_Z)
    DV2model_pvalue <- coef(summary(OLS_modelDV2))["lm_coding", "Pr(>|t|)"]
    
    pvalues_Z <- c(pvalues_Z, DV1model_pvalue, DV2model_pvalue)
    
    #Minimum p-values
    min_pvalueZ <- min(pvalues_Z)
    mInt_pvalueZ <- min(int_pvaluesZ)
    
    if (mInt_pvalueZ <= 0.05) {
      pvalues_scenarioZ[i] <- mInt_pvalueZ
      zvalue_Z <- pval_converter(mInt_pvalueZ)
      zscores_Z[Z] <- zvalue_Z
      Z <- Z + 1
    } else if (min_pvalueZ <= 0.05){
      #"if the effect of condition was significant in any of these analyses"
      pvalues_scenarioZ[i] <- min_pvalueZ
      zvalue_Z <- pval_converter(min_pvalueZ)
      zscores_Z[Z] <- zvalue_Z
      Z <- Z + 1
    } else {
      #Adding 10 observations per cell from Situation B if non significant
      extradata_Z <- rnorm_multi(
        extra_n*2, vars = 2, mu, sd, r,
        varnames = c("DV1","DV2"))
      
      #Column of Sit D "conditions" variables added
      extradata_Z$conditions <- sample(
        rep(c("low", "medium", "high"), length.out = extra_n*2))
      
      #Adding Sit C "gender" variable for each observation
      gender <- rbinom(n = extra_n*2, size = 1, p = 0.5)
      gender <- as.factor(
        ifelse(gender == 1, "Female", "Male"))
      extradata_Z <- extradata_Z %>% mutate(gender = gender)
      
      #Combining additional observations to current data set
      combined_dataZ <- rbind(data_Z, extradata_Z)
      
      #T-tests for all 3 combinations of conditions (Sit D) for DV1 and DV2 (Sit A)
      extra_pvalsZ <- sitZ_ttests(combined_dataZ)
      
      #Gender main effect ANCOVA for DV1 (Sit C)
      dv1_model <- lm(DV1 ~ conditions + gender, data = combined_dataZ)
      dv1_fitP <- summary(dv1_model)$coefficients[
        grep("^conditions", rownames(summary(dv1_model)$coefficients)),
        "Pr(>|t|)"]
      extra_pvalsZ <- c(extra_pvalsZ, dv1_fitP)
      
      #Gender main effect ANCOVA for DV2 (Sit C)
      dv2_model <- lm(DV2 ~ conditions + gender, data = combined_dataZ)
      dv2_fitP <- summary(dv2_model)$coefficients[
        grep("^conditions", rownames(summary(dv2_model)$coefficients)),
        "Pr(>|t|)"]
      extra_pvalsZ <- c(extra_pvalsZ, dv2_fitP)
      
      #Gender*group interaction effect ANCOVA for DV1 (Sit C)
      dv1_intmodel <- lm(DV1 ~ conditions*gender, data = combined_dataZ)
      dv1_intfitP <- summary(
        dv1_intmodel)$coefficients[grep("conditions.*:genderMale",
                                        rownames(summary(
                                          dv1_intmodel)$coefficients)),
                                   "Pr(>|t|)"]
      
      #Gender*group interaction effect ANCOVA for DV2 (Sit C)
      dv2_intmodel <- lm(DV2 ~ conditions*gender, data = combined_dataZ)
      dv2_intfitP <- summary(
        dv2_intmodel)$coefficients[grep("conditions.*:genderMale",
                                        rownames(summary(
                                          dv2_intmodel)$coefficients)),
                                   "Pr(>|t|)"]
      
      #Interaction term ANCOVA model p-values
      extraInt_pvalZ <- c(dv1_intfitP, dv2_intfitP)
      
      #OLS regression for the linear trend of conditions for DV1 (Sit D)
      lm_coding <- ifelse(combined_dataZ$conditions == "low", -1,
                          ifelse(combined_dataZ$conditions == "medium", 0, 1))
      
      OLS_modelDV1 <- lm(DV1~lm_coding, data = combined_dataZ)
      DV1model_pvalue <- coef(summary(OLS_modelDV1))["lm_coding", "Pr(>|t|)"]
      
      OLS_modelDV2 <- lm(DV2~lm_coding, data = combined_dataZ)
      DV2model_pvalue <- coef(summary(OLS_modelDV2))["lm_coding", "Pr(>|t|)"]
      
      extra_pvalsZ <- c(extra_pvalsZ, DV1model_pvalue, DV2model_pvalue)
      
      #Minimum p-values
      ExMin_pvalueZ <- min(extra_pvalsZ)
      ExMInt_pvalueZ <- min(extraInt_pvalZ)
      
      if (ExMInt_pvalueZ <= 0.05) {
        pvalues_scenarioZ[i] <- ExMInt_pvalueZ
        zvalue_Z <- pval_converter(ExMInt_pvalueZ)
        zscores_Z[Z] <- zvalue_Z
        Z <- Z + 1
      } else if (ExMin_pvalueZ <= 0.05) {
        #"if the effect of condition was significant in any of these analyses"
        pvalues_scenarioZ[i] <- ExMin_pvalueZ
        zvalue_Z <- pval_converter(ExMin_pvalueZ)
        zscores_Z[Z] <- zvalue_Z
        Z <- Z + 1
      } else {
        #If non are significant, add lowest p-value to pvalues_ScenarioZ vector
        all_pvalZ <- c(ExMin_pvalueZ, ExMInt_pvalueZ)
        pvalues_scenarioZ[i] <- min(all_pvalZ)
      }
    }
  }
  
  zscores_Z <- zscores_Z[1:(Z - 1)]
  
  fit_Z <- zcurve(zscores_Z, control = list(parallel = TRUE))
  plot(fit_Z, CI = TRUE, annotation = TRUE, main = "Scenario A+B+C+D")
  Z_plot <- recordPlot()
  #Note that the proportion of p-values align with Simmons et al., 2011
  proportions_Z <- sig_pvalues(pvalues_scenarioZ)
  
  Z_list <- list(fit_Z = fit_Z,
                 Z_plot = Z_plot,
                 proportions_Z = proportions_Z)
  
  return(Z_list)
}

#-------------------------------------------------------------------------------
