library(zcurve)
library(faux)
library(truncnorm)
library(tidyverse)
library(ggplot2)
source("helper_functions.R")
set.seed(666)

#----------------------------------------------
#Situation A,B combined, fixed number of significant z-scores
chi_sim <- function(k_sig, n = 20, extra_n = 10, r = 0.5,  
                    control_mu = c(0,0), exp_mu = c(0,0), sd = c(1,1)) {
  zscores_chi <- numeric(k_sig)
  chi <- 1 
  
  while (chi <= k_sig) {
    control_chi <- rnorm_multi(
      n, vars = 2, control_mu, sd, r,
      varnames = c("Control1","Control2"))
    exp_chi <- rnorm_multi(
      n, vars = 2, exp_mu, sd, r,
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
        extra_n, vars = 2, control_mu, sd, r,
        varnames = c("Control1","Control2"))
      extraexp_chi <- rnorm_multi(
        extra_n, vars = 2, exp_mu, sd, r,
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
  
  fit_chi <- zcurve(zscores_chi, control = list(parallel = TRUE))
  plot(fit_chi, CI = TRUE, annotation = TRUE, main = "Scenario Chi")
  chi_plot <- recordPlot()
  
  chi_list <- list(fit_chi = fit_chi,
                    chi_plot = chi_plot)
  
  return(chi_list)
}

chi_500 <- chi_sim(500)
summary(chi_500$fit_chi)

chi_alt <- chi_sim(500, exp_mu = c(0.2, 0.2))
summary(chi_alt$fit_chi)

chi_alt_strong <- chi_sim(500, exp_mu = c(0.8, 0.8))
summary(chi_alt_strong$fit_chi)
#-------------------------------------------------------------------------------
#Situation A,B,C combined, fixed number of significant z-scores
psi_sim <- function(k_sig, n = 20, extra_n = 10, r = 0.5,
                    control_mu = c(0,0), exp_mu = c(0,0), sd = c(1,1)) {
  zscores_psi <- numeric(k_sig)
  psi <- 1
  
  while (psi <= k_sig) {
    control_psi <- rnorm_multi(
      n, vars = 2, control_mu, sd, r,
      varnames = c("DV1","DV2"))
    exp_psi <- rnorm_multi(
      n, vars = 2, exp_mu, sd, r,
      varnames = c("DV1","DV2"))
    
    control_psi$group <- "Control"
    exp_psi$group <- "Experimental"
    
    data_psi <- rbind(control_psi, exp_psi)
    
    gender <- rbinom(n*2, size = 1, p = 0.5)
    gender <- as.factor(
      ifelse(gender == 1, "Female", "Male"))
    data_psi <- data_psi %>% mutate(gender = gender)
    
    pvalues_psi <- sitA_ttests(control_psi$DV1, control_psi$DV2,
                               exp_psi$DV1, exp_psi$DV2)
    
    dv1_model <- lm(DV1 ~ group + gender, data = data_psi)
    dv1_fitP <- summary(dv1_model)$coefficients[
      "groupExperimental","Pr(>|t|)"]
    pvalues_psi <- c(pvalues_psi, dv1_fitP)
    
    dv2_model <- lm(DV2 ~ group + gender, data = data_psi)
    dv2_fitP <- summary(dv2_model)$coefficients[
      "groupExperimental","Pr(>|t|)"]
    pvalues_psi <- c(pvalues_psi, dv2_fitP)
    
    dv1_intmodel <- lm(DV1 ~ group*gender, data = data_psi)
    dv1_intfitP <- summary(dv1_intmodel)$coefficients[
      "groupExperimental:genderMale","Pr(>|t|)"]
    
    dv2_intmodel <- lm(DV2 ~ group*gender, data = data_psi)
    dv2_intfitP <- summary(dv2_intmodel)$coefficients[
      "groupExperimental:genderMale","Pr(>|t|)"]
    int_pvaluespsi <- c(dv1_intfitP, dv2_intfitP)
    
    minpval_psi <- min(pvalues_psi)
    mInt.pval_psi <- min(int_pvaluespsi) 
    
    if (mInt.pval_psi <= 0.05) {
      zvalue_psi <- pval_converter(mInt.pval_psi)
      zscores_psi[psi] <- zvalue_psi
      psi <- psi + 1
    } else if (minpval_psi <= 0.05){
      zvalue_psi <- pval_converter(minpval_psi)
      zscores_psi[psi] <- zvalue_psi
      psi <- psi + 1
    } else {
      extracontrol_psi <- rnorm_multi(
        extra_n, vars = 2, control_mu, sd, r,
        varnames = c("DV1","DV2"))
      extraexp_psi <- rnorm_multi(
        extra_n, vars = 2, exp_mu, sd, r,
        varnames = c("DV1","DV2"))
      
      extracontrol_psi$group <- "Control"
      extraexp_psi$group <- "Experimental"
      
      extra_datapsi <- rbind(extracontrol_psi, extraexp_psi)
      gender <- rbinom(n = extra_n*2, size = 1, p = 0.5)
      gender <- as.factor(
        ifelse(gender == 1, "Female", "Male"))
      extra_datapsi <- extra_datapsi %>% mutate(gender = gender)
      
      combined_datapsi <- rbind(data_psi, extra_datapsi)
      
      extra_pvalpsi <- sitA_ttests(
        combined_datapsi$DV1[combined_datapsi$group=="Control"],
        combined_datapsi$DV2[combined_datapsi$group=="Control"],
        combined_datapsi$DV1[combined_datapsi$group=="Experimental"],
        combined_datapsi$DV2[combined_datapsi$group=="Experimental"]
      )
      
      Exdv1_model <- lm(DV1 ~ group + gender, data = combined_datapsi)
      Exdv1_fitP <- summary(Exdv1_model)$coefficients[
        "groupExperimental","Pr(>|t|)"]
      extra_pvalpsi <- c(extra_pvalpsi, Exdv1_fitP)
      Exdv2_model <- lm(DV2 ~ group + gender, data = combined_datapsi)
      Exdv2_fitP <- summary(Exdv2_model)$coefficients[
        "groupExperimental","Pr(>|t|)"]
      extra_pvalpsi <- c(extra_pvalpsi, Exdv2_fitP)
      
      Exdv1_intmodel <- lm(DV1 ~ group*gender, data = combined_datapsi)
      Exdv1_intfitP <- summary(Exdv1_intmodel)$coefficients[
        "groupExperimental:genderMale","Pr(>|t|)"]
      Exdv2_intmodel <- lm(DV2 ~ group*gender, data = combined_datapsi)
      Exdv2_intfitP <- summary(Exdv2_intmodel)$coefficients[
        "groupExperimental:genderMale","Pr(>|t|)"]
      Extraint_pvaLpsi <- c(Exdv1_intfitP, Exdv2_intfitP)
      
      Exmin_pvalpsi <- min(extra_pvalpsi)
      ExmInt_pvalpsi <- min(Extraint_pvaLpsi)
      
      if (ExmInt_pvalpsi <= 0.05) {
        extrazvalue_psi <- pval_converter(ExmInt_pvalpsi)
        zscores_psi[psi] <- extrazvalue_psi
        psi <- psi+1
      } else if (Exmin_pvalpsi <= 0.05) {
        extrazvalue_psi <- pval_converter(Exmin_pvalpsi)
        zscores_psi[psi] <- extrazvalue_psi
        psi <- psi+1
      }
      
    }}
  
  fit_psi <- zcurve(zscores_psi, control = list(parallel = TRUE))
  plot(fit_psi, CI = TRUE, annotation = TRUE, main = "Scenario Psi")
  psi_plot <- recordPlot()
  
  psi_list <- list(fit_psi = fit_psi,
                    psi_plot = psi_plot)
  
  return(psi_list)
}

psi_500 <- psi_sim(500)
summary(psi_500$fit_psi)

psi_alt <- psi_sim(500, exp_mu = c(0.2, 0.2))
summary(psi_alt$fit_psi)

psi_alt_strong <- psi_sim(500, exp_mu = c(0.8, 0.8))
summary(psi_alt_strong$fit_psi)

#-------------------------------------------------------------------------------
#Situation A,B,C,D combined, fixed number of significant z-scores
zeta_sim <- function(k_sig, n = 20, extra_n = 10, r = 0.5,
                     mu = 0, sd = c(1,1)) {
  zscores_zeta <- numeric(k_sig)
  zeta <- 1
  
  mu_conditions <- list(low = c(-mu, -mu),
                        medium = c(0,0),
                        high = c(mu, mu))
  
  while (zeta <= k_sig) {
    #Generating data frame of Control group DV1, DV2 outcome values for Sit A
    conditions <- sample(rep(c("low", "medium", "high"), length.out = n*2))
    #Column of Sit D "conditions" variables added
    dv_list <- lapply(conditions, zeta_mu_assign,
                 mu_conditions = mu_conditions,
                 sd = sd, r = r)
    dv <- do.call(rbind, dv_list)
    data_zeta <- data.frame(conditions, dv)
    
    #Adding Sit C "gender" variable for each observation
    gender <- rbinom(nrow(data_zeta), size = 1, p = 0.5)
    gender <- as.factor(
      ifelse(gender == 1, "Female", "Male"))
    data_zeta <- data_zeta %>% mutate(gender = gender)
    
    #T-tests for all 3 combinations of conditions (Sit D) for DV1 and DV2 (Sit A)
    pvalues_zeta <- sitZ_ttests(data_zeta)
    
    #Gender main effect ANCOVA for DV1 (Sit C)
    dv1_model <- lm(DV1 ~ conditions + gender, data = data_zeta)
    dv1_fitP <- summary(dv1_model)$coefficients[
      grep("^conditions", rownames(summary(dv1_model)$coefficients)),
      "Pr(>|t|)"]
    pvalues_zeta <- c(pvalues_zeta, dv1_fitP)
    
    #Gender main effect ANCOVA for DV2 (Sit C)
    dv2_model <- lm(DV2 ~ conditions + gender, data = data_zeta)
    dv2_fitP <- summary(dv2_model)$coefficients[
      grep("^conditions", rownames(summary(dv2_model)$coefficients)),
      "Pr(>|t|)"]
    pvalues_zeta <- c(pvalues_zeta, dv2_fitP)
    
    #Gender*group interaction effect ANCOVA for DV1 (Sit C)
    dv1_intmodel <- lm(DV1 ~ conditions*gender, data = data_zeta)
    dv1_intfitP <- summary(
      dv1_intmodel)$coefficients[grep("conditions.*:genderMale",
                                      rownames(summary(
                                        dv1_intmodel)$coefficients)), "Pr(>|t|)"]
    
    #Gender*group interaction effect ANCOVA for DV2 (Sit C)
    dv2_intmodel <- lm(DV2 ~ conditions*gender, data = data_zeta)
    dv2_intfitP <- summary(
      dv2_intmodel)$coefficients[grep("conditions.*:genderMale",
                                      rownames(summary(
                                        dv2_intmodel)$coefficients)), "Pr(>|t|)"]
    
    #Interaction term ANCOVA model p-values
    int_pvalueszeta <- c(dv1_intfitP, dv2_intfitP)
    
    #OLS regression for the linear trend of conditions for DV1 (Sit D)
    lm_coding <- ifelse(data_zeta$conditions == "low", -1,
                        ifelse(data_zeta$conditions == "medium", 0, 1))
    
    OLS_modelDV1 <- lm(DV1~lm_coding, data = data_zeta)
    DV1model_pvalue <- coef(summary(OLS_modelDV1))["lm_coding", "Pr(>|t|)"]
    
    OLS_modelDV2 <- lm(DV2~lm_coding, data = data_zeta)
    DV2model_pvalue <- coef(summary(OLS_modelDV2))["lm_coding", "Pr(>|t|)"]
    
    pvalues_zeta <- c(pvalues_zeta, DV1model_pvalue, DV2model_pvalue)
    
    #Minimum p-values
    min_pvaluezeta <- min(pvalues_zeta)
    mInt_pvaluezeta <- min(int_pvalueszeta)
    
    if (mInt_pvaluezeta <= 0.05) {
      zvalue_zeta <- pval_converter(mInt_pvaluezeta)
      zscores_zeta[zeta] <- zvalue_zeta
      zeta <- zeta + 1
    } else if (min_pvaluezeta <= 0.05){
      #"if the effect of condition was significant in any of these analyses"
      zvalue_zeta <- pval_converter(min_pvaluezeta)
      zscores_zeta[zeta] <- zvalue_zeta
      zeta <- zeta + 1
    } else {
      extra_conditions <- sample(rep(c("low", "medium", "high"), 
                                     length.out = extra_n*2))
      #Column of Sit D "conditions" variables added
      extraDV_list <- lapply(extra_conditions, zeta_mu_assign,
                             mu_conditions = mu_conditions,
                             sd = sd, r = r)
      extra_dv <- do.call(rbind, extraDV_list)
      extradata_zeta <- data.frame(conditions = extra_conditions, extra_dv)
      
      
      #Adding Sit C "gender" variable for each observation
      gender <- rbinom(nrow(extradata_zeta), size = 1, p = 0.5)
      gender <- as.factor(
        ifelse(gender == 1, "Female", "Male"))
      extradata_zeta <- extradata_zeta %>% mutate(gender = gender)
      
      #Combining additional observations to current data set
      combined_datazeta <- rbind(data_zeta, extradata_zeta)
      
      #T-tests for all 3 combos of conditions (Sit D) for DV1 and DV2 (Sit A)
      extra_pvalszeta <- sitZ_ttests(combined_datazeta)
      
      #Gender main effect ANCOVA for DV1 (Sit C)
      dv1_model <- lm(DV1 ~ conditions + gender, data = combined_datazeta)
      dv1_fitP <- summary(dv1_model)$coefficients[
        grep("^conditions", rownames(summary(dv1_model)$coefficients)),
        "Pr(>|t|)"]
      extra_pvalszeta <- c(extra_pvalszeta, dv1_fitP)
      
      #Gender main effect ANCOVA for DV2 (Sit C)
      dv2_model <- lm(DV2 ~ conditions + gender, data = combined_datazeta)
      dv2_fitP <- summary(dv2_model)$coefficients[
        grep("^conditions", rownames(summary(dv2_model)$coefficients)),
        "Pr(>|t|)"]
      extra_pvalszeta <- c(extra_pvalszeta, dv2_fitP)
      
      #Gender*group interaction effect ANCOVA for DV1 (Sit C)
      dv1_intmodel <- lm(DV1 ~ conditions*gender, data = combined_datazeta)
      dv1_intfitP <- summary(
        dv1_intmodel)$coefficients[grep("conditions.*:genderMale",
                                        rownames(summary(
                                          dv1_intmodel)$coefficients)),
                                   "Pr(>|t|)"]
      
      #Gender*group interaction effect ANCOVA for DV2 (Sit C)
      dv2_intmodel <- lm(DV2 ~ conditions*gender, data = combined_datazeta)
      dv2_intfitP <- summary(
        dv2_intmodel)$coefficients[grep("conditions.*:genderMale",
                                        rownames(summary(
                                          dv2_intmodel)$coefficients)),
                                   "Pr(>|t|)"]
      
      #Interaction term ANCOVA model p-values
      extraInt_pvalzeta <- c(dv1_intfitP, dv2_intfitP)
      
      #OLS regression for the linear trend of conditions for DV1 (Sit D)
      lm_coding <- ifelse(combined_datazeta$conditions == "low", -1,
                          ifelse(combined_datazeta$conditions == "medium", 0, 1))
      
      OLS_modelDV1 <- lm(DV1~lm_coding, data = combined_datazeta)
      DV1model_pvalue <- coef(summary(OLS_modelDV1))["lm_coding", "Pr(>|t|)"]
      
      OLS_modelDV2 <- lm(DV2~lm_coding, data = combined_datazeta)
      DV2model_pvalue <- coef(summary(OLS_modelDV2))["lm_coding", "Pr(>|t|)"]
      
      extra_pvalszeta <- c(extra_pvalszeta, DV1model_pvalue, DV2model_pvalue)
      
      #Minimum p-values
      ExMin_pvaluezeta <- min(extra_pvalszeta)
      ExMInt_pvaluezeta <- min(extraInt_pvalzeta)
      
      if (ExMInt_pvaluezeta <= 0.05) {
        zvalue_zeta <- pval_converter(ExMInt_pvaluezeta)
        zscores_zeta[zeta] <- zvalue_zeta
        zeta <- zeta + 1
      } else if (ExMin_pvaluezeta <= 0.05) {
        #"if the effect of condition was significant in any of these analyses"
        zvalue_zeta <- pval_converter(ExMin_pvaluezeta)
        zscores_zeta[zeta] <- zvalue_zeta
        zeta <- zeta + 1
      }
    }
  }
  
  fit_zeta <- zcurve(zscores_zeta, control = list(parallel = TRUE))
  plot(fit_zeta, CI = TRUE, annotation = TRUE, main = "Scenario Zeta")
  zeta_plot <- recordPlot()
  
  zeta_list <- list(fit_zeta = fit_zeta, 
                    zeta_plot = zeta_plot)
  return(zeta_list)
}

zeta_500 <- zeta_sim(500)
summary(zeta_500$fit_zeta)

zeta_alt <- zeta_sim(500, mu = 0.2)
summary(zeta_alt$fit_zeta)

zeta_alt_strong <- zeta_sim(500, mu = 0.8)
summary(zeta_alt_strong$fit_zeta)








