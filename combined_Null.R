source("pHacking_Null.R")
library(tidyverse)

#SITUATIONS A,B COMBINED
zscores_X <- numeric(15000)
X <- 1 

pvalues_scenarioX <- numeric(15000)

for (i in 1:15000) {
  control_X <- rnorm_multi(
    n = 20, vars = 2, mu = c(0,0),sd = c(1,1), r = 0.5,
    varnames = c("Control1","Control2"))
  exp_X <- rnorm_multi(
    n = 20, vars = 2, mu = c(0,0), sd = c(1,1), r = 0.5,
    varnames = c("Dependent1","Dependent2"))
  
  pvalues_X <- sitA_ttests(control_X$Control1, control_X$Control2,
                           exp_X$Dependent1, exp_X$Dependent2)
  min_pvalueX <- min(pvalues_X)
  
  if (min_pvalueX <= 0.05) {
    pvalues_scenarioX[i] <- min_pvalueX
    zvalue_X <- abs(qnorm(min_pvalueX/2,
                          lower.tail = FALSE))
    zscores_X[X] <- zvalue_X
    X <- X + 1
  } else {
    extracontrol_X <- rnorm_multi(
      n = 10, vars = 2, mu = c(0,0),sd = c(1,1), r = 0.5,
      varnames = c("Control1","Control2"))
    extraexp_X <- rnorm_multi(
      n = 10, vars = 2, mu = c(0,0),sd = c(1,1), r = 0.5,
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
      extrazvalue_X <- abs(qnorm(min_ExtraPvalue/2, lower.tail = FALSE))
      zscores_X[X] <- extrazvalue_X
      X <- X+1
    } else {
      pvalues_scenarioX[i] <- min_ExtraPvalue
    }
  }
}


zscores_X <- zscores_X[1:(X - 1)]

fit_X <- zcurve(zscores_X)

X_plot <- plot(fit_X, CI = TRUE, annotation = TRUE, main = "Scenario A+B")

#Note that the proportion of p-values align with Simmons et al., 2011
proportions_X <- sig_pvalues(pvalues_scenarioX)

#SITUATIONS A,B,C COMBINED
zscores_Y <- numeric(15000)
Y <- 1 

pvalues_scenarioY <- numeric(15000)

for (i in 1:15000) {
  control_Y <- rnorm_multi(
    n = 20, vars = 2, mu = c(0,0),sd = c(1,1), r = 0.5,
    varnames = c("DV1","DV2"))
  exp_Y <- rnorm_multi(
    n = 20, vars = 2, mu = c(0,0), sd = c(1,1), r = 0.5,
    varnames = c("DV1","DV2"))
  
  control_Y$group <- "Control"
  exp_Y$group <- "Experimental"
  
  data_Y <- rbind(control_Y, exp_Y)
  
  gender <- rbinom(n = 40, size = 1, p = 0.5)
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
    zvalue_Y <- abs(qnorm(mInt_pvalueY/2,
                          lower.tail = FALSE))
    zscores_Y[Y] <- zvalue_Y
    Y <- Y + 1
  } else if (min_pvalueY <= 0.05){
    pvalues_scenarioY[i] <- min_pvalueY
    zvalue_Y <- abs(qnorm(min_pvalueY/2,
                          lower.tail = FALSE))
    zscores_Y[Y] <- zvalue_Y
    Y <- Y + 1
    } else {
    extracontrol_Y <- rnorm_multi(
      n = 10, vars = 2, mu = c(0,0),sd = c(1,1), r = 0.5,
      varnames = c("DV1","DV2"))
    extraexp_Y <- rnorm_multi(
      n = 10, vars = 2, mu = c(0,0),sd = c(1,1), r = 0.5,
      varnames = c("DV1","DV2"))
    
    extracontrol_Y$group <- "Control"
    extraexp_Y$group <- "Experimental"
    
    extra_dataY <- rbind(extracontrol_Y, extraexp_Y)
    gender <- rbinom(n = 20, size = 1, p = 0.5)
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
      extrazvalue_Y <- abs(qnorm(ExmInt_pvalY/2, lower.tail = FALSE))
      zscores_Y[Y] <- extrazvalue_Y
      Y <- Y+1
    } else if (Exmin_pvalY <= 0.05) {
      pvalues_scenarioY[i] <- Exmin_pvalY
      extrazvalue_Y <- abs(qnorm(Exmin_pvalY/2, lower.tail = FALSE))
      zscores_Y[Y] <- extrazvalue_Y
      Y <- Y+1
    } else {
      all_pvalY <- c(extra_pvalY, Extraint_pvaLY)
      pvalues_scenarioY[i] <- min(all_pvalY)
    }
  }
}

zscores_Y <- zscores_Y[1:(Y - 1)]

fit_Y <- zcurve(zscores_Y)

Y_plot <- plot(fit_Y, CI = TRUE, annotation = TRUE, main = "Scenario A+B+C")

#Note that the proportion of p-values align with Simmons et al., 2011
proportions_Y <- sig_pvalues(pvalues_scenarioY)












