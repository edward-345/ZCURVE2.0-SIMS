library(zcurve)
library(faux)
library(truncnorm)
library(tidyverse)
library(ggplot2)
source("helper_functions.R")
set.seed(666)

#General case of Z-Curve under null hypothesis
z_scores <- numeric(15000)
j <- 1

for (i in 1:15000) {
  control_group <- rnorm(n = 20, mean = 0, sd = 1)
  exp_group <- rnorm(n = 20, mean = 0, sd = 1)
  
  result <- t.test(control_group, exp_group, var.equal = TRUE)
  p_value <- result$p.value
  z_value <- abs(qnorm(p_value/2, lower.tail = FALSE))
  
  if (abs(z_value) >= qnorm(0.975)) {
    z_scores[j] <- z_value
    j <- j + 1
  }
}

z_scores <- z_scores[1:(j - 1)]

fit <- zcurve(z_scores)

sims_plot <- plot(
  fit, CI = TRUE, annotation = TRUE, main = "Simulation under Null")

#Situation A: Two DVs for each observation
zscores_A <- numeric(15000)
A <- 1

pvalues_scenarioA <- numeric(15000)

for (i in 1:15000) {
  control_A <- rnorm_multi(
    n = 20, vars = 2, mu = c(0,0),sd = c(1,1), r = 0.5,
    varnames = c("Control1","Control2"))
  exp_A <- rnorm_multi(
    n = 20, vars = 2, mu = c(0,0), sd = c(1,1), r = 0.5,
    varnames = c("Dependent1","Dependent2"))
  
  pvalue_A <- sitA_ttests(control_A$Control1, control_A$Control2,
              exp_A$Dependent1, exp_A$Dependent2)
  min_pvalue <- min(pvalue_A)
  
  pvalues_scenarioA[i] <- min_pvalue
  
  if (min_pvalue <= 0.05) {
    zvalue_A <- abs(qnorm(min_pvalue/2,
                          lower.tail = FALSE))
    zscores_A[A] <- zvalue_A
    A <- A + 1
  }
}

zscores_A <- zscores_A[1:(A - 1)]

fit_A <- zcurve(zscores_A)

A_plot <- plot(fit_A, CI = TRUE, annotation = TRUE, main = "Scenario A")

#Note that the proportion of p-values align with Simmons et al., 2011
proportions_A <- sig_pvalues(pvalues_scenarioA)

#Situation B: Optional Stopping
zscores_B <- numeric(15000)
B <- 1 

pvalues_scenarioB <- numeric(15000)

for (i in 1:15000) {
  control_B <- rnorm(n = 20, mean = 0, sd = 1)
  exp_B <- rnorm(n = 20, mean = 0, sd = 1)
  result_B <- t.test(control_B, exp_B, var.equal = TRUE)
  pvalue_B <- result_B$p.value
  
  if (pvalue_B <= 0.05) {
    pvalues_scenarioB[i] <- pvalue_B
    zvalue_B <- abs(qnorm(pvalue_B/2, lower.tail = FALSE))
    zscores_B[B] <- zvalue_B
    B <- B+1
  } else {
    extracontrol_B <- rnorm(n = 10, mean = 0, sd = 1)
    extraexp_B <- rnorm(n = 10, mean = 0, sd = 1)
    extraresult_B <- t.test(c(control_B, extracontrol_B),
                            c(exp_B, extraexp_B))
    extrapvalue_B <- extraresult_B$p.value
    
    if (extrapvalue_B <= 0.05) {
      pvalues_scenarioB[i] <- extrapvalue_B
      extrazvalue_B <- abs(qnorm(extrapvalue_B/2, lower.tail = FALSE))
      zscores_B[B] <- extrazvalue_B
      B <- B+1
    } else {
      pvalues_scenarioB[i] <- extrapvalue_B
    }
  }
}

zscores_B <- zscores_B[1:(B - 1)]

fit_B <- zcurve(zscores_B)

B_plot <- plot(fit_B, CI = TRUE, annotation = TRUE, main = "Scenario B")

#Note that the proportion of p-values align with Simmons et al., 2011
proportions_B <- sig_pvalues(pvalues_scenarioB)

#Situation C: Main effect or interaction term ANCOVAs
zscores_C <- numeric(15000)
C <- 1 

pvalues_scenarioC <- numeric(15000)

for (i in 1:15000) {
  groups <- sample(
    rep(c("control", "experimental"), each = 20))
  dv <- rnorm(n = 40, mean = 0, sd = 1)
  gender <- rbinom(n = 40, size = 1, p = 0.5)
  gender <- as.factor(
    ifelse(gender == 1, "female", "male"))
  data_C <- data.frame(dv, gender, groups)
  
  ttest_result <- t.test(dv ~ groups, data = data_C, var.equal = TRUE)
  ttest_pvalue <- ttest_result$p.value
  
  main_model <- lm(dv ~ groups + gender, data = data_C)
  main_pvalue <- summary(
    main_model)$coefficients["groupsexperimental", "Pr(>|t|)"]
  
  int_model <- lm(dv ~ groups*gender, data = data_C)
  coefs <- coef(summary(int_model))
  int_term_name <- grep("group.*:gender.*", rownames(coefs), value = TRUE)
  if (length(int_term_name) == 1) {
    int_pvalue <- coefs[int_term_name, "Pr(>|t|)"]
  } else {
    int_pvalue <- 1
  }
  
  pvalues_C <- c(ttest=ttest_pvalue,
                 main=main_pvalue,
                 int=int_pvalue)
  
  if (pvalues_C["int"] <= 0.05) {
    pvalues_scenarioC[i] <- pvalues_C["int"]
    zvalue_C <- abs(qnorm(pvalues_C["int"]/2, lower.tail = FALSE))
    zscores_C[C] <- zvalue_C
    C <- C + 1
  } else {
    if (min(pvalues_C) <= 0.05) {
      pvalues_scenarioC[i] <- min(pvalues_C)
      zvalue_C <- abs(qnorm(min(pvalues_C)/2, lower.tail = FALSE))
      zscores_C[C] <- zvalue_C
      C <- C + 1
    } else {
      pvalues_scenarioC[i] <- min(pvalues_C)
    }
  }
  
}

zscores_C <- zscores_C[1:(C - 1)]

fit_C <- zcurve(zscores_C)

C_plot <- plot(fit_C, CI = TRUE, annotation = TRUE, main = "Scenario C")

#Note that the proportion of p-values align with Simmons et al., 2011
proportions_C <- sig_pvalues(pvalues_scenarioC)

#Situation D: Ordinal test condition
zscores_D <- numeric(15000)
D <- 1 

pvalues_scenarioD <- numeric(15000)

for (i in 1:15000) {
  conditions <- sample(
    rep(c("low", "medium", "high"), length.out = 40))
  dv <- rnorm(n = 40, mean = 0, sd = 1)
  data_D <- data.frame(conditions, dv)
  
  LowMed_pvalue <- t.test(
    dv ~ conditions,
    data = subset(data_D, conditions %in% c("low","medium")))$p.value
  LowHigh_pvalue <-t.test(
    dv ~ conditions,
    data = subset(data_D, conditions %in% c("low","high")))$p.value
  MedHigh_pvalue <- t.test(
    dv ~ conditions,
    data = subset(data_D, conditions %in% c("medium","high")))$p.value
  
  lm_coding <- ifelse(data_D$conditions == "low", -1,
                      ifelse(data_D$conditions == "medium", 0, 1))
  model_D <- lm(dv~lm_coding)
  model_pvalue <- coef(summary(model_D))["lm_coding", "Pr(>|t|)"]
  
  pvalues_D <- c(LowMed_pvalue, LowHigh_pvalue, MedHigh_pvalue,
                 model_pvalue)
  
  if (min(pvalues_D) <= 0.05) {
    pvalues_scenarioD[i] <- min(pvalues_D)
    zvalue_D <- abs(qnorm(min(pvalues_D)/2, lower.tail = FALSE))
    zscores_D[D] <- zvalue_D
    D <- D + 1
  } else {
    pvalues_scenarioD[i] <- min(pvalues_D)
  }
  
}

zscores_D <- zscores_D[1:(D - 1)]

fit_D <- zcurve(zscores_D)

D_plot <- plot(fit_D, CI = TRUE, annotation = TRUE, main = "Scenario D")

#Note that the proportion of p-values align with Simmons et al., 2011
proportions_D <- sig_pvalues(pvalues_scenarioD)













