library(zcurve)
library(faux)
library(truncnorm)
library(tidyverse)
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

#General case of Z-Curve under alternate hypothesis
zscores_alt <- numeric(15000)
L <- 1

for (i in 1:15000) {
  control_alt <- rnorm(n = 20, mean = 0, sd = 1)
  exp_alt <- rnorm(n = 20, mean = 0.5, sd = 1)
  
  result_alt <- t.test(control_alt, exp_alt, var.equal = TRUE)
  pvalue_alt <- result_alt$p.value
  zvalue_alt <- abs(qnorm(pvalue_alt/2, lower.tail = FALSE))
  
  if (abs(zvalue_alt) >= qnorm(0.975)) {
    zscores_alt[L] <- zvalue_alt
    L <- L + 1
  }
}

zscores_alt <- zscores_alt[1:(L - 1)]

fit_alt <- zcurve(zscores_alt)

alt_plot <- plot(
  fit_alt, CI = TRUE, annotation = TRUE, main = "Simulation under Alternate")

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
  
  result1 <- t.test(control_A$Control1, exp_A$Dependent1,
                    var.equal = TRUE)
  pvalue1 <- result1$p.value
  result2 <- t.test(control_A$Control2, exp_A$Dependent2,
                    var.equal = TRUE)
  pvalue2 <- result2$p.value
  result3 <- t.test(rowMeans(cbind(control_A$Control1, control_A$Control2)),
                    rowMeans(cbind(exp_A$Dependent1, exp_A$Dependent2)),
                    var.equal = TRUE)
  pvalue3 <- result3$p.value
  
  pvalues_scenarioA[i] <- min(
    c(result1$p.value, result2$p.value, result3$p.value))
  pvalues_A <- c(pvalue1, pvalue2, pvalue3)
  zvalue_A <- abs(qnorm(min(pvalues_A)/2,
                       lower.tail = FALSE))
  
  if (abs(zvalue_A) >= qnorm(0.975)) {
    zscores_A[A] <- zvalue_A
    A <- A + 1
  }
}

zscores_A <- zscores_A[1:(A - 1)]

fit_A <- zcurve(zscores_A)

A_plot <- plot(fit_A, CI = TRUE, annotation = TRUE, main = "Scenario A")

#Note that the proportion of p-values align with Simmons et al., 2011
pvalues1_A <- pvalues_scenarioA[pvalues_scenarioA < 0.1]
prop1_A <- (length(pvalues1_A)/15000)*100

pvalues05_A <- pvalues_scenarioA[pvalues_scenarioA < 0.05]
prop05_A <- (length(pvalues05_A)/15000)*100

pvalues01_A <- pvalues_scenarioA[pvalues_scenarioA < 0.01]
prop01_A <- (length(pvalues01_A)/15000)*100

proportions_A <- c(prop1_A, prop05_A, prop01_A)

#Situation B: Optional Stopping
zscores_B <- numeric(15000)
B <- 1 

pvalues_scenarioB <- numeric(15000)

for (i in 1:15000) {
  control_B <- rnorm(n = 10, mean = 0, sd = 1)
  exp_B <- rnorm(n = 10, mean = 0, sd = 1)
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
pvalues1_B <- pvalues_scenarioB[pvalues_scenarioB < 0.1]
prop1_B <- (length(pvalues1_B)/15000)*100

pvalues05_B <- pvalues_scenarioB[pvalues_scenarioB < 0.05]
prop05_B <- (length(pvalues05_B)/15000)*100

pvalues01_B <- pvalues_scenarioB[pvalues_scenarioB < 0.01]
prop01_B <- (length(pvalues01_B)/15000)*100

proportions_B <- c(prop1_B, prop05_B, prop01_B)

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
pvalues1_C <- pvalues_scenarioC[pvalues_scenarioC < 0.1]
prop1_C <- (length(pvalues1_C)/15000)*100

pvalues05_C <- pvalues_scenarioC[pvalues_scenarioC < 0.05]
prop05_C <- (length(pvalues05_C)/15000)*100

pvalues01_C <- pvalues_scenarioC[pvalues_scenarioC < 0.01]
prop01_C <- (length(pvalues01_C)/15000)*100

proportions_C <- c(prop1_C, prop05_C, prop01_C)

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
pvalues1_D <- pvalues_scenarioD[pvalues_scenarioD < 0.1]
prop1_D <- (length(pvalues1_D)/15000)*100

pvalues05_D <- pvalues_scenarioD[pvalues_scenarioD < 0.05]
prop05_D <- (length(pvalues05_D)/15000)*100

pvalues01_D <- pvalues_scenarioD[pvalues_scenarioD < 0.01]
prop01_D <- (length(pvalues01_D)/15000)*100

proportions_D <- c(prop1_D, prop05_D, prop01_D)








