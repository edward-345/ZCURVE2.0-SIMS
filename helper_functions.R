library(zcurve)
library(faux)
library(truncnorm)
library(tidyverse)
library(ggplot2)
set.seed(666)

#ZCURVE3.0 Imports
zcurve3 <- "https://raw.githubusercontent.com/UlrichSchimmack/zcurve3.0/refs/heads/main/Zing.25.07.11.test.R"
source(zcurve3)

#----------------------------------------------
#Run 3 t tests
sitA_ttests <- function(ctrl1, ctrl2, exp1, exp2) {
  p1 <- t.test(ctrl1, exp1, var.equal = TRUE)$p.value
  p2 <- t.test(ctrl2, exp2, var.equal = TRUE)$p.value
  avg_ctrl <- rowMeans(cbind(ctrl1, ctrl2))
  avg_exp  <- rowMeans(cbind(exp1, exp2))
  p3 <- t.test(avg_ctrl, avg_exp, var.equal = TRUE)$p.value
  return(c(p1, p2, p3))
}

#T-Tests for multiple covariates
multivar_ttests <- function(control, exp) {
  
  pvals <- numeric(ncol(control))
  
  for (i in 1:ncol(control)) {
    pvals[i] <- t.test(control[,i], exp[,i], var.equal = TRUE)$p.value
  }
  avg_pval <- t.test(rowMeans(control), rowMeans(exp), var.equal = TRUE)$p.value
  
  pvals <- c(pvals, avg_pval)
  
  return(pvals)
}


#Distribution of significant p-values
sig_pvalues <- function(sig_p) {
  pvalues1 <- sig_p[sig_p <= 0.1]
  prop1 <- (length(pvalues1)/length(sig_p))*100
  pvalues05 <- sig_p[sig_p <= 0.05]
  prop05 <- (length(pvalues05)/length(sig_p))*100
  pvalues01 <- sig_p[sig_p <= 0.01]
  prop01 <- (length(pvalues01)/length(sig_p))*100
  return(c(prop1, prop05, prop01))
}

sitZ_ttests <- function(dataset) {
  LowMed_DV1 <- dataset %>%
    filter(conditions %in% c("low", "medium")) %>%
    t.test(DV1 ~ conditions, data = .) %>%
    .$p.value
  LowHigh_DV1 <- dataset %>%
    filter(conditions %in% c("low", "high")) %>%
    t.test(DV1 ~ conditions, data = .) %>%
    .$p.value
  MedHigh_DV1 <- dataset %>%
    filter(conditions %in% c("medium", "high")) %>%
    t.test(DV1 ~ conditions, data = .) %>%
    .$p.value
  
  LowMed_DV2 <- dataset %>%
    filter(conditions %in% c("low", "medium")) %>%
    t.test(DV2 ~ conditions, data = .) %>%
    .$p.value
  LowHigh_DV2 <- dataset %>%
    filter(conditions %in% c("low", "high")) %>%
    t.test(DV2 ~ conditions, data = .) %>%
    .$p.value
  MedHigh_DV2 <- dataset %>%
    filter(conditions %in% c("medium", "high")) %>%
    t.test(DV2 ~ conditions, data = .) %>%
    .$p.value
  
  return(c(LowMed_DV1, LowHigh_DV1, MedHigh_DV1,
           LowMed_DV2, LowHigh_DV2, MedHigh_DV2))
}

#Convert p-value to z-score
pval_converter <- function(p_val) {
  return(abs(qnorm(p_val/2, lower.tail = FALSE)))
}

#Situation C t-test and ANCOVA model
sitC_tests <- function(dataset) {
  #Results for Situation C were obtained by conducting a t-test...
  ttest_result <- t.test(dv ~ groups, data = dataset, var.equal = TRUE)
  ttest_pvalue <- ttest_result$p.value
  
  #...an analysis of covariance with a gender main effect..
  main_model <- lm(dv ~ groups + gender, data = dataset)
  main_pvalue <- summary(
    main_model)$coefficients["groupsexperimental", "Pr(>|t|)"]
  
  #...and an analysis of covariance with a gender interaction.
  int_model <- lm(dv ~ groups*gender, data = dataset)
  coefs <- coef(summary(int_model))
  int_term_name <- grep("group.*:gender.*", rownames(coefs), value = TRUE)
  
  #Edge case guard of sample consisting entirely of one gender for int_pvalue
  if (length(int_term_name) == 1) {
    int_pvalue <- coefs[int_term_name, "Pr(>|t|)"]
  } else {
    int_pvalue <- 1
  }
  
  pvalues_C <- c(ttest_pvalue, main_pvalue)
  return(list(pvalues_C = pvalues_C,
              int_pvalue = int_pvalue))
}

zeta_mu_assign <- function(g, mu_conditions, sd, r) {
  rnorm_multi(1, vars = 2,
              mu = mu_conditions[[g]], sd, r,
              varnames=c("DV1","DV2"))
}

# True EDR and ERR
true_parameters <- function(n, alpha = .05,
                            pop.es=numeric(5), wgt=c(0.2,0.2,0.2,0.2,0.2)) {
  N <- 2*n
  se <- 2/sqrt(N)
  
  nct = pop.es/se
  p = 1-(alpha/2)
  
  pow.pos <- pt(qt(p,N-2),N-2,nct,lower.tail=FALSE)
  pow.neg <- pt(-qt(p,N-2),N-2,nct,lower.tail=TRUE)
  pow <- pow.pos + pow.neg
  
  true.edr <- sum(pow*wgt)
  
  w.sig <- wgt*pow
  w.sig <- w.sig/sum(w.sig)
  
  true.err <- sum(pow*w.sig)
  
  true_parameters <- c(true.edr, true.err)
  
  return(true_parameters)
}

true_edr <- function(n, alpha = .05,
                     pop.es=numeric(5), wgt=c(0.2,0.2,0.2,0.2,0.2)) {
  N <- 2*n
  se <- 2/sqrt(N)
  
  nct = pop.es/se
  p = 1-(alpha/2)
  
  pow.pos <- pt(qt(p,N-2),N-2,nct,lower.tail=FALSE)
  pow.neg <- pt(-qt(p,N-2),N-2,nct,lower.tail=TRUE)
  pow <- pow.pos + pow.neg
  
  true.edr <- sum(pow*wgt)

  return(true.edr)
}

true_err <- function(n, alpha = .05,
                     pop.es=numeric(5), wgt=c(1,0,0,0,0)) {
  N <- 2*n
  se <- 2/sqrt(N)
  
  nct = pop.es/se
  p = 1-(alpha/2)
  
  pow.pos <- pt(qt(p,N-2),N-2,nct,lower.tail=FALSE)
  pow.neg <- pt(-qt(p,N-2),N-2,nct,lower.tail=TRUE)
  pow <- pow.pos + pow.neg
  
  w.sig <- wgt*pow
  w.sig <- w.sig/sum(w.sig)
  
  true.err <- sum(pow*w.sig)
  
  return(true.err)
}

# ZCURVE 3 EDR extract
z3_EDR <- function(models) {
  EDR_vals <- numeric(length(models))
  i <- 1
  for (i in 1:length(models)) {
    EDR_vals[i] <- models[[i]]$res['EDR']
  }
  return(EDR_vals)
}

# ZCURVE 3 ERR extract
z3_ERR <- function(models) {
  ERR_vals <- numeric(length(models))
  i <- 1
  for (i in 1:length(models)) {
    ERR_vals[i] <- models[[i]]$res['ERR']
  }
  return(ERR_vals)
}


