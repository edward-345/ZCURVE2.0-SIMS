#Run 3 t tests
sitA_ttests <- function(ctrl1, ctrl2, exp1, exp2) {
  p1 <- t.test(ctrl1, exp1, var.equal = TRUE)$p.value
  p2 <- t.test(ctrl2, exp2, var.equal = TRUE)$p.value
  avg_ctrl <- rowMeans(cbind(ctrl1, ctrl2))
  avg_exp  <- rowMeans(cbind(exp1, exp2))
  p3 <- t.test(avg_ctrl, avg_exp, var.equal = TRUE)$p.value
  return(c(p1, p2, p3))
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

#Convert p-values to z-scores
pvect_zvect <- function(pvals) {
  if (any(pvals < 0 | pvals > 1, na.rm = TRUE)) {
    stop("All p-values must be between 0 and 1.")
  }
  qnorm(1 - pvals/2)
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







