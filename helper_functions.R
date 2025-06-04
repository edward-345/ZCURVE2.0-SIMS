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
