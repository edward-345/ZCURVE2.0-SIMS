source("pHacking_Null.R")

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
zscores_alpha <- numeric(15000)
alpha <- 1

pvalues_scenarioAlpha <- numeric(15000)

for (i in 1:15000) {
  control_alpha <- rnorm_multi(
    n = 20, vars = 2, mu = c(0,0),sd = c(1,1), r = 0.5,
    varnames = c("Control1","Control2"))
  exp_alpha <- rnorm_multi(
    n = 20, vars = 2, mu = c(0.5,0.5), sd = c(1,1), r = 0.5,
    varnames = c("Dependent1","Dependent2"))
  
  result1 <- t.test(control_alpha$Control1, exp_alpha$Dependent1,
                    var.equal = TRUE)
  pvalue1 <- result1$p.value
  result2 <- t.test(control_alpha$Control2, exp_alpha$Dependent2,
                    var.equal = TRUE)
  pvalue2 <- result2$p.value
  result3 <- t.test(rowMeans(cbind(control_alpha$Control1, control_alpha$Control2)),
                    rowMeans(cbind(exp_alpha$Dependent1, exp_alpha$Dependent2)),
                    var.equal = TRUE)
  pvalue3 <- result3$p.value
  
  pvalues_scenarioAlpha[i] <- min(
    c(result1$p.value, result2$p.value, result3$p.value))
  pvalues_alpha <- c(pvalue1, pvalue2, pvalue3)
  zvalue_alpha <- abs(qnorm(min(pvalues_alpha)/2,
                        lower.tail = FALSE))
  
  if (abs(zvalue_alpha) >= qnorm(0.975)) {
    zscores_alpha[alpha] <- zvalue_alpha
    alpha <- alpha + 1
  }
}

zscores_alpha <- zscores_alpha[1:(alpha - 1)]

fit_alpha <- zcurve(zscores_alpha)

alpha_plot <- plot(
  fit_alpha, CI = TRUE, annotation = TRUE, main = "Scenario A under Alternate")


