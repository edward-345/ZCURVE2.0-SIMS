#Situations A and B combined
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
  
  result1 <- t.test(control_X$Control1, exp_X$Dependent1,
                    var.equal = TRUE)
  pvalue1 <- result1$p.value
  result2 <- t.test(control_X$Control2, exp_X$Dependent2,
                    var.equal = TRUE)
  pvalue2 <- result2$p.value
  result3 <- t.test(rowMeans(cbind(control_X$Control1, control_X$Control2)),
                    rowMeans(cbind(exp_X$Dependent1, exp_X$Dependent2)),
                    var.equal = TRUE)
  pvalue3 <- result3$p.value
  
  pvalues_X <- c(pvalue1, pvalue2, pvalue3)
  
  if (min(pvalues_X) <= 0.05) {
    pvalues_scenarioX[i] <- min(pvalues_X)
    zvalue_X <- abs(qnorm(min(pvalues_X)/2,
                          lower.tail = FALSE))
    zscores_X[X] <- zvalue_X
    X <- X + 1
  } else {
    extracontrol_X <- rnorm_multi(
      n = 10, vars = 2, mu = c(0,0),sd = c(1,1), r = 0.5,
      varnames = c("ExCon1","ExCon2"))
    extraexp_X <- rnorm_multi(
      n = 10, vars = 2, mu = c(0,0),sd = c(1,1), r = 0.5,
      varnames = c("ExDep1","ExDep2"))
    
    combined_Control1 <- c(control_X$Control1, extracontrol_X$ExCon1)
    combined_Dep1 <- c(exp_X$Dependent1, extraexp_X$ExDep1)
    combined_Control2 <- c(control_X$Control2, extracontrol_X$ExCon2)
    combined_Dep2 <- c(exp_X$Dependent2, extraexp_X$ExDep2) 
    
    combined_control_avg <- rowMeans(
      cbind(combined_Control1, combined_Control2))
    combined_exp_avg     <- rowMeans(cbind(combined_Dep1, combined_Dep2))
    
    ExtraResult1 <- t.test(combined_Control1,
                           combined_Dep1, var.equal = TRUE)
    ExtraPvalue1 <- ExtraResult1$p.value
    
    ExtraResult2 <- t.test(combined_Control2,
                           combined_Dep2, var.equal = TRUE)
    ExtraPvalue2 <- ExtraResult2$p.value
    
    ExtraResult3 <- t.test(combined_control_avg,
                           combined_exp_avg, var.equal = TRUE)
    ExtraPvalue3 <- ExtraResult3$p.value
    
    ExtraPvalues_X <- c(ExtraPvalue1, ExtraPvalue2, ExtraPvalue3)
    
    if (min(ExtraPvalues_X) <= 0.05) {
      pvalues_scenarioX[i] <- min(ExtraPvalues_X)
      extrazvalue_X <- abs(qnorm(min(ExtraPvalues_X)/2, lower.tail = FALSE))
      zscores_X[X] <- extrazvalue_X
      X <- X+1
    } else {
      pvalues_scenarioX[i] <- min(ExtraPvalues_X)
    }
  }
}


zscores_X <- zscores_X[1:(X - 1)]

fit_X <- zcurve(zscores_X)

X_plot <- plot(fit_D, CI = TRUE, annotation = TRUE, main = "Scenario A+B")

#Note that the proportion of p-values align with Simmons et al., 2011
proportions_X <- sig_pvalues(pvalues_scenarioX)

