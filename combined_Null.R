source("pHacking_Null.R")
library(tidyverse)

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

X_plot <- plot(fit_D, CI = TRUE, annotation = TRUE, main = "Scenario A+B")

#Note that the proportion of p-values align with Simmons et al., 2011
proportions_X <- sig_pvalues(pvalues_scenarioX)

