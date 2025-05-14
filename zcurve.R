library(zcurve)
library(faux)
library(truncnorm)
set.seed(666)

#example use of zcurve
#fit <- zcurve(OSC.z)

#zcurve_plot <- plot(
#  fit,CI=TRUE,annotation=TRUE,main="OSC 2015")

#N(0,1)
#n=15,000 simulates studies
#each has n=40, half in control and half in exp group

#General case of Z-Curve
z_scores <- numeric(15000)
j <- 1

for (i in 1:15000) {
  control_group <- rnorm(n=20,mean=0,sd=1)
  exp_group <- rnorm(n=20,mean=0,sd=1)
  result <- t.test(control_group, exp_group,
                   var.equal=TRUE)
  p_value <- result$p.value
  z_value <- abs(qnorm(p_value/2,
                       lower.tail=FALSE))
  if (abs(z_value) >= qnorm(0.975)) {
    z_scores[j] <- z_value
    j <- j+1
  }
}

z_scores <- z_scores[1:(j - 1)]

fit <- zcurve(z_scores)

sims_plot <- plot(fit,CI=TRUE,annotation=TRUE,main="Simulations")

#Two DVs for each observation
zscores_A <- numeric(15000)
A <- 1

pvalues_scenarioA <- numeric(15000)

for (i in 1:15000) {
  control_A <- rnorm_multi(
    n=20,vars=2,mu=c(0,0),sd=c(1,1),
    r=0.5,varnames=c("Control1","Control2"))
  exp_A <- rnorm_multi(
    n=20,vars=2,mu=c(0,0),sd=c(1,1),
    r=0.5,varnames=c("Dependent1","Dependent2"))
  result1 <- t.test(control_A$Control1, exp_A$Dependent1,
                   var.equal=TRUE)
  pvalue1 <- result1$p.value
  result2 <- t.test(control_A$Control2, exp_A$Dependent2,
                    var.equal=TRUE)
  pvalue2 <- result2$p.value
  result3 <- t.test(rowMeans(cbind(control_A$Control1,control_A$Control2)),
                    rowMeans(cbind(exp_A$Dependent1,exp_A$Dependent2)),
                    var.equal = TRUE)
  pvalue3 <- result3$p.value
  pvalues_scenarioA[i] <- min(
    c(result1$p.value,result2$p.value,result3$p.value))
  pvalues_A <- c(pvalue1, pvalue2, pvalue3)
  zvalue_A <- abs(qnorm(min(pvalues_A)/2,
                       lower.tail=FALSE))
  if (abs(zvalue_A) >= qnorm(0.975)) {
    zscores_A[A] <- zvalue_A
    A <- A+1
  }
}

zscores_A <- zscores_A[1:(A - 1)]

fit_A <- zcurve(zscores_A)

A_plot <- plot(fit_A,CI=TRUE,annotation=TRUE,main="Scenario A")

#Note that the proportion of p-values align with simmons et al., 2011
pvalues1_A <- pvalues_scenarioA[pvalues_scenarioA < 0.1]
prop1_A <- (length(pvalues1_A)/15000)*100

pvalues05_A <- pvalues_scenarioA[pvalues_scenarioA < 0.05]
prop05_A <- (length(pvalues05_A)/15000)*100

pvalues01_A <- pvalues_scenarioA[pvalues_scenarioA < 0.01]
prop01_A <- (length(pvalues01_A)/15000)*100

#Optional Stopping
zscores_B <- numeric(15000)
B <- 1 

pvalues_scenarioB <- numeric(15000)

for (i in 1:15000) {
  control_B <- rnorm(n=10,mean=0,sd=1)
  exp_B <- rnorm(n=10,mean=0,sd=1)
  result_B <- t.test(control_B, exp_B,
                   var.equal=TRUE)
  pvalue_B <- result_B$p.value
  if (pvalue_B <= 0.05) {
    pvalues_scenarioB[i] <- pvalue_B
    zvalue_B <- abs(qnorm(pvalue_B/2,lower.tail=FALSE))
    zscores_B[B] <- zvalue_B
    B <- B+1
  } else {
    extracontrol_B <- rnorm(n=10,mean=0,sd=1)
    extraexp_B <- rnorm(n=10,mean=0,sd=1)
    extraresult_B <- t.test(c(control_B, extracontrol_B),
                            c(exp_B, extraexp_B))
    extrapvalue_B <- extraresult_B$p.value
    if (extrapvalue_B <= 0.05) {
      pvalues_scenarioB[i] <- extrapvalue_B
      extrazvalue_B <- abs(qnorm(extrapvalue_B/2,lower.tail=FALSE))
      zscores_B[B] <- extrazvalue_B
      B <- B+1
    } else {
      pvalues_scenarioB[i] <- extrapvalue_B
    }
  }
}

zscores_B <- zscores_B[1:(B - 1)]

fit_B <- zcurve(zscores_B)

B_plot <- plot(fit_B,CI=TRUE,annotation=TRUE,main="Scenario B")

#Note that the proportion of p-values align with simmons et al., 2011
pvalues1_B <- pvalues_scenarioB[pvalues_scenarioB < 0.1]
prop1_B <- (length(pvalues1_B)/15000)*100

pvalues05_B <- pvalues_scenarioB[pvalues_scenarioB < 0.05]
prop05_B <- (length(pvalues05_B)/15000)*100

pvalues01_B <- pvalues_scenarioB[pvalues_scenarioB < 0.01]
prop01_B <- (length(pvalues01_B)/15000)*100
