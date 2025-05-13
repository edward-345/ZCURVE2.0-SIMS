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






