library(zcurve)
library(faux)
library(truncnorm)
set.seed(666)

fit <- zcurve(OSC.z)

zcurve_plot <- plot(
  fit,CI=TRUE,annotation=TRUE,main="OSC 2015")

#N(0,1)
#n=15,000 simulates studies
#each has n=40, half in control and half in exp group

n <- 15000
z_scores <- numeric(n)


for (i in 1:15000) {
  control_group <- rnorm(n=20,mean=0,sd=1)
  exp_group <- rnorm(n=20,mean=0,sd=1)
  result <- t.test(control_group, exp_group)
  p_value <- result$p.value
  z_scores[i] <- p_value
}




