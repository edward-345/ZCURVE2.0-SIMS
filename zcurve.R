library(zcurve)
library(faux)
library(truncnorm)
set.seed(666)

#example use of zcurve
fit <- zcurve(OSC.z)

zcurve_plot <- plot(
  fit,CI=TRUE,annotation=TRUE,main="OSC 2015")

#N(0,1)
#n=15,000 simulates studies
#each has n=40, half in control and half in exp group


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
  if (abs(z_value) > qnorm(0.975)) {
    z_scores[j] <- z_value
    j <- j+1
  }
}

z_scores <- z_scores[1:(j - 1)]





