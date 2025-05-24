#Run 3 t tests
sitA_ttests <- function(ctrl1, ctrl2, exp1, exp2) {
  p1 <- t.test(ctrl1, exp1, var.equal = TRUE)$p.value
  p2 <- t.test(ctrl2, exp2, var.equal = TRUE)$p.value
  avg_ctrl <- rowMeans(cbind(ctrl1, ctrl2))
  avg_exp  <- rowMeans(cbind(exp1, exp2))
  p3 <- t.test(avg_ctrl, avg_exp, var.equal = TRUE)$p.value
  return(c(p1, p2, p3))
}
