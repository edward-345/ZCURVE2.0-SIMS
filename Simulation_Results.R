# SIMULATION RESULTS
source("pHacking_Null.R")
source("combined_Null.R")
source("fixed_sigZ.R")
source("combined_fixed.R")

#-------------------------------------------------------------------------------
# 500 studies under null hypothesis, both groups come from N(0,1)
alpha_500 <- alpha_sim(500)
summary(alpha_500$fit_alpha)
alpha_500.plot <- plot(alpha_500$fit_alpha,
                       CI = TRUE, annotation = TRUE, main = "Scenario Alpha")

alpha_500.pvalModel <- zcurve(p = alpha_500$pval_list,
                              control = list(parallel = TRUE))
alpha_500.pvalPlot <- plot(alpha_500.pvalModel,
                           CI = TRUE, annotation = TRUE, main = "Scenario Alpha")
alpha_500.pvals <- hist(alpha_500$pval_list)

# 500 studies under null hypothesis, both groups come from N(0,1)
beta_500 <- beta_sim(500)
summary(beta_500$fit_beta)
beta_500.plot <- plot(beta_500$fit_beta,
                      CI = TRUE, annotation = TRUE, main = "Scenario Beta")

beta_500.pvalModel <- zcurve(p = beta_500$pval_list,
                             control = list(parallel = TRUE))
beta_500.pvalPlot <- plot(beta_500.pvalModel,
                          CI = TRUE, annotation = TRUE, main = "Scenario Beta")
beta_500.pvals <- hist(beta_500$pval_list)

#Null hypothesis
gamma_500 <- gamma_sim(500)
summary(gamma_500$fit_gamma)
gamma_500.plot <- plot(gamma_500$fit_gamma,
                       CI = TRUE, annotation = TRUE, main = "Scenario Gamma")

gamma_500.pvalModel <- zcurve(p = gamma_500$pval_list,
                              control = list(parallel = TRUE))
gamma_500.pvalPlot <- plot(gamma_500.pvalModel,
                           CI = TRUE, annotation = TRUE, main = "Scenario Gamma")
gamma_500.pvals <- hist(gamma_500$pval_list)

delta_500 <- delta_sim(500)
summary(delta_500$fit_delta)
delta_500.plot <- plot(delta_500$fit_delta,
                       CI = TRUE, annotation = TRUE, main = "Scenario Delta")

delta_500.pvalModel <- zcurve(p = delta_500$pval_list,
                              control = list(parallel = TRUE))
delta_500.pvalPlot <- plot(delta_500.pvalModel,
                           CI = TRUE, annotation = TRUE, main = "Scenario Delta")
delta_500.pvals <- hist(delta_500$pval_list)

#-------------------------------------------------------------------------------

chi_500 <- chi_sim(500)
summary(chi_500$fit_chi)
chi_500.plot <- plot(chi_500$fit_chi,
                     CI = TRUE, annotation = TRUE, main = "Scenario Chi")

chi_500.pvalModel <- zcurve(p = chi_500$pval_list,
                            control = list(parallel = TRUE))
chi_500.pvalPlot <- plot(chi_500.pvalModel,
                         CI = TRUE, annotation = TRUE, main = "Scenario Chi")
chi_500.pvals <- hist(chi_500$pval_list)

psi_500 <- psi_sim(500)
summary(psi_500$fit_psi)
psi_500.plot <- plot(psi_500$fit_psi,
                     CI = TRUE, annotation = TRUE, main = "Scenario Psi")

psi_500.pvalModel <- zcurve(p = psi_500$pval_list,
                            control = list(parallel = TRUE))
psi_500.pvalPlot <- plot(psi_500.pvalModel,
                         CI = TRUE, annotation = TRUE, main = "Scenario Psi")
psi_500.pvals <- hist(psi_500$pval_list)

zeta_500 <- zeta_sim(500)
summary(zeta_500$fit_zeta)
zeta_500.plot <- plot(zeta_500$fit_zeta,
                      CI = TRUE, annotation = TRUE, main = "Scenario Zeta")

zeta_500.pvalModel <- zcurve(p = zeta_500$pval_list,
                             control = list(parallel = TRUE))
zeta_500.pvalPlot <- plot(zeta_500.pvalModel,
                          CI = TRUE, annotation = TRUE, main = "Scenario Zeta")
zeta_500.pvals <- hist(zeta_500$pval_list)

#-------------------------------------------------------------------------------

