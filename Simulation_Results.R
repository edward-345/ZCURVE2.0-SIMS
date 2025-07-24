# SIMULATION RESULTS
source("pHacking_Null.R")
source("combined_Null.R")
source("fixed_sigZ.R")
source("combined_fixed.R")
#-------------------------------------------------------------------------------
# FIXED NUMBER OF TESTS
#-------------------------------------------------------------------------------

A_null <- A_sim(5000)
summary(A_null$fit_A)
A_null.plot <- plot(A_null$fit_A,
                       CI = TRUE, annotation = TRUE,
               main = "Scenario A under Null Hypothesis")

A_null.pvals <- zcurve(p = A_null$pvalues_scenarioA,
                       control = list(parallel = TRUE))
A_null.pvals.plot <- plot(A_null.pvals, ymax = 10,
                          CI = TRUE, annotation = TRUE,
                          main = "Scenario A P-vals under Null Hypothesis")
# ZCURVE 3.0
ymax <- 0.8
TEST4HETEROGENEITY <- 0
TEST4BIAS <- TRUE
A_null.3.0 <- Zing(pval_converter(A_null$pvalues_scenarioA))

#----------------------------------------------

B_null <- B_sim(5000)
summary(B_null$fit_B)
B_null.plot <- plot(B_null$fit_B,
                    CI = TRUE, annotation = TRUE,
                    main = "Scenario B under Null Hypothesis")

B_null.pvals <- zcurve(p = B_null$pvalues_scenarioB,
                       control = list(parallel = TRUE))
B_null.pvals.plot <- plot(B_null.pvals, ymax = 10,
                          CI = TRUE, annotation = TRUE,
                          main = "Scenario B P-vals under Null Hypothesis")
# ZCURVE 3.0
ymax <- 0.8
TEST4HETEROGENEITY <- 0
TEST4BIAS <- TRUE
B_null.3.0 <- Zing(pval_converter(B_null$pvalues_scenarioB))

#----------------------------------------------

C_null <- C_sim(5000)
summary(C_null$fit_C)
C_null.plot <- plot(C_null$fit_C,
                    CI = TRUE, annotation = TRUE,
                    main = "Scenario C under Null Hypothesis")

C_null.pvals <- zcurve(p = C_null$pvalues_scenarioC,
                       control = list(parallel = TRUE))
C_null.pvals.plot <- plot(C_null.pvals, ymax = 10,
                          CI = TRUE, annotation = TRUE,
                          main = "Scenario C P-vals under Null Hypothesis")
# ZCURVE 3.0
ymax <- 0.8
TEST4HETEROGENEITY <- 0
TEST4BIAS <- TRUE
C_null.3.0 <- Zing(pval_converter(C_null$pvalues_scenarioC))

#----------------------------------------------

D_null <- D_sim(5000)
summary(D_null$fit_D)
D_null.plot <- plot(D_null$fit_D,
                    CI = TRUE, annotation = TRUE,
                    main = "Scenario D under Null Hypothesis")

D_null.pvals <- zcurve(p = D_null$pvalues_scenarioD,
                       control = list(parallel = TRUE))
D_null.pvals.plot <- plot(D_null.pvals, ymax = 10,
                          CI = TRUE, annotation = TRUE,
                          main = "Scenario D P-vals under Null Hypothesis")
# ZCURVE 3.0
ymax <- 0.8
TEST4HETEROGENEITY <- 0
TEST4BIAS <- TRUE
D_null.3.0 <- Zing(pval_converter(D_null$pvalues_scenarioD))

#-------------------------------------------------------------------------------
# COMBINED SCENARIOS WITH FIXED NUMBER OF SIGNIFICANT Z-SCORES
#-------------------------------------------------------------------------------

X_null <- X_sim(5000)
summary(X_null$fit_X)
X_null.plot <- plot(X_null$fit_X,
                    CI = TRUE, annotation = TRUE,
                    main = "Scenario X under Null Hypothesis")

X_null.pvals <- zcurve(p = X_null$pvalues_scenarioX,
                       control = list(parallel = TRUE))
X_null.pvals.plot <- plot(X_null.pvals, ymax = 10,
                          CI = TRUE, annotation = TRUE,
                          main = "Scenario X P-vals under Null Hypothesis")
# ZCURVE 3.0
ymax <- 0.8
TEST4HETEROGENEITY <- 0
TEST4BIAS <- TRUE
#Int.Beg <- z.crit + just
X_null.3.0 <- Zing(pval_converter(X_null$pvalues_scenarioX))

#----------------------------------------------

Y_null <- Y_sim(5000)
summary(Y_null$fit_Y)
Y_null.plot <- plot(Y_null$fit_Y,
                    CI = TRUE, annotation = TRUE,
                    main = "Scenario Y under Null Hypothesis")

Y_null.pvals <- zcurve(p = Y_null$pvalues_scenarioY,
                       control = list(parallel = TRUE))
Y_null.pvals.plot <- plot(Y_null.pvals, ymax = 10,
                          CI = TRUE, annotation = TRUE,
                          main = "Scenario Y P-vals under Null Hypothesis")
# ZCURVE 3.0
ymax <- 0.8
TEST4HETEROGENEITY <- 0
TEST4BIAS <- TRUE
Y_null.3.0 <- Zing(pval_converter(Y_null$pvalues_scenarioY))

#----------------------------------------------

Z_null <- Z_sim(5000)
summary(Z_null$fit_Z)
Z_null.plot <- plot(Z_null$fit_Z,
                    CI = TRUE, annotation = TRUE,
                    main = "Scenario Z under Null Hypothesis")

Z_null.pvals <- zcurve(p = Z_null$pvalues_scenarioZ,
                       control = list(parallel = TRUE))
Z_null.pvals.plot <- plot(Z_null.pvals, ymax = 10,
                          CI = TRUE, annotation = TRUE,
                          main = "Scenario Z P-vals under Null Hypothesis")
# ZCURVE 3.0
ymax <- 1.0
TEST4HETEROGENEITY <- 0
TEST4BIAS <- TRUE
Z_null.3.0 <- Zing(pval_converter(Z_null$pvalues_scenarioZ))


#-------------------------------------------------------------------------------
# FIXED NUMBER OF SIGNIFICANT Z-SCORES
#-------------------------------------------------------------------------------

alpha_null <- alpha_sim(5000)
summary(alpha_null$fit_alpha)
alpha_null.plot <- plot(alpha_null$fit_alpha,
                    CI = TRUE, annotation = TRUE,
                    main = "Scenario Alpha under Null Hypothesis")

alpha_null.pvals <- zcurve(p = alpha_null$pval_list,
                       control = list(parallel = TRUE))
alpha_null.pvals.plot <- plot(alpha_null.pvals, ymax = 10,
                          CI = TRUE, annotation = TRUE,
                          main = "Scenario alpha P-vals under Null Hypothesis")
# ZCURVE 3.0
ymax <- 0.8
TEST4HETEROGENEITY <- 0
TEST4BIAS <- TRUE
alpha_null.3.0 <- Zing(pval_converter(alpha_null$pval_list))

#----------------------------------------------

beta_null <- beta_sim(5000)
summary(beta_null$fit_beta)
beta_null.plot <- plot(beta_null$fit_beta,
                        CI = TRUE, annotation = TRUE,
                        main = "Scenario Beta under Null Hypothesis")

beta_null.pvals <- zcurve(p = beta_null$pval_list,
                           control = list(parallel = TRUE))
beta_null.pvals.plot <- plot(beta_null.pvals, ymax = 10,
                              CI = TRUE, annotation = TRUE,
                              main = "Scenario Beta P-vals under Null Hypothesis")
# ZCURVE 3.0
ymax <- 0.8
TEST4HETEROGENEITY <- 0
TEST4BIAS <- TRUE
beta_null.3.0 <- Zing(pval_converter(beta_null$pval_list))

#----------------------------------------------

gamma_null <- gamma_sim(5000)
summary(gamma_null$fit_gamma)
gamma_null.plot <- plot(gamma_null$fit_gamma,
                       CI = TRUE, annotation = TRUE,
                       main = "Scenario Gamma under Null Hypothesis")

gamma_null.pvals <- zcurve(p = gamma_null$pval_list,
                          control = list(parallel = TRUE))
gamma_null.pvals.plot <- plot(gamma_null.pvals, ymax = 10,
                             CI = TRUE, annotation = TRUE,
                             main = "Scenario Gamma P-vals under Null Hypothesis")
# ZCURVE 3.0
ymax <- 0.8
TEST4HETEROGENEITY <- 0
TEST4BIAS <- TRUE
gamma_null.3.0 <- Zing(pval_converter(gamma_null$pval_list))

#----------------------------------------------

delta_null <- delta_sim(5000)
summary(delta_null$fit_delta)
delta_null.plot <- plot(delta_null$fit_delta,
                       CI = TRUE, annotation = TRUE,
                       main = "Scenario Delta under Null Hypothesis")

delta_null.pvals <- zcurve(p = delta_null$pval_list,
                          control = list(parallel = TRUE))
delta_null.pvals.plot <- plot(delta_null.pvals, ymax = 10,
                             CI = TRUE, annotation = TRUE,
                             main = "Scenario Delta P-vals under Null Hypothesis")
# ZCURVE 3.0
ymax <- 0.8
TEST4HETEROGENEITY <- 0
TEST4BIAS <- TRUE
delta_null.3.0 <- Zing(pval_converter(delta_null$pval_list))

#-------------------------------------------------------------------------------
#COMBINED SCENARIOS WITH FIXED NUMBER OF SIGNIFICANT ZSCORES
#-------------------------------------------------------------------------------

chi_null <- chi_sim(5000)
summary(chi_null$fit_chi)
chi_null.plot <- plot(chi_null$fit_chi,
                       CI = TRUE, annotation = TRUE,
                       main = "Scenario Chi under Null Hypothesis")

chi_null.pvals <- zcurve(p = chi_null$pval_list,
                          control = list(parallel = TRUE))
chi_null.pvals.plot <- plot(chi_null.pvals, ymax = 10,
                             CI = TRUE, annotation = TRUE,
                             main = "Scenario Chi P-vals under Null Hypothesis")
# ZCURVE 3.0
ymax <- 0.8
TEST4HETEROGENEITY <- 0
TEST4BIAS <- TRUE
chi_null.3.0 <- Zing(pval_converter(chi_null$pval_list))

#----------------------------------------------

psi_null <- psi_sim(5000)
summary(psi_null$fit_psi)
psi_null.plot <- plot(psi_null$fit_psi,
                       CI = TRUE, annotation = TRUE,
                       main = "Scenario Psi under Null Hypothesis")

psi_null.pvals <- zcurve(p = psi_null$pval_list,
                          control = list(parallel = TRUE))
psi_null.pvals.plot <- plot(psi_null.pvals, ymax = 10,
                             CI = TRUE, annotation = TRUE,
                             main = "Scenario Psi P-vals under Null Hypothesis")
# ZCURVE 3.0
ymax <- 0.8
TEST4HETEROGENEITY <- 0
TEST4BIAS <- TRUE
psi_null.3.0 <- Zing(pval_converter(psi_null$pval_list))

#----------------------------------------------

zeta_null <- zeta_sim(5000)
summary(zeta_null$fit_zeta)
zeta_null.plot <- plot(zeta_null$fit_zeta,
                       CI = TRUE, annotation = TRUE,
                       main = "Scenario Zeta under Null Hypothesis")

zeta_null.pvals <- zcurve(p = zeta_null$pval_list,
                          control = list(parallel = TRUE))
zeta_null.pvals.plot <- plot(zeta_null.pvals, ymax = 10,
                             CI = TRUE, annotation = TRUE,
                             main = "Scenario Zeta P-vals under Null Hypothesis")
# ZCURVE 3.0
ymax <- 1.2
TEST4HETEROGENEITY <- 0
TEST4BIAS <- TRUE
zeta_null.3.0 <- Zing(pval_converter(zeta_null$pval_list))

#-------------------------------------------------------------------------------

