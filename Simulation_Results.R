# SIMULATION RESULTS
source("pHacking_Null.R")
source("combined_Null.R")
source("fixed_sigZ.R")
source("combined_fixed.R")
source("MultiVar_sims.R")
#-------------------------------------------------------------------------------
# FIXED NUMBER OF TESTS
#-------------------------------------------------------------------------------

#Default setting (no p-hacking)
sim_null <- simulation(5000)
summary(sim_null$fit)
sim_null.plot <- plot(sim_null$fit,
                    CI = TRUE, annotation = TRUE,
                    main = "Default under Null Hypothesis")

sim_null.pvals <- zcurve(p = sim_null$pvals,
                       control = list(parallel = TRUE))
sim_null.pvals.plot <- plot(sim_null.pvals, ymax = 10,
                          CI = TRUE, annotation = TRUE,
                          main = "Default P-vals under Null Hypothesis")
# ZCURVE 3.0
ymax <- 0.8
TEST4HETEROGENEITY <- 0
TEST4BIAS <- TRUE
sim_null.3.0 <- Zing(pval_converter(sim_null$pvals))

#----------------------------------------------

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
#MULTIVAR SIMULATIONS
#-------------------------------------------------------------------------------

# 3 Dependent Variables under null
multi3 <- multivar_sim(5000, control_mu = c(0,0,0), exp_mu = c(0,0,0))
summary(multi3$fit_A)
multi3.plot <- plot(multi3$fit_A,
                    CI = TRUE, annotation = TRUE,
                    main = "3 DVs under Null Hypothesis")

multi3.pvals <- zcurve(p = multi3$pvalues_scenarioA,
                       control = list(parallel = TRUE))
summary(multi3.pvals)
multi3.pvals.plot <- plot(multi3.pvals, ymax = 10,
                          CI = TRUE, annotation = TRUE,
                          main = "3DVs P-vals under Null Hypothesis")
# ZCURVE 3.0
source(zcurve3)
ymax <- 1.2
TEST4HETEROGENEITY <- 0
TEST4BIAS <- TRUE
multi3_3.0 <- Zing(pval_converter(multi3$pvalues_scenarioA))

#----------------------------------------------
# 4 Dependent Variables
multi4 <- multivar_sim(5000, control_mu = c(0,0,0,0), exp_mu = c(0,0,0,0))
summary(multi4$fit_A)
multi4.plot <- plot(multi4$fit_A,
                    CI = TRUE, annotation = TRUE,
                    main = "4 DVs under Null Hypothesis")

multi4.pvals <- zcurve(p = multi4$pvalues_scenarioA,
                       control = list(parallel = TRUE))
summary(multi4.pvals)
multi4.pvals.plot <- plot(multi4.pvals, ymax = 10,
                          CI = TRUE, annotation = TRUE,
                          main = "4DVs P-vals under Null Hypothesis")
# ZCURVE 3.0
source(zcurve3)
ymax <- 1.2
TEST4HETEROGENEITY <- 0
TEST4BIAS <- TRUE
multi4_3.0 <- Zing(pval_converter(multi4$pvalues_scenarioA))

#----------------------------------------------
# 5 Dependent Variables
multi5 <- multivar_sim(5000, control_mu = c(0,0,0,0,0), exp_mu = c(0,0,0,0,0))
summary(multi5$fit_A)
multi5.plot <- plot(multi5$fit_A,
                    CI = TRUE, annotation = TRUE,
                    main = "5 DVs under Null Hypothesis")

multi5.pvals <- zcurve(p = multi5$pvalues_scenarioA,
                       control = list(parallel = TRUE))
summary(multi5.pvals)
multi5.pvals.plot <- plot(multi5.pvals, ymax = 10,
                          CI = TRUE, annotation = TRUE,
                          main = "5DVs P-vals under Null Hypothesis")
# ZCURVE 3.0
source(zcurve3)
ymax <- 1.2
TEST4HETEROGENEITY <- 0
TEST4BIAS <- TRUE
multi5_3.0 <- Zing(pval_converter(multi5$pvalues_scenarioA))

#----------------------------------------------
# 5 Dependent Variables with Fixed
multi5_fixed <- multivar_fixed(5000,
                               control_mu = c(0,0,0,0,0),
                               exp_mu = c(0,0,0,0,0))
summary(multi5_fixed$fit_A)
multi5.plot <- plot(multi5_fixed$fit_A,
                    CI = TRUE, annotation = TRUE,
                    main = "5 DVs under Null Hypothesis (Fixed)")

multi5.pvals <- zcurve(p = multi5_fixed$pvals,
                       control = list(parallel = TRUE))
summary(multi5.pvals)
multi5.pvals.plot <- plot(multi5.pvals, ymax = 10,
                          CI = TRUE, annotation = TRUE,
                          main = "5DVs P-vals under Null Hypothesis")
# ZCURVE 3.0
source(zcurve3)
ymax <- 1.2
TEST4HETEROGENEITY <- 0
TEST4BIAS <- TRUE
multi5_fixed.3.0 <- Zing(pval_converter(multi5_fixed$pvals))

#----------------------------------------------
# 5 Dependent Variables True Effect
multi5_med <- multivar_sim(5000, control_mu = c(0,0,0,0,0),
                           exp_mu = c(0.5,0.5,0.5,0.5,0.5))
summary(multi5_med$fit_A)
multi5_med.plot <- plot(multi5_med$fit_A,
                        CI = TRUE, annotation = TRUE,
                        main = "5 DVs under Null Hypothesis")

multi5_med.pvals <- zcurve(p = multi5_med$pvalues_scenarioA,
                           control = list(parallel = TRUE))
summary(multi5_med.pvals)
multi5_med.pvals.plot <- plot(multi5_med.pvals, ymax = 10,
                              CI = TRUE, annotation = TRUE,
                              main = "5DVs P-vals under Null Hypothesis")
# ZCURVE 3.0
source(zcurve3)
ymax <- 1.2
TEST4HETEROGENEITY <- 0
TEST4BIAS <- TRUE
multi5_3.0 <- Zing(pval_converter(multi5_med$pvalues_scenarioA))

#----------------------------------------------
tester_5000 <- multivar_sim(1000, n = 40,
                            control_mu = rep(0, 5), exp_mu = rep(0.5, 5))
summary(tester_5000$fit_A)
tester_5000.plot <- plot(tester_5000$fit_A,
                         CI = TRUE, annotation = TRUE,
                         main = "5 DVs under Null Hypothesis")

tester_5000.pvals <- zcurve(p = tester_5000$total_pvals,
                            control = list(parallel = TRUE))
summary(tester_5000.pvals)
tester_5000.pvals.plot <- plot(tester_5000.pvals, ymax = 10,
                               CI = TRUE, annotation = TRUE,
                               main = "Total Pvals with 5DVunder Null Hypothesis")
# ZCURVE 3.0
source(zcurve3)
ymax <- 0.8
TEST4HETEROGENEITY <- 0
#TEST4BIAS <- TRUE
Title <- paste("EDR: ", round(true.edr*100), " ERR: ", round(true.err*100))
sim.z <- pval_converter(tester_5000$total_pvals)
results_tester1000 <- Zing(sim.z)

#----------------------------------------------

corr_test <- multivar_sim(3000, n = 40,
                          control_mu = rep(0, 5), exp_mu = rep(0.5, 5), r=0.7)
summary(corr_test$fit_A)
corr_test.plot <- plot(corr_test$fit_A,
                       CI = TRUE, annotation = TRUE,
                       main = "5DVs correlated by 0.9")

corr_test.pvals <- zcurve(p = corr_test$pvalues_scenarioA,
                          control = list(parallel = TRUE))
summary(corr_test.pvals)
corr_test.pvals.plot <- plot(corr_test.pvals, ymax = 10,
                             CI = TRUE, annotation = TRUE,
                             main = "5DVs correlated by 0.9 (pvals)")

corr_test.total <- zcurve(p = corr_test$total_pvals,
                          control = list(parallel = TRUE))
summary(corr_test.total)
corr_test.total.plot <- plot(corr_test.pvals, ymax = 10,
                             CI = TRUE, annotation = TRUE,
                             main = "5DVs correlated by 0.9 (total)")

# ZCURVE 3.0
source(zcurve3)
ymax <- 1.2
TEST4HETEROGENEITY <- 0
#TEST4BIAS <- TRUE
#Title <- paste("EDR: ", round(true.edr*100), " ERR: ", round(true.err*100))
sim.z <- pval_converter(corr_test$total_pvals)
#zsds = rep(0.7, 7)
results_corr_test <- Zing(sim.z)

source(zcurve3)
ymax <- 1.2
TEST4HETEROGENEITY <- 0
#TEST4BIAS <- TRUE
#Title <- paste("EDR: ", round(true.edr*100), " ERR: ", round(true.err*100))
sim.z <- pval_converter(corr_test$total_pvals)
Est.Method = "EXT"
ncz = 2
zsds = rep(5, 1)
Int.Beg = 0
new_corr <- Zing(sim.z)

