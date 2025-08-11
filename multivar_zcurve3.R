source('helper_functions.R')
source('MultiVar_sims.R')
source(zcurve3)
library(caret)
#-------------------------------------------------------------------------------
# EDR AND ERR FUNCTIONS USING ZCURVE3
#-------------------------------------------------------------------------------
# EDR of MultiVar sim
edr_3 <- function(k_sims, n = 20, r = 0.5, 
                      control_mu = c(0,0), exp_mu = c(0,0), sd = 1) {
  zscores_A <- numeric(k_sims)
  A <- 1
  
  for (i in 1:k_sims) {
    st.d <- rep(sd, times = length(control_mu))
    
    control_A <- rnorm_multi(
      n, vars = length(control_mu), control_mu, st.d, r)
    exp_A <- rnorm_multi(
      n, vars = length(exp_mu), exp_mu, st.d, r)
    
    pvalue_A <- multivar_ttests(control_A, exp_A)
    min_pvalue <- min(pvalue_A)
    
    #If the smallest p-value of the three t-tests is significant at .05, convert 
    #to z-score and add to zscores_A vector
    if (min_pvalue <= 0.05) {
      zvalue_A <- pval_converter(min_pvalue)
      zscores_A[A] <- zvalue_A
      A <- A + 1
    }
  }
  
  #Trim unused cells
  zscores_A <- zscores_A[1:(A - 1)]
  
  fit_A <- Zing(zscores_A)
  
  edr <- fit_A$res['EDR']
  
  return(edr)
}

#-------------------------------------------------------------------------------
#ERR of MultiVar sim
err_3 <- function(k_sims, n = 20, r = 0.5, 
                      control_mu = c(0,0), exp_mu = c(0,0), sd = 1) {
  zscores_A <- numeric(k_sims)
  A <- 1
  
  for (i in 1:k_sims) {
    st.d <- rep(sd, times = length(control_mu))
    
    control_A <- rnorm_multi(
      n, vars = length(control_mu), control_mu, st.d, r)
    exp_A <- rnorm_multi(
      n, vars = length(exp_mu), exp_mu, st.d, r)
    
    pvalue_A <- multivar_ttests(control_A, exp_A)
    min_pvalue <- min(pvalue_A)
    
    #If the smallest p-value of the three t-tests is significant at .05, convert 
    #to z-score and add to zscores_A vector
    if (min_pvalue <= 0.05) {
      zvalue_A <- pval_converter(min_pvalue)
      zscores_A[A] <- zvalue_A
      A <- A + 1
    }
  }
  
  #Trim unused cells
  zscores_A <- zscores_A[1:(A - 1)]
  
  fit_A <- Zing(zscores_A)
  
  err <- coef(fit_A)["ERR"] 
  
  return(err)
}



#-------------------------------------------------------------------------------
# ZCURVE3 GLOBAL PARAMETERS
#-------------------------------------------------------------------------------
ymax <- 1.3
# TEST4HETEROGENEITY <- 0
# TEST4BIAS <- TRUE
# Title <- past("title string")
# Est.method <- "EXT"
# boot.iter <- 500
# ncz <- 0:6                          # Component locations (z-values at which densities are centered)
# components <- length(ncz)           # Number of components
# zsd <- 1                            # SD of standard normal z-distribution
# zsds = rep(zsd,components)          # one SD for each component
# Int.Beg <- 1.96
# Int.End <- 6
#-------------------------------------------------------------------------------
# BEHAVIOR AS n INCREASES
n_vals <- seq(20, 100, by = 10)

het_truEDR.n <- sapply(n_vals, function(n)
  true_edr(
    n = n,
    alpha = .05,
    pop.es = c(0, 0.2,0.4,0.6,0.8),
    wgt = rep(0.2, 5)))

null_truEDR.n <- sapply(n_vals, function(n)
  true_edr(
    n = n,
    alpha = .05,
    pop.es = rep(0,5),
    wgt = c(1,0,0,0,0)))

two_truEDR.n <- sapply(n_vals, function(n)
  true_edr(
    n = n,
    alpha = .05,
    pop.es = c(0.2, 0.8, 0.2, 0.8, 0.2),
    wgt = c(0.333,0,0.333,0,0.333)))

hom_truEDR.n <- sapply(n_vals, function(n)
  true_edr(
    n = n,
    alpha = .05,
    pop.es = rep(0.6, 5),
    wgt = c(1,0,0,0,0)))

het_truERR.n <- sapply(n_vals, function(n)
  true_err(
    n = n,
    alpha = .05,
    pop.es = c(0, 0.2,0.4,0.6,0.8),
    wgt = rep(0.2, 5)))

null_truERR.n <- sapply(n_vals, function(n)
  true_err(
    n = n,
    alpha = .05,
    pop.es = rep(0,5),
    wgt = c(1,0,0,0,0)))

two_truERR.n <- sapply(n_vals, function(n)
  true_err(
    n = n,
    alpha = .05,
    pop.es = c(0.2, 0.8, 0.2, 0.8, 0.2),
    wgt = c(0.333,0,0.333,0,0.333)))

hom_truERR.n <- sapply(n_vals, function(n)
  true_err(
    n = n,
    alpha = .05,
    pop.es = rep(0.6, 5),
    wgt = c(1,0,0,0,0)))
#-------------------------------------------------------------------------------

# HETEROGENEOUS ES = [0, .2, .4, .6, .8]
het_sim.n_20 <- multivar_sim(k_sims = 1000, n = 20,
                             control_mu = rep(0, 5),
                             exp_mu = c(0, .2, .4, .6, .8))
het_model.n_20 <- Zing(het_sim.n_20$zscores_A)
het_tru.model.n_20 <- Zing(het_sim.n_20$total_zscores)

het_sim.n_30 <- multivar_sim(k_sims = 1000, n = 30,
                             control_mu = rep(0, 5),
                             exp_mu = c(0, .2, .4, .6, .8))
het_model.n_30 <- Zing(het_sim.n_30$zscores_A)
het_tru.model.n_30 <- Zing(het_sim.n_30$total_zscores)

het_sim.n_40 <- multivar_sim(k_sims = 1000, n = 40,
                             control_mu = rep(0, 5),
                             exp_mu = c(0, .2, .4, .6, .8))
het_model.n_40 <- Zing(het_sim.n_40$zscores_A)
het_tru.model.n_40 <- Zing(het_sim.n_40$total_zscores)

het_sim.n_50 <- multivar_sim(k_sims = 1000, n = 50,
                             control_mu = rep(0, 5),
                             exp_mu = c(0, .2, .4, .6, .8))
het_model.n_50 <- Zing(het_sim.n_50$zscores_A)
het_tru.model.n_50 <- Zing(het_sim.n_50$total_zscores)

het_sim.n_60 <- multivar_sim(k_sims = 1000, n = 60,
                             control_mu = rep(0, 5),
                             exp_mu = c(0, .2, .4, .6, .8))
het_model.n_60 <- Zing(het_sim.n_60$zscores_A)
het_tru.model.n_60 <- Zing(het_sim.n_60$total_zscores)

het_sim.n_70 <- multivar_sim(k_sims = 1000, n = 70,
                             control_mu = rep(0, 5),
                             exp_mu = c(0, .2, .4, .6, .8))
het_model.n_70 <- Zing(het_sim.n_70$zscores_A)
het_tru.model.n_70 <- Zing(het_sim.n_70$total_zscores)

het_sim.n_80 <- multivar_sim(k_sims = 1000, n = 80,
                             control_mu = rep(0, 5),
                             exp_mu = c(0, .2, .4, .6, .8))
het_model.n_80 <- Zing(het_sim.n_80$zscores_A)
het_tru.model.n_80 <- Zing(het_sim.n_80$total_zscores)

het_sim.n_90 <- multivar_sim(k_sims = 1000, n = 90,
                             control_mu = rep(0, 5),
                             exp_mu = c(0, .2, .4, .6, .8))
het_model.n_90 <- Zing(het_sim.n_90$zscores_A)
het_tru.model.n_90 <- Zing(het_sim.n_90$total_zscores)

het_sim.n_100 <- multivar_sim(k_sims = 1000, n = 100,
                             control_mu = rep(0, 5),
                             exp_mu = c(0, .2, .4, .6, .8))
het_model.n_100 <- Zing(het_sim.n_100$zscores_A)
het_tru.model.n_100 <- Zing(het_sim.n_100$total_zscores)

het_EDR.n <- z3_EDR(list(het_model.n_20,
                      het_model.n_30,
                      het_model.n_40,
                      het_model.n_50,
                      het_model.n_60,
                      het_model.n_70,
                      het_model.n_80,
                      het_model.n_90,
                      het_model.n_100))
het_ERR.n <- z3_ERR(list(het_model.n_20,
                      het_model.n_30,
                      het_model.n_40,
                      het_model.n_50,
                      het_model.n_60,
                      het_model.n_70,
                      het_model.n_80,
                      het_model.n_90,
                      het_model.n_100))

het_EDR.n_df <- data.frame(n = n_vals,
                           Simulated = het_EDR.n,
                           True = het_truEDR.n)

het_EDR.n_long <- het_EDR.n_df %>%
  pivot_longer(cols = c(Simulated, True),
               names_to = "Type",
               values_to = "EDR")

het_EDR.n_plot <- ggplot(het_EDR.n_long, aes(x = n, y = EDR, group = Type)) +
  geom_line(aes(linetype = Type)) +
  geom_point(aes(shape = Type)) +
  scale_x_continuous(breaks = n_vals) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "EDR vs Sample Size (ES = [0, 0.8])",
       x = "Sample Size (n)",
       y = "Estimated Discovery Rate",
       linetype = NULL, shape = NULL) +
  theme_minimal()

het_EDR.n_plot

het_ERR.n_df <- data.frame(n = n_vals,
                           Simulated = het_ERR.n,
                           True = het_truERR.n)

het_ERR.n_long <- het_ERR.n_df %>%
  pivot_longer(cols = c(Simulated, True),
               names_to = "Type",
               values_to = "ERR")

het_ERR.n_plot <- ggplot(het_ERR.n_long, aes(x = n, y = ERR, group = Type)) +
  geom_line(aes(linetype = Type)) +
  geom_point(aes(shape = Type)) +
  scale_x_continuous(breaks = n_vals) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "ERR vs Sample Size (ES = [0, 0.8])",
       x = "Sample Size (n)",
       y = "Estimated Replicability Rate",
       linetype = NULL, shape = NULL) +
  theme_minimal()

het_ERR.n_plot


#-------------------------------------------------------------------------------
# TRUE NULL ES = [0, 0, 0, 0, 0]
null_sim.n_20 <- multivar_sim(k_sims = 1000, n = 20,
                             control_mu = rep(0, 5),
                             exp_mu = rep(0, 5))
null_model.n_20 <- Zing(null_sim.n_20$zscores_A)
null_tru.model.n_20 <- Zing(null_sim.n_20$total_zscores)

null_sim.n_30 <- multivar_sim(k_sims = 1000, n = 30,
                             control_mu = rep(0, 5),
                             exp_mu = rep(0, 5))
null_model.n_30 <- Zing(null_sim.n_30$zscores_A)
null_tru.model.n_30 <- Zing(null_sim.n_30$total_zscores)

null_sim.n_40 <- multivar_sim(k_sims = 1000, n = 40,
                             control_mu = rep(0, 5),
                             exp_mu = rep(0, 5))
null_model.n_40 <- Zing(null_sim.n_40$zscores_A)
null_tru.model.n_40 <- Zing(null_sim.n_40$total_zscores)

null_sim.n_50 <- multivar_sim(k_sims = 1000, n = 50,
                             control_mu = rep(0, 5),
                             exp_mu = rep(0, 5))
null_model.n_50 <- Zing(null_sim.n_50$zscores_A)
null_tru.model.n_50 <- Zing(null_sim.n_50$total_zscores)

null_sim.n_60 <- multivar_sim(k_sims = 1000, n = 60,
                             control_mu = rep(0, 5),
                             exp_mu = rep(0, 5))
null_model.n_60 <- Zing(null_sim.n_60$zscores_A)
null_tru.model.n_60 <- Zing(null_sim.n_60$total_zscores)

null_sim.n_70 <- multivar_sim(k_sims = 1000, n = 70,
                             control_mu = rep(0, 5),
                             exp_mu = rep(0, 5))
null_model.n_70 <- Zing(null_sim.n_70$zscores_A)
null_tru.model.n_70 <- Zing(null_sim.n_70$total_zscores)

null_sim.n_80 <- multivar_sim(k_sims = 1000, n = 80,
                             control_mu = rep(0, 5),
                             exp_mu = rep(0, 5))
null_model.n_80 <- Zing(null_sim.n_80$zscores_A)
null_tru.model.n_80 <- Zing(null_sim.n_80$total_zscores)

null_sim.n_90 <- multivar_sim(k_sims = 1000, n = 90,
                             control_mu = rep(0, 5),
                             exp_mu = rep(0, 5))
null_model.n_90 <- Zing(null_sim.n_90$zscores_A)
null_tru.model.n_90 <- Zing(null_sim.n_90$total_zscores)

null_sim.n_100 <- multivar_sim(k_sims = 1000, n = 100,
                              control_mu = rep(0, 5),
                              exp_mu = rep(0, 5))
null_model.n_100 <- Zing(null_sim.n_100$zscores_A)
null_tru.model.n_100 <- Zing(null_sim.n_100$total_zscores)

null_EDR.n <- z3_EDR(list(null_model.n_20,
                         null_model.n_30,
                         null_model.n_40,
                         null_model.n_50,
                         null_model.n_60,
                         null_model.n_70,
                         null_model.n_80,
                         null_model.n_90,
                         null_model.n_100))
null_ERR.n <- z3_ERR(list(null_model.n_20,
                         null_model.n_30,
                         null_model.n_40,
                         null_model.n_50,
                         null_model.n_60,
                         null_model.n_70,
                         null_model.n_80,
                         null_model.n_90,
                         null_model.n_100))

null_EDR.n_df <- data.frame(n = n_vals,
                           Simulated = null_EDR.n,
                           True = null_truEDR.n)

null_EDR.n_long <- null_EDR.n_df %>%
  pivot_longer(cols = c(Simulated, True),
               names_to = "Type",
               values_to = "EDR")

null_EDR.n_plot <- ggplot(null_EDR.n_long, aes(x = n, y = EDR, group = Type)) +
  geom_line(aes(linetype = Type)) +
  geom_point(aes(shape = Type)) +
  scale_x_continuous(breaks = n_vals) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "EDR vs Sample Size (True Null)",
       x = "Sample Size (n)",
       y = "Estimated Discovery Rate",
       linetype = NULL, shape = NULL) +
  theme_minimal()

null_EDR.n_plot

null_ERR.n_df <- data.frame(n = n_vals,
                           Simulated = null_ERR.n,
                           True = null_truERR.n)

null_ERR.n_long <- null_ERR.n_df %>%
  pivot_longer(cols = c(Simulated, True),
               names_to = "Type",
               values_to = "ERR")

null_ERR.n_plot <- ggplot(null_ERR.n_long, aes(x = n, y = ERR, group = Type)) +
  geom_line(aes(linetype = Type)) +
  geom_point(aes(shape = Type)) +
  scale_x_continuous(breaks = n_vals) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "ERR vs Sample Size (Heterogeneous ES)",
       x = "Sample Size (n)",
       y = "Estimated Replicability Rate",
       linetype = NULL, shape = NULL) +
  theme_minimal()

null_ERR.n_plot

#-------------------------------------------------------------------------------
# TWO EFFECT SIZES ES = [0.2, 0.8, 0.2, 0.8, 0.2]
two_sim.n_20 <- multivar_sim(k_sims = 1000, n = 20,
                              control_mu = rep(0, 5),
                              exp_mu = c(0.2, 0.8, 0.2, 0.8, 0.2))
two_model.n_20 <- Zing(two_sim.n_20$zscores_A)
two_tru.model.n_20 <- Zing(two_sim.n_20$total_zscores)

two_sim.n_30 <- multivar_sim(k_sims = 1000, n = 30,
                              control_mu = rep(0, 5),
                              exp_mu = c(0.2, 0.8, 0.2, 0.8, 0.2))
two_model.n_30 <- Zing(two_sim.n_30$zscores_A)
two_tru.model.n_30 <- Zing(two_sim.n_30$total_zscores)


two_sim.n_40 <- multivar_sim(k_sims = 1000, n = 40,
                              control_mu = rep(0, 5),
                              exp_mu = c(0.2, 0.8, 0.2, 0.8, 0.2))
two_model.n_40 <- Zing(two_sim.n_40$zscores_A)
two_tru.model.n_40 <- Zing(two_sim.n_40$total_zscores)

two_sim.n_50 <- multivar_sim(k_sims = 1000, n = 50,
                              control_mu = rep(0, 5),
                              exp_mu = c(0.2, 0.8, 0.2, 0.8, 0.2))
two_model.n_50 <- Zing(two_sim.n_50$zscores_A)
two_tru.model.n_50 <- Zing(two_sim.n_50$total_zscores)

two_sim.n_60 <- multivar_sim(k_sims = 1000, n = 60,
                              control_mu = rep(0, 5),
                              exp_mu = c(0.2, 0.8, 0.2, 0.8, 0.2))
two_model.n_60 <- Zing(two_sim.n_60$zscores_A)
two_tru.model.n_60 <- Zing(two_sim.n_60$total_zscores)

two_sim.n_70 <- multivar_sim(k_sims = 1000, n = 70,
                              control_mu = rep(0, 5),
                              exp_mu = c(0.2, 0.8, 0.2, 0.8, 0.2))
two_model.n_70 <- Zing(two_sim.n_70$zscores_A)
two_tru.model.n_70 <- Zing(two_sim.n_70$total_zscores)

two_sim.n_80 <- multivar_sim(k_sims = 1000, n = 80,
                              control_mu = rep(0, 5),
                              exp_mu = c(0.2, 0.8, 0.2, 0.8, 0.2))
two_model.n_80 <- Zing(two_sim.n_80$zscores_A)
two_tru.model.n_80 <- Zing(two_sim.n_80$total_zscores)

two_sim.n_90 <- multivar_sim(k_sims = 1000, n = 90,
                              control_mu = rep(0, 5),
                              exp_mu = c(0.2, 0.8, 0.2, 0.8, 0.2))
two_model.n_90 <- Zing(two_sim.n_90$zscores_A)
two_tru.model.n_90 <- Zing(two_sim.n_90$total_zscores)

two_sim.n_100 <- multivar_sim(k_sims = 1000, n = 100,
                               control_mu = rep(0, 5),
                               exp_mu = c(0.2, 0.8, 0.2, 0.8, 0.2))
two_model.n_100 <- Zing(two_sim.n_100$zscores_A)
two_tru.model.n_100 <- Zing(two_sim.n_100$total_zscores)

two_EDR.n <- z3_EDR(list(two_model.n_20,
                         two_model.n_30,
                         two_model.n_40,
                         two_model.n_50,
                         two_model.n_60,
                         two_model.n_70,
                         two_model.n_80,
                         two_model.n_90,
                         two_model.n_100))
two_ERR.n <- z3_ERR(list(two_model.n_20,
                         two_model.n_30,
                         two_model.n_40,
                         two_model.n_50,
                         two_model.n_60,
                         two_model.n_70,
                         two_model.n_80,
                         two_model.n_90,
                         two_model.n_100))

two_EDR.n_df <- data.frame(n = n_vals,
                            Simulated = two_EDR.n,
                            True = two_truEDR.n)

two_EDR.n_long <- two_EDR.n_df %>%
  pivot_longer(cols = c(Simulated, True),
               names_to = "Type",
               values_to = "EDR")

two_EDR.n_plot <- ggplot(two_EDR.n_long, aes(x = n, y = EDR, group = Type)) +
  geom_line(aes(linetype = Type)) +
  geom_point(aes(shape = Type)) +
  scale_x_continuous(breaks = n_vals) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "EDR vs Sample Size (ES = 0.2 and 0.8)",
       x = "Sample Size (n)",
       y = "Estimated Discovery Rate",
       linetype = NULL, shape = NULL) +
  theme_minimal()

two_EDR.n_plot

two_ERR.n_df <- data.frame(n = n_vals,
                            Simulated = two_ERR.n,
                            True = two_truERR.n)

two_ERR.n_long <- two_ERR.n_df %>%
  pivot_longer(cols = c(Simulated, True),
               names_to = "Type",
               values_to = "ERR")

two_ERR.n_plot <- ggplot(two_ERR.n_long, aes(x = n, y = ERR, group = Type)) +
  geom_line(aes(linetype = Type)) +
  geom_point(aes(shape = Type)) +
  scale_x_continuous(breaks = n_vals) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "ERR vs Sample Size (ES = 0.2 and 0.8)",
       x = "Sample Size (n)",
       y = "Estimated Replicability Rate",
       linetype = NULL, shape = NULL) +
  theme_minimal()

two_ERR.n_plot

#-------------------------------------------------------------------------------
# HOMOGENOUS MED ES = [0.6, 0.6, 0.6, 0.6, 0.6]
hom_sim.n_20 <- multivar_sim(k_sims = 1000, n = 20,
                             control_mu = rep(0, 5),
                             exp_mu = rep(0.6, 5))
hom_model.n_20 <- Zing(hom_sim.n_20$zscores_A)
hom_tru.model.n_20 <- Zing(hom_sim.n_20$total_zscores)

hom_sim.n_30 <- multivar_sim(k_sims = 1000, n = 30,
                             control_mu = rep(0, 5),
                             exp_mu = rep(0.6, 5))
hom_model.n_30 <- Zing(hom_sim.n_30$zscores_A)
hom_tru.model.n_30 <- Zing(hom_sim.n_30$total_zscores)

hom_sim.n_40 <- multivar_sim(k_sims = 1000, n = 40,
                             control_mu = rep(0, 5),
                             exp_mu = rep(0.6, 5))
hom_model.n_40 <- Zing(hom_sim.n_40$zscores_A)
hom_tru.model.n_40 <- Zing(hom_sim.n_40$total_zscores)

hom_sim.n_50 <- multivar_sim(k_sims = 1000, n = 50,
                             control_mu = rep(0, 5),
                             exp_mu = rep(0.6, 5))
hom_model.n_50 <- Zing(hom_sim.n_50$zscores_A)
hom_tru.model.n_50 <- Zing(hom_sim.n_50$total_zscores)

hom_sim.n_60 <- multivar_sim(k_sims = 1000, n = 60,
                             control_mu = rep(0, 5),
                             exp_mu = rep(0.6, 5))
hom_model.n_60 <- Zing(hom_sim.n_60$zscores_A)
hom_tru.model.n_60 <- Zing(hom_sim.n_60$total_zscores)

hom_sim.n_70 <- multivar_sim(k_sims = 1000, n = 70,
                             control_mu = rep(0, 5),
                             exp_mu = rep(0.6, 5))
hom_model.n_70 <- Zing(hom_sim.n_70$zscores_A)
hom_tru.model.n_70 <- Zing(hom_sim.n_70$total_zscores)

hom_sim.n_80 <- multivar_sim(k_sims = 1000, n = 80,
                             control_mu = rep(0, 5),
                             exp_mu = rep(0.6, 5))
hom_model.n_80 <- Zing(hom_sim.n_80$zscores_A)
hom_tru.model.n_80 <- Zing(hom_sim.n_80$total_zscores)

hom_sim.n_90 <- multivar_sim(k_sims = 1000, n = 90,
                             control_mu = rep(0, 5),
                             exp_mu = rep(0.6, 5))
hom_model.n_90 <- Zing(hom_sim.n_90$zscores_A)
hom_tru.model.n_90 <- Zing(hom_sim.n_90$total_zscores)

hom_sim.n_100 <- multivar_sim(k_sims = 1000, n = 100,
                              control_mu = rep(0, 5),
                              exp_mu = rep(0.6, 5))
hom_model.n_100 <- Zing(hom_sim.n_100$zscores_A)
hom_tru.model.n_100 <- Zing(hom_sim.n_100$total_zscores)

hom_EDR.n <- z3_EDR(list(hom_model.n_20,
                         hom_model.n_30,
                         hom_model.n_40,
                         hom_model.n_50,
                         hom_model.n_60,
                         hom_model.n_70,
                         hom_model.n_80,
                         hom_model.n_90,
                         hom_model.n_100))
hom_ERR.n <- z3_ERR(list(hom_model.n_20,
                         hom_model.n_30,
                         hom_model.n_40,
                         hom_model.n_50,
                         hom_model.n_60,
                         hom_model.n_70,
                         hom_model.n_80,
                         hom_model.n_90,
                         hom_model.n_100))

hom_EDR.n_df <- data.frame(n = n_vals,
                            Simulated = hom_EDR.n,
                            True = hom_truEDR.n)

hom_EDR.n_long <- hom_EDR.n_df %>%
  pivot_longer(cols = c(Simulated, True),
               names_to = "Type",
               values_to = "EDR")

hom_EDR.n_plot <- ggplot(hom_EDR.n_long, aes(x = n, y = EDR, group = Type)) +
  geom_line(aes(linetype = Type)) +
  geom_point(aes(shape = Type)) +
  scale_x_continuous(breaks = n_vals) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "EDR vs Sample Size (Homogenous ES = 0.6)",
       x = "Sample Size (n)",
       y = "Estimated Discovery Rate",
       linetype = NULL, shape = NULL) +
  theme_minimal()

hom_EDR.n_plot

hom_ERR.n_df <- data.frame(n = n_vals,
                            Simulated = hom_ERR.n,
                            True = hom_truERR.n)

hom_ERR.n_long <- hom_ERR.n_df %>%
  pivot_longer(cols = c(Simulated, True),
               names_to = "Type",
               values_to = "ERR")

hom_ERR.n_plot <- ggplot(hom_ERR.n_long, aes(x = n, y = ERR, group = Type)) +
  geom_line(aes(linetype = Type)) +
  geom_point(aes(shape = Type)) +
  scale_x_continuous(breaks = n_vals) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "ERR vs Sample Size (Homogenous ES = 0.6)",
       x = "Sample Size (n)",
       y = "Estimated Replicability Rate",
       linetype = NULL, shape = NULL) +
  theme_minimal()

hom_ERR.n_plot
#-------------------------------------------------------------------------------
# BEHAVIOR AS r INCREASES
#-------------------------------------------------------------------------------
# HETEROGENEOUS ES = [0, .2, .4, .6, .8]
het_sim.r.0 <- multivar_sim(k_sims = 1000, r = 0,
                             control_mu = rep(0, 5),
                             exp_mu = c(0, .2, .4, .6, .8))
het_model.r.0 <- Zing(het_sim.r.0$zscores_A)
het_tru.model.r.0 <- Zing(het_sim.r.0$total_zscores)

het_sim.r.1 <- multivar_sim(k_sims = 1000, r = 0.1,
                            control_mu = rep(0, 5),
                            exp_mu = c(0, .2, .4, .6, .8))
het_model.r.1 <- Zing(het_sim.r.1$zscores_A)
het_tru.model.r.1 <- Zing(het_sim.r.1$total_zscores)

het_sim.r.2 <- multivar_sim(k_sims = 1000, r = 0.2,
                            control_mu = rep(0, 5),
                            exp_mu = c(0, .2, .4, .6, .8))
het_model.r.2 <- Zing(het_sim.r.2$zscores_A)
het_tru.model.r.2 <- Zing(het_sim.r.2$total_zscores)

het_sim.r.3 <- multivar_sim(k_sims = 1000, r = 0.3,
                            control_mu = rep(0, 5),
                            exp_mu = c(0, .2, .4, .6, .8))
het_model.r.3 <- Zing(het_sim.r.3$zscores_A)
het_tru.model.r.3 <- Zing(het_sim.r.3$total_zscores)

het_sim.r.4 <- multivar_sim(k_sims = 1000, r = 0.4,
                            control_mu = rep(0, 5),
                            exp_mu = c(0, .2, .4, .6, .8))
het_model.r.4 <- Zing(het_sim.r.4$zscores_A)
het_tru.model.r.4 <- Zing(het_sim.r.4$total_zscores)

het_sim.r.5 <- multivar_sim(k_sims = 1000, r = 0.5,
                            control_mu = rep(0, 5),
                            exp_mu = c(0, .2, .4, .6, .8))
het_model.r.5 <- Zing(het_sim.r.5$zscores_A)
het_tru.model.r.5 <- Zing(het_sim.r.5$total_zscores)

het_sim.r.6 <- multivar_sim(k_sims = 1000, r = 0.6,
                            control_mu = rep(0, 5),
                            exp_mu = c(0, .2, .4, .6, .8))
het_model.r.6 <- Zing(het_sim.r.6$zscores_A)
het_tru.model.r.6 <- Zing(het_sim.r.6$total_zscores)

het_sim.r.7 <- multivar_sim(k_sims = 1000, r = 0.7,
                            control_mu = rep(0, 5),
                            exp_mu = c(0, .2, .4, .6, .8))
het_model.r.7 <- Zing(het_sim.r.7$zscores_A)
het_tru.model.r.7 <- Zing(het_sim.r.7$total_zscores)

het_sim.r.8 <- multivar_sim(k_sims = 1000, r = 0.8,
                            control_mu = rep(0, 5),
                            exp_mu = c(0, .2, .4, .6, .8))
het_model.r.8 <- Zing(het_sim.r.8$zscores_A)
het_tru.model.r.8 <- Zing(het_sim.r.8$total_zscores)

het_sim.r.9 <- multivar_sim(k_sims = 1000, r = 0.9,
                            control_mu = rep(0, 5),
                            exp_mu = c(0, .2, .4, .6, .8))
het_model.r.9 <- Zing(het_sim.r.9$zscores_A)
het_tru.model.r.9 <- Zing(het_sim.r.9$total_zscores)

het_EDR.r <- z3_EDR(list(het_model.r.0,
                         het_model.r.1,
                         het_model.r.2,
                         het_model.r.3,
                         het_model.r.4,
                         het_model.r.5,
                         het_model.r.6,
                         het_model.r.7,
                         het_model.r.8,
                         het_model.r.9))
het_ERR.r <- z3_ERR(list(het_model.r.0,
                         het_model.r.1,
                         het_model.r.2,
                         het_model.r.3,
                         het_model.r.4,
                         het_model.r.5,
                         het_model.r.6,
                         het_model.r.7,
                         het_model.r.8,
                         het_model.r.9))

#-------------------------------------------------------------------------------
# TRUE NULL ES = [0, 0, 0, 0, 0]
null_sim.r.0 <- multivar_sim(k_sims = 1000, r = 0,
                            control_mu = rep(0, 5),
                            exp_mu = rep(0, 5))
null_model.r.0 <- Zing(null_sim.r.0$zscores_A)
null_tru.model.r.0 <- Zing(null_sim.r.0$total_zscores)

null_sim.r.1 <- multivar_sim(k_sims = 1000, r = 0.1,
                            control_mu = rep(0, 5),
                            exp_mu = rep(0, 5))
null_model.r.1 <- Zing(null_sim.r.1$zscores_A)
null_tru.model.r.1 <- Zing(null_sim.r.1$total_zscores)

null_sim.r.2 <- multivar_sim(k_sims = 1000, r = 0.2,
                            control_mu = rep(0, 5),
                            exp_mu = rep(0, 5))
null_model.r.2 <- Zing(null_sim.r.2$zscores_A)
null_tru.model.r.2 <- Zing(null_sim.r.2$total_zscores)

null_sim.r.3 <- multivar_sim(k_sims = 1000, r = 0.3,
                            control_mu = rep(0, 5),
                            exp_mu = rep(0, 5))
null_model.r.3 <- Zing(null_sim.r.3$zscores_A)
null_tru.model.r.3 <- Zing(null_sim.r.3$total_zscores)

null_sim.r.4 <- multivar_sim(k_sims = 1000, r = 0.4,
                            control_mu = rep(0, 5),
                            exp_mu = rep(0, 5))
null_model.r.4 <- Zing(null_sim.r.4$zscores_A)
null_tru.model.r.4 <- Zing(null_sim.r.4$total_zscores)

null_sim.r.5 <- multivar_sim(k_sims = 1000, r = 0.5,
                            control_mu = rep(0, 5),
                            exp_mu = rep(0, 5))
null_model.r.5 <- Zing(null_sim.r.5$zscores_A)
null_tru.model.r.5 <- Zing(null_sim.r.5$total_zscores)

null_sim.r.6 <- multivar_sim(k_sims = 1000, r = 0.6,
                            control_mu = rep(0, 5),
                            exp_mu = rep(0, 5))
null_model.r.6 <- Zing(null_sim.r.6$zscores_A)
null_tru.model.r.6 <- Zing(null_sim.r.6$total_zscores)

null_sim.r.7 <- multivar_sim(k_sims = 1000, r = 0.7,
                            control_mu = rep(0, 5),
                            exp_mu = rep(0, 5))
null_model.r.7 <- Zing(null_sim.r.7$zscores_A)
null_tru.model.r.7 <- Zing(null_sim.r.7$total_zscores)

null_sim.r.8 <- multivar_sim(k_sims = 1000, r = 0.8,
                            control_mu = rep(0, 5),
                            exp_mu = rep(0, 5))
null_model.r.8 <- Zing(null_sim.r.8$zscores_A)
null_tru.model.r.8 <- Zing(null_sim.r.8$total_zscores)

null_sim.r.9 <- multivar_sim(k_sims = 1000, r = 0.9,
                            control_mu = rep(0, 5),
                            exp_mu = rep(0,5))
null_model.r.9 <- Zing(null_sim.r.9$zscores_A)
null_tru.model.r.9 <- Zing(null_sim.r.9$total_zscores)

null_EDR.r <- z3_EDR(list(null_model.r.0,
                         null_model.r.1,
                         null_model.r.2,
                         null_model.r.3,
                         null_model.r.4,
                         null_model.r.5,
                         null_model.r.6,
                         null_model.r.7,
                         null_model.r.8,
                         null_model.r.9))
null_ERR.r <- z3_ERR(list(null_model.r.0,
                         null_model.r.1,
                         null_model.r.2,
                         null_model.r.3,
                         null_model.r.4,
                         null_model.r.5,
                         null_model.r.6,
                         null_model.r.7,
                         null_model.r.8,
                         null_model.r.9))

#-------------------------------------------------------------------------------
# TWO EFFECT SIZES ES = [0.2, 0.8, 0.2, 0.8, 0.2]
duo_sim.r.0 <- multivar_sim(k_sims = 1000, r = 0,
                            control_mu = rep(0, 5),
                            exp_mu = c(0.2, 0.8, 0.2, 0.8, 0.2))
duo_model.r.0 <- Zing(duo_sim.r.0$zscores_A)
duo_tru.model.r.0 <- Zing(duo_sim.r.0$total_zscores)

duo_sim.r.1 <- multivar_sim(k_sims = 1000, r = 0.1,
                            control_mu = rep(0, 5),
                            exp_mu = c(0.2, 0.8, 0.2, 0.8, 0.2))
duo_model.r.1 <- Zing(duo_sim.r.1$zscores_A)
duo_tru.model.r.1 <- Zing(duo_sim.r.1$total_zscores)

duo_sim.r.2 <- multivar_sim(k_sims = 1000, r = 0.2,
                            control_mu = rep(0, 5),
                            exp_mu = c(0.2, 0.8, 0.2, 0.8, 0.2))
duo_model.r.2 <- Zing(duo_sim.r.2$zscores_A)
duo_tru.model.r.2 <- Zing(duo_sim.r.2$total_zscores)

duo_sim.r.3 <- multivar_sim(k_sims = 1000, r = 0.3,
                            control_mu = rep(0, 5),
                            exp_mu = c(0.2, 0.8, 0.2, 0.8, 0.2))
duo_model.r.3 <- Zing(duo_sim.r.3$zscores_A)
duo_tru.model.r.3 <- Zing(duo_sim.r.3$total_zscores)

duo_sim.r.4 <- multivar_sim(k_sims = 1000, r = 0.4,
                            control_mu = rep(0, 5),
                            exp_mu = c(0.2, 0.8, 0.2, 0.8, 0.2))
duo_model.r.4 <- Zing(duo_sim.r.4$zscores_A)
duo_tru.model.r.4 <- Zing(duo_sim.r.4$total_zscores)

duo_sim.r.5 <- multivar_sim(k_sims = 1000, r = 0.5,
                            control_mu = rep(0, 5),
                            exp_mu = c(0.2, 0.8, 0.2, 0.8, 0.2))
duo_model.r.5 <- Zing(duo_sim.r.5$zscores_A)
duo_tru.model.r.5 <- Zing(duo_sim.r.5$total_zscores)

duo_sim.r.6 <- multivar_sim(k_sims = 1000, r = 0.6,
                            control_mu = rep(0, 5),
                            exp_mu = c(0.2, 0.8, 0.2, 0.8, 0.2))
duo_model.r.6 <- Zing(duo_sim.r.6$zscores_A)
duo_tru.model.r.6 <- Zing(duo_sim.r.6$total_zscores)

duo_sim.r.7 <- multivar_sim(k_sims = 1000, r = 0.7,
                            control_mu = rep(0, 5),
                            exp_mu = c(0.2, 0.8, 0.2, 0.8, 0.2))
duo_model.r.7 <- Zing(duo_sim.r.7$zscores_A)
duo_tru.model.r.7 <- Zing(duo_sim.r.7$total_zscores)

duo_sim.r.8 <- multivar_sim(k_sims = 1000, r = 0.8,
                            control_mu = rep(0, 5),
                            exp_mu = c(0.2, 0.8, 0.2, 0.8, 0.2))
duo_model.r.8 <- Zing(duo_sim.r.8$zscores_A)
duo_tru.model.r.8 <- Zing(duo_sim.r.8$total_zscores)

duo_sim.r.9 <- multivar_sim(k_sims = 1000, r = 0.9,
                            control_mu = rep(0, 5),
                            exp_mu = c(0.2, 0.8, 0.2, 0.8, 0.2))
duo_model.r.9 <- Zing(duo_sim.r.9$zscores_A)
duo_tru.model.r.9 <- Zing(duo_sim.r.9$total_zscores)

duo_EDR.r <- z3_EDR(list(duo_model.r.0,
                         duo_model.r.1,
                         duo_model.r.2,
                         duo_model.r.3,
                         duo_model.r.4,
                         duo_model.r.5,
                         duo_model.r.6,
                         duo_model.r.7,
                         duo_model.r.8,
                         duo_model.r.9))
duo_ERR.r <- z3_ERR(list(duo_model.r.0,
                         duo_model.r.1,
                         duo_model.r.2,
                         duo_model.r.3,
                         duo_model.r.4,
                         duo_model.r.5,
                         duo_model.r.6,
                         duo_model.r.7,
                         duo_model.r.8,
                         duo_model.r.9))

#-------------------------------------------------------------------------------
# HOMOGENOUS MED ES = [0.6, 0.6, 0.6, 0.6, 0.6]
hom_sim.r.0 <- multivar_sim(k_sims = 1000, r = 0,
                            control_mu = rep(0, 5),
                            exp_mu = rep(0.6, 5))
hom_model.r.0 <- Zing(hom_sim.r.0$zscores_A)
hom_tru.model.r.0 <- Zing(hom_sim.r.0$total_zscores)

hom_sim.r.1 <- multivar_sim(k_sims = 1000, r = 0.1,
                            control_mu = rep(0, 5),
                            exp_mu = rep(0.6, 5))
hom_model.r.1 <- Zing(hom_sim.r.1$zscores_A)
hom_tru.model.r.1 <- Zing(hom_sim.r.1$total_zscores)

hom_sim.r.2 <- multivar_sim(k_sims = 1000, r = 0.2,
                            control_mu = rep(0, 5),
                            exp_mu = rep(0.6, 5))
hom_model.r.2 <- Zing(hom_sim.r.2$zscores_A)
hom_tru.model.r.2 <- Zing(hom_sim.r.2$total_zscores)

hom_sim.r.3 <- multivar_sim(k_sims = 1000, r = 0.3,
                            control_mu = rep(0, 5),
                            exp_mu = rep(0.6, 5))
hom_model.r.3 <- Zing(hom_sim.r.3$zscores_A)
hom_tru.model.r.3 <- Zing(hom_sim.r.3$total_zscores)

hom_sim.r.4 <- multivar_sim(k_sims = 1000, r = 0.4,
                            control_mu = rep(0, 5),
                            exp_mu = rep(0.6, 5))
hom_model.r.4 <- Zing(hom_sim.r.4$zscores_A)
hom_tru.model.r.4 <- Zing(hom_sim.r.4$total_zscores)

hom_sim.r.5 <- multivar_sim(k_sims = 1000, r = 0.5,
                            control_mu = rep(0, 5),
                            exp_mu = rep(0.6, 5))
hom_model.r.5 <- Zing(hom_sim.r.5$zscores_A)
hom_tru.model.r.5 <- Zing(hom_sim.r.5$total_zscores)

hom_sim.r.6 <- multivar_sim(k_sims = 1000, r = 0.6,
                            control_mu = rep(0, 5),
                            exp_mu = rep(0.6, 5))
hom_model.r.6 <- Zing(hom_sim.r.6$zscores_A)
hom_tru.model.r.6 <- Zing(hom_sim.r.6$total_zscores)

hom_sim.r.7 <- multivar_sim(k_sims = 1000, r = 0.7,
                            control_mu = rep(0, 5),
                            exp_mu = rep(0.6, 5))
hom_model.r.7 <- Zing(hom_sim.r.7$zscores_A)
hom_tru.model.r.7 <- Zing(hom_sim.r.7$total_zscores)

hom_sim.r.8 <- multivar_sim(k_sims = 1000, r = 0.8,
                            control_mu = rep(0, 5),
                            exp_mu = rep(0.6, 5))
hom_model.r.8 <- Zing(hom_sim.r.8$zscores_A)
hom_tru.model.r.8 <- Zing(hom_sim.r.8$total_zscores)

hom_sim.r.9 <- multivar_sim(k_sims = 1000, r = 0.9,
                            control_mu = rep(0, 5),
                            exp_mu = rep(0.6, 5))
hom_model.r.9 <- Zing(hom_sim.r.9$zscores_A)
hom_tru.model.r.9 <- Zing(hom_sim.r.9$total_zscores)

hom_EDR.r <- z3_EDR(list(hom_model.r.0,
                         hom_model.r.1,
                         hom_model.r.2,
                         hom_model.r.3,
                         hom_model.r.4,
                         hom_model.r.5,
                         hom_model.r.6,
                         hom_model.r.7,
                         hom_model.r.8,
                         hom_model.r.9))
hom_ERR.r <- z3_ERR(list(hom_model.r.0,
                         hom_model.r.1,
                         hom_model.r.2,
                         hom_model.r.3,
                         hom_model.r.4,
                         hom_model.r.5,
                         hom_model.r.6,
                         hom_model.r.7,
                         hom_model.r.8,
                         hom_model.r.9))

#-------------------------------------------------------------------------------