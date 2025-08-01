source("helper_functions.R")
set.seed(666)

#ZCURVE3.0 Imports
zcurve3 <- "https://raw.githubusercontent.com/UlrichSchimmack/zcurve3.0/refs/heads/main/Zing.25.07.11.test.R"
source(zcurve3)
#-------------------------------------------------------------------------------

random_mean <- function(x) {
  if (x == 1) {exp_mean = 0}
  else if (x == 2) {exp_mean = 0.2}
  else if (x == 3) {exp_mean = 0.4}
  else if (x == 4) {exp_mean = 0.6}
  else {exp_mean = 0.8}
  
  return(exp_mean)
}

# Function for returning heterogeneity estimates from TEST4HETEROGENEITY
het_estimates <- function(x) {
  het_test_results <- list(het_25 = x[2, 4],
                           het_median = x[7, 4],
                           het_95 = x[12, 4])
  return(het_test_results)
}

#-------------------------------------------------------------------------------

het_sim <- function(k_sims, n = 20) {
  z_scores <- numeric(k_sims)
  j <- 1
  pvals <- numeric(k_sims)
  for (i in 1:k_sims) {
    exp_mean <- random_mean(round(runif(1, min = 1, max = 4), 0))
    
    control_group <- rnorm(n, mean = 0, sd = 1)
    exp_group <- rnorm(n, mean = exp_mean, sd = 1)
    
    result <- t.test(control_group, exp_group, var.equal = TRUE)
    p_value <- result$p.value
    z_value <- pval_converter(p_value)
    
    if (abs(z_value) >= qnorm(0.975)) {
      pvals[i] <- p_value
      z_scores[j] <- z_value
      j <- j + 1
    } else {
      pvals[i] <- p_value
    }
  }
  
  z_scores <- z_scores[1:(j - 1)]
  fit <- zcurve(z_scores, control = list(parallel = TRUE))
  
  sim_list <- list(fit = fit,
                   z_scores = z_scores,
                   pvals = pvals)
  return(sim_list)
}

het_trial <- het_sim(5000)
summary(het_trial$fit)
het_trial.plot <- plot(het_trial$fit,
                      CI = TRUE, annotation = TRUE,
                      main = "Heterogenous under Null Hypothesis")

het_trial.pvals <- zcurve(p = het_trial$pvals,
                         control = list(parallel = TRUE))
het_trial.pvals.plot <- plot(het_trial.pvals, ymax = 10,
                            CI = TRUE, annotation = TRUE,
                            main = "Heterogenous P-vals under Null Hypothesis")
# ZCURVE 3.0
ymax <- 0.8
TEST4HETEROGENEITY <- 500
TEST4BIAS <- TRUE
het_trial.3.0 <- Zing(pval_converter(het_trial$pvals))
het_estimates(het_trial.3.0$fit.comp)

#-------------------------------------------------------------------------------
#SITUATION A: Two DVs for each observation
A_het <- function(k_sims, n = 20, r = 0.5, control_mu = c(0,0), sd = c(1,1)) {
  zscores_A <- numeric(k_sims)
  A <- 1
  
  #Vector of all p-values generated
  pvals <- numeric(k_sims)
  
  for (i in 1:k_sims) {
    exp_mu <- c(random_mean(round(runif(1, min = 1, max = 4), 0)),
                random_mean(round(runif(1, min = 1, max = 4), 0)))
    
    control_A <- rnorm_multi(
      n, vars = 2, control_mu, sd = c(1,1), r,
      varnames = c("Control1","Control2"))
    exp_A <- rnorm_multi(
      n, vars = 2, exp_mu, sd = c(1,1), r,
      varnames = c("Dependent1","Dependent2"))
    
    #Helper function conducts 3 t-tests, one on each of two dependent variables 
    #and a third on the average of these two variables
    pval <- sitA_ttests(control_A$Control1, control_A$Control2,
                            exp_A$Dependent1, exp_A$Dependent2)
    min_pval <- min(pval)
    
    #Add p-value to pvalues_scenarioA regardless of significance
    pvals[i] <- min_pval
    
    #If the smallest p-value of the three t-tests is significant at .05, convert 
    #to z-score and add to zscores_A vector
    if (min_pval <= 0.05) {
      zvalue_A <- pval_converter(min_pval)
      zscores_A[A] <- zvalue_A
      A <- A + 1
    }
  }
  
  #Trim unused cells
  zscores_A <- zscores_A[1:(A - 1)]
  
  fit <- zcurve(zscores_A, control = list(parallel = TRUE))
  
  A_list <- list(fit = fit,
                 zscores_A = zscores_A,
                 pvals = pvals)
  
  return(A_list)
}


A_hetero <- A_het(5000)
summary(A_hetero$fit)
A_hetero.plot <- plot(A_hetero $fit,
                       CI = TRUE, annotation = TRUE,
                       main = "Heterogenous under Null Hypothesis (A)")

A_hetero.pvals <- zcurve(p = A_hetero$pvals,
                          control = list(parallel = TRUE))
A_hetero.pvals.plot <- plot(
  A_hetero.pvals, ymax = 10,
  CI = TRUE, annotation = TRUE,
  main = "Heterogenous P-vals under Null Hypothesis (A)")

# ZCURVE 3.0
source(zcurve3)
ymax <- 0.8
TEST4HETEROGENEITY <- 500
TEST4BIAS <- TRUE
A_hetero.3.0 <- Zing(pval_converter(A_hetero$pvals))
het_estimates(A_hetero.3.0$fit.comp)

#-------------------------------------------------------------------------------

#SITUATION B: Optional Stopping
B_het <- function(k_sims, n = 20, extra_n = 10,
                  control_mu = 0, sd = 1) {
  zscores_B <- numeric(k_sims)
  B <- 1 
  
  pvals <- numeric(k_sims)
  
  for (i in 1:k_sims) {
    exp_mu <- random_mean(round(runif(1, min = 1, max = 4), 0))
    
    control_B <- rnorm(n, control_mu, sd)
    exp_B <- rnorm(n, exp_mu, sd)
    result_B <- t.test(control_B, exp_B, var.equal = TRUE)
    pvalue_B <- result_B$p.value
    
    if (pvalue_B <= 0.05) {
      #If the result is significant, the researcher stops collecting data and 
      #reports the result
      pvals[i] <- pvalue_B
      zvalue_B <- pval_converter(pvalue_B)
      zscores_B[B] <- zvalue_B
      B <- B + 1
    } else {
      #If the result is non significant, the researcher collects 10 additional 
      #observations per condition
      extracontrol_B <- rnorm(extra_n, control_mu, sd)
      extraexp_B <- rnorm(extra_n, exp_mu, sd)
      extraresult_B <- t.test(c(control_B, extracontrol_B),
                              c(exp_B, extraexp_B), var.equal = TRUE)
      extrapvalue_B <- extraresult_B$p.value
      
      if (extrapvalue_B <= 0.05) {
        #then again tests for significance
        pvals[i] <- extrapvalue_B
        extrazvalue_B <- pval_converter(extrapvalue_B)
        zscores_B[B] <- extrazvalue_B
        B <- B+1
      } else {
        pvals[i] <- extrapvalue_B
      }
    }
  }
  
  zscores_B <- zscores_B[1:(B - 1)]
  
  fit_B <- zcurve(zscores_B, control = list(parallel = TRUE))
  
  B_list <- list(fit_B = fit_B,
                 zscores_B = zscores_B,
                 pvals = pvals)
  
  return(B_list)
}

B_hetero <- B_het(5000)
summary(B_hetero$fit)
B_hetero.plot <- plot(B_hetero $fit,
                      CI = TRUE, annotation = TRUE,
                      main = "Heterogenous under Null Hypothesis (B)")

B_hetero.pvals <- zcurve(p = B_hetero$pvals,
                         control = list(parallel = TRUE))
B_hetero.pvals.plot <- plot(
  B_hetero.pvals, ymax = 10,
  CI = TRUE, annotation = TRUE,
  main = "Heterogenous P-vals under Null Hypothesis (B)")

# ZCURVE 3.0
source(zcurve3)
ymax <- 0.8
TEST4HETEROGENEITY <- 100
TEST4BIAS <- TRUE
B_hetero.3.0 <- Zing(pval_converter(B_hetero$pvals))
het_estimates(B_hetero.3.0$fit.comp)

#-------------------------------------------------------------------------------
#SITUATION C: Main effect or interaction term ANCOVAs
C_het <- function(k_sims, n = 20, control_mu = 0, sd = 1) {
  zscores_C <- numeric(k_sims)
  C <- 1 
  
  pvals <- numeric(k_sims)
  
  for (i in 1:k_sims) {
    exp_mu <- random_mean(round(runif(1, min = 1, max = 4), 0))
    
    groups <- sample(
      rep(c("control", "experimental"), each = n))
    dv <- rnorm(n*2, 
                mean = ifelse(groups=="control", control_mu, exp_mu), sd)
    #Each observation was assigned a 50% probability of being female
    gender <- rbinom(n*2, size = 1, p = 0.5)
    gender <- as.factor(
      ifelse(gender == 1, "female", "male"))
    data_C <- data.frame(dv, gender, groups)
    
    #Results for Situation C were obtained by conducting a t-test...
    ttest_result <- t.test(dv ~ groups, data = data_C, var.equal = TRUE)
    ttest_pvalue <- ttest_result$p.value
    
    #...an analysis of covariance with a gender main effect..
    main_model <- lm(dv ~ groups + gender, data = data_C)
    main_pvalue <- summary(
      main_model)$coefficients["groupsexperimental", "Pr(>|t|)"]
    
    #...and an analysis of covariance with a gender interaction.
    int_model <- lm(dv ~ groups*gender, data = data_C)
    coefs <- coef(summary(int_model))
    int_term_name <- grep("group.*:gender.*", rownames(coefs), value = TRUE)
    
    #Edge case guard of sample consisting entirely of one gender for int_pvalue
    if (length(int_term_name) == 1) {
      int_pvalue <- coefs[int_term_name, "Pr(>|t|)"]
    } else {
      int_pvalue <- 1
    }
    
    pvalues_C <- c(ttest_pvalue, main_pvalue)
    
    min_pvalueC <- min(pvalues_C)
    
    #We report a significant effect if the effect of condition was significant in 
    #any of these analyses or if the Gender×Condition interaction was significant.
    if (int_pvalue <= 0.05) {
      pvals[i] <- int_pvalue
      zvalue_C <- pval_converter(int_pvalue)
      zscores_C[C] <- zvalue_C
      C <- C + 1
    } else {
      if (min_pvalueC <= 0.05) {
        pvals[i] <- min_pvalueC
        zvalue_C <- pval_converter(min_pvalueC)
        zscores_C[C] <- zvalue_C
        C <- C + 1
      } else {
        pvals[i] <- min_pvalueC
      }
    }
    
  }
  
  zscores_C <- zscores_C[1:(C - 1)]
  
  fit <- zcurve(zscores_C, control = list(parallel = TRUE))
  
  C_list <- list(fit = fit,
                 zscores_C = zscores_C,
                 pvals = pvals)
  
  return(C_list)
}

C_hetero <- C_het(5000)
summary(C_hetero$fit)
C_hetero.plot <- plot(C_hetero $fit,
                      CI = TRUE, annotation = TRUE,
                      main = "Heterogenous under Null Hypothesis (C)")

C_hetero.pvals <- zcurve(p = C_hetero$pvals,
                         control = list(parallel = TRUE))
C_hetero.pvals.plot <- plot(
  C_hetero.pvals, ymax = 10,
  CI = TRUE, annotation = TRUE,
  main = "Heterogenous P-vals under Null Hypothesis (C)")

# ZCURVE 3.0
source(zcurve3)
ymax <- 0.8
TEST4HETEROGENEITY <- 100
TEST4BIAS <- TRUE
C_hetero.3.0 <- Zing(pval_converter(C_hetero$pvals))
het_estimates(C_hetero.3.0$fit.comp)

#-------------------------------------------------------------------------------
