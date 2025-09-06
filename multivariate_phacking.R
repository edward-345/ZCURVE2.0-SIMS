#### Co-Authord with Edward Lee 2025/08/18
###############################################################################
# MULTI-VARIATE P-HACKING SIMULATION
################################################################################
#### SET-UP --------------------------------------------------------------------
options(scipen = 999)
set.seed(67)

# Load z-curve and p-curve functions directly from GitHub
zcurve3 <- "https://raw.githubusercontent.com/UlrichSchimmack/zcurve3.0/refs/heads/main/Zing.25.07.11.test.R"
source(zcurve3)
pcurve <- "https://github.com/UlrichSchimmack/zcurve3.0/raw/refs/heads/main/Pcurve.Function.R"
source(pcurve)

# Install required packages (run once only)
library(zcurve)
library(KernSmooth)
library(stringr)
library(parallel)
library(pwr)
library(tidyverse)
library(faux)
library(caret)

#T-Tests for multiple covariates
multi_ttests <- function(exp, control) {
  
  ttest_res <- list()
  
  for (i in seq_len(ncol(control))) {
    nm <- paste0("dv", i, "_ttest")
    
    ttest_res[[nm]] <- t.test(exp[, i],
                              control[, i],
                              var.equal = TRUE)
  }
  ttest_res$avg_ttest <- t.test(rowMeans(exp),
                                rowMeans(control),
                                var.equal = TRUE)
  return(ttest_res)
}

#### SIMULATION CONDITIONS -----------------------------------------------------
# EMPTY DATAFRAME TO STORE MODEL SUMMARY VALUES
model.values <- data.frame(tot.mdl.edr = numeric(),
                           tot.mdl.err = numeric(), 
                           
                           sel.mdl.edr = numeric(), 
                           sel.mdl.err = numeric(), 
                           
                           tot.mdl_bias.edr = numeric(),
                           tot.mdl_bias.err = numeric(),
                           
                           tot.mdl_bias.OJS = numeric(),
                           tot.mdl_bias.EJS = numeric(),
                           tot.mdl_bias.P = numeric(),
                           
                           sel.mdl_bias.edr = numeric(), 
                           sel.mdl_bias.err = numeric(),
                           
                           sel.mdl_bias.OJS = numeric(),
                           sel.mdl_bias.EJS = numeric(),
                           sel.mdl_bias.P = numeric())

k.sig <- 1000             # how many significant results?

es.mean <- seq(0, 1, .25)    # mean effect size, Cohen's d
n.obs <- c(15, 20, 30)   # sample size per cell
r.vars <- c(0, .25, .50, .75)      # CORRELATION R BETWEEN COVARIATES
n.vars <- c(2, 5, 10) # NUMBER OF COVARIATES 

# MATRIX OF ALL POSSIBLE SIMULATION CONDITION COMBINATIONS 1296 X 5
simul.cons <- c()
simul.cons.i <- 0
for (i in 1:length(es.mean) ) {
  for (j in 1:length(r.vars)) {
    for (k in 1:length(n.obs) ) {
      for (l in 1:length(n.vars) ) {
        simul.cons.i = simul.cons.i + 1
        simul.cons = rbind(simul.cons,c(simul.cons.i,
                            es.mean[i],
                            r.vars[j],
                            n.obs[k],
                            n.vars[l]) )
      }}}}
simul.cons <- data.frame(simul.cons)
colnames(simul.cons) <- c("run","es.mean","r.vars","n.obs", "n.vars")
dim(simul.cons)

#POWER OF EACH SIMULATION CONDITION COMBO
simul.cons$upper.pwr <- pt(qt(.975, (2*simul.cons$n.obs - 2)),
                           df = (2*simul.cons$n.obs - 2),
                           ncp = simul.cons$es.mean*sqrt(simul.cons$n.obs/2),
                           lower.tail = FALSE)
simul.cons$lower.pwr <- pt(-qt(.975, (2*simul.cons$n.obs - 2)),
                           df = (2*simul.cons$n.obs - 2),
                           ncp = simul.cons$es.mean*sqrt(simul.cons$n.obs/2),
                           lower.tail = TRUE)
simul.cons$obs.power <- simul.cons$upper.pwr + simul.cons$lower.pwr

e <- nrow(simul.cons)

#### SIMULATION/DATA GENERATION ------------------------------------------------
#### WE ARE TESTING FOR ANY SIGNIFICANT DIFFERENCE, USING TWO-SIDED T-TESTS

for (run.i in 1:e) {
  
  #Run iteration tracker
  print(paste("------------------------------------------------ RUN: ",run.i))
  simul.cons[run.i,]
  
  #SIGNIFCANT RESULT COUNTER
  count.sig <- 0
  #NUMBER OF RUNS COUNTER
  count.run <- 0
  
  # EMPTY DATAFRAME TO STORE ALL RESULTS IN EACH ITERATION 
  all.results <- data.frame(tval = numeric(),
                            df = numeric(),
                            N = numeric(),
                            es = numeric(),
                            delta = numeric(),
                            pval = numeric(),
                            z.score = numeric())
  
  # EMPTY DATAFRAME TO STORE P-HACKED SELECTED RESULTS IN EACH ITERATION 
  selected.results <- data.frame(tval = numeric(),
                                 df = numeric(),
                                 N = numeric(),
                                 es = numeric(),
                                 delta = numeric(),
                                 pval = numeric(),
                                 z.score = numeric())
  
  # GENERATING DATA -------------------------------------------------
  while (count.sig < k.sig) {
    
    #GENERATE EXP AND CONTROL SAMPLES AND THEIR Y VALUES
    #NOTE THIS USES faux::rnorm_multi()
    n_vars <- simul.cons$n.vars[run.i]
    exp_group <- rnorm_multi(n = simul.cons$n.obs[run.i],
                             mu = rep(simul.cons$es.mean[run.i], n_vars),
                             sd = 1,
                             r = simul.cons$r.vars[run.i])
    control_group <- rnorm_multi(n = simul.cons$n.obs[run.i],
                                 mu = rep(0, n_vars),
                                 sd = 1,
                                 r = simul.cons$r.vars[run.i])
    
    #T-TESTS OF EACH COVARIATE AND THE AVERAGE OF ALL COVARIATES
    ttest_res <- multi_ttests(exp_group, control_group) #note this is a list
    
    #EXTRACTING VALUES FROM EVERY T-TEST IN THE CURRENT ITERATION
    all.res <- data.frame()
    for (i in seq_along(ttest_res)) {
      tval <- unname(ttest_res[[i]]$statistic)
      df <- unname(ttest_res[[i]]$parameter)
      N <- unname(ttest_res[[i]]$parameter) + 2
      es <- unname(ttest_res[[i]]$estimate[1]) #technically isnt used
      delta <- simul.cons$es.mean[run.i]      #technically isnt used
      pval <- unname(ttest_res[[i]]$p.value)
      z.score <- abs(qnorm(ttest_res[[i]]$p.value/2, lower.tail = FALSE))
      values <- c(tval, df, N, es, delta, pval, z.score)
      
      all.res <- rbind(all.res, values)
    }
    
    colnames(all.res) <- c("tval", "df", "N", "es", "delta", "pval", "z.score")

    #ADDING RESULTS TO TOTAL RESULTS
    all.results <- rbind(all.results, all.res)
    
    
    #FINDING THE INDEX OF THE T-TEST WITH LOWEST P-VALUE
    i_min <- which.min(vapply(ttest_res, function(m) m$p.value, numeric(1)))
    
    #T-TEST WITH THE LOWEST PVALUE
    sig_ttest <- ttest_res[[i_min]]
    min.pvalue <- sig_ttest$p.value
    
    #ADD ONLY MOST SIGNIFCANT TEST RESULTS TO selected.results
    m.tval <- sig_ttest$statistic
    m.df <- sig_ttest$parameter
    m.N <- sig_ttest$parameter + 2
    m.es <- unname(sig_ttest$estimate[1]) 
    m.delta <- simul.cons$es.mean[run.i]
    m.pval <- sig_ttest$p.value
    m.zscore <- abs(qnorm(sig_ttest$p.value/2, lower.tail = FALSE))
    m.values <- c(m.tval, m.df, m.N, m.es, m.delta, m.pval, m.zscore)
    
    #ADDING RESULTS TO SELECTED RESULTS
    selected.results <- rbind(selected.results, m.values)
    colnames(selected.results) <- c("tval", "df", "N", "es",
                                    "delta", "pval", "z.score")
    
    #COUNTER FOR SPECIFIED NUMBER OF SIGNIFCANT P-VALUES
    if (min.pvalue < 0.05) {
      count.sig <- count.sig + 1
    }
  } # END OF WHILE LOOP/DATA GENERATOR -----------------------------------
  
  # DATA FRAMES WE WILL WORK WITH 
  head(all.results)
  head(selected.results)
  
  
  #FITTING ZCURVE MODELS----------------------------------
  #EMPTY DATA FRAME FOR STORING MODEL VALUES OF CURRENT ITERATION (SINGLE ROW)
  model.val <- data.frame()
  
  # TOTAL T-TESTS
  source(zcurve3)
  tot.mdl <- Zing(all.results$z.score)
  tot.mdl.edr <- unname(tot.mdl$res[3])
  tot.mdl.err <- unname(tot.mdl$res[2])
  
  # SELECTED T-TESTS
  source(zcurve3)
  sel.mdl <- Zing(selected.results$z.score)
  sel.mdl.edr <- unname(sel.mdl$res[3]) 
  sel.mdl.err <- unname(sel.mdl$res[2])
  
  
  # TOTAL T-TESTS WITH BIAS TEST
  source(zcurve3)
  TEST4BIAS <- TRUE
  tot.mdl_bias <- Zing(all.results$z.score)
  tot.mdl_bias.edr <- unname(tot.mdl_bias$res[3])
  tot.mdl_bias.err <- unname(tot.mdl_bias$res[2])
  
  tot.mdl_bias.OJS <- unname(tot.mdl_bias$bias[1])
  tot.mdl_bias.EJS <- unname(tot.mdl_bias$bias[2])
  tot.mdl_bias.P <- unname(tot.mdl_bias$bias[3])
  
  # SELECTED T-TESTS WITH BIAS TEST
  source(zcurve3)
  TEST4BIAS <- TRUE
  sel.mdl_bias <- Zing(selected.results$z.score)
  sel.mdl_bias.edr <- unname(sel.mdl_bias$res[3])
  sel.mdl_bias.err <- unname(sel.mdl_bias$res[2])
  
  sel.mdl_bias.OJS <- unname(sel.mdl_bias$bias[1])
  sel.mdl_bias.EJS <- unname(sel.mdl_bias$bias[2])
  sel.mdl_bias.P <- unname(sel.mdl_bias$bias[3])
  
  # TRUE EDR AND ERR !! NOTE EDR = ERR = POWER SINCE EACH ITERATION IS HOMOGENOUS
  upper_pwr <- pt(qt(.975, (2*simul.cons$n.obs[run.i] - 2)),
                             df = (2*simul.cons$n.obs[run.i] - 2),
                             ncp = simul.cons$es.mean[run.i]*sqrt(simul.cons$n.obs[run.i]/2),
                             lower.tail = FALSE)
  lower_pwr <- pt(-qt(.975, (2*simul.cons$n.obs[run.i] - 2)),
                             df = (2*simul.cons$n.obs[run.i] - 2),
                             ncp = simul.cons$es.mean[run.i]*sqrt(simul.cons$n.obs[run.i]/2),
                             lower.tail = TRUE)
  true.edr <- upper_pwr + lower_pwr
  
  true.err <- ((nrow(all.results)*true.edr)*upper_pwr)/(nrow(all.results)*true.edr)
  
  
  mdl.vals <- c(true.edr,
                true.err,
                tot.mdl.edr,
                tot.mdl.err,
                sel.mdl.edr,
                sel.mdl.err,
                tot.mdl_bias.edr,
                tot.mdl_bias.err,
                tot.mdl_bias.OJS,
                tot.mdl_bias.EJS,
                tot.mdl_bias.P,
                sel.mdl_bias.edr, 
                sel.mdl_bias.err,
                sel.mdl_bias.OJS,
                sel.mdl_bias.EJS,
                sel.mdl_bias.P)
  
  #cant add to simul.cons yet- will need to after so that col lens equal
  model.val <- rbind(model.val, mdl.vals)
  colnames(model.val) <- c("true.edr",
                           "true.err",
                           "tot.mdl.edr",
                           "tot.mdl.err",
                           "sel.mdl.edr",
                           "sel.mdl.err",
                           "tot.mdl_bias.edr",
                           "tot.mdl_bias.err",
                           "tot.mdl_bias.OJS",
                           "tot.mdl_bias.EJS",
                           "tot.mdl_bias.P",
                           "sel.mdl_bias.edr", 
                           "sel.mdl_bias.err",
                           "sel.mdl_bias.OJS",
                           "sel.mdl_bias.EJS",
                           "sel.mdl_bias.P")
  
  model.values <- rbind(model.values, model.val)
  
} # END OF LOOP- WE NOW HAVE simul.cons and model.values

simul.data <- cbind(simul.cons, model.values)




