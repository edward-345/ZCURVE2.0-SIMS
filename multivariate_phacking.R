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
    ttest_res[[nm]] <- t.test(exp[, i], control[, i], var.equal = TRUE)
  }
  ttest_res$avg_ttest <- t.test(rowMeans(exp), rowMeans(control),
                                var.equal = TRUE)
  return(ttest_res)
}

#### SIMULATION CONDITIONS -----------------------------------------------------

runs <- 100               # repeat each simulation 100 times
k.sig <- 10000             # how many significant results?

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

# BUILDING RES COLUMN NAMES
res_names <- c("all",
               "selected",
               "bias.all",
               "bias.selected")


RES_COLS <- c("run",
              "es.mean",
              "r.vars",
              "n.obs",
              "n.vars",
              "true.err",
              paste0("err.", res_names),
              "true.edr",
              paste0("edr.", res_names)
              )

res <- data.frame()


e <- nrow(simul.cons)
# for testing, assigns e as 1 so loops only once
# e <- 1

#### SIMULATION/DATA GENERATION ------------------------------------------------

for (run.i in 1:e) {
  
  #Run iteration tracker
  print(paste("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RUN: ",run.i))
  sims[run.i,]
  
  
  
  
  
  
  
  
  
}

