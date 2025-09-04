count.sig <- 0
k.sig <- 10

all.results <- data.frame(tval = numeric(),
                           df = numeric(),
                           N = numeric(),
                           es = numeric(),
                           delta = numeric(),
                           pval = numeric())
selected.results <- data.frame(tval = numeric(),
                                df = numeric(),
                                N = numeric(),
                                es = numeric(),
                                delta = numeric(),
                                pval = numeric())
while (count.sig < k.sig) {
  
  #GENERATE EXP AND CONTROL SAMPLES AND THEIR Y VALUES
  #NOTE THIS USES faux::rnorm_multi()
  exp_group <- rnorm_multi(n = 20,
                           mu = c(1, 2, 3),
                           sd = 1,
                           r = .5)
  control_group <- rnorm_multi(n = 20,
                               mu = c(0,1,0),
                               sd = 1,
                               r = .5)
  
  #T-TESTS OF EACH COVARIATE AND THE AVERAGE OF ALL COVARIATES
  ttest_res <- multi_ttests(exp_group, control_group) #note this is a list
  
  #EXTRACTING VALUES FROM EVERY T-TEST IN THE CURRENT ITERATION
  all.res <- data.frame()
  for (i in seq_along(ttest_res)) {
    tval <- ttest_res[[i]]$statistic
    df <- ttest_res[[i]]$parameter
    N <- ttest_res[[i]]$parameter + 2
    es <- unname(ttest_res[[i]]$estimate[1]) 
    delta <- 0
    pval <- ttest_res[[i]]$p.value
    values <- c(tval, df, N, es, delta, pval)
    
    all.res <- rbind(all.res, values)
  }
  
  colnames(all.res) <- c("tval", "df", "N", "es", "delta", "pval")
  
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
  m.delta <- 0
  m.pval <- sig_ttest$p.value
  m.values <- c(m.tval, m.df, m.N, m.es, m.delta, m.pval)
  
  #ADDING RESULTS TO SELECTED RESULTS
  selected.results <- rbind(selected.results, m.values)
  colnames(selected.results) <- c("tval", "df", "N", "es", "delta", "pval")
  
  if (min.pvalue < 0.5) {
    count.sig <- count.sig + 1
  }
  
}
