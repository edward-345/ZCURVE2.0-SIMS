

#### Co-Authord with Edward Lee 2025/08/18

###############################################################################
# Chapter 9 – P-Hacking Simulation 1: Pilot Dropping / Sample Patchwork
###############################################################################

### 0.0 – SETUP

# Disable scientific notation
options(scipen = 999)

# Set work directory
setwd("C:/Users/ulric/Dropbox/PHPCurves/DOCS/z-curve/Tutorial")

# Install required packages (run once only)
#install.packages(c("zcurve", "KernSmooth", "stringr", "parallel","pwr"))

# Load z-curve and P-curve functions from a local file
setwd("C:/Users/ulric/Dropbox/PHPCurves/DOCS/z-curve/Tutorial")
zcurve3 <- "Zing.25.07.11.test.R"
source(zcurve3)
pcurve <- "Pcurve.Function.R"
source(pcurve)


# Alternatively, load z-curve and p-curve functions directly from GitHub
zcurve3b <- "https://raw.githubusercontent.com/UlrichSchimmack/zcurve3.0/refs/heads/main/Zing.25.07.11.test.R"
source(zcurve3b)
pcurve <- "https://github.com/UlrichSchimmack/zcurve3.0/raw/refs/heads/main/Pcurve.Function.R"
source(pcurve)

#T-Tests for multiple covariates
multi_ttests <- function(control, exp) {
  
  ttest_res <- list()
  
  for (i in seq_len(ncol(control))) {
    nm <- paste0("dv", i, "_ttest")
    ttest_res[[nm]] <- t.test(control[, i], exp[, i], var.equal = TRUE)
  }
  ttest_res$avg_ttest <- t.test(rowMeans(control), rowMeans(exp_grp),
                                var.equal = TRUE)
  return(ttest_res)
  }

################################################################################
### 1.0 – ILLUSTRATION WITH EXAMPLES


################################################################################
### 2.0 – SIMULATION 

sim.runs <- 10                 # repeat each simulation 100 times
sim.k.sig <- 1000             # how many significant results?

sim.es.mean <- seq(0,1,.2)    # mean effect size, Cohen's d
sim.n.obs <- c(10, 20, 40, 60, 80, 100)   # sample size per cell
sim.r.var <- seq(0,1,.2)      # CORRELATION R BETWEEN COVARIATES
sim.n.vars <- c(2, 4, 6, 8, 10, 12) # NUMBER OF COVARIATES 
#sim.es.sd <- seq(0,.6,.2) # heterogeneity of es (tau) 

# MATRIX OF ALL POSSIBLE SIMULATION CONDITION COMBINATIONS 1296 X 5
sims <- c()
sim.i <- 0
for (i in 1:length(sim.es.mean) ) {
	for (j in 1:length(sim.r.var)) {
		for (k in 1:length(sim.n.obs) ) {
		  for (l in 1:length(sim.n.vars) ) {
			sim.i = sim.i + 1
			sims = rbind(sims,c(sim.i,
			                    sim.es.mean[i],
			                    sim.r.var[j],
			                    sim.n.obs[k],
			                    sim.n.vars[l]) )
}}}}
sims <- data.frame(sims)
colnames(sims) <- c("run","es.mean","r.var","n.obs", "n.vars")
dim(sims)

run.i <- 66       # used for testing without loop
b <- 1            # form row b
e <- 5            # to row e, used for testing
e <- nrow(sims);e # to last row, run all conditions

res = c()        # collect the results 

################################################################################
################################################################################
# CURRENTLY for(66 in 1:5) {}
for (run.i in b:e ) { # for loop to run the simulations
  
  #Run iteration tracker
  print(paste("Run: ",run.i))
  sims[run.i,]
  
  #EMPTY LIST OF TTEST SUMMARIES FROM SIMULATIONS
  results <- data.frame()
  #SIGNIFCANT RESULT COUNTER
  count.sig <- 0
  #NUMBER OF RUNS COUNTER
  count.run <- 0
  
  #WHILE LOOP GENERATING SAMPLES AND DATA
  while (count.sig < sim.k.sig) { 
    
    #LOOP ITERATION COUNTER
    count.run <- count.run + 1
    
    #GENERATE EXP AND CONTROL SAMPLES AND THEIR Y VALUES
    exp_group <- rnorm_multi(n = sims$n.obs[run.i],
                             vars = rep(sims$es.mean[run.i], sims$n.vars[run.i]),
                             sd = 1,
                             r = sim$r.var[run.i])
    control_group <- rnorm_multi(n = sims$n.obs[run.i],
                                 vars = rep(0, sims$n.vars[run.i]),
                                 sd = 1,
                                 r = 0)
    
    #T-TESTS OF EACH COVARIATE AND THE AVERAGE OF ALL COVARIATES
    ttest_res <- multi_ttests(exp_group, control_group)
    
    #EXTRACT PVALUES FOR COMPARISON
    pvalues <- sapply(ttest_res, function(x) x$p.value)
    
    #MINIMUM PVALUE OF THE T-TESTS
    min_pvalue <- min(pvalues)
    
    #EXTRACT VALUES FROM THE MOST SIG TEST FOR RESULTS DATA FRAME
    i_min <- which.min(vapply(ttest_res, function(m) m$p.value, numeric(1)))
    sig_ttest <- ttest_res[[i_min]]
    
    res.run <- data.frame()
    for (i in seq_len(length(ttest_res))) {
      tval <- ttest_res[[i]]$statistic
      dfs <- ttest_res[[i]]$parameter
      se <- ttest_res[[i]]$stderr
      es <- ttest_res[[i]]$estimate
      delta <- es[1]-es[2]
      pval <- ttest_res[[i]]$p.value
      
      res.run <- rbind(res.run,
                       c(count.run,
                         sims$es.mean[run.i],
                         delta,
                         se,
                         tval,
                         dfs,
                         pval))
    }
    
    results <- rbind(results, res.run)
    
    #IF THE LOWEST PVALUE IS SIGNIFCANT ADD TO SIGNIFCANT PVAL COUNTER
    if (min_pvalue < .05) {
      count.sig = count.sig + 1
      }
    
  } #END OF WHILE LOOP  
  
  results$N <- results$df+2
  results$se <- 2/sqrt(results$N)
  results$nct <- abs(results[,2]/results$se)  
  results$abs.t <- abs(results$t)
  results$z <- qnorm(pt(results$abs.t,results$N-2,log.p=TRUE),log.p=TRUE)
  results$pow.dir <- pt(qt(.975,results$df),results$df,results$nct,lower.tail=FALSE)
  results$pow.sign.err <- pt(-qt(.975,results$df),results$df,results$nct,lower.tail=TRUE)
  results$pow <- results$pow.dir + results$pow.sign.err
  
  summary(results)
  
  plot(results$abs.t,results$z)
  cor(results$abs.t,results$z)
  
  col1b = rgb(0, 100, 100, 255, max = 255,alpha = 50)
  col2b = rgb(100,0,100, max = 255, alpha = 50)
  
  cor(results$N,results$pow)
  cor(abs(results[,2]),results$pow) #  only if effect sizes are heterogeneous
  
  true.edr = mean(results$pow);true.edr
  true.err = sum(results$pow*results$pow.dir)/sum(results$pow);true.err
  
  tab = table(results[,7] < .05);tab/sum(tab) #ODR
  
  true.edr
  true.err
  
  sim.z.all = results$z
  sim.t.all = results$abs.t
  
  
  sim.min.p = apply(matrix(results$V7,7,),2,function(x) as.numeric(x == min(x)))
  results$min.p = c(sim.min.p)
  table(results$min.p)
  table(sim.z.all[results$min.p == TRUE] > 1.96)
  
  sim.z.max = sim.z.all[sim.min.p == TRUE]
  sim.t.max = sim.t.all[sim.min.p == TRUE]
  sim.df.max = results$df[sim.min.p == TRUE]
  
  
  source(zcurve3)
  Title = paste(round(true.edr*100),"  ",round(true.err*100))
  res.1 = Zing(sim.z.all);res.1
  res.1 = res.1$res[2:3];res.1
  
  source(zcurve3)
  Title = paste(round(true.edr*100),"  ",round(true.err*100))
  res.2 = Zing(sim.t.all);res.2
  res.2 = res.2$res[2:3];res.2
  
  source(zcurve3)
  Title = paste(round(true.edr*100),"  ",round(true.err*100))
  Est.Method = "DF"
  res.3 = Zing(sim.t.all,df=median(results$df));res.3
  res.3 = res.3$res[2:3];res.3
  
  source(zcurve3)
  Title = paste(round(true.edr*100),"  ",round(true.err*100))
  res.4 = Zing(sim.z.max);res.4
  res.4 = res.4$res[2:3];res.4
  
  source(zcurve3)
  Title = paste(round(true.edr*100),"  ",round(true.err*100))
  TEST4BIAS = TRUE
  just = 1
  Int.Beg = 2 + just
  res.5 = Zing(sim.z.max);res.5
  res.5 = res.5$res[2:3];res.5
  
  source(zcurve3)
  Title = paste(round(true.edr*100),"  ",round(true.err*100))
  TEST4BIAS = TRUE
  just = 1
  Int.Beg = 2 + just
  res.6 = Zing(sim.t.max,median(sim.df.max));res.6
  res.6 = res.6$res[2:3];res.6
  
  pcurve.input = paste0("t(",results$df,") = ",sim.t.all)
  res.pcurve.1 = pcurve_app(pcurve.input,SHOW.PLOT=TRUE)
  res.pcurve.1 
  if (res.pcurve.1[2] > res.pcurve.1[1]) res.pcurve.1[1] = res.pcurve.1[2]
  
  
  pcurve.input = paste0("t(",sim.df.max,") = ",sim.t.max)
  res.pcurve.2 = pcurve_app(pcurve.input,SHOW.PLOT=TRUE)
  res.pcurve.2 
  true.err
  if (res.pcurve.2[2] > res.pcurve.2[1]) res.pcurve.2[1] = res.pcurve.2[2]
  
  res.run = rbind(
    res.1,
    res.2,
    res.3,
    res.4,
    res.5,
    res.6,
    c(res.pcurve.1[1],NA),
    c(res.pcurve.2[1],NA)
  )
  
  rownames(res.run) = c("z.all","t.all","t.curve.all","z.sig","z.sig.sel",
                        "t.sig.sel","pcurve.all","pcurve.sig")
  
  round(rbind(c(true.err,true.edr),res.run),2)
  
  res.run = c(unlist(sims[run.i,]),true.err,res.run[,1],true.edr,res.run[,2])
  print(res.run)
  
  res = rbind(res,res.run)
  
  write.table(res,"sim.phack1.dat")    # write results each trial
  # can resume if stopped 
  # by loading the completed results

} # End of for loop

################################################################################
################################################################################

dim(res)  # 

round(res,2)

#write the completed results / overwrites the data from the for loop
#write.table(res,"sim.phack1.dat.dat")


########################################################
### GET STORED RESULTS
########################################################

# load the saved data
setwd("C:/Users/ulric/Dropbox/PHPCurves/DOCS/z-curve/Tutorial")
check = read.table("sim.phack1.dat",row.names=NULL)
check = data.frame(check)
check = check[,2:dim(check)[2]]
dim(check)
check[,1]

# careful, if you run this code any results in the res file will be replaced 
# with the saved data. 
# table(round(res[,c(1:14)],5) == round(check[,c(1:14)],5) )
# round(cbind(res[,7],check[,7]),5)

#res = check  # use only if sure to replace res; run analysis with "res" matrix 

dim(res)

#res = res[67:132,]

summary(res)

res[,1:4]

round(res[2,],2)

########################################################################
### Evaluation of ERR estimates: z-values versus t(28)->p->z values
########################################################################

# columns
# -  5: true ERR
# -  6: res.1
# -  7: res.2 
# -  8: res.3
# -  9: res.4 z.sig
# - 10: res.5 
# - 11: res.6 
# - 12: pcurve.1
# - 13: pcurve.2  t.sig

summary(res[,5:13])

cor(res[,5:13])[,1]

rmse.err.1 = sqrt(mean((res[,5]-res[,6])^2));rmse.err.1
rmse.err.2 = sqrt(mean((res[,5]-res[,7])^2));rmse.err.2
rmse.err.3 = sqrt(mean((res[,5]-res[,8])^2,na.rm=TRUE));rmse.err.3
rmse.err.4 = sqrt(mean((res[,5]-res[,9])^2,na.rm=TRUE));rmse.err.4
rmse.err.5 = sqrt(mean((res[,5]-res[,10])^2,na.rm=TRUE));rmse.err.5
rmse.err.6 = sqrt(mean((res[,5]-res[,11])^2,na.rm=TRUE));rmse.err.6
rmse.pcu.1 = sqrt(mean((res[,5]-res[,12])^2,na.rm=TRUE));rmse.pcu.1
rmse.pcu.2 = sqrt(mean((res[,5]-res[,13])^2,na.rm=TRUE));rmse.pcu.2

### all sig
rmse.err.1
rmse.err.2
rmse.err.3
rmse.pcu.1

### min.p
rmse.err.4
rmse.err.5
rmse.err.6
rmse.pcu.2


res = res[order(res[,5]),]

lwd.ci = .5

### all significant results

graphics.off()
plot(res[,5],res[,12],pch=15,cex=1,xlim=c(0,1),ylim=c(0,1),col="purple3",ylab="",xlab="")  # plot OF and ER estimates
par(new=TRUE)
plot(res[,5],res[,6],xlim=c(0,1),ylim=c(0,1),pch=16,cex=1,col="forestgreen",ylab="Estimated ERR",
	xlab = "Simulated True ERR")  # plot OF and ER estimates
abline(a = 0, b = 1,lty=2)
abline(v = .5,lty=2)
abline(h = .5,lty=2)
legend(.1,.9,legend=c("z-curve", "p-curve"),lwd=2,col=c("forestgreen","purple3"))

### lowest p-value per study 

graphics.off()
plot(res[,5],res[,13],pch=15,cex=1,xlim=c(0,1),ylim=c(0,1),col="purple3",ylab="",xlab="")  # plot OF and ER estimates
par(new=TRUE)
plot(res[,5],res[,9],xlim=c(0,1),ylim=c(0,1),pch=16,cex=1,col="forestgreen",ylab="Estimated ERR",
	xlab = "Simulated True ERR")  # plot OF and ER estimates
abline(a = 0, b = 1,lty=2)
abline(v = .5,lty=2)
abline(h = .5,lty=2)
legend(.1,.9,legend=c("z-curve", "p-curve"),lwd=2,col=c("forestgreen","purple3"))

dir.bias.z = cbind(res[,9]-res[,5],res[,c(1:5,9)])
round(dir.bias.z[abs(dir.bias.z[,1]) > .15,],2)
summary(lm(dir.bias.z[,1] ~ res[,2] + res[,3] + res[,4]))
tapply(dir.bias.z[,1],res[,3],mean)

dir.bias.p = cbind(res[,13]-res[,5],res[,c(1:5,13)])
dir.bias.p[dir.bias.p < -.50] = NA
round(dir.bias.p[abs(dir.bias.p[,1]) > .15,],2)
summary(lm(dir.bias.p[,1] ~ res[,2] + res[,3] + res[,4] + res[,5]))

cor.lw = function(x) cor(x,use="complete.obs")
cor.lw(cbind(dir.bias.z[,1],dir.bias.p[,1]))

###

plot(res[,5],res[,9],pch=15,cex=1,xlim=c(0,1),ylim=c(0,1),col="purple3",ylab="",xlab="")  # plot OF and ER estimates
par(new=TRUE)
plot(res[,5],res[,6],xlim=c(0,1),ylim=c(0,1),pch=16,cex=1,col="forestgreen",ylab="Estimated ERR",
	xlab = "Simulated True ERR")  # plot OF and ER estimates
abline(a = 0, b = 1,lty=2)
legend(.1,.9,legend=c("z-values", "t->power->z"),lwd=2,col=c("forestgreen","purple3"))

###

plot(res[,5],res[,10],pch=15,cex=1,xlim=c(0,1),ylim=c(0,1),col="purple3",ylab="",xlab="")  # plot OF and ER estimates
par(new=TRUE)
plot(res[,5],res[,6],xlim=c(0,1),ylim=c(0,1),pch=16,cex=1,col="forestgreen",ylab="Estimated ERR",
	xlab = "Simulated True ERR")  # plot OF and ER estimates
abline(a = 0, b = 1,lty=2)
legend(.1,.9,legend=c("z-values", "t-curve"),lwd=2,col=c("forestgreen","purple3"))


########################################################################
### Evaluation of EDR estimates: z-values versus t(28)->p->z values
########################################################################

# columns
# - 14: true EDR
# - 15: res.1
# - 16: res.2 
# - 17: res.3
# - 18: res.4
# - 19: res.5
# - 20: res.6

rmse.edr.1 = sqrt(mean((res[,14]-res[,15])^2));rmse.edr.1
rmse.edr.2 = sqrt(mean((res[,14]-res[,16])^2));rmse.edr.2
rmse.edr.3 = sqrt(mean((res[,14]-res[,17])^2));rmse.edr.3
rmse.edr.4 = sqrt(mean((res[,14]-res[,18])^2));rmse.edr.4
rmse.edr.5 = sqrt(mean((res[,14]-res[,19])^2));rmse.edr.5
rmse.edr.6 = sqrt(mean((res[,14]-res[,20])^2));rmse.edr.6

rmse.edr.1
rmse.edr.2
rmse.edr.3
rmse.edr.4
rmse.edr.5
rmse.edr.6

###

res = res[order(res[,11]),]

lwd.ci = .5

### all significant results

graphics.off()
plot(res[,14],res[,18],pch=15,cex=1,xlim=c(0,1),ylim=c(0,1),col="purple3",ylab="",xlab="")  # plot OF and ER estimates
par(new=TRUE)
plot(res[,14],res[,15],xlim=c(0,1),ylim=c(0,1),pch=16,cex=1,col="forestgreen",
	ylab="Estimated EDR",
	xlab = "Simulated True EDR")  
abline(a = 0, b = 1,lty=2)
abline(v = .5,lty=2)
abline(h = .5,lty=2)
legend(.1,.9,legend=c("all sig", "max sig"),lwd=2,col=c("forestgreen","purple3"))

