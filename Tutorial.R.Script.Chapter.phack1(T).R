

#### Co-Authord with Edward Lee 2025/08/25

###############################################################################
# Chapter 9 – P-Hacking Simulation 1: Sample Patchwork
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
zcurve3 <- "https://raw.githubusercontent.com/UlrichSchimmack/zcurve3.0/refs/heads/main/Zing.25.07.11.test.R"
source(zcurve3)
pcurve <- "https://github.com/UlrichSchimmack/zcurve3.0/raw/refs/heads/main/Pcurve.Function.R"
source(pcurve)


### 1.0 – ILLUSTRATION WITH EXAMPLES












### 2.0 – SIMULATION 

sim.runs = 1                 # repeat each simulation 100 times
sim.k.sig = 10000            # how many significant results?

sim.es.mean = seq(0,1,.2)    # mean effect size, Cohen's d
sim.es.sd   = seq(0,.6,.2)   # heterogeneity of es (tau)
sim.n.obs   = c(15,20,30)    # sample size per cell 

# create the simulation conditions

sims = c()
sim.i = 0
for (i in 1:length(sim.es.mean) ) {
	for (j in 1:length(sim.es.sd)) {
		for (k in 1:length(sim.n.obs) ) {
			sim.i = sim.i + 1
			sims = rbind(sims,c(sim.i,sim.es.mean[i],sim.es.sd[j],sim.n.obs[k]) )
}}}
sims = data.frame(sims)
colnames(sims) = c("run","es.mean","es.sd","n.obs")
dim(sims)

sims

b = 1            # form row b
e = 5            # to row e, used for testing
e = nrow(sims);e # to last row, run all conditions

res = c()        # collect the results 

###

############################
### Start Simulation Loop
############################

for (run.i in b:e ) { # for loop to run the simulations

#run.i = 15       # used for testing without loop
print(paste("Run: ",run.i))
sims[run.i,]


sample1 = 1:sims$n.obs[run.i]
sample2 = (sims$n.obs[run.i]+1):(2*sims$n.obs[run.i])
sample3 = (2*sims$n.obs[run.i]+1):(3*sims$n.obs[run.i])
samples <- list(sample1 = sample1,
                    sample2 = sample2,
                    sample3 = sample3,
                    sample12 = c(sample1, sample2),
                    sample13 = c(sample1, sample3),
                    sample23 = c(sample2, sample3),
                    sample123 = c(sample1, sample2, sample3))

results = c()
count.sig = 0
count.run = 0

set.seed(run.i)

while (count.sig < sim.k.sig) {

count.run = count.run + 1

pop.es = rnorm(1,sims$es.mean[run.i],sims$es.sd[run.i]);pop.es
exp_group <- rnorm(n = 3*sims$n.obs[run.i], mean = pop.es, sd = 1)
summary(exp_group)

control <- rnorm(n = 3*sims$n.obs[run.i], mean = 0, sd = 1)
summary(control)

i.patch = 1   
res.run = c()
for (i.patch in 1:7) {
      sub_exp <- exp_group[unlist(samples[i.patch])]
      sub_control <- control[unlist(samples[i.patch])] 
      length(sub_exp)
      length(sub_control)
      tval <- t.test(sub_exp, sub_control, paird=FALSE,var.equal=TRUE)$statistic
      dfs <- t.test(sub_exp, sub_control, var.equal=TRUE)$parameter
      pval <- t.test(sub_exp, sub_control, var.equal=TRUE)$p.value
      se <- t.test(sub_exp, sub_control, var.equal=TRUE)$stderr
      es <- t.test(sub_exp, sub_control, var.equal=TRUE)$estimate
      res.run = rbind(res.run,c(count.run,pop.es,es[1]-es[2],se,tval,dfs,pval)	)
    }

res.run

results = rbind(results,res.run)

if (min(res.run[,7]) < .05) count.sig = count.sig + 1

print(count.run)
print(count.sig)


} # EOF while loop 

count.run
count.sig

dim(results)

results = data.frame(results)
results$N = results$df+2
results$se = 2/sqrt(results$N)
results$nct = abs(results[,2]/results$se)  
results$abs.t = abs(results$t)
results$z = qnorm(pt(results$abs.t,results$N-2,log.p=TRUE),log.p=TRUE)
results$pow.dir = pt(qt(.975,results$df),results$df,results$nct,lower.tail=FALSE)
results$pow.sign.err = pt(-qt(.975,results$df),results$df,results$nct,lower.tail=TRUE)
results$pow = results$pow.dir + results$pow.sign.err

summary(results)

results$rn = rep(1:7,nrow(results)/7)

tapply(results$df,results$rn,mean)
tapply(results$pow,results$rn,mean)


true.edr = tapply(results$pow,results$rn,mean)[c(1,4,7)];true.edr

true.err.1 = tapply(results$pow*results$pow.dir,results$rn,sum)
true.err.2 = tapply(results$pow,results$rn,sum)
true.err = true.err.1/true.err.2
true.err = true.err[c(1,4,7)];true.err 

true.edr
true.err

########################################
### Start P-Hacking 
########################################

phack = function(p) {
  sel = rep(0,7)
  if (p[1] < .05) sel[1] = 1	
  if (p[2] < .05) sel[2] = 1	
  if (p[3] < .05) sel[3] = 1	
  if(sum(as.numeric(p[c(1,2)] < .05)) == 0 & p[4] < .05) sel[4] = 1
  if(sum(as.numeric(p[c(1,3,4)] < .05)) == 0 & p[5] < .05) sel[5] = 1
  if(sum(as.numeric(p[c(2,3,4,5)] < .05)) == 0 & p[6] < .05) sel[6] = 1
  if(sum(as.numeric(p[1:6] < .05)) == 0 & p[7] < .05) sel[7] = 1
return(sel)
}
round(cbind(results$sim.sel.p,results$V7)[1:21,],2)


sim.sel.p = apply(matrix(results$V7,7,),2,function(x) phack(x) )
table(sim.sel.p)
results$sim.sel.p = c(sim.sel.p)
table(results$sim.sel.p,results$df)

attempts1  = nrow(results)/7*3;attempts1
success1 = sum(as.numeric(results$V7[results$rn %in% 1:3] < .05)) / attempts1;success1

attempts2  = nrow(results)/7;attempts2
success2 = sim.k.sig / attempts2;success2

table(results$df,results$rn)
tapply(results$pow,results$rn,mean)
tapply(results$pow,list(sim.sel.p,results$df),mean)



#################################################
### Analyze the Simulated Data 
#################################################


source(zcurve3)
Title = paste(round(true.edr[1]*100),"  ",round(true.err[1]*100))
res.1 = Zing(results$z[results$rn %in% 1:3]);res.1
bw.est = .02
res.1 = res.1$res
res.1

source(zcurve3)
Title = paste(round(true.edr[1]*100),"  ",round(true.err[1]*100))
Est.Method = "EXT"
FIXED.ZSDS = FALSE
ncp = 3
zsds = 2
bw.est = .02
res.2 = Zing(results$z[results$rn %in% 1:3]);res.2
res.2 = res.2$res
res.2

source(zcurve3)
Title = paste(round(true.edr[1]*100),"  ",round(true.err[1]*100))
CURVE.TYPE = "t"
Int.Beg = 2.1
bw.est = .02
res.3 = Zing(results$abs.t[results$rn %in% 1:3],df=median(results$df[results$rn %in% 1:3]));res.3
res.3 = res.3$res
res.3


source(zcurve3)
Title = paste(round(true.edr[1]*100),"  ",round(true.err[1]*100))
CURVE.TYPE = "t"
Est.Method = "EXT"
ncp = c(0,1,3)
zsds = c(1,1,1)
Int.Beg = 2.1
bw.est = .02
Augment = FALSE
res.4 = Zing(results$abs.t[results$rn %in% 1:3],df=median(results$df[results$rn %in% 1:3]));res.4
res.4 = res.4$res
res.4

pcurve.input = paste0("t(",df=results$df[results$rn %in% 1:3],") = ",
  results$abs.t[results$rn %in% 1:3])
pcurve.input = pcurve.input[results$V7[results$rn %in% 1:3] < .05]
res.pcurve.1 = pcurve_app(pcurve.input,SHOW.PLOT=TRUE)
if (res.pcurve.1[2] > res.pcurve.1[1]) res.pcurve.1[1] = res.pcurve.1[2]
res.pcurve.1 

###########################
### patchwork 
###########################

table(results$sim.sel.p,results$df)
tapply(results$pow[results$sim.sel.p == 1],results$df[results$sim.sel.p == 1],mean)

true.edr.sel = mean(results$pow[results$sim.sel.p == 1]);true.edr.sel

true.err.sel = sum(results$pow[results$sim.sel.p == 1]*results$pow.dir[results$sim.sel.p == 1]) /
  sum(results$pow[results$sim.sel.p == 1]);true.err.sel

source(zcurve3)
Title = paste(round(true.edr.sel[1]*100),"  ",round(true.err.sel[1]*100))
ymax = 2
ymin = -.08
bw.est = .02
res.5 = Zing(results$z[results$sim.sel.p == 1]);res.5
res.5 = res.5$res
res.5

source(zcurve3)
Title = paste(round(true.edr.sel[1]*100),"  ",round(true.err.sel[1]*100))
Est.Method = "EXT"
FIXED.ZSDS = FALSE
ncp = 3
zsds = 2
Int.Beg = 2
bw.est = .02
res.6 = Zing(results$z[results$sim.sel.p == 1]);res.6
res.6 = res.6$res
res.6

source(zcurve3)
Title = paste(round(true.edr.sel[1]*100),"  ",round(true.err.sel[1]*100))
CURVE.TYPE = "t"
res.7 = Zing(results$t[results$sim.sel.p == 1],df=median(results$df[results$sim.sel.p == 1]));res.3
res.7 = res.7$res
res.7

source(zcurve3)
Title = paste(round(true.edr.sel[1]*100),"  ",round(true.err.sel[1]*100));Title
CURVE.TYPE = "t"
Est.Method = "EXT"
ncp = c(0,1,3)
zsds = c(1,1,1)
Int.Beg = 2.1
bw.est = .02
Augment = FALSE
res.8 = Zing(results$t[results$sim.sel.p == 1],df=median(results$df[results$sim.sel.p == 1]))
res.8 = res.8$res
res.8

source(pcurve)
Title = paste(round(true.edr.sel[1]*100),"  ",round(true.err.sel[1]*100));Title
pcurve.input = paste0("t(",df=results$df[results$sim.sel.p == 1],") = ",
  results$abs.t[results$sim.sel.p == 1])
pcurve.input = pcurve.input[results$V7[results$sim.sel.p == 1] < .05]
res.pcurve.2 = pcurve_app(pcurve.input,SHOW.PLOT=TRUE)
if (res.pcurve.2[2] > res.pcurve.2[1]) res.pcurve.2[1] = res.pcurve.2[2]
res.pcurve.2 

res.run = rbind(
	res.1,
	res.2,
	res.3,
	res.4,
	c(NA,res.pcurve.1[1],NA,NA,NA),
	res.5,
	res.6,
	res.7,
	res.8,
	c(NA,res.pcurve.2[1],NA,NA,NA)
)

res.run = res.run[,2:3]


rownames(res.run) = c(
  "z.curve.of.no.bias","z.curve.ext.no.bias",
  "t.curve.of.no.bias","t.curve.ext.no.bias","pcurve.no.bias",
  "z.curve.of.bias","z.curve.ext.bias",
  "t.curve.of.bias","t.curve.ext.bias","pcurve.bias")

r1 = round(c(true.err[1],res.run[1:5,1]),3)
r2 = round(c(true.err.sel[1],res.run[6:10,1]),3)
r3 = round(c(true.edr[1],res.run[1:5,2]),3)
r4 = round(c(true.edr.sel[1],res.run[6:10,2]),3)

res.run = c(unlist(sims[run.i,]),
	success1,r1,r2,
	success2,r3,r4
)

print(res.run)

res = rbind(res,res.run)

write.table(res,"sim.patchwork.dat")    # write results each trial
                                   # can resume if stopped 
                                   # by loading the completed results

} # End of for loop

dim(res)  # 

#write the completed results / overwrites the data from the for loop
#write.table(res,"sim.patchwork.dat.dat")


########################################################
### GET STORED RESULTS
########################################################

# load the saved data
setwd("C:/Users/ulric/Dropbox/PHPCurves/DOCS/z-curve/Tutorial")
check = read.table("sim.patchwork.dat",row.names=NULL)
check = data.frame(check)
check = check[,2:dim(check)[2]]
dim(check)
table(round(res[,c(1:14)],5) == round(check[,c(1:14)],5) )

check[,c(1:6,13)]

# careful, if you run this code any results in the res file will be replaced 
# with the saved data. 

#res = check  # use only if sure to replace res; run analysis with "res" matrix 

res = data.frame(res)
dim(res)

summary(res)



#############################################
### Inflation (False Positive Rate for EDR = 5%
#############################################

res$inflation = res[,5]-res[,13]
summary(res$inflation)

H0true = as.numeric(res[,2] == 0 & res[,3] == 0)
table(H0true)
tapply(res$inflation,H0true,mean)
summary(lm(res$inflation ~ H0true + res[,2] + res[,3] + res[,13]))

round(cbind(res[,c(1:5,13)],res$inflation),3)

summary(res[,c(5,7)])
plot(res[,7],res[,5],ylim=c(0,1),xlim=c(0,1))
abline(a = 0,b = 1)



########################################################################
### Evaluation of ERR estimates (no bias) 
########################################################################

# columns
# -  6: true ERR
# -  7: z.1
# -  8: z.2 
# -  9: t.1
# - 10: t.2
# - 11: pcurve

summary(res[,6:11])

cor(res[,6:11])[,1]

rmse.err.1 = sqrt(mean((res[,6]-res[,7])^2));rmse.err.1
rmse.err.2 = sqrt(mean((res[,6]-res[,8])^2));rmse.err.2
rmse.err.3 = sqrt(mean((res[,6]-res[,9])^2,na.rm=TRUE));rmse.err.3
rmse.err.4 = sqrt(mean((res[,6]-res[,10])^2,na.rm=TRUE));rmse.err.4
rmse.err.5 = sqrt(mean((res[,6]-res[,11])^2,na.rm=TRUE));rmse.err.5

### all sig
rmse.err.1
rmse.err.2
rmse.err.3
rmse.err.4
rmse.err.5


########################################################################
### Evaluation of ERR estimates (patchwork p-hacked) 
########################################################################

# columns
# - 12: true ERR
# - 13: z.1
# - 14: z.2 
# - 15: t.1
# - 16: t.2
# - 17: pcurve

summary(res[,12:17])

cor(res[,12:17])[,1]

rmse.err.sel.1 = sqrt(mean((res[,12]-res[,13])^2));rmse.err.sel.1
rmse.err.sel.2 = sqrt(mean((res[,12]-res[,14])^2));rmse.err.sel.2
rmse.err.sel.3 = sqrt(mean((res[,12]-res[,15])^2,na.rm=TRUE));rmse.err.sel.3
rmse.err.sel.4 = sqrt(mean((res[,12]-res[,16])^2,na.rm=TRUE));rmse.err.sel.4
rmse.err.sel.5 = sqrt(mean((res[,12]-res[,17])^2,na.rm=TRUE));rmse.err.sel.5

### all sig
rmse.err.sel.1
rmse.err.sel.2
rmse.err.sel.3
rmse.err.sel.4
rmse.err.sel.5


res = res[order(res[,5]),]

lwd.ci = .5

### all significant results

graphics.off()
plot(res[,12],res[,17],pch=15,cex=1,xlim=c(0,1),ylim=c(0,1),col="purple3",ylab="",xlab="")  # plot OF and ER estimates
par(new=TRUE)
plot(res[,12],res[,13],xlim=c(0,1),ylim=c(0,1),pch=16,cex=1,col="forestgreen",ylab="Estimated ERR",
	xlab = "Simulated True ERR")  # plot OF and ER estimates
par(new=TRUE)
plot(res[,12],res[,15],xlim=c(0,1),ylim=c(0,1),pch=16,cex=1,col="firebrick3",ylab="Estimated ERR",
	xlab = "Simulated True ERR")  # plot OF and ER estimates
abline(a = 0, b = 1,lty=2)
abline(v = .5,lty=2)
abline(h = .5,lty=2)
legend(.1,.9,legend=c("z-curve", "t-curve","p-curve"),lwd=2,col=c("forestgreen","firebrick3","purple3"))






### sel significant results

graphics.off()
plot(res[,6],res[,11],pch=15,cex=1,xlim=c(0,1),ylim=c(0,1),col="purple3",ylab="",xlab="")  # plot OF and ER estimates
par(new=TRUE)
plot(res[,6],res[,7],xlim=c(0,1),ylim=c(0,1),pch=16,cex=1,col="forestgreen",ylab="Estimated ERR",
	xlab = "Simulated True ERR")  # plot OF and ER estimates
abline(a = 0, b = 1,lty=2)
abline(v = .5,lty=2)
abline(h = .5,lty=2)
legend(.1,.9,legend=c("z-curve", "p-curve"),lwd=2,col=c("forestgreen","purple3"))


graphics.off()
plot(res[,6],res[,12],pch=15,cex=1,xlim=c(0,1),ylim=c(0,1),col="purple3",ylab="",xlab="")  # plot OF and ER estimates
par(new=TRUE)
plot(res[,6],res[,10],xlim=c(0,1),ylim=c(0,1),pch=16,cex=1,col="forestgreen",ylab="Estimated ERR",
	xlab = "Simulated True ERR")  # plot OF and ER estimates
abline(a = 0, b = 1,lty=2)
abline(v = .5,lty=2)
abline(h = .5,lty=2)
legend(.1,.9,legend=c("z-curve", "p-curve"),lwd=2,col=c("forestgreen","purple3"))

res = res[order(res[,1]),]

res[,1:6]

dir.bias.z = cbind(res[,1],res[,10]-res[,6],res[,c(2:5,10,6)])
round(dir.bias.z[abs(dir.bias.z[,2]) > .30,],2)
summary(lm(dir.bias.z[,1] ~ res[,2] + res[,3] + res[,4]))
tapply(dir.bias.z[,1],res[,3],mean)

dir.bias.p = cbind(res[,13]-res[,5],res[,c(1:5,13)])
dir.bias.p[dir.bias.p < -.25] = NA
round(dir.bias.p[abs(dir.bias.p[,1]) > .25,],2)
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



df = 125
t = 7.11 
se = 1/sqrt(df+1)

d = t*se
d
ci.lb = d - 2*se
d + 2*se

ci.lb
N = 60
se = 1/sqrt(N)
ci.lb/se

pt(qt(.975,59),59,5.53,lower.tail=FALSE)





