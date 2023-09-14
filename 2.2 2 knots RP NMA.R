#~~~~~~~~~~~~~~~~~~~~~~2.2 RP NMA, 2 knots~~~~~~~~~~~~~~~~~~~~~#
# Royston-Parmar model in R - fixed effect
#time unit: wk, mo
# Empty environment
library("dplyr")
install.packages("tidyverse")

rm(list = ls()) 
# Load libraries
library(survival)
library(foreach)
library(R2WinBUGS)

# Working diectory
path <- "your working directory/R"
setwd(path)

# Load IPD data
data <- read.csv("NMA_data.csv")#month
data <- read.csv("NMA_data.csv") %>%mutate(time=time*(365.25/7)/12) #mo->wk
# Set the location for WinBUGS
bugs.directory <- "your winbugs location/winbugs143_unrestricted/winbugs14_full_patched/WinBUGS14"

# WinBUGS burn-in & simulation size
num.sims <- 10000
burn.in <- 10000

# Number of studies
num_studies <- length(unique(data$studyCode))

# Sort data into trialid order
data <- data[order(data$studyCode),]

# Number of patients in each study
# Note: if data is re-ordered will now to re-order this vector
num_patients <- c(712,157,847)

# Maximum time
maxt <-2080#40 yrs
maxt <-480
# Number of treatments
nt <- length(unique(data$txCode))

# Create treatment indicator variables - caution needed to ensure consistency equations hold
data$trt2 <- 0
data$trt2[data$txCode==2] <- 1

data$trt3 <- 0
data$trt3[data$txCode==3] <- 1

data$trt4 <- 0
data$trt4[data$txCode==4] <- 1


# Add  trialid variable - easier than changing trialid to studyCode throughout this file!
data$trialid <- data$studyCode

# Put eventtime onto the ln scale
data$eventtime <- data$time
data$lnt <- log(data$eventtime)

# Create trt*lnt variables
data$trt2lnt <- data$trt2*data$lnt
data$trt3lnt <- data$trt3*data$lnt
data$trt4lnt <- data$trt4*data$lnt


# Total number of patients
pts <- nrow(data)

# Set location of knots for each trial - 33rd and 67th percentiles of uncensored survival times
#two knots, k1 & k4 boundary knots (Kmin and Kmax), k2=first knot at 33rd pctl, k3=second knot at 67 pctl
knots <- data.frame(trialid=NA, k1=NA, k2=NA, k3=NA, k4=NA) 
foreach(i=1:num_studies) %do% {
  k <- quantile(data$lnt[data$status==1 & data$trialid==i], c(0, 0.33, 0.67, 1))
  knots <- rbind(knots, data.frame(trialid=i, k1=k[1], k2=k[2], k3=k[3], k4=k[4]))
}
# Remove the empty top row
knots <- knots[-1,]

# Add knot values to data #for each trial, you will have 4 survival time points (2 knots)
foreach(i=1:num_studies) %do% {
  data$k1[data$trialid==i] <- knots$k1[i]
  data$k2[data$trialid==i] <- knots$k2[i]
  data$k3[data$trialid==i] <- knots$k3[i]
  data$k4[data$trialid==i] <- knots$k4[i]
}  


#calculate basis function first
data$t1 <- data$lnt-data$k2
data$t2 <- data$lnt-data$k1  #min
data$t3 <- data$lnt-data$k4  #max
data$t4 <- data$lnt-data$k3
#t2(kmin) t3(kmax) is used in every  v
foreach(i=1:pts) %do% {
  data$a[i] <- max(c(data$t1[i], 0))
  data$b[i] <- max(c(data$t2[i], 0))
  data$c[i] <- max(c(data$t3[i], 0))
  data$e[i] <- max(c(data$t4[i], 0))
  #2 knots, 2 phi
  data$phi1[i] <- (data$k4[i]-data$k2[i])/(data$k4[i]-data$k1[i]) #first knot kmax-k1st knot/kmax-kmin
  data$phi2[i] <- (data$k4[i]-data$k3[i])/(data$k4[i]-data$k1[i]) #second knot:kmax-k2nd knot/kmax-kmin
  #2 knots, 3 v
  data$v0[i] <- data$lnt[i]
  data$v1[i] <- data$a[i]^3 - (data$phi1[i]*(data$b[i]^3)) - ((1-data$phi1[i])*(data$c[i]^3))  #1st knot
  data$v2[i] <- data$e[i]^3 - (data$phi2[i]*(data$b[i]^3)) - ((1-data$phi2[i])*(data$c[i]^3))  #2nd knot
}
# Orthogonalisation: trransfer v to mu
#1. v0->u0
#first calculate the gathered value, 
df1 <- data.frame(trialid=NA, mean=NA, sd=NA)
foreach(i=1:num_studies) %do% {
  mv0 <- mean(data$v0[data$trialid==i])
  sdv0 <- sd(data$v0[data$trialid==i])
  df1 <- rbind(df1, data.frame(trialid=i, mean=mv0, sd=sdv0))
}
df1 <- df1[-1,] #each trial have the mean and sd of v0

# then Add values back to data
foreach(i=1:num_studies) %do% {
  data$mv0[data$trialid==i] <- df1$mean[i]
  data$sdv0[data$trialid==i] <- df1$sd[i]
} 
#then calculate the transformed v (which is u) for each person each trial
foreach(i=1:pts) %do% {
  data$u0[i] <- (data$v0[i]-data$mv0[i])/data$sdv0[i] #normalising v0 to u0
  data$f[i] <- data$v1[i]*data$u0[i]#for following calc. of v1
  data$g[i] <- data$u0[i]*data$u0[i]#for following calc.of v1
}
#2.calculate gathered value, v1->u1
df2 <- data.frame(trialid=NA, v1u0=NA, u0u0=NA)
foreach(i=1:num_studies) %do% {
  v1u0 <- sum(data$f[data$trialid==i])#inner product, v1(pt1)*u0(pt2)+v1(pt2)*u0(pt2)....
  u0u0 <- sum(data$g[data$trialid==i])
  df2 <- rbind(df2, data.frame(trialid=i, v1u0=v1u0, u0u0=u0u0))
}
df2 <- df2[-1,]

# Add values back to data: these are sums
foreach(i=1:num_studies) %do% {
  data$v1u0[data$trialid==i] <- df2$v1u0[i]
  data$u0u0[data$trialid==i] <- df2$u0u0[i]
}
#calculate un1
foreach(i=1:pts) %do% {
  data$u1[i] <- data$v1[i] - ((data$v1u0[i]/data$u0u0[i])*data$u0[i])
}
#calculate the mean and sd for u1, normalising un1->u1
df3 <- data.frame(trialid=NA, mean=NA, sd=NA)
foreach(i=1:num_studies) %do% {
  mu1 <- mean(data$u1[data$trialid==i])
  sdu1 <- sd(data$u1[data$trialid==i])
  df3 <- rbind(df3, data.frame(trialid=i, mean=mu1, sd=sdu1))
}
df3 <- df3[-1,]

# Add values back to data
foreach(i=1:num_studies) %do% {
  data$mu1[data$trialid==i] <- df3$mean[i]
  data$sdu1[data$trialid==i] <- df3$sd[i]
} 
#normalising u1, prepare to use it to calculate v2->u2
foreach(i=1:pts) %do% {
  data$u1n[i] <- (data$u1[i]-data$mu1[i])/data$sdu1[i]
  data$h[i] <- data$v2[i]*data$u0[i]
  data$aa[i] <- data$v2[i]*data$u1n[i]
  data$j[i] <- data$u1n[i]*data$u1n[i]
}
#3. v2->u2
df4 <- data.frame(trialid=NA, v2u0=NA, v2u1n=NA, u1nu1n=NA)
foreach(i=1:num_studies) %do% {
  v2u0 <- sum(data$h[data$trialid==i])
  v2u1n <- sum(data$aa[data$trialid==i])
  u1nu1n <- sum(data$j[data$trialid==i])
  df4 <- rbind(df4, data.frame(trialid=i, v2u0=v2u0, v2u1n=v2u1n, u1nu1n=u1nu1n))
}
df4 <- df4[-1,]

# Add values back to data
foreach(i=1:num_studies) %do% {
  data$v2u0[data$trialid==i] <- df4$v2u0[i]
  data$v2u1n[data$trialid==i] <- df4$v2u1n[i]
  data$u1nu1n[data$trialid==i] <- df4$u1nu1n[i]
}
#calculate un2
foreach(i=1:pts) %do% {
  data$u2[i] <- data$v2[i] - ((data$v2u0[i]/data$u0u0[i])*data$u0[i]) - ((data$v2u1n[i]/data$u1nu1n[i])*data$u1n[i])
}
#normalising un2->u2
df5 <- data.frame(trialid=NA, mean=NA, sd=NA)
foreach(i=1:num_studies) %do% {
  mu2 <- mean(data$u2[data$trialid==i])
  sdu2 <- sd(data$u2[data$trialid==i])
  df5 <- rbind(df5, data.frame(trialid=i, mean=mu2, sd=sdu2))
}
df5 <- df5[-1,]

#Add values back to data
foreach(i=1:num_studies) %do% {
  data$mu2[data$trialid==i] <- df5$mean[i]
  data$sdu2[data$trialid==i] <- df5$sd[i]
} 


foreach(i=1:pts) %do% {
  data$u2n[i] <- (data$u2[i]-data$mu2[i])/data$sdu2[i]#normalising u2
  #these are for likelihood estimation
  data$du0[i] <- 1/data$sdv0[i]
  data$k[i] <- (3*(data$a[i]^2)) -((3*data$phi1[i])*(data$b[i]^2))- (3*(1-data$phi1[i])*(data$c[i]^2))
  data$l[i] <- data$sdu1[i]
  data$du1n[i] <- (1/data$l[i])*data$k[i] - (((data$v1u0[i]/data$u0u0[i])/data$l[i])*data$du0[i])
  data$m[i] <- (3*(data$e[i]^2)) - ((3*data$phi2[i])*(data$b[i]^2))- (3*(1-data$phi2[i])*(data$c[i]^2))
  data$n[i] <- data$sdu2[i]
  data$du2n[i] <- ((1/data$n[i])*data$m[i]) - (((data$v2u0[i]/data$u0u0[i])/data$n[i])*data$du0[i]) - (((data$v2u1n[i]/data$u1nu1n[[i]])/data$n[i])*data$du1n[i])
}
#cumulative pts number (0, num in s1, num in s1+s2, num in s1+s2+s3... )
offset <- 0
foreach(i=1:num_studies) %do% {
  offset[i+1] <- offset[i] + num_patients[i]
}


# Prepare to fit model in WinBUGS

bugs_data <- list(d=data$status, trt2=data$trt2, trt3=data$trt3, trt4=data$trt4, 
                  trt2lnt=data$trt2lnt, trt3lnt=data$trt3lnt, trt4lnt=data$trt4lnt, 
                  u0=data$u0, u1n=data$u1n, u2n=data$u2n, du0=data$du0,
                  du1n=data$du1n, du2n=data$du2n, Ntrials=num_studies, offset=offset,
                  k0=knots$k1[1], k1=knots$k2[1], k2=knots$k3[1], k3=knots$k4[1],#reference study:study 1
                  meanv0=df1$mean[1], sdv0=df1$sd[1], v1w0=df2$v1u0[1], w0w0=df2$u0u0[1],
                  meanw1=df3$mean[1], sdw1=df3$sd[1], v2w0=df4$v2u0[1], v2w1n=df4$v2u1n[1],
                  w1nw1n=df4$u1nu1n[1], meanw2=df5$mean[1], sdw2=df5$sd[1],
                  nt=nt,maxt=maxt)

# Create initial values for model:3 chains
beta1 <- c(rep(0.1, 6))
beta2 <- c(rep(0.2, 6))
beta3 <- c(rep(0.3, 6))

gamma1 <- array(rep(c(0.2), 4*num_studies), dim=c(4,num_studies))
gamma2 <- array(rep(c(0.4), 4*num_studies), dim=c(4,num_studies))
gamma3 <- array(rep(c(0.1), 4*num_studies), dim=c(4,num_studies))

inits <- list(list(beta=beta1, gamma=gamma1), 
              list(beta=beta2, gamma=gamma2),
              list(beta=beta3, gamma=gamma3))

beta1 <- c(rep(0.1, 6))
beta2 <- c(rep(0.2, 6))
beta3 <- c(rep(0.3, 6))

gamma1 <- array(rep(c(0.2), 4*num_studies), dim=c(4,num_studies))
gamma2 <- array(rep(c(0.4), 4*num_studies), dim=c(4,num_studies))
gamma3 <- array(rep(c(0.1), 4*num_studies), dim=c(4,num_studies))

inits <- list(list(beta=beta1, gamma=gamma1), 
              list(beta=beta2, gamma=gamma2),
              list(beta=beta3, gamma=gamma3))
# Fit model in WinBUGS
start_time <- system.time(
bugs.object<-bugs(data=bugs_data, inits=inits, 
                  parameters.to.save=c("surv", "gamma","beta"), 
                  model.file="RP model - 2 knots.txt", clearWD=F, 
                  summary.only=FALSE, n.iter=(num.sims+burn.in), n.sims=num.sims,
                  n.burnin=burn.in, n.chains=3, 
                  bugs.seed=3, bugs.directory=bugs.directory, 
                  debug=F, working.directory=path)
)
elapsed_time <- as.numeric(start_time["elapsed"])
print(paste("Elapsed time:", elapsed_time, "seconds"))
results <- bugs.object$summary #get the survival function at each month for each trt
#summary(results)
# Save results as a csv file
write.csv(results, "RP_2k_results_wk_40.csv")#fully extrapolated for 40 yrs
write.csv(results, "RP_2k_results_mo_40.csv")#fully extrapolated for 40 yrs

