
#~~~~~~~~~~~~~~~~~~~~3.2 2ND order FP NMA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

####################
# Start with an empty environment
rm(list = ls())

# Fractional Polynomial model - fixed effect
path <- "your working directory/R"
setwd(path)

# Load libraries
library(coda)
library(survival)
library(R2WinBUGS)
library("dplyr")
# Set the working directory


# Load the data
data <- read.csv("aRCC_aggregate_data.csv")#0,3,6,12 mo

# Create a treatment code variable as an integer
data$txCode[data$treatment=="sun"] <- 1
data$txCode[data$treatment=="pem+len"] <- 2
data$txCode[data$treatment=="cab"] <- 3
data$txCode[data$treatment=="niv+ipi"] <- 4

# Order data
data <- data[order(data$trialid, data$txCode, data$spgrp),]#ordered by trialid->txcode->spgrp

# Set the location for WinBUGS
bugs.directory <- "your winbugs location/winbugs143_unrestricted/winbugs14_full_patched/WinBUGS14"

# WinBUGS burn-in & simulation size
######change the setting here

num.sims <- 100000
burn.in <- 100000


# Fractional polynomial powers
P1<--0.5
#change P2
P2<--2
P2<--1
P2<--0.5
P2<-0
P2<-1
P2<-2
P2<-3
#-----------------------------------------------------------------------------
# Data formatting
#-----------------------------------------------------------------------------
#re-run each time when change parameters

# Need to number the treatment arms within each trial
data$arm[data$treatment=="sun"] <- 1
data$arm[data$treatment=="pem+len" | data$treatment=="cab" | data$treatment=="niv+ipi"] <- 2

# Check all arms coded
data$arm==1 | data$arm==2 

# Length of time intervals
data$length <- data$time-data$start

# Number of treatments
nt <- length(unique(data$treatment))

# Number of studies
ns <- length(unique(data$trialid)) 

# Number of rows in dataset
N <- nrow(data)

# Mean & precision :3 parameters
mean <- c(0,0,0)
prec <- array(c(0.0001, 0, 0, 0,0.0001,0,0,0,0.0001), dim=c(3,3))#precision is used in JAGS and WinBugs instead of sd or variance

# Number of treatment arms for each trial
na <- c(2, 2, 2)

# Treatment in each trial arm - This fills the down the columns first
t <- array(data=c(1,1,1,#first column #do i need to change this?
                  2,3,4),#second column
           dim=c(3, 2))#3 rows, 2 columns

bugs_data <- list(s=data$trialid, r=data$nevents, z=data$natrisk, a=data$arm, time=data$time,
                  dt=data$length, P1=P1,P2=P2, N=N, nt=nt, ns=ns, mean=mean, prec=prec,
                  t=t,  na=na)


#-----------------------------------------------------------------------------
# Model
#-----------------------------------------------------------------------------
#re-run each time when change parameters

# Create initial values for model: 3 chains, 3 parameters-->change to three, is it correct?
d1 <- array(c(NA, rep(0.1, 3), NA, rep(0.2, 3)), dim=c(nt,3))#4x3 in total
d2 <- array(c(NA, rep(0.2, 3), NA, rep(-0.1, 3)), dim=c(nt,3))
d3 <- array(c(NA, rep(-0.1, 3), NA, rep(0.1, 3)), dim=c(nt,3))
#nt was defined-->num of tx #what does the numbers stand for?
mu1 <- array(rep(c(0.4, 0.5, 0.6, 0.2, 0.3, 0.1, 0.1),2), dim=c(ns,3)) 
mu2 <- array(rep(c(0.3, 0.4, 0.5, 0.6, 0.7, -0.1, -0.2),2), dim=c(ns,3))
mu3 <- array(rep(c(0.5, 0.6, 0.7, 0.1, -0.1,  0.2, 0.2),2), dim=c(ns,3))
#ns was defined-->num of studies #what does the numbers stand for?do I need to change them?

inits <- list(list(d=d1, mu=mu1), 
              list(d=d2, mu=mu2),
              list(d=d3, mu=mu3))#I think this is 3 chains

bugs.object <- bugs(data=bugs_data, inits=inits, 
                    parameters.to.save=c( "alpha0", "alpha1", "alpha2","d","mu"), 
                    model.file="FP_2nd_order_model.txt", clearWD=F, 
                    summary.only=FALSE, n.iter=(num.sims+burn.in), 
                    n.sims=num.sims, n.burnin=burn.in, n.chains=3, 
                    bugs.seed=385916, bugs.directory=bugs.directory, 
                    debug=F, working.directory=path)
results <- bugs.object$summary

save.image(file = "P2.m2.b10i10.mo.RData")
save.image(file = "P2.m1.b10i10.mo.RData")
save.image(file = "P2.m0.5.b10i10.mo.RData")
save.image(file = "P2.0.b10i10.mo.RData")
# Fractional Polynomial model survival graph
write.csv(results,file="2FP_p2_m0.5_parameters.csv")
write.csv(results,file="2FP_p2_m1_parameters.csv")
write.csv(results,file="2FP_p2_0_parameters.csv")
write.csv(results,file="2FP_p2_m2_parameters.csv")
