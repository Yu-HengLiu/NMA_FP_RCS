
#~~~~~~~~~~~~~~~~~~~~3.1st order FP NMA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Start with an empty environment
rm(list = ls())

# Fractional Polynomial model - 1st Order fixed effect
path <- "your working directory/R"
setwd(path)

# Load libraries
library(coda)
library(survival)
library(R2WinBUGS)
library("dplyr")
# Set the working directory


# Load the data
data <- read.csv("aRCC_aggregate_data.csv") #mo

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
num.sims <- 30000
burn.in <- 50000

# Fractional polynomial powers
######change the setting here
P1<--2
P1<--1
P1<--0.5

P1 <- 0

P1<-0.5
P1 <- 1
P1<-2
P1<-3

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

# Mean & precision
mean <- c(0,0)
prec <- array(c(0.0001, 0, 0, 0.0001), dim=c(2,2))#precision is used in JAGS and WinBugs instead of sd or variance

# Number of treatment arms for each trial
na <- c(2, 2, 2)

# Treatment in each trial arm - This fills the down the columns first
t <- array(data=c(1,1,1,#first column
                  2,3,4),#second column
           dim=c(3, 2))#3 rows, 2 columns

#Generate survival from Winbugs
#Generate parameters only from winbugs
bugs_data <- list(s=data$trialid, r=data$nevents, z=data$natrisk, a=data$arm, time=data$time,
                  dt=data$length, P1=P1, N=N, nt=nt, ns=ns,  mean=mean, prec=prec,
                  t=t,  na=na)


#-----------------------------------------------------------------------------
# Model
#-----------------------------------------------------------------------------
#re-run each time when change parameters

# Create initial values for model
d1 <- array(c(NA, rep(0.1, 3), NA, rep(0.2, 3)), dim=c(nt,2))#4x2 in total
d2 <- array(c(NA, rep(0.2, 3), NA, rep(-0.1, 3)), dim=c(nt,2))
d3 <- array(c(NA, rep(-0.1, 3), NA, rep(0.1, 3)), dim=c(nt,2))
#nt was defined-->num of tx #what does the numbers stand for?do I need to change them?why is there 3 parameters?
mu1 <- array(rep(c(0.4, 0.5, 0.6, 0.2, 0.3, 0.1, 0.1),2), dim=c(ns,2)) 
mu2 <- array(rep(c(0.3, 0.4, 0.5, 0.6, 0.7, -0.1, -0.2),2), dim=c(ns,2))
mu3 <- array(rep(c(0.5, 0.6, 0.7, 0.1, -0.1,  0.2, 0.2),2), dim=c(ns,2))
#ns was defined-->num of studies #what does the numbers stand for?do I need to change them?

inits <- list(list(d=d1, mu=mu1), 
              list(d=d2, mu=mu2),
              list(d=d3, mu=mu3))#I think this is 3 chains


#only generate parameters
bugs.object <- bugs(data=bugs_data, inits=inits, 
                    parameters.to.save=c( "alpha0", "alpha1", "d","mu"), 
                    model.file="FP_1st_order_model.txt", clearWD=F, 
                    summary.only=FALSE, n.iter=(num.sims+burn.in), 
                    n.sims=num.sims, n.burnin=burn.in, n.chains=3, 
                    bugs.seed=385916, bugs.directory=bugs.directory, 
                    debug=F, working.directory=path)
#alpha0, alpha1 can be used to calculate extrapolated estimation in excel
results <- bugs.object$summary

# Save results in csv file
save.image(file = "P1.m2.b5i3.mo.RData")
save.image(file = "P1.m1.b5i3.mo.RData")
save.image(file = "P1.m0.5.b5i3.mo.RData")
save.image(file = "P1.0.b5i3.mo.RData")

write.csv(results,file="1FP_p1_m0.5_parameters.csv")
write.csv(results,file="1FP_p1_m1_parameters.csv")
write.csv(results,file="1FP_p1_0_parameters.csv")
