
#~~~~~~~~~~~~~~3.3 output generation for 1st order FP NMA~~~~~~~~~~~~~~~#
################################################
#method 1: use parameters estimated in month
################################################
rm(list = ls())
#1. load data first
results <- bugs.object$summary #this is the estimation using month as time scale
alpha1<-bugs.object$sims.list$alpha1
alpha0 <- bugs.object$sims.list$alpha0

#2.time1 calc.
  #check: if monitor the prob monthly: estimation looks the same as prob generated from winbugs
  maxt<-480#480 mo, 40 yrs
  time1 <- numeric(maxt)
  for (t in 1:maxt) {
    if (P1 == 0) {
      time1[t] <- log(t)
    } else {
      time1[t] <- (1 - (P1 == 0)) * t ^ P1
    }
  }

#if monitor prob weekly: estimation doesn't look correct??????
  maxt<-2080 #week as a cycle
  time1 <- numeric(maxt)
  converted_time <- numeric(maxt)
  for (t in 1:maxt) {
    converted_time[t] <- t *(12*7/365.25) #week as time unit, 1 wk=0.2308 months
  #  converted_time <- t*(365.25/7)/12
    # Calculate time1 based on conditional operations
    if (P1 == 0) {
      time1[t] <- log(converted_time[t])
    } else {
      time1[t] <- (1 - (P1 == 0)) * converted_time[t] ^ P1
    }
  }
  #replace the time1 estimates with NA within 1 mo
  converted_time[1:4] <- c(0, 0, 0, 0)
  time1[1:4] <- c(NA, NA, NA, NA)
#3. Initialize matrices/arrays for hazards, cumulative hazards, mortality, and survival
#mo: array 30000 x 480 x 4
#wk: array 30000 x 2080 x 4
HAZARD <- array(NA, dim = c(nrow(alpha1), maxt, nt))
CUM_H <- array(NA, dim = c(nrow(alpha1), maxt, nt))
T <- array(NA, dim = c(nrow(alpha1), maxt, nt))
S <- array(NA, dim = c(nrow(alpha1), maxt, nt))
  #dim(HAZARD)

#4. generate prob in each sample
#month
for (i in 1:nrow(alpha1)) {
  for (k in 1:nt) {
    for (t in 1:maxt) {
      # Calculate hazard
        HAZARD[i, t, k] <- exp(alpha0[i, k] + alpha1[i, k] * time1[t])
      
      # Calculate cumulative hazard
      if (t == 1) {
        CUM_H[i, t, k] <- HAZARD[i, t, k]  # For the first time point, cumulative hazard is equal to hazard
      } else {
        CUM_H[i, t, k] <- CUM_H[i, t - 1, k] + HAZARD[i, t, k]
      }
      
      # Calculate mortality over time
      T[i, t, k] <- 1 - exp(-CUM_H[i, t, k])
      
      # Calculate survival over time
      S[i, t, k] <- 1 - T[i, t, k]
    }
  }
}
#week
for (i in 1:nrow(alpha1)) {
  for (k in 1:nt) {
    for (t in 1:maxt) {
      # Calculate hazard
      if (t<5){
        HAZARD[i, t, k] <-0
      }else{
        HAZARD[i, t, k] <- exp(alpha0[i, k] + alpha1[i, k] * time1[t])*(converted_time[t]-converted_time[t-1])
      }
      
      # Calculate cumulative hazard
      if (t == 1) {
        CUM_H[i, t, k] <- HAZARD[i, t, k]  # For the first time point, cumulative hazard is equal to hazard
      } else {
        CUM_H[i, t, k] <- CUM_H[i, t - 1, k] + HAZARD[i, t, k]
      }
      
      # Calculate mortality over time
      T[i, t, k] <- 1 - exp(-CUM_H[i, t, k])
      
      # Calculate survival over time
      S[i, t, k] <- 1 - T[i, t, k]
    }
  }
}

#5.calc mean prob.
# Calculate the mean value for HAZARD across dimensions m and k
mean_HAZARD <- apply(HAZARD, c(2, 3), mean)
# Calculate the mean value for CUM_H across dimensions m and k
mean_CUM_H <- apply(CUM_H, c(2, 3), mean)
# Calculate the mean value for T across dimensions m and k
mean_T <- apply(T, c(2, 3), mean)
# Calculate the mean value for S across dimensions m and k
mean_S <- apply(S, c(2, 3), mean)
Survival<-data.frame(HAZARD=mean_HAZARD, CUM_H=mean_CUM_H,T=mean_T, S=mean_S)

#calculate HR to prepare for plotting
HR_PL_sun<-mean_HAZARD[,2]/mean_HAZARD[,1]
HR_Cab_sun<-mean_HAZARD[,3]/mean_HAZARD[,1]
HR_NI_sun<-mean_HAZARD[,4]/mean_HAZARD[,1]
HR<-data.frame(HR_PL_sun=HR_PL_sun,HR_Cab_sun=HR_Cab_sun,HR_NI_sun=HR_NI_sun)
HR <- HR[-1,]
# Save results in csv file
write.csv(Survival, "P1.m2.b5i3.mo_40_S.csv")#fully extrapolated for 40 yrs
write.csv(Survival, "P1.m1.b5i3.mo_40_S.csv")#fully extrapolated for 40 yrs
write.csv(Survival, "P1.m0.5.b5i3.mo_40_S.csv")#fully extrapolated for 40 yrs
write.csv(Survival, "P1.0.b5i3.mo_40_S.csv")#fully extrapolated for 40 yrs

#write.csv(Survival, "P1.m2.b5i3.wk_40_S.csv")#fully extrapolated for 40 yrs
write.csv(Survival, "P1.m1.b5i3.wk_40_S.csv")#fully extrapolated for 40 yrs
write.csv(Survival, "P1.m0.5.b5i3.wk_40_S.csv")#fully extrapolated for 40 yrs
write.csv(Survival, "P1.0.b5i3.wk_40_S.csv")#fully extrapolated for 40 yrs
