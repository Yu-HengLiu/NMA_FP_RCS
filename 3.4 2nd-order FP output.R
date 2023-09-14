
#~~~~~~~~~~~~~~3.4 output generation for 2nd order FP NMA~~~~~~~~~~~~~~~#
################################################
#method 1: use parameters estimated in month
################################################
rm(list = ls())
#1. load the data first
results <- bugs.object$summary 
alpha2 <- bugs.object$sims.list$alpha2
alpha1<-bugs.object$sims.list$alpha1
alpha0 <- bugs.object$sims.list$alpha0

#week: restrict samples to 30000 samples (same as 1FP)
alpha0<-alpha0[1:(nrow(alpha0)-70002),]
alpha1<-alpha1[1:(nrow(alpha1)-70002),]
alpha2<-alpha2[1:(nrow(alpha2)-70002),]
#2. time1 calc.
#check: if monitor the prob monthly
maxt<-480#480 mo, 40 yrs
t1 <- numeric(maxt)
t2 <- numeric(maxt)
for (m in 1:maxt) {
  if (P1 == 0) {
    t1[m] <- log(m)
  } else {
    t1[m] <- m^P1
  }
  
  if (P2 == P1) {
    if (P2 == 0) {
      t2[m] <- log(m)^2
    } else {
      t2[m] <- m^P2 * log(m)
    }
  } else {
    if (P2 == 0) {
      t2[m] <- log(m)
    } else {
      t2[m] <- (1 - (P2 == P1)) * (P2 == 0) * log(m) + (1 - (P2 == P1)) * (1 - (P2 == 0)) * m^P2
    }
  }
}

#if monitor prob weekly
maxt<-2080 #week as a cycle
t1 <- numeric(maxt)
t2 <- numeric(maxt)
converted_time <- numeric(maxt)
for (t in 1:maxt) {
  converted_time[t] <- t *(12*7/365.25) #week as time unit, 1 wk=0.2308 months
  
  if (P1 == 0) {
    t1[t] <- log(converted_time[t])
  } else {
    t1[t] <- converted_time[t]^P1
  }
  
  if (P2 == P1) {
    if (P2 == 0) {
      t2[t] <- log(converted_time[t])^2
    } else {
      t2[t] <- converted_time[t]^P2 * log(converted_time[t])
    }
  } 
  else {
    if (P2 == 0) {
      t2[t] <- log(converted_time[t])
    } else {
      t2[t] <- (1 - (P2 == P1)) * (P2 == 0) * log(converted_time[t]) + (1 - (P2 == P1)) * (1 - (P2 == 0)) * converted_time[t]^P2
    }
  }
}

#replace the time1 estimates with NA within 1 mo
converted_time[1:4] <- c(0, 0, 0, 0)
t1[1:4] <- c(NA, NA, NA, NA)
t2[1:4] <- c(NA, NA, NA, NA)
#3. Initialize matrices/arrays for hazards, cumulative hazards, mortality, and survival
#mo: array 10000 x 480 x 4
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
      HAZARD[i, t, k] <- exp(alpha0[i, k] + alpha1[i, k] * t1[t]+ alpha2[i, k] * t2[t])
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
      if (t<5){
        HAZARD[i, t, k] <-0
      }else{
        HAZARD[i, t, k] <- exp(alpha0[i, k] + alpha1[i, k] * t1[t]+alpha2[i, k] * t2[t])*(converted_time[t]-converted_time[t-1])
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
write.csv(Survival, "P2.m2.b10i10.mo_40_S.csv")#fully extrapolated for 40 yrs
write.csv(Survival, "P2.m1.b10i10.mo_40_S.csv")#fully extrapolated for 40 yrs
write.csv(Survival, "P2.m0.5.b10i10.mo_40_S.csv")#fully extrapolated for 40 yrs
write.csv(Survival, "P2.0.b10i10.mo_40_S.csv")#fully extrapolated for 40 yrs

write.csv(Survival, "P2.m2.b10i10.wk_40_S.csv")#fully extrapolated for 40 yrs
write.csv(Survival, "P2.m1.b10i10.wk_40_S.csv")#fully extrapolated for 40 yrs
write.csv(Survival, "P2.m0.5.b10i10.wk_40_S.csv")#fully extrapolated for 40 yrs
write.csv(Survival, "P2.0.b10i10.wk_40_S.csv")#fully extrapolated for 40 yrs

rm()