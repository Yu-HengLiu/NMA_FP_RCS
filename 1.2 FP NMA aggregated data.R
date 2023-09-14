
#~~~~~~~~~~~~~~~~~~~~~~1.3 aggregated IPD for FP NMA~~~~~~~~~~~~~~~~~~~~~#
rm(list = ls()) 
setwd("your working directoryy/R")


# Create function that takes a generated dataset and formats the data ready to apply the anova
# parameterisation
#timepoints, timepoints2, ref.study=1, df will be defined
anova_data <- function(timepoints, timepoints2, ref.study=1, df){
  
  # Split the data at timepoints
  df2 <- survSplit(Surv(time, status) ~., data=df,
                   cut=timepoints, episode ="timegroup")
                  #cut = timepoints is where you provide the vector of timepoints that will be used to define the intervals for splitting the survival data.
                  #episode = "timegroup" indicates that you want to use the "timegroup" variable in your data frame to define episodes for grouping the split survival data.
  
  # Calculate offset
  df2$y <- df2$time - df2$tstart   #character string with the name of a start time variable (will be created if needed)
  
  # Add a variable that equals one for all patients - this is so the number at risk
  # can be calculated when we collapse the data
  df2$n <- 1
  
  # Collapse data
  df3 <- summaryBy(y + status + n ~ timegroup + treatment + studyCode, FUN=c(sum, max), data=df2)
  df3 <- subset(df3, select=-c(status.max, n.max))
  names(df3) <- c("spgrp", "treatment", "trialid", "y", "nevents", "natrisk", "y.max")
  #spgrp:group depending on split time
  
  # Add in a start time variable
  df3$start <- NA
  for(i in unique(df3$spgrp)){
    df3$start[df3$spgrp==i] <- timepoints2[i]
  }
  
  # Add in a time variable (i.e. how long since time 0 to max value of y for each row)
  df3$time <- df3$start + df3$y.max
  
  # Return the formatted dataset
  return(df3)
  
}


# Format data for use with fractional polynomial models

library(survival)
library(doBy)
NMA_data <- read.csv("NMA_data.csv")
# Start with IPD data and convert to ANOVA format by choosing time points for aggregating over

# Add treatment as a factor variable
NMA_data$treatment[NMA_data$txCode==1] <- "sun"
NMA_data$treatment[NMA_data$txCode==2] <- "pem+len"
NMA_data$treatment[NMA_data$txCode==3] <- "cab"
NMA_data$treatment[NMA_data$txCode==4] <- "niv+ipi"

NMA_data$treatment <- as.factor(NMA_data$treatment)
# Need to create treatment as a factor variable before applying this function

# Select time points for aggregating data
timepoints=c(3,6, 12, 120)#months
# Time points including zero
timepoints2=c(0, 3, 6, 12,120)
# Apply function
anova <- data.frame(spgrp=NA, treatment=NA, trialid=NA, y=NA, nevents=NA,
                    natrisk=NA, y.max=NA, start=NA, time=NA, trt=NA, treatnumf=NA,
                    studynumf=NA)

anova <- anova_data(timepoints=timepoints, timepoints2=timepoints2, ref.study=1, 
                    df=NMA_data)

# Save as csv file
write.csv(anova, "aRCC_aggregate_data.csv")
