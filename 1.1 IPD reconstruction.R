
#~~~~~~~~~~~~~~~~~~~~~~1-1.IPD reconstruction~~~~~~~~~~~~~~~~~~~~~#
#PFS
setwd("your working directory/R")

install.packages("IPDfromKM")
library(IPDfromKM)

################################################
#Extract data points from published K-M curves.
################################################
#PFS
  #pem plus len vs. sun
  predata_PL_PFS <- getpoints("your working directory/PFS_PL.jpg",x1=0,x2=40,y1=0,y2=1)
  predata_sun_PL_PFS <- getpoints("your working directory/PFS_PL.jpg",x1=0,x2=40,y1=0,y2=1)
  colnames(predata_PL_PFS) <- c("survTime","survProb")
  colnames(predata_sun_PL_PFS) <- c("survTime","survProb")
  
  #cab vs. sun
  predata_cab_PFS <- getpoints("your working directory/PFS_C.jpg",x1=0,x2=30,y1=0,y2=1)
  predata_sun_cab_PFS <- getpoints("your working directory/PFS_C.jpg",x1=0,x2=30,y1=0,y2=1)
  colnames(predata_cab_PFS) <- c("survTime","survProb")
  colnames(predata_sun_cab_PFS) <- c("survTime","survProb")
 
  # niv plus ipi vs. sun
  predata_NI_PFS <- getpoints("your working directory/PFS_NI.jpg",x1=0,x2=51,y1=0,y2=1)
  predata_sun_NI_PFS <- getpoints("your working directory/PFS_NI.jpg",x1=0,x2=51,y1=0,y2=1)
  colnames(predata_NI_PFS) <- c("survTime","survProb")
  colnames(predata_sun_NI_PFS) <- c("survTime","survProb")
  
##############################
#Process extracted data points.
##############################
  #pem plus len vs. sun
  trisk <- seq(from = 0, to = 40, by = 2)
  nrisk<-c(355,321,300,276,259,235,213,186,160,136,126,106,80,56,30,14,6,3,1,1,0)
  process_PL_PFS <-preprocess (dat=predata_PL_PFS, trisk=trisk, nrisk=nrisk, maxy=1)
  nrisk<-c(357,262,218,145,124,107,85,69,62,49,42,32,25,16,9,3,2,1,0)
  process_sun_PL_PFS <-preprocess (dat=predata_sun_PL_PFS, trisk=trisk, nrisk=nrisk, maxy=1)
  
  #cab vs. sun
  trisk <- seq(from = 0, to = 30, by = 3)
  nrisk<-c(79,51,37,24,22,18,12,5,2,1,0)
  process_cab_PFS <-preprocess (dat=predata_cab_PFS, trisk=trisk, nrisk=nrisk, maxy=1)
  nrisk<-c(78,36,21,12,9,5,3,2,1,0,0)
  process_sun_cab_PFS <-preprocess (dat=predata_sun_cab_PFS, trisk=trisk, nrisk=nrisk, maxy=1)
  
  # niv plus ipi vs. sun
  trisk <- seq(from = 0, to = 51, by = 3)
  nrisk<-c(425,304,233,187,164,149,129,116,99,96,91,83,76,71,56,34,13,2)
  process_NI_PFS <-preprocess (dat=predata_NI_PFS, trisk=trisk, nrisk=nrisk, maxy=1)
  nrisk<-c(422,281,189,139,107,90,75,61,47,37,31,26,23,18,11,8,3,0)
  process_sun_NI_PFS <-preprocess (dat=predata_sun_NI_PFS, trisk=trisk, nrisk=nrisk, maxy=1)
  
  
##################################################
#Reconstruct IPD using the modified-iKM algorithm.
##################################################
  #pem plus len vs. sun
  est_PL_PFS <- getIPD(prep= process_PL_PFS,armID=1,tot.events=NULL)
  est_sun_PL_PFS <- getIPD(prep= process_sun_PL_PFS,armID=0,tot.events=NULL)
  #cab vs. sun
  est_cab_PFS <- getIPD(prep= process_cab_PFS,armID=1,tot.events=43)
  est_sun_cab_PFS <- getIPD(prep= process_sun_cab_PFS,armID=0,tot.events=49)
  # niv plus ipi vs. sun
  est_NI_PFS <- getIPD(prep= process_NI_PFS,armID=1,tot.events=NULL)
  est_sun_NI_PFS <- getIPD(prep= process_sun_NI_PFS,armID=0,tot.events=NULL)
  
###########################################
#Assess the accuracy of the reconstruction.
###########################################
# the summary() produce statistics to assess accuracy
  # what shows the reconstructed IPD is accurate
    #  The small values for RMSE (< 0.05)
    # mean absolute error (< 0.02)
    # max absolute error (< 0.05)
    # large p-value of the Kolmogorov-Smirnov test 
#The output using the plot() shows
  #(1) K-M curves based on reconstructed IPD and read-in data points, respectively
  #(2) number of patients at risk using the reconstructed IPD versus reported
  #(3) difference in survival probability at different time points for the reconstructed IPD and the read-in data.

#PFS
  #pem plus len vs. sun
  summary(est_PL_PFS)
  plot(est_PL_PFS)
  summary(est_sun_PL_PFS)
  plot(est_sun_PL_PFS)
 
  #cab vs. sun
  summary(est_cab_PFS)
  plot(est_cab_PFS)
  summary(est_sun_cab_PFS)
  plot(est_sun_cab_PFS)
 
  # niv plus ipi vs. sun
  summary(est_NI_PFS)
  plot(est_NI_PFS)
  summary(est_sun_NI_PFS)
  plot(est_sun_NI_PFS)

##################
#Secondary analysis
##################
  #we have separated IPD dataset for intervention(arm2) and comparator(arm1)
  #can use the function "survreport" to generate survival probability plot and cumulative hazard plot
  #can print the result to generate point estimate and 95% CI for survival time (e.g. median survival time at which 50% of pts survive)
  #or can use the IPD data to do other analysis (if perform other analysis might need to stack IPD1 and IPD2)
  
  #pem plus len vs. sun
  report_PL_sun_PFS <- survreport(ipd1=est_sun_PL_PFS $IPD,ipd2=est_PL_PFS$IPD, arms=2,interval=6,s=c(0.50,0.75,0.95))
  print(report_PL_sun_PFS$arm1)
  #cab vs. sun
  report_cab_sun_PFS <- survreport(ipd1=est_sun_cab_PFS $IPD,ipd2=est_cab_PFS$IPD, arms=2,interval=6,s=c(0.50,0.75,0.95))
  print(report_cab_sun_PFS$arm1)
  print(report_cab_sun_PFS$arm2)
  # niv plus ipi vs. sun
  report_NI_sun_PFS <- survreport(ipd1=est_sun_NI_PFS $IPD,ipd2=est_NI_PFS$IPD, arms=2,interval=6,s=c(0.50,0.75,0.95))
  print(report_NI_sun_PFS$arm1)

###################
#Generate IPD list
##################

#merge the dataset first: include each arms and each clinical trial in one dataset
CLEAR_PL_PFS<-est_PL_PFS$IPD
CLEAR_sun_PFS<-est_sun_PL_PFS$IPD
CABOSUN_cab_PFS<-est_cab_PFS$IPD
CABOSUN_sun_PFS<-est_sun_cab_PFS$IPD
CheckMate214_NI_PFS<-est_NI_PFS$IPD
CheckMate214_sun_PFS<-est_sun_NI_PFS$IPD
#Generate Row number
CLEAR_PL_PFS$patid <- seq.int(nrow(CLEAR_PL_PFS)) 
CLEAR_sun_PFS$patid <- seq.int(nrow(CLEAR_sun_PFS)) 
CABOSUN_cab_PFS$patid <- seq.int(nrow(CABOSUN_cab_PFS)) 
CABOSUN_sun_PFS$patid <- seq.int(nrow(CABOSUN_sun_PFS)) 
CheckMate214_NI_PFS$patid <- seq.int(nrow(CheckMate214_NI_PFS)) 
CheckMate214_sun_PFS$patid <- seq.int(nrow(CheckMate214_sun_PFS)) 
#generate study var
CLEAR_PL_PFS$study <- "CLEAR"
CLEAR_sun_PFS$study <- "CLEAR"
CABOSUN_cab_PFS$study <- "CABOSUN"
CABOSUN_sun_PFS$study <- "CABOSUN"
CheckMate214_NI_PFS$study <- "CheckMate214"
CheckMate214_sun_PFS$study <- "CheckMate214"
#generate study code
CLEAR_PL_PFS$studyCode <-1
CLEAR_sun_PFS$studyCode <- 1
CABOSUN_cab_PFS$studyCode <- 2
CABOSUN_sun_PFS$studyCode <- 2
CheckMate214_NI_PFS$studyCode <- 3
CheckMate214_sun_PFS$studyCode <- 3
#generate tx Code
CLEAR_PL_PFS$txCode <-2
CLEAR_sun_PFS$txCode <-1
CABOSUN_cab_PFS$txCode <- 3
CABOSUN_sun_PFS$txCode <- 1
CheckMate214_NI_PFS$txCode <- 4
CheckMate214_sun_PFS$txCode <- 1

#combine all dataset
NMA_data<-rbind(CLEAR_PL_PFS,CLEAR_sun_PFS,CABOSUN_cab_PFS,CABOSUN_sun_PFS,CheckMate214_NI_PFS,CheckMate214_sun_PFS)
write.csv(NMA_data, "NMA_data.csv", row.names = TRUE)
write.csv(CLEAR_sun_PFS, "sun_CLEAR.csv", row.names = TRUE)
