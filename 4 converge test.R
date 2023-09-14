#~~~~~~~~~~~~~~4 Convergence check~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##########
#FP#
#########
#see the convergence of d, alpha parameters
results2 <- bugs.object$sims.matrix[,grep("d",rownames(bugs.object$summary))]
results2 <- bugs.object$sims.matrix[,grep("alpha0",rownames(bugs.object$summary))]
results2 <- bugs.object$sims.matrix[,grep("alpha1",rownames(bugs.object$summary))]
results2 <- bugs.object$sims.matrix[,grep("alpha2",rownames(bugs.object$summary))]
#####
#RP#
#####
#check the convergence of beta and gamma
results2 <- bugs.object$sims.matrix[,grep("beta",rownames(bugs.object$summary))]
results2 <- bugs.object$sims.matrix[,grep("gamma",rownames(bugs.object$summary))]


#save file
results2 <- cbind(rep(0,dim(results2)[1]),results2)
summary(results2)
results_mcmc<-mcmc(results2)
#summary(results_mcmc)
#plot(results_mcmc)

#for d:2-7
#for beta:2-7
#for gamma:2-13
dev.off()
##################################
#1.Rhat
##################################
summary(bugs.object$summary)
##################################
#2. Check autocorrelation
##################################
#1st FP d
autocorr.plot(results_mcmc[,2:7])
#2nd FP d,
autocorr.plot(results_mcmc[,2:10])
#FP alpha
autocorr.plot(results_mcmc[,2:5])
#RP beta
autocorr.plot(results_mcmc[,2:7])
#RP gamma 1 knot
autocorr.plot(results_mcmc[,2:10])
#RP gamma 2 knot
autocorr.plot(results_mcmc[,2:13])
#RP gamma 3 knot
autocorr.plot(results_mcmc[,2:16])

#If the autocorrelation coefficients drop quickly to zero or stay within the confidence interval, it suggests good mixing and low autocorrelation, which are favorable signs of convergence.
##################################
# 3.Check trace for convergence
##################################
par(mfrow=c(2,2))#row x col
par(mfrow=c(2,3))#row x col
par(mfrow=c(3,3))#row x col
par(mfrow=c(4,3))
par(mfrow=c(5,3))

traceplot(results_mcmc[,2])
traceplot(results_mcmc[,3])
traceplot(results_mcmc[,4])
traceplot(results_mcmc[,5])
traceplot(results_mcmc[,6])
traceplot(results_mcmc[,7])

traceplot(results_mcmc[,8])
traceplot(results_mcmc[,9])
traceplot(results_mcmc[,10])
traceplot(results_mcmc[,11])
traceplot(results_mcmc[,12])
traceplot(results_mcmc[,13])

traceplot(results_mcmc[,14])
traceplot(results_mcmc[,15])
traceplot(results_mcmc[,16])

# Histograms of posterior distributions
densplot(results_mcmc[,2])
densplot(results_mcmc[,3])
densplot(results_mcmc[,4])
densplot(results_mcmc[,5])
densplot(results_mcmc[,6])
densplot(results_mcmc[,7])

densplot(results_mcmc[,8])
densplot(results_mcmc[,9])
densplot(results_mcmc[,10])
densplot(results_mcmc[,11])
densplot(results_mcmc[,12])
densplot(results_mcmc[,13])

densplot(results_mcmc[,14])
densplot(results_mcmc[,15])
densplot(results_mcmc[,16])