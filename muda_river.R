## Read in the data (Table 1 of Julien et al, 2010)
muda=read.table('muda_river.txt')

# Extract the ranked discharge
Q_muda = c(muda[,3], muda[,6])

# Number of data points
n = length(Q_muda)

# Compute the ranks
Q_rank = seq(1,n)

# Compute the estimate of the AEP for each discharge
alpha = 0.4 # Chosen based on Kuczera and Franks, Draft ARR
Q_AEP_est  = (Q_rank - alpha)/(n + 1 - 2*alpha)

# Compute the estimated ARI
Q_ARI_est = -1/(log(1-Q_AEP_est))

# Make a probability plot, then save it to a pdf
plot(Q_AEP_est, Q_muda, xlab='AEP estimate', ylab='Discharge (m^3/s)', main = 'Muda River Probability Plot', pch=19, log='x')
dev.copy(pdf, file='muda_discharge_aep.pdf')
dev.off()

# Make a recurrence interval plot, then save it to a pdf
plot(Q_ARI_est, Q_muda, xlab='ARI estimate', ylab='Discharge (m^3/s)', main = 'Muda River Recurrence Interval Plot', pch=19, log='x')
dev.copy(pdf, file='muda_discharge_ari.pdf')
dev.off()

###########################################################################
#
# Try fitting with maximum likelihood
#
# Note that the confidence intervals in the 'ismev' package are simple
# asymptotic approximations, which are probably not very accurate (and from
# experience, are probably too narrow).
#
###########################################################################

# Import R library for fitting GEV, Gumbel using Maximum Liklihood methods
library(ismev)

# Fit the gumbel distribution
muda_gum=gum.fit(Q_muda) # muda_gum describes the fitted distribution

# Look at diagnostic plot.
gum.diag(muda_gum)

# Save the diagnostic plot using a more complex (and more flexible) method than
# we used above
pdf(file='muda_gum_mle_diagnostic.pdf', width=8,height=6)
gum.diag(muda_gum)
dev.off()

## Note -- to fit a GEV, use similar commands, with 'gev' instead of 'gum'
muda_gev=gev.fit(Q_muda) # muda_gum describes the fitted distribution

# Look at diagnostic plot.
gev.diag(muda_gev)

gev.prof(muda_gev, 50, xlow=800, xup=3000)
# And so on ...

############################################################################
library(MCMCpack)
LogLik<-function(par) {
    # Return positive log likelihood
    # Get rid of NA's, since MCMCpack doesn't work with them 
    out=sum(log(dgev(Q_muda, par[1],par[2], par[3])))+ # Add uniform prior
        log(par[1]> -10 & par[1] < 10) + 
        log(par[2]<10000 & par[2] > 0) + 
        log(par[3] <10000 & par[3] > 0)
    if(is.na(out)) out=-Inf
    return(out)
}
bayesFit=MCMCmetrop1R(LogLik, theta.init=c(0.1, 450, 160), 
                      thin=10, mcmc=1000000, tune=2)
# Run again for Gelman diagnostic
bayesFit2=MCMCmetrop1R(LogLik, theta.init=c(0.1, 451, 161), 
                      thin=10, mcmc=1000000, tune=2,seed=1)
bayesFit3=MCMCmetrop1R(LogLik, theta.init=c(0.12, 453, 163), 
                      thin=10, mcmc=1000000, tune=2)
bdraws=mcmc.list(list(bayesFit,bayesFit2, bayesFit3))
# Check convergence
library(coda)
geweke.diag(bayesFit) # z score comparing means in first 10% and last 50% Shouldn't be a large z score
cumuplot(bayesFit) # Plot of running mean + spread of each variable
autocorr.plot(bayesFit)
gelman.diag(bdraws) # Estimate scale reduction factor
heidel.diag(bayesFit) # Null hypothesis of stationary chains
# Compute log likelihood
ll=apply(as.matrix(bayesFit), 1,LogLik)
# Get q50 for each MCMC sample
q50=qgev(0.98,bayesFit[,1], bayesFit[,2],bayesFit[,3])
quantile(c(q50), c(0.025, 0.975))
# What about Maximum likelihood principles?
mlConfR=which(ll>max(ll)-qchisq(0.95,1)/2)
range(q50[mlConfR])
## Try highest density interval
HPDinterval(q50)
HPDinterval(qgev(0.98, bayesFit2[,1], bayesFit2[,2], bayesFit2[,3]))
HPDinterval(qgev(0.98, bayesFit3[,1], bayesFit3[,2], bayesFit3[,3]))
q100=qgev(0.99,bayesFit[,1], bayesFit[,2],bayesFit[,3])
HPDinterval(q100)
HPDinterval(qgev(0.99, bayesFit2[,1], bayesFit2[,2], bayesFit2[,3]))
HPDinterval(qgev(0.99, bayesFit3[,1], bayesFit3[,2], bayesFit3[,3]))
# We can see the need for a longer MCMC run here

########################################################################
# The 'fExtremes' package provides the option for a more realistic 'profile
# likelihood' confidence interval plot for a given return level.
library(fExtremes)
muda_gev2 = gevFit(Q_muda)

par(mfrow=c(2,2))
# Diagnostic plot
summary(muda_gev2)

# We can do a profile-likelihood plot for a single return level
gevrlevelPlot(muda_gev2,100,ci=0.9)

# It's nice to have the complete return level plot -- one way to do this is
# with a loop to calculate the return levels for enough points to do the plot. 
# A better way is coded in the prof_likelihood.R script. The following code
# has been used to test the latter script.
#
#return_levels=seq(2.0,100,by=4)
#store_CIs = matrix(NA,ncol=3,nrow=length(return_levels)) # Use this to store output
#for(i in 1:length(return_levels)){
#    # Compute confidence intervals
#    x = gevrlevelPlot(muda_gev2,return_levels[i], ci=0.9) # 90% Confidence intervals -- there can be problems with these
#    # Store the confidence intervals and parameter estimate for the return_levels[i] event
#    store_CIs[i,]=as.numeric(x)[1:3]
#}
## Plot the result
#plot(return_levels,store_CIs[,2],ylim=range(store_CIs[,c(1,3)]),t='o',
#     xlab='AEP of 1/Y',ylab='Discharge', log='x')
#points(return_levels,store_CIs[,1],col=2,t='l',lty='dotted')
#points(return_levels,store_CIs[,3],col=2,t='l',lty='dotted')
##points(Q_ARI_est,Q_muda,col='blue',pch=19)
#points(1/Q_AEP_est,Q_muda,col='blue',pch=19)

###########################################
#
# Try L-moments for the Gumbel
#
###########################################
write.table(Q_muda,file='Muda_discharge.txt', row.names=F, col.names=F)
source('river_freq_lmom.R') #This script does L-moments.

# To fit a gev, change the 'distribution' parameter to 'gev' in the
# 'river_freq_lmom.R' file, and then source('river_freq_lmom.R')


## FURTHER EXPERIMENTATION BELOW


# Import R library for fitting distributions with L-moments
#library(lmomco) 
#
#lmoms = lmom.ub(Q_muda) # Compute L-moments of the data
#
## Compute the parameters of the Gumbel distribution that have the same
## L-moments as the data
#gum_consts_lmom = pargum(lmoms) 
## Compare with Maximum-liklihood fit
#print(muda_gum$mle)
#print(gum_consts_lmom)
#
## Comparitive plot of the two gumbel distributions and the data.
#
#temp_points = seq(min(Q_muda)/3, max(Q_muda)*1.2,len=100) # We use these 'discharge' points to plot the curve
#
## Calculate aep at points x. We do this via the cumulative distribution
## function, noting that AEP = 1-cdf
#gum_aep_lmom = 1-cdfgum(temp_points, gum_consts_lmom) 
#
#
## To calculate the aep for the maximum liklihood parameters, we define an
## object similar to gum_consts_lmom, and change the parameter values to the
## maximum liklihood values
#gum_consts_mle = gum_consts_lmom
#gum_consts_mle$para = muda_gum$mle
#
#gum_aep_mle = 1-cdfgum(temp_points, gum_consts_mle) # Calculate cumulative distribution function at points x. Note: AEP = 1-cdf
#
## Convert AEP to ARI
#gum_ari_lmom = -1/log(1-gum_aep_lmom)
#gum_ari_mle = -1/log(1-gum_aep_mle)
#
## Plot the L-moments fit and the mle fit, and the data
#pdf('muda_gumbel_lmom_mle.pdf', width=8,height=6)
#plot(gum_ari_lmom,temp_points, t='l', log='x', xlab='Average Recurrence Interval (yrs)',ylab='Discharge (m^3/s)')
#points(gum_ari_mle,temp_points, t='l',col=2)
#points(Q_ARI_est,Q_muda, pch=19)
#legend('topleft', 
#       c('Gumbel (L-moments)', 'Gumbel (Maximum Likelihood)', 'Muda River Data'), 
#       lty=c(1,1, 0), col=c(1,2,1), pch=c(NA, NA, 19))
#dev.off()
#
################################################################
#
# Try fitting Generalised Extreme Value with Maximum likelihood.
#
################################################################
#muda_gev=gev.fit(Q_muda)
#gev.diag(muda_gev)




