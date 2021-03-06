## Example using profile likelihood to compute the confidence limits for a
## return level plot. Presently supports gumbel, gev, and log-pearson 3 distribution

# This script has been tested by comparing its confidence intervals with the
# profile likelihood CIs that can be computed in the fExtremes package. See the
# muda_river.R script for more info. 

#######################################################
# PARAMETERS THAT THE USER WILL PROABALY WANT TO ADJUST
#######################################################
river_site='Muda River'
input_data=scan('Muda_discharge.txt') # Input data = vector of discharges
distribution='gev' #'gev'  # Set to 'gum', 'gev', or 'lp3'
alpha=0.4 # Adjustment factor in empirical AEPs. See Kuczera and Franks, Draft ARR
cilevel = 0.90 # Level of confidence intervals
flood_return=c(1.1,1.3,1.5,1.8,2,3,5,7,10,seq(15,100,by=5)) # Return levels


######################################################
# MAIN PART OF THE CODE
######################################################
Q_muda=sort(input_data, decreasing=T)


# Guess start parameters using L-moments
library(lmomco)
if(distribution=='lp3'){
    Muda_lmoms = lmom.ub(log(Q_muda))
    lmomfit=lmom2par(Muda_lmoms, 'pe3')
    # The lp3 is parameterised differently in lmomco and FAdist. This converts
    # between them
    tmp=as.numeric(lmomfit$para)
    startpars =list(x1=4/tmp[3]^2,
                    x2=1/(0.5*tmp[2]*abs(tmp[3])), 
                    x3=-(tmp[1]-2*tmp[2]/tmp[3]))
    
}else{
    Muda_lmoms = lmom.ub(Q_muda)
    lmomfit = lmom2par(Muda_lmoms, distribution)
    tmp=as.numeric(lmomfit$para)
    if(distribution=='gev'){
        # The gev is parameterised slightly differently in lmomco and FAdist
        startpars=list(x1=-tmp[3],x2=tmp[2],x3=tmp[1])
    }else{
        startpars=list(x1=tmp[2],x2=tmp[1])
    }
}


distribution_df = length(startpars) # Number of free parameters for the distribution
# Sanity check
if(distribution %in% c('gev','lp3')){
    if(distribution_df!=3){
        stop(paste('Error: Wrong number of starting parameters for', distribution))
    }
}else if(distribution %in% c('gum')){
    if(distribution_df!=2){
        stop(paste('Error: Wrong number of starting parameters for', distribution))
    }

}

# Number of data points
n = length(Q_muda)

# Compute the ranks
Q_rank = seq(1,n)

# Compute the estimate of the AEP for each discharge
Q_AEP_est  = (Q_rank - alpha)/(n + 1 - 2*alpha)


### STEP 1: Define a function which calculates the negative log likelihood of
### the distribution we want to fit, and the related quantile and probability
### functions
library(stats4)
library(FAdist)
source('lp3.R')
source('dgev.R')
source('profile_function.R')

gev_negloglik<-function(x1, x2, x3=NaN){
    # Compute the negative log likelihood of a range of distributions for the
    # Q_muda data. The meaning of x1, x2 and x3 depends on the particular
    # distribution chosen. However, two parameter distributions should only
    # refer to x1 and x2, and use the default x3 value.

    # 'distribution' is defined at the top of the script
    switch(distribution,
        gev= {
            -sum((FAdist::dgev(Q_muda, x1, x2, x3, log=TRUE)))
            #-sum((dgev(Q_muda, x1, x2, x3, log=TRUE)))
        },
        gum= {
            -sum((FAdist::dgumbel(Q_muda, x1, x2, log=TRUE)))
        },
        lp3 = {
            -sum((dlp3(Q_muda, x1, x2, x3, log=TRUE)))
        }
        )
}

gev_quantile<-function(data,x1,x2,x3=NaN){
    # Compute the quantile function for a range of distributions for the data
    # The meaning of x1, x2 and x3 depends on the particular
    # distribution chosen. However, two parameter distributions should only
    # refer to x1 and x2, and use the default x3 value.
    switch(distribution,
        gev= {
            FAdist::qgev(data, x1, x2, x3)
        },
        gum= {
            FAdist::qgumbel(data, x1, x2)
        },
        lp3 = {
            qlp3(data, x1, x2, x3)
        }
        )

}

gev_probability<-function(data,x1,x2,x3=NaN){
    # Compute the quantile function for a range of distributions for the data
    # The meaning of x1, x2 and x3 depends on the particular
    # distribution chosen. However, two parameter distributions should only
    # refer to x1 and x2, and use the default x3 value.
    switch(distribution,
        gev= {
            FAdist::pgev(data, x1, x2, x3)
        },
        gum= {
            FAdist::pgumbel(data, x1, x2)
        },
        lp3 = {
            plp3(data, x1, x2, x3)
        }
        )

}

### STEP 2: Compute the maximum likelihood estimate, and confidence limits on
### the parameters

# Set the names of the startpars as required by 'mle'
muda_startpars=startpars
if(length(muda_startpars)==3){
    names(muda_startpars)=c('x1','x2','x3')
    }else if(length(muda_startpars)==2){
        names(muda_startpars)=c('x1','x2')
        }

mleFit = mle(gev_negloglik, start=muda_startpars, nobs=n, method='Nelder-Mead')

### Compute confidence limits (also for return periods)
#if(FALSE){
#    # Make positive log likelihood function and quantile function with right interface
#    # , so we can use profile_function
#    if(distribution_df==3){
#         log_lik<-function(x) -gev_negloglik(x[1], x[2],x[3])
#         gev_qf<-function(x, p) gev_quantile(p,x[1],x[2],x[3])
#    }else if(distribution_df==2){
#         log_lik<-function(x) -gev_negloglik(x[1], x[2])
#         gev_qf<-function(x, p) gev_quantile(p,x[1],x[2])
#    }
#
#    conf_limits=matrix(NA,ncol=2,nrow=length(flood_return))
#    for(i in 1:length(flood_return)){
#        conf_limits[i,]=profile_function(var_to_profile<-function(x) gev_qf(x,p=1-1/(flood_return[i])),
#                                         log_lik=log_lik,
#                                         maxlikpar=mleFit@coef,
#                                         level=cilevel,
#                                         method='Nelder-Mead')
#    }
#}else{

    ### Compute confidence limits
    ## Use the asymptotic ci's from mle
    tmp = coef(summary(mleFit))
    ## According to Bolker, the parameter estimates are asymptotically normally distributed
    ## Here, we take confidence limits that are deliberately too large. Later
    ## we will search through this region for parameter values which are within
    ## the asymptotic profile likelihood confidence limits
    zvalue = qnorm(1-(1-cilevel)/5)  
    ## Note -- in theory, zvalue=qnorm(1-(1-cilevel)/2) should be enough. But as
    ## the latter result is asymptotic, I am trying to be conservative by
    ## dividing by 3 instead
    ## FIXME: Should put a check in the code to ensure that the boundary of this
    ## search region really is outside the profile likelihood confidence
    ## intervals
    ci.x = cbind(tmp[,1]-zvalue*tmp[,2], tmp[,1]+zvalue*tmp[,2])
    #
    #
    #
    #### Step 3: Numerically compute the cilevel confidence intervals for a range
    #### of flood return periods, using a 'brute-force' method.
    #### METHOD: Take a box around the parameter values contained in ci.x, and
    #### search through this. Not all points in this 'box' will be
    #### inside the cilevel confidence interval. We record the likelihood value of all
    #### points that we search, and later determine confidence intervals only for
    #### parameter values inside the cilevel confidence limits 
    #
    nn=60 # We divide the ci 'box' into n^3 values for searching
    storeme = matrix(NA,ncol=length(flood_return)+1,nrow=nn^distribution_df) # Store the confidence limits
    countme=0 # Used for counting in the loop
    ijk=matrix(NA,ncol=3,nrow=nn^distribution_df)
    #
    ## Begin search
    if(distribution_df==3){

        for(i in seq(ci.x[1,1], ci.x[1,2],len=nn)){
            for(j in seq(ci.x[2,1], ci.x[2,2],len=nn)){
                for(k in seq(ci.x[3,1], ci.x[3,2],len=nn)){
                    countme=countme+1
                    ijk[countme,]=c(i,j,k)  
                    storeme[countme,1:length(flood_return)] = gev_quantile(1-1/flood_return, i, j,k)
                    storeme[countme,length(flood_return)+1] = gev_negloglik(i,j,k)
                }
            }
        }
    }else if(distribution_df==2){

        for(i in seq(ci.x[1,1], ci.x[1,2],len=nn)){
            for(j in seq(ci.x[2,1], ci.x[2,2],len=nn)){
                countme=countme+1  
                ijk[countme,]=c(i,j)  
                storeme[countme,1:length(flood_return)] = gev_quantile(1-1/flood_return, i, j)
                storeme[countme,length(flood_return)+1] = gev_negloglik(i,j)
            }
        }



    }
    # End search

    # Examples of computing confidence intervals for the return period can be found
    # in gevrlevelPlot
    # Compute indicies of points which are inside the confidence limit
    tmp = which(storeme[,length(flood_return)+1]< 
                gev_negloglik(coef(mleFit)[1], coef(mleFit)[2], coef(mleFit)[3]) + qchisq(cilevel,1)/2 )
    conf_limits1= apply(storeme[tmp,1:length(flood_return)], 2, range) # The max and min = confidence limit
    conf_limits1=t(conf_limits1)
    # Find parameter values associated with these confidence limits
    conf_par_min=matrix(NA,ncol=distribution_df,nrow=length(flood_return))
    conf_par_max=conf_par_min
    for(i in 1:length(flood_return)){
        conf_par_min[i,]=ijk[tmp[which.min(storeme[tmp,i])],]
        conf_par_max[i,]=ijk[tmp[which.max(storeme[tmp,i])],]
    }

    # Make positive log likelihood function and quantile function with right interface
    # , so we can use profile_function
    if(distribution_df==3){
         log_lik<-function(x) -gev_negloglik(x[1], x[2],x[3])
         gev_qf<-function(x, p) gev_quantile(p,x[1],x[2],x[3])
    }else if(distribution_df==2){
         log_lik<-function(x) -gev_negloglik(x[1], x[2])
         gev_qf<-function(x, p) gev_quantile(p,x[1],x[2])
    }

    conf_limits=matrix(NA,ncol=2,nrow=length(flood_return))
    for(i in 1:length(flood_return)){
        conf_limits[i,]=profile_function(var_to_profile<-function(x) gev_qf(x,p=1-1/(flood_return[i])),
                                         log_lik=log_lik,
                                         maxlikpar=mleFit@coef,
                                         searchStart=rbind(conf_par_min[i,],conf_par_max[i,]),
                                         level=cilevel,
                                         method='BFGS')
    }



# Make the plot
pdf(file=paste('Flood_frequency_plot_', distribution,'_maxLike.pdf',sep=""), width=7,height=5)

fitted_model = gev_quantile(1-1/flood_return,coef(mleFit)[1],coef(mleFit)[2],coef(mleFit)[3])
plot(flood_return, fitted_model ,
     log='x',t='l', ylim=c(min(conf_limits),max(conf_limits)), xlab='AEP of 1/Y Years', ylab='Discharge (m^3/s)',
     main=river_site,cex.main=1.5)
points(flood_return,conf_limits[,1],t='l',col=2,lty='dashed')
points(flood_return,conf_limits[,2],t='l',col=2,lty='dashed')
points(1/Q_AEP_est, Q_muda,col='steelblue',pch=19)
grid(nx=10,ny=10)
legend('topleft', c(paste('Fitted curve (', distribution, ', maxlikelihood)', sep=""), 
       paste(cilevel*100, '% Confidence Limits (Profile Likelihood)',sep=""), 'Data'), 
       lty=c('solid','dashed',NA),col=c('black', 'red', 'steelblue'), pch=c(NA,NA,19), 
       bty='o',bg='white')
dev.off()

## Write data for later plotting / investigation
write.table(cbind(flood_return,fitted_model, conf_limits[,1], conf_limits[,2]), 
            file=paste('Fitted_',distribution,'_', river_site, '_proflike.txt', sep=''), 
            col.names=c('AEP of 1/Y', 'Model (Max Likelihood)', paste('lower ci',cilevel), 
            paste('upper ci',cilevel)), row.names=FALSE, sep="," )

# Quantile-quantile plot
theoretical_quantiles=gev_quantile(1-Q_AEP_est, coef(mleFit)[1], coef(mleFit)[2], coef(mleFit)[3])
pdf(file=paste('Quantile_plot_',distribution,'_maxLike.pdf',sep=""),width=8,height=6)
plot(Q_muda, theoretical_quantiles, xlab='Measured Discharge (m^3/s)', ylab='Theoretical Discharge (m^3/s)', main=river_site)
abline(0,1)
dev.off()

# Probability plot
theoretical_probs=gev_probability(Q_muda, coef(mleFit)[1], coef(mleFit)[2], coef(mleFit)[3])
pdf(file=paste('Probability_plot_',distribution,'_maxLike.pdf',sep=""),width=8,height=6)
plot(Q_AEP_est, 1-theoretical_probs, xlab='Empirical AEP', ylab='Theoretical AEP', main=river_site)
abline(0,1)
dev.off()


###############################
# Try with MCMC, but still using profile-likelihood reasoning
log_lik_prt<-function(x){
     out=log_lik(x)
     if(is.nan(out)|!is.finite(out)) out=-Inf
     return(out)
 }
library(MCMCpack)
xxx=MCMCmetrop1R(log_lik_prt,mleFit@coef,burnin=10000,mcmc=1e+06)
# Compute log-likelihood for each
xxxLL=apply(as.matrix(xxx),1,log_lik_prt)
# Compute q100 stats
xxQ100=apply(as.matrix(xxx),1, f<-function(x) gev_qf(x,p=0.99))
#xxQ100=qgev(0.99,xxx[,1],xxx[,2],xxx[,3])
# Find where LL is in the acceptable zone
xxTmp=which(xxxLL>log_lik_prt(mleFit@coef)-qchisq(cilevel,1)/2)
summary(xxQ100[xxTmp])
maxPars=xxx[xxTmp[which.max(xxQ100[xxTmp])],]
minPars=xxx[xxTmp[which.min(xxQ100[xxTmp])],]
# This may suggest that the profile likelihood search tends to miss the most
# extreme regions that are still within the 'acceptable zone'

# Try to use these as starting values for prof likelihood -- seems a legit method, could work well
i=length(flood_return)
profile_function(var_to_profile<-function(x) gev_qf(x,p=1-1/(flood_return[i])),
                                         log_lik=log_lik,
                                         maxlikpar=mleFit@coef,
                                         searchStart=rbind(minPars,maxPars),
                                         level=cilevel,
                                         method='Nelder-Mead')

