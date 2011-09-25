## Example using profile likelihood to compute the confidence limits for a
## return level plot.

# This script has been tested by comparing its confidence intervals with the
# profile likelihood CIs that can be computed in the fExtremes package. See the
# muda_river.R script for more info. Note that fExtremes assumes that the
# likelihood is chisquare distributed with one degree of freedom, whereas in
# this script, we follow Bolker (20XX, Ecological models and Data in R) and
# assume that the likelihood is chi-square distributed with degrees of freedom
# = number of free parameters in the model (e.g. 3 for a gev, 2 for a normal,
# etc).


#######################################################
# PARAMETERS THAT THE USER WILL PROABALY WANT TO ADJUST
#######################################################
river_site='Muda River'
input_data=scan('Muda_discharge.txt') # Input data = vector of discharges
distribution='lp3'
alpha=0.4 # Adjustment factor in empirical AEPs. See Kuczera and Franks, Draft ARR
cilevel = 0.90 # Level of confidence intervals


profile_cis=FALSE
# profile_cis: Flag to use profile likelihood for parameter confidence limits
# (otherwise just invert hessian). This should not matter unless there is a
# problem with the confidence limits. However, TRUE is more sensitive to
# numerical problems, so FALSE should probably be used always

startpars=list(9.0,7.0, 1.0) # First guess of parameters for distribution
#startpars=list(9.0,7.0) # First guess of parameters for distribution



######################################################
# MAIN PART OF THE CODE
######################################################
Q_muda=sort(input_data, decreasing=T)
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
### the distribution we want to fit
#library(fExtremes) # This has a dgev function (gev density)
library(stats4)
library(FAdist)


gev_negloglik<-function(x1, x2, x3=NaN){
    # Compute the negative log likelihood of a range of distributions for the
    # Q_muda data. The meaning of x1, x2 and x3 depends on the particular
    # distribution chosen. However, two parameter distributions should only
    # refer to x1 and x2, and use the default x3 value.

    # 'distribution' is defined at the top of the script
    switch(distribution,
        gev= {
            -sum((dgev(Q_muda, x1, x2, x3, log=TRUE)))
        },
        gum= {
            -sum((dgumbel(Q_muda, x1, x2, log=TRUE)))
        },
        lp3 = {
            -sum((dlgamma3(Q_muda, x1, x2, x3, log=TRUE)))
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
            qgev(data, x1, x2, x3)
        },
        gum= {
            qgumbel(data, x1, x2)
        },
        lp3 = {
            qlgamma3(data, x1, x2, x3)
        }
        )

}

### Experimental effort to generalise this function
# library(lmomco)
# distribution='gev'
# gev_negloglik<-function(...){
#    # Note that I am taking the log myself, because the 'log=TRUE' argument
#    # seems to fail in this function
#    xxx$type=distribution
#    xxx$para=c(...)
#    -sum(log(dlmomco(Q_muda, xxx)))
# }
# muda_lmoms=lmom.ub(Q_muda)
# source('mle3.R')
# x = mle3(gev_negloglik, 
#          start=as.list(lmom2par(muda_lmoms,distribution)$para),
#          nobs=n,
#          method="Nelder-Mead")

### STEP 2: Compute the maximum likelihood estimate, and confidence limits on
### the parameters

#startpars=list(xi=0., mu=400,beta=150)
#startpars=list(x1=0., x2=400,x3=150)
#muda_lmoms=lmom.ub(Q_muda)
#muda_startpars = as.list(lmom2par(muda_lmoms, distribution)$para)
# Set the names of the startpars as required by 'mle'
muda_startpars=startpars
if(length(muda_startpars)==3){
    names(muda_startpars)=c('x1','x2','x3')
    }else if(length(muda_startpars)==2){
        names(muda_startpars)=c('x1','x2')
        }

x = mle(gev_negloglik, start=muda_startpars, nobs=n, method='Nelder-Mead')


## Compute confidence limits
if(profile_cis){
    # Compute profile likelihood confidence limits. This can be numerically
    # hard for some cases, in which case errors will occur
    ci.x = confint(x, level=cilevel)
}else{
    # Use the asymptotic ci's from mle
    tmp = coef(summary(x))
    # According to Bolker, the parameter estimates are normally distributed
    # Here, we take confidence limits that are deliberately too large. Later
    # we will search through this region for parameter values which are within
    # the asymptotic profile likelihood confidence limits
    zvalue = qnorm(1-(1-cilevel)/3)  
    # Note -- in theory, zvalue=qnorm(1-(1-cilevel)/2) should be enough. But as
    # the latter result is asymptotic, I am trying to be conservative by
    # dividing by 3 instead
    ci.x = cbind(tmp[,1]-zvalue*tmp[,2], tmp[,1]+zvalue*tmp[,2])
}


### Step 3: Numerically compute the cilevel confidence intervals for a range
### of flood return periods, using a 'brute-force' method.
### METHOD: Take a box around the parameter values contained in ci.x, and
### search through this. Not all points in this 'box' will be
### inside the cilevel confidence interval. We record the likelihood value of all
### points that we search, and later determine confidence intervals only for
### parameter values inside the cilevel confidence limits 

flood_return=c(1.1,1.3,1.5,1.8,2,5,7,10,15,20,30, 50, 70, 90,100)
nn=60 # We divide the ci 'box' into n^3 values for searching
storeme = matrix(NA,ncol=length(flood_return)+1,nrow=nn^distribution_df) # Store the confidence limits
countme=0 # Used for counting in the loop
#loglik_thresh=gev_negloglik(coef(x)[1],coef(x)[2],coef(x)[3])

if(distribution_df==3){

    for(i in seq(ci.x[1,1], ci.x[1,2],len=nn)){
        for(j in seq(ci.x[2,1], ci.x[2,2],len=nn)){
            for(k in seq(ci.x[3,1], ci.x[3,2],len=nn)){
                countme=countme+1  
                storeme[countme,1:length(flood_return)] = gev_quantile(1-1/flood_return, i, j,k)
                storeme[countme,length(flood_return)+1] = gev_negloglik(i,j,k)
            }
        }
    }
}else if(distribution_df==2){

    for(i in seq(ci.x[1,1], ci.x[1,2],len=nn)){
        for(j in seq(ci.x[2,1], ci.x[2,2],len=nn)){
            countme=countme+1  
            storeme[countme,1:length(flood_return)] = gev_quantile(1-1/flood_return, i, j)
            storeme[countme,length(flood_return)+1] = gev_negloglik(i,j)
        }
    }



}
# Examples of computing confidence intervals for the return period can be found
# in gevrlevelPlot
# Note that according to Bolker, the -2loglikelihood is approximately chisq
# distributed with df = number of free parameters (3 in the case of a gev).
# This means that we search for all likelihoods within qchisq(cilevel,num_free_parameters)/2 of
# the maximum likelihood. On the other hand, the fExtremes package always uses df=1

# Compute indicies of points which are inside the confidence limit
tmp = which(storeme[,length(flood_return)+1]< 
            gev_negloglik(coef(x)[1], coef(x)[2], coef(x)[3]) + qchisq(cilevel,distribution_df)/2 )
conf_limits= apply(storeme[tmp,1:length(flood_return)], 2, range) # The max and min = confidence limit

# Make the plot
pdf(file=paste('Flood_frequency_plot_', distribution,'_maxLike.pdf',sep=""), width=7,height=5)
plot(flood_return,gev_quantile(1-1/flood_return,coef(x)[1],coef(x)[2],coef(x)[3]) ,
     log='x',t='o', ylim=c(0,max(conf_limits)), xlab='AEP of 1/Y Years', ylab='Discharge (m^3/s)',
     main=river_site,cex.main=1.5)
points(flood_return,conf_limits[1,],t='l',col=2,lty='dashed')
points(flood_return,conf_limits[2,],t='l',col=2,lty='dashed')
points(1/Q_AEP_est, Q_muda,col='steelblue',pch=19)
grid(nx=10,ny=10)
legend('topleft', c(paste('Fitted curve (', distribution, ')', sep=""), paste(cilevel*100, '% Confidence Limits (Profile Likelihood)',sep=""), 'Data'), lty=c('solid','dashed',NA),col=c('black', 'red', 'steelblue'), pch=c(NA,NA,19) ,bty='o',bg='white')
dev.off()
