## Example using profile likelihood to compute the confidence limits for a
## return level plot

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


### STEP 1: Define a function which calculates the negative log likelihood of
### the distribution we want to fit

library(fExtremes) # This has a dgev function (gev density)
gev_negloglik<-function(xi, mu, beta){
    # Note that I am taking the log myself, because the 'log=TRUE' argument
    # seems to fail in this function
    -sum(log(dgev(Q_muda, xi, mu, beta, log=FALSE)))
}

# library(lmomco)
#gev_negloglik<-function(Q_lmoms){
#    # Note that I am taking the log myself, because the 'log=TRUE' argument
#    # seems to fail in this function
#    -sum(log(dlmomco(Q_muda, Q_lmoms)))
#}

### STEP 2: Compute the maximum likelihood estimate, and confidence limits on
### the parameters

startpars=list(xi=0., mu=400,beta=150)
x = mle(gev_negloglik, start=startpars, nobs=n, method='Nelder-Mead')

cilevel = 0.9
ci.x = confint(x, level=cilevel)


### Step 3: Numerically compute the 90% confidence intervals for a range
### of flood return periods
### Do this by taking a box around the parameter values contained in ci.x, and
### searching through this, ignoring values for which the negative log-liklihood
### makes that parameter set outside of the real confidence interval. 

flood_return=c(1.1,1.3,1.5,1.8,2,5,7,10,15,20,30, 50, 70, 90,100)
nn=20 # We divide the ci 'box' into n^3 values for searching
storeme = matrix(NA,ncol=length(flood_return)+1,nrow=nn^3) # Store the confidence limits
countme=0 # Used for counting in the loop
loglik_thresh=gev_negloglik(coef(x)[1],coef(x)[2],coef(x)[3])

for(i in seq(ci.x[1,1], ci.x[1,2],len=nn)){
    for(j in seq(ci.x[2,1], ci.x[2,2],len=nn)){
        for(k in seq(ci.x[3,1], ci.x[3,2],len=nn)){
        countme=countme+1  
        storeme[countme,1:length(flood_return)] = qgev(1-1/flood_return, xi=i, mu=j,beta=k)
        storeme[countme,length(flood_return)+1] = gev_negloglik(i,j,k)
        }
    }
}
# Examples of computing confidence intervals for the return period can be found
# in gevrlevelPlot
# Note that according to Bolker, the -2loglikelihood is approximately chisq
# distributed with df = number of free parameters (3 in the case of a gev).
# This means that we search for all likelihoods within qchisq(cilevel,3)/2 of
# the maximum likelihood
gev_df=3 # Note that gevrlevelPlot in fExtremes always uses gev_df=1
# Compute indicies of points which are inside the confidence limit
tmp = which(storeme[,length(flood_return)+1]< 
            gev_negloglik(coef(x)[1], coef(x)[2], coef(x)[3]) + qchisq(cilevel,gev_df)/2 )
conf_limits= apply(storeme[tmp,1:length(flood_return)], 2, range) # The max and min = confidence limit

# Make the plot
plot(flood_return,qgev(1-1/flood_return,coef(x)[1],coef(x)[2],coef(x)[3]) , log='x',t='o', ylim=c(0,max(conf_limits)), xlab='AEP of 1/Y Years', ylab='Discharge')
points(flood_return,conf_limits[1,],t='l',col=2,lty='dashed')
points(flood_return,conf_limits[2,],t='l',col=2,lty='dashed')
points(1/Q_AEP_est, Q_muda,col='steelblue',pch=19)
grid(nx=10,ny=10)
legend('topleft', c('Fitted curve', paste(cilevel*100, '% Confidence Limits (Profile likelihood)',sep=""), 'Data'), lty=c('solid','dashed',NA),col=c('black', 'red', 'steelblue'), pch=c(NA,NA,19) ,bty='o',bg='white')
dev.copy(pdf,'Flood_frequency_plot_gev_profLike.pdf')
dev.off()
