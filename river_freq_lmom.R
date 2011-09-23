# Code to analyse the styx river data, from Example 7 in the ARR Flood
# frequency analysis draft update chapter (Kuczera and Franks)

# This uses L-moments to fit a distribution supported in the L-moments package,
# and a parametric bootstrap to compute quantile confidence limits.


# The code was initially written for the styx river, using a gev distribution. 
# Thus, 

# User modified parameters
site_name = 'Muda River'
input_data = scan('Muda_discharge.txt')
probability_distribution= 'gev' #'gev' # 'gum', 'pe3', 'lp3', other choices from lmomco
alpha =0.4 # Parameter to adjust empirical AEP estimates (Kuczera and Franks)
nrand = 25000 # Number of bootstrap samples for confidence interval calculation -- 5000 was not really enough. 




############################################################
### The user may less often want to modify things below here
############################################################

# Adjust script to treat the log pearson 3 distribution
if(probability_distribution == 'lp3'){
    # Treat log-pearson type 3 using the pearson type 3, by log transforming
    # the data, then transforming back later
    distribution = 'pe3'
    input_data = log(input_data)
}else{
    distribution = probability_distribution
}

styx_Q = sort(input_data, decreasing=T)
n = length(styx_Q)

# Compute empirical AEPs for interest / later plotting
# Note the data is already sorted
styx_AEPs = (seq(1,n)-alpha)/(n+1-2*alpha)

# Import lmoments library
library(lmomco)

# Find Lmoments of sample
styx_lmoms = lmom.ub(styx_Q)

# Fit gev
styx_gev = lmom2par(styx_lmoms, distribution) #pargev(styx_lmoms)

# Parameters for bootstrap
#nrand = 5000 # Number of samples -- 5000 was not really enough
storeAEPs = 1/c(100, 90, 80, 70, 60, 50, 40, 30, 20, 15, 10, 7, 5, 1.1) # We will store simulated quantiles at these AEP values

# Storage 
store_gev = matrix(NA, nrow=nrand,ncol=length(styx_gev$para))
store_aep = matrix(NA, nrow=nrand, ncol=length(storeAEPs))

# Bootstrap sampling loop
for(i in 1:nrand){
    # Take a random sample and compute lmoments and gev
    rand_data = rlmomco(n, styx_gev)
    # rand_data = sample(styx_Q,n, replace=TRUE) # In the styx river example,
    # this gave smaller CIs, consistent with the literature

    rand_lmoms = lmom.ub(rand_data)

    rand_gev = lmom2par(rand_lmoms, distribution) #pargev(rand_lmoms)


    # Store the gev parameters
    store_gev[i,] = rand_gev$para

    # Store the quantiles at the desired AEP values
    # Note that the input is probability values = 1-AEP
    store_aep[i,] = qlmomco(1-storeAEPs, rand_gev)
        
}

# Compute 5% and 95% confidence limits. Note that according to Kuczera and
# Franks, these should underestimate the true 90% confidence limits, because
# the sampling assumed that the fitted parameters were the true parameters
# This seems to fit with my experience
five_quant = storeAEPs*NA
ninetyfive_quant = storeAEPs*NA
for(i in 1:length(storeAEPs)){
    tmp = quantile(store_aep[,i], probs = c(0.05, 0.95))
    five_quant[i] = tmp[1]
    ninetyfive_quant[i] = tmp[2]
}

# Compute quantities that will be useful for later plotting
if(probability_distribution == 'lp3'){
    # Take exponential to undo the logarithm, which we used in the calculation 
    five_quant = exp(five_quant)
    ninetyfive_quant = exp(ninetyfive_quant)
    theoretical_probabilities = plmomco(styx_Q, styx_gev)

    styx_Q = exp(styx_Q)
    fitted = exp(qlmomco(1-storeAEPs, styx_gev))
    theoretical_quantiles =  exp(qlmomco(1-styx_AEPs, styx_gev))
}else{
    fitted = qlmomco(1-storeAEPs, styx_gev)
    theoretical_quantiles =  qlmomco(1-styx_AEPs, styx_gev)
    theoretical_probabilities = plmomco(styx_Q, styx_gev)
}

#####################################
#
# Plot it all like Kuczera and Franks
#
#####################################
pdf(file=paste('Flood_frequency_plot_',probability_distribution,'_lmom.pdf',sep=""),width=8,height=6)

plot(1/storeAEPs, fitted, log='x', xlab='AEP of 1/Y years', ylab='Discharge (m^3/s)', ylim=c(min(five_quant), max(ninetyfive_quant)), t='l', main=site_name)
points(1/storeAEPs, five_quant,t='l',col='red',lty='dashed')
points(1/storeAEPs, ninetyfive_quant,t='l',col='red',lty='dashed')
points(1/styx_AEPs, styx_Q,col='blue',pch=19)
# add some grid lines
grid(nx=10,ny=10)
legend('topleft', c(paste('Fitted model (', probability_distribution, ')'), 'Approximate 90% Confidence Intervals', 'Data'), 
        lty=c('solid', 'dashed', NA),
        col=c('black', 'red', 'blue'),
        pch=c(NA, NA,19), 
        bg='white'
        )

dev.off()

#####################################
#
# Quantile-Quantile plot
#
#####################################
pdf(file=paste('Quantile_plot_',probability_distribution,'_lmom.pdf',sep=""),width=8,height=6)
plot(styx_Q, theoretical_quantiles, xlab='Measured Discharge (m^3/s)', ylab='Theoretical Discharge (m^3/s)', main=site_name)
abline(0,1)
dev.off()

pdf(file=paste('Probability_plot_',probability_distribution,'_lmom.pdf',sep=""),width=8,height=6)
plot(1-styx_AEPs,theoretical_probabilities , xlab='Empirical Probability', ylab='Theoretical Probability', main=site_name)
abline(0,1)
dev.off()

write.table(cbind(storeAEPs,fitted, five_quant, ninetyfive_quant), file=paste('Fitted_',probability_distribution,'_', site_name, '.txt', sep=''), col.names=c('AEP', 'Model', '5%', '95%'), row.names=FALSE, sep="," )
