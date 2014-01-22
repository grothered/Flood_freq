#' Profile some function of the parameters of the log_likehood function 
#' i.e. find it's max/min over the region where the log-likelihood is within qchisq(level,1) of the maximum
#' @param var_to_profile -- function which takes as input a vector of the model
#'                    parameters (e.g. from coef(mle(..)), and returns the variable we want to profile
#' @param log_lik -- log likelihood function [NOT NEGATIVE LOG_LIKELIHOOD!]
#' @param maxlikpar -- the maximum likelihood values
#' @param level -- the confidence level
#' @param details -- TRUE/FALSE -- return a list with the optim results from the lower/upper CIs, or just the values
#' @param optim.method -- parameter for optim. Probably should only be 'Nelder-Mead', or 'SANN' as a last resort
#' @param maxit -- parameters to optim
#
# optim.method might be Nelder-Mead or SANN
#
profile_function<-function(var_to_profile, log_lik, maxlikpar, 
                           searchStart=rbind(maxlikpar,maxlikpar),
                           level=0.95, details=FALSE,
                           optim.method='Nelder-Mead', maxit=2e+06, ...){
    f<-function(x, lowerCI=TRUE){
        # Function used with optim to get range of a function within the acceptable likelihood region
        if(lowerCI){
            # Minimised to get lower CI
            out=var_to_profile(x) + ( (log_lik(x) < log_lik(maxlikpar)-qchisq(level,1)/2))*.Machine$double.xmax
        }else{
            # Maximised to get upper CI
            out=var_to_profile(x) - ( (log_lik(x) < log_lik(maxlikpar)-qchisq(level,1)/2))*.Machine$double.xmax
        }
        return(out)
    }
    # Find the min/max value of var_to_profile within the acceptable likelihood region
    lowerP=optim(searchStart[1,], f , method=optim.method, control=list(maxit=maxit)) 
    lowerP=optim(lowerP$par, f , method=optim.method, control=list(maxit=maxit)) 
    upperP=optim(searchStart[2,],f,lowerCI=FALSE, method=optim.method, control=list(fnscale=-1,maxit=maxit))
    upperP=optim(upperP$par,f,lowerCI=FALSE, method=optim.method, control=list(fnscale=-1,maxit=maxit))

    if(lowerP$convergence!=0 | upperP$convergence!=0) warning('profile_function did not have convergence==0')

    if(!details){
        return(c(lowerP$value, upperP$value))
    }else{
        return(list(lowerCI=lowerP, upperCI=upperP))
    }
}    
