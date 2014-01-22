##############################################################################

dlp3<-function(x,shape, scale, thresh=0,log=FALSE){
    # Here is a log pearson type 3 implementation which permits scale to be
    # negative
    #
    # dens = abs(scale)/(x*Gamma(shape))*[scale*(log(x)-thresh)]^{shape-1} x exp( -scale*(log(x)-thresh))
    #
    fx=dgamma(scale*(log(x)-thresh), shape, 1, log=T)+log(abs(scale))-log(x) 
    if(log==FALSE){ fx=exp(fx)}
    return(fx)
}

plp3<-function(x, shape, scale, thresh=0,lower.tail=TRUE, log.p=FALSE){
    # CDF for log pearson type 3 distribution
    #
    # Idea: sign(scale)*[log(x)-thresh] has a gamma distribution with parameters shape,abs(scale)
    #
    # Idea: pgamma( q, shape, scale) = pgamma(scale*q, shape, 1)

    # Ensure shape, scale,thresh have length=length(x)
    scale=x*0+scale
    shape=x*0+shape
    thresh=x*0+thresh

    # Find where scale >0
    mm=(scale>=0)
    log_p=numeric(length(scale))

    # Cases with positive scale parameter
    if(sum(mm)>0){
        indz=which(mm)
        # Compute log(p)
        log_p[indz]=pgamma( (log(x[indz]) - thresh[indz]), shape[indz], scale[indz], lower.tail=lower.tail, log.p=TRUE)
    
        if(log.p==FALSE) log_p[indz]=exp(log_p)
    }

    # Cases with negative scale parameter
    if(sum(mm)<length(scale)){
        indz=which(!mm)
        # Need to convert p to 1-p, due to the reversal of sign caused by the
        # sign(scale) operation, which reverses the limits in the cdf
        # integral
        # Since we work directly with log(p), use expm1 to get 1-p
        #
        # NOTE: expm1 = exp(x) - 1, and is accurate for x close to 0
        #log_p[indz]=-expm1(
        #                  pgamma(sign(scale[indz])*(log(x[indz])-thresh[indz]), 
        #                              shape[indz],abs(scale[indz]),
        #                               lower.tail=lower.tail,log.p=TRUE)
        #                  )
        #if(log.p==TRUE) log_p[indz]=log(log_p)

        log_p[indz]=
                     1- pgamma(sign(scale[indz])*(log(x[indz])-thresh[indz]), 
                               shape[indz],abs(scale[indz]),
                               lower.tail=lower.tail,log.p=FALSE)
        if(log.p==TRUE) log_p[indz]=log(log_p)
    }

    return(log_p)
}

qlp3<-function(p,shape,scale,thresh=0, lower.tail=TRUE,log.p=FALSE){
    # Inverse CDF for log-pearson type 3 distribution
    #
    # Idea: sign(scale)*[log(x)-thresh] has a gamma distribution with parameters shape,abs(scale)
    #
    # Because of the sign(scale) operation, we must replace p with 1-p for
    # negative scale
    #
    
    # Ensure shape, scale,thresh have length=length(p)
    scale=p*0+scale
    shape=p*0+shape
    thresh=p*0+thresh

    # Find where scale >0
    mm=(scale>=0)
    log_p=numeric(length(scale))
  
    qua_nolog=numeric(length(scale)) 
    if(sum(mm)>0){ 
        # Cases with positive scale
        indz=which(mm)
        qua_nolog[indz]=qgamma(p[indz], shape[indz],abs(scale[indz]),lower.tail=lower.tail,log.p=log.p)
    }

    if(sum(mm)<length(scale)){
        # Cases with negative scale
        # WARNING: Can be subject to round-off error -- could potentially be improved by comparing with plp3?
        #          using optimize to minimise abs(plp3(y) - x), with dlp3(x) = y
        indz=which(!mm)
        if(log.p==FALSE){
             #qua_nolog[indz]=qgamma(1-p[indz], shape[indz],abs(scale[indz]),lower.tail=lower.tail,log.p=log.p)
             qua_nolog[indz]=qgamma(p[indz], shape[indz],abs(scale[indz]),lower.tail=!lower.tail,log.p=log.p)
        }else{
             #qua_nolog[indz]=qgamma(log1p(-exp(p[indz])), shape[indz],abs(scale[indz]),lower.tail=lower.tail,log.p=TRUE)
             #qua_nolog[indz]=qgamma(1-exp(p[indz]), shape[indz],abs(scale[indz]),lower.tail=lower.tail,log.p=FALSE)
             qua_nolog[indz]=qgamma(exp(p[indz]), shape[indz],abs(scale[indz]),lower.tail=!lower.tail,log.p=FALSE)
        }

    }
    
    qua_nolog=(qua_nolog/sign(scale)+thresh)

    return(exp(qua_nolog))
}

rlp3<-function(n, shape, scale, thresh=0){
    # Random numbers for log-pearson type 3 distribution
    #
    # Idea: sign(scale)*[log(x)-thresh] has a gamma distribution with parameters shape,abs(scale)
    #
    myrand=exp((rgamma(n, shape, abs(scale))/sign(scale) + thresh))
    return(myrand)
}

test_lp3<-function(){
    # Test against FAdist where that was applicable
    require(FAdist)

    shape=187
    scale=9.6
    thresh=-13

    # quantile
    p=0.5
    checkQ=all.equal(qlgamma3(p,shape,scale,thresh) , qlp3(p,shape,scale,thresh))
    if(!checkQ){
        stop('quantile function failed')
    }else{
        print('pass')
    }

    # Probability
    q=100
    checkp=all.equal(plgamma3(q,shape,scale,thresh) , plp3(q,shape,scale,thresh))
    if(!checkp){
        stop('quantile function failed')
    }else{
        print('pass')
    }

    # Density
    x=seq(0,500,len=500)
    checkd=all.equal(dlgamma3(x,shape,scale,thresh) , dlp3(x,shape,scale,thresh))
    if(!checkd){
        stop('density function failed')
    }else{
        print('pass')
    }

    # Random numbers
    print('Making random numbers ...')
    N=10000000
    m1=rlp3(N, shape, scale, thresh)
    m1_dens=density(m1,n=N/4)
    plotshift=min(m1_dens$x)+1

    par(mfrow=c(2,2))
    plot(m1_dens$x+plotshift,m1_dens$y,t='l',log='xy',main='Empirical density from rlp3 vs analytical density, positive scale')
    points(m1_dens$x+plotshift, dlp3(m1_dens$x,shape,scale,thresh),t='l',col=2)
    legend('bottomleft',c('Empirical Density from random deviates','Density'),lwd=c(1,1),col=c(1,2))

    ################################################################################
    #
    # Check for accuracy / consistency with negative scale
    shape=187
    scale=-9.6
    thresh=13

    # Check plp3 is the integral of dlp3
    xx=seq(1.0e-6,10,len=200)
    pxx=plp3(xx,shape,scale,thresh)
    passTest=rep(0,length(xx))
    for(i in 1:length(pxx)){
        myInt=integrate(dlp3,0,xx[i], shape=shape,scale=scale,thresh=thresh,rel.tol=1.0e-12)
        passTest[i]=all.equal(myInt$value,pxx[i])    
    }
    
    if(any(passTest==0)){
        stop('CDF fails with negative scale parameter')
    }else{
        print('pass')
    }

    # Check that exp( plp3(..., log.p=TRUE) ) = plp3(...)
    pxxLog=plp3(xx,shape,scale,thresh,log.p=TRUE)
    if(!all.equal(exp(pxxLog),pxx)){
        stop('CDF with log.p=TRUE not consistent with log.p=FALSE')
    }else{
        print('pass')
    }
    
    # Check that 1-plp3(..., lower.tail=FALSE) ) = plp3(...)
    pxxUpper=plp3(xx,shape,scale,thresh,lower.tail=FALSE)
    if(!all.equal(1-pxxUpper,pxx)){
        stop('CDF with lower.tail=FALSE not consistent with default')
    }else{
        print('pass')
    }
    
    # Check that 1-exp(plp3(..., lower.tail=FALSE,log.p=TRUE) )) = plp3(...)
    pxxUpperLog=plp3(xx,shape,scale,thresh,lower.tail=FALSE,log.p=TRUE)
    if(!all.equal(1-exp(pxxUpperLog),pxx)){
        stop('CDF with lower.tail=FALSE and log.p=TRUE not consistent with default')
    }else{
        print('pass')
    }

    # Use qlp3 to get quantile values at pxx, assuming it has passed
    quant_pxx=qlp3(pxxLog, shape,scale,thresh,log.p=TRUE)
    quant_pxx2=qlp3(pxx, shape,scale,thresh,log.p=FALSE)
    quant_pxx3=qlp3(pxxUpper, shape,scale,thresh,lower.tail=FALSE,log.p=FALSE)
    quant_pxx4=qlp3(pxxUpperLog, shape,scale,thresh,lower.tail=FALSE,log.p=TRUE)
    if((!all.equal(quant_pxx,quant_pxx2) ) |
       (!all.equal(quant_pxx,quant_pxx3) ) |
       (!all.equal(quant_pxx,quant_pxx4) ) 
        ){
        stop('Inv CDF not consistent using some combinations of log.p=TRUE/FALSE and lower.tail=TRUE/FALSE')
    }else{
        print('pass')
    }

    # NOTE: This is subject to round-off error!
    if(!all.equal(quant_pxx,xx, tol=1.0e-03)){
        stop('Quantile function fails with negative scale parameter')
    }else{
        print('pass -- but quantile function is subject to round-off error near p=1.0')
    }

    plot(quant_pxx,xx, main='x vs qlp3(plp3(x))')
    abline(0,1,col=2)
    plot(quant_pxx-xx, main='x - qlp3(plp3(x))  (ROUND OFF ERRORS)')
    abline(h=0,col=2)
    
    print('Making random numbers ...')
    N=10000000
    m1=rlp3(N, shape, scale, thresh)
    m1_dens=density(m1,n=N/4)
    plotshift=min(m1_dens$x)+1

    plot(m1_dens$x+plotshift,m1_dens$y,t='l',log='xy',main='Empirical density from rlp3 vs analytical density, negative scale')
    points(m1_dens$x+plotshift, dlp3(m1_dens$x,shape,scale,thresh),t='l',col=2)
    legend('bottomleft',c('Empirical Density from random deviates','Density'),lwd=c(1,1),col=c(1,2))

}
