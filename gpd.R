## All this code is slightly adapted from fExtremes
## so I can distribute gpd code without the
## many dependencies of fExtremes

dgpd<-function (x, location = 0, scale = 1, shape = 0, log = FALSE) 
{
    stopifnot(min(scale) > 0)
    stopifnot(length(shape) == 1)
    d <- (x - location)/scale
    nn <- length(d)
    scale <- rep(scale, length.out = nn)
    index <- (d > 0 & ((1 + shape * d) > 0)) | is.na(d)
    if (shape == 0) {
        d[index] <- log(1/scale[index]) - d[index]
        d[!index] <- -Inf
    }
    else {
        d[index] <- log(1/scale[index]) - (1/shape + 1) * log(1 + 
            shape * d[index])
        d[!index] <- -Inf
    }
    if (!log) 
        d <- exp(d)
    #attr(d, "control") = data.frame(location = location[1], scale = scale[1], 
    #    shape = shape[1], log = log, row.names = "")
    return(d)
}


qgpd<-function (p, location = 0, scale = 1, shape = 0, lower.tail = TRUE){
    stopifnot(min(scale) > 0)
    stopifnot(length(shape) == 1)
    stopifnot(min(p, na.rm = TRUE) >= 0)
    stopifnot(max(p, na.rm = TRUE) <= 1)
    if (lower.tail) 
        p <- 1 - p
    if (shape == 0) {
        q = location - scale * log(p)
    }
    else {
        q = location + scale * (p^(-shape) - 1)/shape
    }
    #attr(q, "control") = data.frame(location = location[1], scale = scale[1], 
    #    shape = shape[1], lower.tail = lower.tail, row.names = "")
    return(q)
}


pgpd<-function (q, location = 0, scale = 1, shape = 0, lower.tail = TRUE) 
{
    stopifnot(min(scale) > 0)
    stopifnot(length(shape) == 1)
    q <- pmax(q - location, 0)/scale
    if (shape == 0) 
        p <- 1 - exp(-q)
    else {
        p <- pmax(1 + shape * q, 0)
        p <- 1 - p^(-1/shape)
    }
    if (!lower.tail) 
        p <- 1 - p
    #attr(p, "control") = data.frame(location = location[1], scale = scale[1], 
    #    shape = shape[1], lower.tail = lower.tail, row.names = "")
    return(p)
}


rgpd<-function (n, location = 0, scale = 1, shape = 0) 
{
    stopifnot(min(scale) > 0)
    stopifnot(length(shape) == 1)
    if (shape == 0) {
        r = location + scale * rexp(n)
    }
    else {
        r = location + scale * (runif(n)^(-shape) - 1)/shape
    }
    attr(r, "control") = data.frame(location = location[1], scale = scale[1], 
        shape = shape[1], row.names = "")
    return(r)
}

#############################
test_gpd_code<-function(){
    x=seq(0,1000,len=100)
    p=seq(0,1,len=101)
    #library(fExtremes) # Should be identical to this
    
    xd=dgpd(x, shape=0.2, location=200,scale=100)
    xq=qgpd(p, shape=0.2, location=200,scale=100)
    xp=pgpd(x, shape=0.2, location=200,scale=100)
    set.seed(1) 
    xr=rgpd(100, shape=0.2, location=200,scale=100)


    library(fExtremes)
    yd=fExtremes::dgpd(x,xi=0.2,mu=200,beta=100)
    
    stopifnot(all(xd==yd))
    print('PASS - same as dgpd fExtremes')

    yq=fExtremes::qgpd(p,xi=0.2,mu=200,beta=100)
    stopifnot(all(xq==yq))
    print('PASS - same as qgpd fExtremes')

    yp=fExtremes::pgpd(x,xi=0.2,mu=200,beta=100)
    stopifnot(all(xp==yp))
    print('PASS - same as pgpd fExtremes')
   
    set.seed(1) 
    yr=fExtremes::rgpd(100,xi=0.2,mu=200,beta=100)
    stopifnot(all(xr==yr))
    print('PASS - same as rgpd fExtremes')

}
