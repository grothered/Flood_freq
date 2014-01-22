
dgev<-function (x, shape = 1, scale = 1, location = 0, log = FALSE) 
{
    #fx <- 1/scale * (1 + shape * ((x - location)/scale))^(-1/shape - 
    #    1) * exp(-(1 + shape * ((x - location)/scale))^(-1/shape))
    fx = -log(scale) + log(1+shape*((x-location)/scale))*(-1/shape-1) + (-(1+shape*((x-location)/scale))^(-1/shape))
    if (!log) 
        return(exp(fx))
    else return(fx)
}
