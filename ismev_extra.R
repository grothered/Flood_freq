########################################################################
# Code to make a nicer Return level plot, and write out fitted values +
# confidence intervals to a file called 'gumbel_model.txt'
# This is a slight modification of the function gum.rl in ismev
########################################################################

gum.rl2<-function (gumbel_model) {
    # This is a modification of the gum.rl function from ismev, designed to
    # make a nicer plot, and write fitted values to a file

    # Input to gum.rl based on gum.diag
    a = c(gumbel_model$mle, 0) 
    mat = gumbel_model$cov
    dat = gumbel_model$data

    # Plotting position adjustment factor
    alpha=0.4
  
    # Core of function 
    eps <- 1e-06
    a1 <- a
    a2 <- a
    a1[1] <- a[1] + eps
    a2[2] <- a[2] + eps

    #f <- c(seq(0.01, 0.09, by = 0.01), 0.1, 0.2, 0.3, 0.4, 0.5, 
    #    0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 0.995)
    plot_upper_bound = max(100,2*length(dat))
    AEPs = 1/c(1.1,1.3,1.5, 1.8, 2,3,5,seq(10,plot_upper_bound,by=10))

    f = 1-AEPs 
    q <- gevq(a, 1 - f)
    d1 <- (gevq(a1, 1 - f) - q)/eps
    d2 <- (gevq(a2, 1 - f) - q)/eps
    d <- cbind(d1, d2)
    v <- apply(d, 1, q.form, m = mat)
    #plot(-1/log(f), q, log = "x", type = "n", xlim = c(0.1, 1000), 
    #    ylim = c(min(dat, q-1.96*sqrt(v)), max(dat, q+1.96*sqrt(v))), xlab = "Return Period", 
    #    ylab = "Return Level")
    
    # Note -- here we plot in terms of AEP of 1/Y
    plot(1/AEPs, q, log = "x", type = "n", xlim = c(1, plot_upper_bound), 
        ylim = c(min(dat, q-1.96*sqrt(v)), max(dat, q+1.96*sqrt(v))), ylab = "Discharge (m^3/s)", 
        xlab = "AEP of 1/Y years")
    #title("Return Level Plot")
    lines(1/AEPs, q)
    lines(1/AEPs, q + 1.96 * sqrt(v), col = 2, lty='dashed')
    lines(1/AEPs, q - 1.96 * sqrt(v), col = 2, lty='dashed')
    #points(-1/log((1:length(dat))/(length(dat) + 1)), sort(dat), col='steelblue', pch=19)
    points(1/(1-(((1:length(dat))-alpha)/(length(dat) + 1-2*alpha))), sort(dat), col='steelblue', pch=19)
    legend('topleft', c('Fitted Gumbel Model', 'Approximate 95% Confidence Intervals', 'Data'), pch=c(NA, NA, 19),
           lwd=c(1,1,0),lty=c('solid', 'dashed', NA), col=c(1, 2, 'steelblue'), bty='n')

    # Write the model fit and confidence intervals to a file
    write.table(file='gumbel_model.txt', cbind(1/(1-f), q, q+1.96*sqrt(v), q-1.96*sqrt(v)), col.names=c("AEP of 1/Y", "Fitted Model", "Upper 95%CI", "Lower 95% CI"), row.names=FALSE)
}

