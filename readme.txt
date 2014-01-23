These scripts can be used to do certain types of flood frequency analysis. 

The script 'river_freq_lmom.R' fits various curves to annual maxima data using L-moments, with parametric bootstrap confidence intervals.

The script 'prof_liklihood.R' fits various curves to annual maxima data using the method of maximum likelihood, with profile likelihood confidence intervals. It also has a 'check' code which uses MCMC to profile the likelihood, to check that it is working. Needs cleaning up / perhaps use MCMC instead of brute-force to create starting values for profiling.

The script 'muda_river.R' illustrates the use of several other R packages for frequency analysis, for example 'ismev' and 'fExtremes'

The script 'ismev_extra.R' contains a function which makes a pretty plot with the outputs of gum.fit in ismev, and also outputs information on the fit to a file.
