######################################################################
###
### Supplementary R script to run analysis presented in paper, using
### Bayesian methods as implemented in R-jags.
###
### Pascal A. Niklaus & Bernhard Schmid
###
### History:
### - 20-04-2016 file created
### - 30-05-2016 version for final test
### - 07-06-2016 minor corrections
### - 15-06-2016 script separated from main script

######################################################################
###
### Notes on formatting
###
### This script is self-contained and can be executed in sequence,
### provided that the required libraries and extra software packages
### are installed (see below), and that the required data files are
### available.
###
### Depending on your software environment you may want to set the
### working directory using 'setwd("/path/to/your/directory")'.

######################################################################
###
### Required R libraries and extra software
###
### Libraries can be installed from the command line with
### 'install.packages(...)'.
###
### Library for analysis with Bayesian model:
### - rjags
###
### Note that jags needs to be installed
### (available at http://mcmc-jags.sourceforge.net).

library(rjags)
options(digits=4)                     # set number of post-decimal digits in output to 4
                                      # increased to show more digits if required
rm(list=ls())                         # clear workspace

d14X.plot <- read.csv("derived_data/pilot_pool_X_plot_ba.csv")

######################################################################
###
### Data analysis with Bayesian statistics using rjags :

## Write model to file;
## here, we use abbreviated variable names:
##   P   == "plot"
##   L   == "light"
##   C   == "com"
##   LxC == "light:com"
## The estimated parameters are:
##   beta_0      == overall mean
##   beta_L[1:2] == effect of light=="c" and "s" (in this order)
##   beta_D[1:3] == effects of diversity (levels 1,2,4, in this order)

m11.bugs <- '
data {
    for(i in 1:n_P) {
        y[i] <- ba[i]
    }    
}
model {
    for(i in 1:n_P) {
        y[i] ~ dnorm( yhat_plot[i], tau ) 
        yhat_plot[i] <- mu_LxC[ LxC[i] ]
    }
    for(j in 1:n_LxC) {
        mu_LxC[j]  ~ dnorm( yhat_LxC[j], tau_LxC )
        yhat_LxC[j] <- mu_C[ C[j] ] + beta_L[ L[j] ]
    }
    for(k in 1:n_C) {
        mu_C[k] ~ dnorm(yhat_C[ k ], tau_C ) 
        yhat_C[k] <- beta_0 + beta_D[ D[k] ]
    }

    ## ----- Non-informative gamma priors for variance components:
    tau_C   ~ dgamma(0.001,0.001)
    tau_LxC ~ dgamma(0.001,0.001)
    tau     ~ dgamma(0.001,0.001)

    ## ----- Diffuse prior for grand mean = terminal node:
    beta_0 ~ dnorm(20, 0.001)       # (mean,1/variance)

    ## ----- Normal, non-informative priors for fixed parameters:
    beta_L[1] <- 0                  # resolve aliasing 
    beta_L[2] ~ dnorm(0, 0.00001)
    beta_D[1] <- 0                  # resolve aliasing
    beta_D[2] ~ dnorm(0, 0.00001)
    beta_D[3] ~ dnorm(0, 0.00001)
 
    ## ----- Compute standard deviations:
    sd     <- 1/sqrt(tau)
    sd_C   <- 1/sqrt(tau_C)
    sd_LxC <- 1/sqrt(tau_LxC)

    ## ----- Create new names for (derived) parameter:
    grand_mean <- beta_0
    shade_eff  <- beta_L[2] 
    div_2      <- beta_D[2]
    div_4      <- beta_D[3]
    VC_C       <- sd_C^2
    VC_LxC     <- sd_LxC^2
    VC_P       <- sd^2
}
'

writeLines(m11.bugs, "m11.bug")

## Prepare list with data for rjags model:

## Create factor levels for interaction light x community composition.
## Also, re-create factor 'community'. This is important since
## otherwise the factor 'remembers' that there once were more levels
## than currently are present in the subset.

tmp.LxC <- factor(paste(d14X.plot$light,
                        d14X.plot$com,
                        sep=":"))
tmp.C   <- factor(as.character(d14X.plot$com))

## Create list with number of levels for the different factors:
m11.d <- list(n_P = nrow(d14X.plot),
              n_C = nlevels(tmp.C),
              n_LxC = nlevels(tmp.LxC))

## Add vectors to map groups to next aggregate level in hierarchy:
## - LxC maps plot to combination of light and com
## - C   maps LxC to com
## - D   maps C to fdiv
m11.d$LxC <- as.numeric( tmp.LxC )
m11.d$C <- aggregate(as.numeric(tmp.C),  
                     by=list( tmp.LxC ),
                     mean )$x
m11.d$L <- aggregate(as.numeric(d14X.plot$light),
                     by=list(tmp.LxC),
                     mean)$x
m11.d$D <- aggregate(as.numeric(d14X.plot$fdiv),
                     by=list(tmp.C),
                     mean)$x

## Check if treatment coding is correct:
stopifnot(all( m11.d$D[ m11.d$C[ m11.d$LxC ]] == as.numeric(d14X.plot$fdiv) ))
stopifnot(all( m11.d$L[ m11.d$LxC ] == as.numeric(d14X.plot$light) ))

## Add dependent variable to analyse:
m11.d$ba <- d14X.plot$ba

## set up jags model...
m11<-jags.model('m11.bug',
                data=m11.d,
                inits=list(tau=10,tau_LxC=5,tau_C=20,beta_0=10),
                n.chains=4,
                n.adapt=1000)

## ... and run Gibb sampler for initial phase:
update(m11,n.iter=100000)

## Extract data from chain by more sampling:
par.vec <- c("grand_mean","shade_eff","div_2","div_4","VC_C","VC_LxC","VC_P")
m11.res <- jags.samples(m11,
                        par.vec,
                        n.iter=100000)

## Extract parameter means and 95% quantiles and store them in data
## frame for better printing:

m11.tab <- data.frame(mean=rep(NA,length(par.vec)),CI.L=NA,CI.U=NA)
rownames(m11.tab) <- par.vec
for(par in par.vec) {
    m11.tab[par,"mean"] <- mean(m11.res[[par]])
    m11.tab[par,"CI.L"] <- quantile(m11.res[[par]],0.025)
    m11.tab[par,"CI.U"] <- quantile(m11.res[[par]],0.975)
}

m11.tab
#|                mean       CI.L    CI.U
#| grand_mean  19.2702   9.124948  29.464
#| shade_eff  -11.3604 -14.184333  -8.536
#| div_2        1.3661 -11.756926  14.425
#| div_4        0.6403 -22.284567  23.394
#| VC_C       100.2471  29.376322 286.512
#| VC_LxC       5.7797   0.002134  25.049
#| VC_P        20.5122  14.236570  29.215

coef(m11)
#| $beta_0
#| [1] 24.86
#| 
#| $beta_D
#| [1]     NA  5.814 -8.786
#| 
#| $beta_L
#| [1]     NA -12.99
#| 
#| $div_2
#| [1] 5.814
#| 
#| $div_4
#| [1] -8.786
#| 
#| $grand_mean
#| [1] 24.86
#| 
#| $mu_C
#|  [1] 13.37 27.82 23.98 12.62 18.64 34.32 26.82 30.75 12.65 16.94 18.72
#| 
#| $mu_LxC
#|  [1] 13.67341 27.48206 23.54847 12.79501 18.59729 34.35002 27.24403 31.24757
#|  [9] 12.78423 16.61199 19.11422  0.09365 14.43465 10.54191  0.02881  5.72966
#| [17] 21.18787 13.74581 17.54326  0.09013  4.27424  5.57829
#| 
#| $shade_eff
#| [1] -12.99
#| 
#| $tau
#| [1] 0.03913
#| 
#| $tau_C
#| [1] 0.009772
#| 
#| $tau_LxC
#| [1] 13.21
#| 
#| $yhat_plot
#|  [1] 13.67341 13.67341 13.67341 13.67341  0.09365  0.09365  0.09365  0.09365
#|  [9] 34.35002 34.35002 34.35002 34.35002 21.18787 21.18787 21.18787 21.18787
#| [17] 27.48206 27.48206 27.48206 27.48206 14.43465 14.43465 14.43465 14.43465
#| [25] 12.78423 12.78423 12.78423 12.78423  0.09013  0.09013  0.09013  0.09013
#| [33] 12.79501 12.79501 12.79501 12.79501  0.02881  0.02881  0.02881  0.02881
#| [41] 27.24403 27.24403 27.24403 27.24403 13.74581 13.74581 13.74581 13.74581
#| [49] 19.11422 19.11422 19.11422 19.11422  5.57829  5.57829  5.57829  5.57829
#| [57] 18.59729 18.59729 18.59729 18.59729  5.72966  5.72966  5.72966  5.72966
#| [65] 31.24757 31.24757 31.24757 31.24757 17.54326 17.54326 17.54326 17.54326
#| [73] 16.61199 16.61199 16.61199 16.61199  4.27424  4.27424  4.27424  4.27424
#| [81] 23.54847 23.54847 23.54847 23.54847 10.54191 10.54191 10.54191 10.54191
#| 

## Calculate residual standard deviation. We use pnorm to get the quantiles
## that correspond to minus and plus 1 s.d.:
m11.limits <- unname(quantile(m11.res$grand_mean, pnorm(c(-1,1))))
(m11.limits[2] - m11.limits[1]) / 2
#| [1] 4.739

###
### End of script.
###
######################################################################
