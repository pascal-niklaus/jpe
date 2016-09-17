######################################################################
###
### Supplementary R script S4. Analysis of growth rates (slopes) as
### described in manuscript.
###
### Pascal A. Niklaus & Bernhard Schmid
###
### History:
### - 20-04-2016 file created
### - 30-05-2016 version for final test
### - 07-06-2016 minor corrections
### - 16-09-2016 code separated from main file

######################################################################
###
### Notes on formatting
###
### This script is self-contained and can be executed in sequence,
### provided that the required libraries and extra software packages
### are installed (see below).
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
### Libraries for analysis by maximum likelihood:
### - nlme, lmer, lmerTest, pbkrtest

library(lmerTest)
library(lme4)                         # not really needed since loaded by lmerTest
options(digits=4)                     # set number of post-decimal digits in output to 4
                                      # increased to show more digits if required
rm(list=ls())                         # clear workspace

for(d in c("derived_data","figures")) # create directories if these don't exist
    if(!dir.exists(d))                # data will be stored there
        dir.create(d)

### Analyse slopes of height ~ time for each individual:

## Create data file with slopes of height ~ time for each individual.
## First, create a data frame with one row per tree. Then, calculate
## the slope of height against time by linear regression for each
## individual separately:

d <- read.csv("derived_data/pilot_bd_all_pools.csv")

d$ftime <- factor(sprintf("t-%02d",d$time))  # time as factor

d.slopes <- aggregate(div ~ block + plot + com + pool + uind + fdiv + light + sp,                         
                      data = d,
                      FUN = mean)

d.slopes$h.slope <- NA

for(i in 1:nrow(d.slopes)) {
    tmp <- subset(d, uind == d.slopes$uind[i] & is.finite(height) )
    if(nrow(tmp) >= 2)
        d.slopes$h.slope[i] <- coef(lm(height~time,data=tmp))["time"]
}

## Inspection of the data shows that there are some really bad outliers,
## which we remove for this example.

idx <- is.finite(d.slopes$h.slope) & (d.slopes$h.slope < -9 | d.slopes$h.slope > 14)
d.slopes$h.slope[idx] <- NA

write.csv(d.slopes,"derived_data/pilot_height_slopes.csv",row.names=FALSE)

head(d.slopes)
#|   block    plot com pool         uind fdiv light sp div h.slope
#| 1    B1 B1.Ya04  cg    Y B1.Ya04|i-01   D1     c cg   1  2.0000
#| 2    B1 B1.Ya04  cg    Y B1.Ya04|i-02   D1     c cg   1  0.9275
#| 3    B1 B1.Ya04  cg    Y B1.Ya04|i-03   D1     c cg   1  4.9142
#| 4    B1 B1.Ya04  cg    Y B1.Ya04|i-04   D1     c cg   1  0.0000
#| 5    B1 B1.Ya04  cg    Y B1.Ya04|i-05   D1     c cg   1  2.3860
#| 6    B1 B1.Ya04  cg    Y B1.Ya04|i-06   D1     c cg   1  2.5453

## Check how many individuals with valid slopes there are:
table(is.finite(d.slopes$h.slope))

## Analyse slope in dependence of the treatments, using a model equivalent
## to the one we would use for a single time point:

mslope <- lmer(h.slope ~ block + light + div + sp + light:div 
               + light:sp + div:sp
               + (1|com) + (1|com:light) + (1|com:sp) + (1|com:light:sp) + (1|plot),
               data=d.slopes)
summary(mslope)
#| Random effects:
#|  Groups       Name        Variance Std.Dev.
#|  plot         (Intercept) 5.27e-01 7.26e-01
#|  com:light:sp (Intercept) 9.70e-02 3.11e-01
#|  com:light    (Intercept) 4.42e-17 6.65e-09
#|  com:sp       (Intercept) 5.73e-02 2.39e-01
#|  com          (Intercept) 0.00e+00 0.00e+00
#|  Residual                 3.65e+00 1.91e+00
#| Number of obs: 3863, groups:  
#| plot, 257; com:light:sp, 120; com:light, 66; com:sp, 60; com, 33

anova(mslope,type=1,ddf="Kenward-Roger")
#| Analysis of Variance Table of type I  with  Kenward-Roger 
#| approximation for degrees of freedom
#|           Sum Sq Mean Sq NumDF DenDF F.value  Pr(>F)    
#| block        119      40     3 190.4    10.9 1.2e-06 ***
#| light        111     111     1  25.0    30.5 9.7e-06 ***
#| div            0       0     1  14.0     0.1    0.74    
#| sp          3814     347    11  29.3    95.0 < 2e-16 ***
#| light:div     10      10     1  20.0     2.7    0.12    
#| light:sp     453      41    11  36.2    11.3 1.1e-08 ***
#| div:sp        37       3    11  35.5     0.9    0.53    

###
### End of script.
###
######################################################################
