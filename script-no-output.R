######################################################################
###
### R script to run analyses presented in paper
###
### Pascal A. Niklaus & Bernhard Schmid
###
### History:
### - 20-04-2016 file created
### - 30-05-2016 version for final test
### - 07-06-2016 minor corrections
### - 16-09-2016 final revision
### - 31-10-2017 comment added with respect to mm6

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
###
### We provide two versions of this script:
###
### script-no-output.R:
###   contains the code only
###
### script.R:
###   containts code amended with important output that is inserted
###   abbreviated to save space. The purpouse of this file is to
###   provide an opportunity to study the code and output without
###   a computer at hand.

######################################################################
###
### Required R libraries and extra software
###
### Libraries can be installed from the command line with
### 'install.packages(...)'.
###
### Libraries for analysis by maximum likelihood:
### - nlme, lmer, lmerTest, pbkrtest
###
### The library 'pascal' provided by one of the authors via github
### contains some convenience functions that are used in this script. 
### It can easily be installed as follows:
###   install.packages("devtools")
###   library(devtools)          
###   install_github("pascal-niklaus/pascal/pascal")

library(lmerTest)
library(lme4)                         # not really needed since loaded by lmerTest
library(nlme)
options(digits=4)                     # set number of post-decimal digits in output to 4
                                      # increased to show more digits if required
rm(list=ls())                         # clear workspace

asr <- require(pascal) &              # load ASReml, if installed, plus convenience functions
       require(asreml)                # if both libraries were loaded, asr == TRUE

for(d in c("derived_data","figures")) # create directories if these don't exist
    if(!dir.exists(d))                # data will be stored there
        dir.create(d)
    

######################################################################
###
### Based on a single original file, the following derived data files
### were prepared:
###
### pilot_bd_all_pools.csv
###    original file with basal diameter and height of trees
###    in all pools, at all times.
###
### pilot_pool_X_plot_ba.csv
###    sums of basal area for pool X and time == 14.
###
### pilot_height_slopes.csv
###    slopes of regression height~time,
###    for individual trees of all pools.

d <- read.csv("derived_data/pilot_bd_all_pools.csv")       # read in data file

## create data set with plot sums, for pool X, time == 14:

d14X <- subset(d, time == 14 & pool == "X" )  # use only single time point

d14X$ba <- ( (d14X$bd/10) / 2)^2 * pi         # calculate basal [cm^2] area from diameter [mm]

d14X.plot <- aggregate(                       # aggregate at plot level (sums)
    ba ~ block + plot + div + fdiv + light + com +
        cg + ch + cl + cm + cs + dh + ed + lg + pm + qs + sm + ss,
    data=d14X,
    FUN=function(x) sum(x,na.rm=TRUE))

write.csv(d14X.plot,                          # save result for reference
          "derived_data/pilot_pool_X_plot_ba.csv",
          row.names=FALSE)

## Structure of data set in 'pilot_pool_X_plot_ba.csv':


######################################################################
###
### Plots of figure 2

d14X.plot <- read.csv("derived_data/pilot_pool_X_plot_ba.csv")

## Fig. 2: Plot basal area by diversity separately for control and
## shade treatment, using different symbols for the plots containing
## Elaeocarpus decipiens.

pdf("figures/fig_2.pdf",width=5,height=5)

## Solution with 'lattice' panel functions...
print(
    xyplot(ba ~ jitter(div) | light,
           xlab="Species richness",
           ylab=expression("Plot basal area ("*cm^2*")"),
           groups = ed,
           par.settings = list(
               superpose.symbol = list(
                   pch=c(1,10),
                   cex=2,
                   col=c('red','blue'))),
           data=d14X.plot))

## ...or with the simple but flexible "base graphics":
par(mfrow=c(1,2), mai=c(1,0,1,0), oma=c(1,6,1,2), cex.lab = 1.5)
yrange <- range(c(0,d14X.plot$ba), na.rm=TRUE)
for(li in c("c","s")) {               # loop over two panels
    tmp <- subset(d14X.plot, light == li)
    plot(ba ~ jitter(as.numeric(fdiv)),
         xaxt = "n", yaxt = "n",      # don't plot axis ticks and labels
         pch = ifelse(tmp$ed > 0, 16, 1),
         main = if(li == "c") "light" else "shade",
         xlab = "",
         ylab = if(li == "c") expression("Plot basal area ("*cm^2*")") else "",
         ylim = yrange,
         xpd = NA,
         cex = 1.5,
         data = tmp)    
    axis(1, at = 1:3, labels=c(1,2,4)) # X-axis
    if(li == "c")                      # Y-axis only for 1st plot
        axis(2, las=1)
}
mtext("Species richness",side=1,line=-2,outer=TRUE,cex=1.5)

dev.off()

######################################################################
###
### Model lm1 (Table 1a)

## Fixed-effects model:
## We use the R-function 'aov' but one could also use 'lm' but then
## would have to call 'summary.aov(m1)' instead.
lm1 <- aov(ba ~ light + fdiv,
           data = d14X.plot)
summary(lm1)


######################################################################
###
### Models lm2/lm3 (Tables 2a/c)

## Error-only model:
## The model formula is wrapped in 'terms( ... , keep.order=TRUE)' to
## prevent plot from being moved in front of the interaction 'light:com'.
lm2 <- aov(terms(ba ~ com + light:com + plot,
                 keep.order = TRUE ),
           data=d14X.plot)
summary(lm2)

## Since 'plot' is equivalent to the residual, one can also omit 'plot':
lm2 <- aov(ba ~ com + light:com,
           data = d14X.plot)
summary(lm2)

## Treatment and error model combined:
lm3 <- aov(ba ~ light + fdiv + com + light:com,
           data = d14X.plot)
summary(lm3)

## Note that in the summary of 'lm3' all statistical tests use the
## residual as error term, which is wrong here since 'fdiv' should be
## tested using 'com' as error term.
## The correct F-test for 'fdiv' can be calculated manually using 'pf':
pf(q=20/650,df1=2,df2=8,lower.tail = FALSE) 

## Note: If you often calculate F-tests manually, then have a look at
## https://github.com/pascal-niklaus/pascal, function aov.ftest:
##
## To apply a test, provide it in the form 'fixed_effect ~ error_term'.
## Multiple tests can be provided as list:
##
## > aov.ftest(lm3,fdiv~com,table=TRUE)

## Model with several error strata:
lm3b <- aov(ba ~ light + fdiv + Error(com/light),
            data=d14X.plot)
summary(lm3b)
## This time, we obtained the correct test for 'fdiv'.

######################################################################
###
### Model mm3 (Table 3a,b)

## ANOVA using 'lme' (library: 'nlme'):
mm3 <- lme(ba ~ light + fdiv,
           random = ~1 | com/light,
           data = d14X.plot)
anova(mm3)
summary(mm3)

## Mixed-model analysis with 'lmer' (library: lme4)
## (the two models below are equivalent):
mm3 <- lmer(ba ~ light + fdiv + (1|com/light),
           data = d14X.plot)
mm3 <- lmer(ba ~ light + fdiv + (1|com) + (1|com:light),
           data = d14X.plot)
anova(mm3, ddf="lme4")
summary(mm3)

## 'lmer' does not provide P-values in the output; for details about
## the reasons see:
## https://stat.ethz.ch/pipermail/r-help/2006-May/094765.html

## To nevertheless obtain P values using standard methods of
## estimating approximate denominator degrees of freedom for F-tests,
## the libraries "lmerTest" and "pbkrtest" need to be installed and
## 'lmerTest' loaded instead of 'lme4'.

## See the following reference for a comparison of methods:
## Schaalje GB, McBride JB, Fellingham GW (2002)
## Adequacy of approximations to distributions of test statistics in
## complex mixed linear models.
## Journal of Agricultural, Biological, and Environmental Statistics, 7:512-524
## http://link.springer.com/article/10.1198/108571102726.

## One should request sequential tests (type=1) and we recommend the
## use of ddf estimated by the Kenward-Roger method:

anova(mm3, type=1, ddf="Kenward-Roger")


######################################################################
###
### Models LM4 and MM4 (Table 4a,b; Fig. 4c)

## This model includes a contrast for the presence of
## Elaeocarpus decipiens (ed):
lm4 <- aov(terms(ba ~ block
                 + light
                 + (ed + div)
                 + com
                 + light:(ed + div)
                 + light:com,
                 keep.order=TRUE),
           data = d14X.plot)
summary(lm4)

aov.ftest(lm4,
          list(block ~ Residuals,
               light ~ light:com,
               ed ~ com,
               div ~ com,
               com ~ light:com,
               light:ed ~ light:com,
               light:div ~ light:com,
               light:com ~ Residuals),
          table=TRUE)

lm4b <- aov(ba ~ block + light * (ed + div) + Error(com + com:light),
              data = d14X.plot)
summary(lm4b)

mm4 <- lme(ba ~ block + light * ( ed + div ),
              random =~ 1 | com/light,
              data= d14X.plot)
anova(mm4)


######################################################################
###
### Individual-based model

### In the following analyses, we use the entire time series. Note
### that in fact serial correlations would have to be considered,
### which we ignore here since they cannot easily be fit in 'lmer' and
### don't affect the statistical tests of interest much as long as
### interactions of time-contrasts with random effect factors are
### properly included as random-effects terms in the model.

d <- read.csv("derived_data/pilot_bd_all_pools.csv")

d$ftime <- factor(sprintf("t-%02d",d$time))  # time as factor

## Fit basic error model with tree as error:

## With lm/aov, the design matrix becomes too large (on most computers):
m12.aov <- aov(height ~ uind*ftime,data=d)

## On a computer with A LOT of memory, we get:

## With lmer, specifying the full error model does not work because
## it does not leave a residual:
m12 <- lmer(height ~ (1|uind) + (1|ftime) + (1|uind:ftime), data=d)

## However, we can omit (1|uind:ftime), a term which corresponds to the
## residual in the model below:
m <- lmer(height ~ (1|uind) + (1|ftime), data=d)
summary(m)

## (Note that here we consider 'ftime' as random term because we
##  use it to specify an error model. In other contexts, it may make
##  more sense to consider it as fixed term because its levels are
##  on an interval scale.)

######################################################################
###
### Fig. 5: Superimpose height vs. time curves for all individuals,
###  separately for each pool:

pdf("figures/fig_5.pdf",width=12,height=8)

colrs <- rainbow(12,alpha=.25)
par(mfrow=c(1,3))
for(p in sort(unique(d$pool))) {
    tmp <- subset(d, pool == p)
    sp.set <- sort(unique(tmp$sp))
    plot(NA,NA,
         xlim=range(d$time),
         ylim=c(0,max(d$height,na.rm=TRUE)),
         xlab="Time (months)",
         ylab="Height (cm)",
         main=sprintf("Pool = %s",p),
         las=1)
    for(i in unique(tmp$uind)) {
        idx <- d$uind == i
        clr <- colrs[d$sp[idx][1]]                     
        lines(tmp$height[idx]  ~ tmp$time[idx], data = d,col=clr)
    }
    legend("topleft",legend=sp.set,col=colrs[sp.set],lty=1)
}

dev.off()

######################################################################
###
### Model mm5

## Now we consider the grouping of trees in plots and mixtures, plus
## their interaction with 'light' and 'time' as additional random effects.
## (This will take a few minutes):
mm5 <- lmer(height ~ (1|com)       + (1|com:light)       + (1|plot)      + (1|uind) +
                     (1|ftime) +
                     (1|com:ftime) + (1|com:light:ftime) + (1|plot:ftime),
            data=d)
summary(mm5)   # this again will take a few minutes

######################################################################
###
### Model mm6 (Table 5a,b)
###
### Important note: For continuous random effects, both a random
### intercept and a random slope are determined. While the slope does
### not depend on the chosen origin of the scale used for the
### continous random effect, the intercept does. Hence, all tests of
### fixed effects that refer to this random intercept will depend
### on the choice of origin made. 
### In mm6 for example, tests of contrasts within com (e.g. div)
### refer to time = 0 (this is the time the random intercept refers to).
### In other words, shifting the origin of our time measurement
### will change significances in the ANOVA table for terms that 
### are related to this random intercept.

## (This will take a few minutes):
mm6 <- lmer(height ~ block
            + light
            + div
            + sp
            + light:div
            + light:sp
            + div:sp
            + time   # note that time will be moved in front of the interactions!
            + light:time
            + div:time
            + sp:time
            + (1|plot)
            + (1|com:light:sp)
            + (time|com) + (time|com:light) + (time|com:sp),
            data=d)  
summary(mm6)

## We use Sattertwaite's method to approximate degrees of freedom instead of the
## often more preferable Kenward-Roger (KR) procedure because it is less resource-
## demanding. In fact, we were not able to get ANOVA results with the KR-method
## for this model. Still, the following command will take a few minutes to run:
anova(mm6,type=1,ddf="Satterthwaite") 


## Fit the same model with ASReml, if available. In this model, we
## allow for negative estimates of variance components. This is
## relatively complicated to do directly in asreml, and we therefore
## provide a more convenient wrapper function 'asreml.nvc' in library
## 'pascal'.
if(asr) {
    mm6asr <- asreml.nvc(height ~ block
                         + light + div + sp
                         + light:div + light:sp
                         + div:sp
                         + time + light:time + div:time + sp:time,
                         random =~ com
                         + com:light
                         + com:sp
                         + light:com:sp
                         + plot 
                         + time:com
                         + time:com:light
                         + time:com:sp,
                         keep.order = TRUE,
                         control=asreml.control(maxiter = 20),
                         data=d)
    
    test.asreml(mm6asr)
}

######################################################################
###
### Non-orthogonality:

## Create data frame with last time point (month 14) only:
d <- read.csv("derived_data/pilot_bd_all_pools.csv")
d14 <- subset(d, time == 14)

## Fit model using aov:
lm7 <- aov(terms(height ~ block
                 + light
                 + ed + dh + sm
                 + div
                 + sp
                 + light:ed + light:dh + light:sm
                 + light:div
                 + light:sp
                 + div:sp
                 + light:sp
                 + light:div:sp
                 + com           # approx. error term for ed, dh, sm, div
                 + light:com     # approx. error term for light:div ... light:sm
                 + plot          # approx. error term for block, light:com
                 + com:sp        # approx. error term for div:sp
                 + light:com:sp, # approx. error term for light:div:sp (but see text)
                 keep.order=TRUE),
           data=d14)
summary(lm7)

## Tests can be calculated by choosing an approximate error term
## Because these F-tests have to be calculated manually, it is best to use
## a convenience function provided by one of the authors.
## This will only work if the library 'pascal' has been loaded.
aov.ftest(lm7,
          list(block ~ plot, light ~ light:com,
               ed ~ com, dh ~ com, sm ~ com, div ~ com,
               sp ~ com:sp,
               light:ed ~ light:com, light:dh ~ light:com, light:sm ~ light:com, light:div ~ light:com,
               light:sp ~ light:com:sp, # see text !
               div:sp ~ com:sp,
               light:sp:div ~ light:sp:com,
               com ~ light:com,
               light:com ~ plot,
               plot ~ Residuals,
               com:sp ~ light:com:sp,
               light:com:sp ~ Residuals),
          table=TRUE)

lm8 <- aov(terms(height ~ block
                 + light
                 + div
                 + ed + dh + sm
                 + sp
                 + light:div
                 + light:ed + light:dh + light:sm
                 + light:sp
                 + div:sp
                 + light:sp
                 + light:div:sp
                 + com           # approx. error term for ed, dh, sm, div
                 + light:com     # approx. error term for light:div ... light:sm
                 + plot          # approx. error term for block, light:com
                 + com:sp        # approx. error term for div:sp
                 + light:com:sp, # approx. error term for light:div:sp (but see text)
                 keep.order=TRUE),
           data=d14)
summary(lm8)

## Tests using manual method outlined above for lm7:
aov.ftest(lm8,
          list(block ~ plot, light ~ light:com,
               ed ~ com, dh ~ com, sm ~ com, div ~ com,
               sp ~ com:sp,
               light:ed ~ light:com, light:dh ~ light:com, light:sm ~ light:com, light:div ~ light:com,
               light:sp ~ light:com:sp, # see text !
               div:sp ~ com:sp,
               light:sp:div ~ light:sp:com,
               com ~ light:com,
               light:com ~ plot,
               plot ~ Residuals,
               com:sp ~ light:com:sp,
               light:com:sp ~ Residuals),
          table=TRUE)

## Inspect the effect of div after adjusting for dominant species:
m <- aov(height ~ block + light + ed + dh + sm + div,
           data = d14)
summary.lm(m)

## Inspect the effect of div before adjusting for dominant species:
m <- aov(height ~ block + light + div,
           data = d14)
summary.lm(m)


###
### End of script.
###
######################################################################
