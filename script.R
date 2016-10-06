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
###
### We provide two versions of this script. This version contains the
### code with important output inserted as comment. The output often
### is abbreviated shortened to save space. Output is marked by lines
### beginning with '#|'. This version of the script will allow studying
### code and results without a computer at hand.

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

for(d in c("figures"))                # create directory if it does not exist
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
#|    block    plot         uind center pool com div fdiv light  ind sp time   bd height cg ch cl cm cs dh ed lg pm qs sm ss     ba
#| 14    B1 B1.Xa01 B1.Xa01|i-01 buffer    X  ss   1   D1     c i-01 ss   14 15.0    104  0  0  0  0  0  0  0  0  0  0  0  1 1.7671
#| 31    B1 B1.Xa01 B1.Xa01|i-02 buffer    X  ss   1   D1     c i-02 ss   14  9.0     74  0  0  0  0  0  0  0  0  0  0  0  1 0.6362
#| 48    B1 B1.Xa01 B1.Xa01|i-03 buffer    X  ss   1   D1     c i-03 ss   14 11.0     76  0  0  0  0  0  0  0  0  0  0  0  1 0.9503
#| 65    B1 B1.Xa01 B1.Xa01|i-04 buffer    X  ss   1   D1     c i-04 ss   14 12.5     85  0  0  0  0  0  0  0  0  0  0  0  1 1.2272
#| 82    B1 B1.Xa01 B1.Xa01|i-05 buffer    X  ss   1   D1     c i-05 ss   14 11.5     91  0  0  0  0  0  0  0  0  0  0  0  1 1.0387
#| 99    B1 B1.Xa01 B1.Xa01|i-06 center    X  ss   1   D1     c i-06 ss   14 10.5    101  0  0  0  0  0  0  0  0  0  0  0  1 0.8659


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
#|             Df Sum Sq Mean Sq F value  Pr(>F)    
#| light        1   2824    2824   34.44 8.5e-08 ***
#| fdiv         2     39      20    0.24    0.79    
#| Residuals   84   6888      82                    


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
#|             Df Sum Sq Mean Sq
#| com         10   5242     524
#| com:light   11   3270     297
#| plot        66   1239      19

## Since 'plot' is equivalent to the residual, one can also omit 'plot':
lm2 <- aov(ba ~ com + light:com,
           data = d14X.plot)
summary(lm2)
#|             Df Sum Sq Mean Sq F value  Pr(>F)    
#| com         10   5242     524    27.9 < 2e-16 ***
#| com:light   11   3270     297    15.8 1.4e-14 ***
#| Residuals   66   1239      19                    

## Treatment and error model combined:
lm3 <- aov(ba ~ light + fdiv + com + light:com,
           data = d14X.plot)
summary(lm3)
#|             Df Sum Sq Mean Sq F value Pr(>F)    
#| light        1   2824    2824  150.44 <2e-16 ***
#| fdiv         2     39      20    1.04  0.358    
#| com          8   5202     650   34.65 <2e-16 ***
#| light:com   10    447      45    2.38  0.018 *  
#| Residuals   66   1239      19                   

## Note that in the summary of 'lm3' all statistical tests use the
## residual as error term, which is wrong here since 'fdiv' should be
## tested using 'com' as error term.
## The correct F-test for 'fdiv' can be calculated manually using 'pf':
pf(q=20/650,df1=2,df2=8,lower.tail = FALSE) 
#| 0.9698 

## Note: If you often calculate F-tests manually, then have a look at
## https://github.com/pascal-niklaus/pascal, function aov.ftest:
##
## To apply a test, provide it in the form 'fixed_effect ~ error_term'.
## Multiple tests can be provided as list:
##
## > aov.ftest(lm3,fdiv~com,table=TRUE)
#|    nom den df ddf    SS    MS    F     P s
#| 1 fdiv com  2   8 39.22 19.61 0.03 0.971  

## Model with several error strata:
lm3b <- aov(ba ~ light + fdiv + Error(com/light),
            data=d14X.plot)
summary(lm3b)
#| library(pascal)
#| 
#| Error: com
#|           Df Sum Sq Mean Sq F value Pr(>F)
#| fdiv       2     39      20    0.03   0.97
#| Residuals  8   5202     650               
#| 
#| Error: com:light
#|           Df Sum Sq Mean Sq F value  Pr(>F)    
#| light      1   2824    2824    63.2 1.2e-05 ***
#| Residuals 10    447      45                    
#| 
#| Error: Within
#|           Df Sum Sq Mean Sq F value Pr(>F)
#| Residuals 66   1239    18.8               
## This time, we obtained the correct test for 'fdiv'.

######################################################################
###
### Model mm3 (Table 3a,b)

## ANOVA using 'lme' (library: 'nlme'):
mm3 <- lme(ba ~ light + fdiv,
           random = ~1 | com/light,
           data = d14X.plot)
anova(mm3)
#|             numDF denDF F-value p-value
#| (Intercept)     1    66   27.95  <.0001
#| light           1    10   63.22  <.0001
#| fdiv            2     8    0.03  0.9704
summary(mm3)
#| Random effects:
#|  Formula: ~1 | com
#|         (Intercept)
#| StdDev:       8.701
#| 
#|  Formula: ~1 | light %in% com
#|         (Intercept) Residual
#| StdDev:       2.544    4.332
#| 

## Mixed-model analysis with 'lmer' (library: lme4)
## (the two models below are equivalent):
mm3 <- lmer(ba ~ light + fdiv + (1|com/light),
           data = d14X.plot)
mm3 <- lmer(ba ~ light + fdiv + (1|com) + (1|com:light),
           data = d14X.plot)
anova(mm3, ddf="lme4")
#| Analysis of Variance Table
#|       Df Sum Sq Mean Sq F value
#| light  1   1187    1187   63.22
#| fdiv   2      1       1    0.03
summary(mm3)
#| Random effects:
#|  Groups    Name        Variance Std.Dev.
#|  com:light (Intercept)  6.47    2.54    
#|  com       (Intercept) 75.71    8.70    
#|  Residual              18.77    4.33    
#| Number of obs: 88, groups:  com:light, 22; com, 11

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
#| Analysis of Variance Table of type I  with  Kenward-Roger 
#| approximation for degrees of freedom
#|       Sum Sq Mean Sq NumDF DenDF F.value   Pr(>F)    
#| light 1186.5  1186.5     1    10   63.22 1.24e-05 ***
#| fdiv     1.1     0.6     2     8    0.03     0.97    


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
#|             Df Sum Sq Mean Sq F value  Pr(>F)    
#| block        3    210      70    4.28  0.0082 ** 
#| light        1   2824    2824  172.86 < 2e-16 ***
#| ed           1   4137    4137  253.28 < 2e-16 ***
#| div          1    681     681   41.71 1.7e-08 ***
#| com          8    423      53    3.24  0.0038 ** 
#| light:ed     1    189     189   11.59  0.0012 ** 
#| light:div    1    140     140    8.59  0.0047 ** 
#| light:com    8    117      15    0.90  0.5261    
#| Residuals   63   1029      16                    

## If you installed library 'pascal', you can compute the F-tests as
## follows:
## aov.ftest(lm4,
##           list(block ~ Residuals,
##                light ~ light:com,
##                ed ~ com,
##                div ~ com,
##                com ~ light:com,
##                light:ed ~ light:com,
##                light:div ~ light:com,
##                light:com ~ Residuals),
##           table=TRUE)
#|         nom       den df ddf      SS      MS      F     P   s
#| 1     block Residuals  3  63  209.66   69.89   4.28 0.009  **
#| 2     light light:com  1   8 2823.64 2823.64 193.12 0.001 ***
#| 3        ed       com  1   8 4137.35 4137.35  78.26 0.001 ***
#| 4       div       com  1   8  681.37  681.37  12.89 0.008  **
#| 5       com light:com  8   8  422.95   52.87   3.62 0.044   *
#| 6  light:ed light:com  1   8  189.40  189.40  12.95 0.007  **
#| 7 light:div light:com  1   8  140.28  140.28   9.59 0.015   *
#| 8 light:com Residuals  8  63  116.97   14.62   0.90 0.527    

lm4b <- aov(ba ~ block + light * (ed + div) + Error(com + com:light),
              data = d14X.plot)
summary(lm4b)
#| Error: com
#|           Df Sum Sq Mean Sq F value  Pr(>F)    
#| ed         1   4137    4137    78.3 2.1e-05 ***
#| div        1    681     681    12.9  0.0071 ** 
#| Residuals  8    423      53                    
#| 
#| Error: com:light
#|           Df Sum Sq Mean Sq F value Pr(>F)    
#| light      1   2824    2824  193.12  7e-07 ***
#| light:ed   1    189     189   12.95  0.007 ** 
#| light:div  1    140     140    9.59  0.015 *  
#| Residuals  8    117      15                   
#| 
#| Error: Within
#|           Df Sum Sq Mean Sq F value Pr(>F)   
#| block      3    210    69.9    4.28 0.0082 **
#| Residuals 63   1029    16.3                  

mm4 <- lme(ba ~ block + light * ( ed + div ),
              random =~ 1 | com/light,
              data= d14X.plot)
anova(mm4)
#|             numDF denDF F-value p-value
#| (Intercept)     1    63 343.772  <.0001
#| block           3    63   4.329  0.0077
#| light           1     8 174.926  <.0001
#| ed              1     8  78.257  <.0001
#| div             1     8  12.888  0.0071
#| light:ed        1     8  11.733  0.0090
#| light:div       1     8   8.690  0.0185


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
#| Error: cannot allocate vector of size 27.9 Gb

## On a computer with A LOT of memory, we get:
#|                Df   Sum Sq Mean Sq 
#| uind         4105 72049602   17552 
#| ftime          16 19799497 1237469 
#| uind:ftime  49462 14706122     297                     

## With lmer, specifying the full error model does not work because
## it does not leave a residual:
m12 <- lmer(height ~ (1|uind) + (1|ftime) + (1|uind:ftime), data=d)
#| Error: number of levels of each grouping factor must be < number of observations

## However, we can omit (1|uind:ftime), a term which corresponds to the
## residual in the model below:
m <- lmer(height ~ (1|uind) + (1|ftime), data=d)
summary(m)
#| Random effects:
#|  Groups   Name        Variance Std.Dev.
#|  uind     (Intercept) 1095     33.1    
#|  ftime    (Intercept)  485     22.0    
#|  Residual              297     17.2    
#| Number of obs: 53584, groups:  uind, 4106; ftime, 17

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
#| Random effects:
#|  Groups          Name        Variance Std.Dev.
#|  plot:ftime      (Intercept)  18.2     4.26   
#|  uind            (Intercept) 600.4    24.50   
#|  com:light:ftime (Intercept)  37.2     6.10   
#|  com:ftime       (Intercept)  79.0     8.89   
#|  plot            (Intercept)  19.8     4.45   
#|  com:light       (Intercept)  15.4     3.93   
#|  com             (Intercept) 568.6    23.85   
#|  ftime           (Intercept) 450.7    21.23   
#|  Residual                    178.6    13.36   
#| Number of obs: 53584, groups:  
#| plot:ftime, 4174; uind, 4106; com:light:ftime, 1122; com:ftime, 561; plot, 257; com:light, 66; com, 33; ftime, 17

######################################################################
###
### Model mm6 (Table 5a,b)

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
#| Random effects:
#|  Groups       Name        Variance Std.Dev. Corr 
#|  plot         (Intercept) 5.56e+01 7.45e+00      
#|  com:light:sp (Intercept) 9.49e+00 3.08e+00      
#|  com:light    (Intercept) 4.52e+00 2.12e+00      
#|               time        3.21e-01 5.67e-01 -1.00
#|  com:sp       (Intercept) 5.70e-01 7.55e-01      
#|               time        1.73e-01 4.16e-01 -1.00
#|  com          (Intercept) 5.98e-08 2.45e-04      
#|               time        3.20e-10 1.79e-05 -1.00
#|  Residual                 4.41e+02 2.10e+01      
#| Number of obs: 53584, groups:  
#| plot, 257; com:light:sp, 120; com:light, 66; com:sp, 60; com, 33

## We use Sattertwaite's method to approximate degrees of freedom instead of the
## often more preferable Kenward-Roger (KR) procedure because it is less resource-
## demanding. In fact, we were not able to get ANOVA results with the KR-method
## for this model. Still, the following command will take a few minutes to run:
anova(mm6,type=1,ddf="Satterthwaite") 
#| Analysis of Variance Table of type I  with  Satterthwaite 
#| approximation for degrees of freedom
#|             Sum Sq Mean Sq NumDF DenDF F.value  Pr(>F)    
#| block        19017    6339     3 228.7      14 1.3e-08 ***
#| light         1060    1060     1  79.9       2   0.125    
#| div            443     443     1  78.6       1   0.320    
#| sp         1285714  116883    11  56.5     265 < 2e-16 ***
#| time        732090  732090     1  72.0    1659 < 2e-16 ***
#| light:div      379     379     1 246.2       1   0.355    
#| light:sp     93917    8538    11  70.7      19 < 2e-16 ***
#| div:sp        9931     903    11 102.7       2   0.031 *  
#| light:time   13739   13739     1  50.7      31 9.4e-07 ***
#| div:time         0       0     1  73.2       0   0.991    
#| sp:time     476704   43337    11  47.8      98 < 2e-16 ***


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
#| ---- Wald tests:
#|             Df denDF F.inc      Pr    
#| (Intercept)  1  17.0  7811 < 2e-16 ***
#| block        3 192.6     2 0.10774    
#| light        1  15.7    34 2.8e-05 ***
#| div          1   9.4   100 2.5e-06 ***
#| sp          11  30.6   330 < 2e-16 ***
#| light:div    1   6.9     1 0.46263    
#| light:sp    11  40.3     9 5.0e-08 ***
#| div:sp      11  32.8     4 0.00079 ***
#| time         1  18.0  4761 < 2e-16 ***
#| light:time   1  32.4    24 2.9e-05 ***
#| div:time     1   7.4     1 0.32499    
#| sp:time     11  44.0   124 < 2e-16 ***
#| 
#| ---- Stratum variances:
#|                       df Variance    com time:com com:light time:com:light       plot     com:sp time:com:sp light:com:sp R!variance
#| com               19.975  1434.23 461.58   760.76 232.60693         380.38  59.208980 265.854050    276.5784   133.324103          1
#| time:com          17.756  3827.40   0.00 55336.41  -0.94167       27657.14  -0.579519   0.818908  24517.9702     0.019734          1
#| com:light         24.047  2178.64   0.00     0.00 178.27330        -159.60  44.715758  -0.074116      7.7433    99.160672          1
#| time:com:light    33.464  5789.40   0.00     0.00   0.00000       13345.99   0.044988  -0.089690     13.7836    -1.390594          1
#| plot             189.429 12635.53   0.00     0.00   0.00000           0.00 199.145517   1.162712    -29.5327     1.275979          1
#| com:sp            14.725   992.98   0.00     0.00   0.00000           0.00   0.000000 217.410175   -193.3783   108.629180          1
#| time:com:sp       29.003  4385.26   0.00     0.00   0.00000           0.00   0.000000   0.000000  19919.0054    -0.121424          1
#| light:com:sp      22.061  7464.10   0.00     0.00   0.00000           0.00   0.000000   0.000000      0.0000   322.834717          1
#| R!variance     53179.540   441.42   0.00     0.00   0.00000           0.00   0.000000   0.000000      0.0000     0.000000          1
#| 
#| ---- Variance components:
#|                               gamma component std.error   z.ratio    constraint
#| com!com.var              0.00313961   1.38589  2.974511   0.46592 Unconstrained
#| time:com!time.var       -0.00051496  -0.22731  0.063268  -3.59289 Unconstrained
#| com:light!com.var       -0.03930374 -17.34957  5.482804  -3.16436 Unconstrained
#| time:com:light!time.var  0.00091187   0.40252  0.106054   3.79543 Unconstrained
#| plot!plot.var            0.13857405  61.16976  6.519580   9.38247 Unconstrained
#| com:sp!com.var          -0.01847604  -8.15575  3.864397  -2.11048 Unconstrained
#| time:com:sp!time.var     0.00044883   0.19812  0.057814   3.42693 Unconstrained
#| light:com:sp!light.var   0.04927928  21.75300  6.961295   3.12485 Unconstrained
#| R!variance               1.00000000 441.42291  2.707058 163.06370      Positive
#| 
#| ---- Dispersion:
#| 21.01 


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
#|                Df  Sum Sq Mean Sq F value  Pr(>F)    
#| block           3    5922    1974    3.30 0.01970 *  
#| light           1   26103   26103   43.58 4.9e-11 ***
#| ed              1  739383  739383 1234.35 < 2e-16 ***
#| dh              1  136628  136628  228.09 < 2e-16 ***
#| sm              1 1223734 1223734 2042.94 < 2e-16 ***
#| div             1  254804  254804  425.38 < 2e-16 ***
#| sp             11 3164059  287642  480.20 < 2e-16 ***
#| light:ed        1   12230   12230   20.42 6.5e-06 ***
#| light:dh        1     312     312    0.52 0.47027    
#| light:sm        1   12554   12554   20.96 4.9e-06 ***
#| light:div       1    1438    1438    2.40 0.12134    
#| light:sp       11  123386   11217   18.73 < 2e-16 ***
#| div:sp         11   39980    3635    6.07 6.8e-10 ***
#| light:div:sp   11   21638    1967    3.28 0.00017 ***
#| com            24   49942    2081    3.47 2.4e-08 ***
#| light:com      24   37143    1548    2.58 3.8e-05 ***
#| plot          188  380816    2026    3.38 < 2e-16 ***
#| sp:com          9   30246    3361    5.61 1.0e-07 ***
#| light:sp:com    9   12404    1378    2.30 0.01426 *  
#| Residuals    2777 1663438     599                    
#| 1024 observations deleted due to missingness

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
#|             nom          den  df  ddf         SS         MS      F     P   s
#| 1         block         plot   3  188    5922.24    1974.08   0.97 0.406    
#| 2         light    light:com   1   24   26102.78   26102.78  16.87 0.001 ***
#| 3            ed          com   1   24  739383.41  739383.41 355.32 0.001 ***
#| 4            dh          com   1   24  136627.78  136627.78  65.66 0.001 ***
#| 5            sm          com   1   24 1223734.32 1223734.32 588.08 0.001 ***
#| 6           div          com   1   24  254803.79  254803.79 122.45 0.001 ***
#| 7            sp       com:sp  11    9 3164058.89  287641.72  85.59 0.001 ***
#| 8      light:ed    light:com   1   24   12230.08   12230.08   7.90 0.010  **
#| 9      light:dh    light:com   1   24     312.37     312.37   0.20 0.658    
#| 10     light:sm    light:com   1   24   12553.87   12553.87   8.11 0.009  **
#| 11    light:div    light:com   1   24    1438.46    1438.46   0.93 0.345    
#| 12     light:sp light:com:sp  11    9  123386.37   11216.94   8.14 0.002  **
#| 13       div:sp       com:sp  11    9   39979.78    3634.53   1.08 0.461    
#| 14 light:sp:div light:sp:com  11    9   21638.11    1967.10   1.43 0.302    
#| 15          com    light:com  24   24   49941.52    2080.90   1.34 0.237    
#| 16    light:com         plot  24  188   37142.88    1547.62   0.76 0.779    
#| 17         plot    Residuals 188 2777  380816.15    2025.62   3.38 0.001 ***
#| 18       com:sp light:com:sp   9    9   30246.10    3360.68   2.44 0.101    
#| 19 light:com:sp    Residuals   9 2777   12404.47    1378.27   2.30 0.015   *
    
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
#|                Df  Sum Sq Mean Sq F value  Pr(>F)    
#| block           3    5922    1974    3.30 0.01970 *  
#| light           1   26103   26103   43.58 4.9e-11 ***
#| div             1    5807    5807    9.69 0.00187 ** 
#| ed              1  738340  738340 1232.61 < 2e-16 ***
#| dh              1  156312  156312  260.95 < 2e-16 ***
#| sm              1 1454090 1454090 2427.51 < 2e-16 ***
#| sp             11 3164059  287642  480.20 < 2e-16 ***
#| light:div       1    6549    6549   10.93 0.00096 ***
#| light:ed        1    9706    9706   16.20 5.8e-05 ***
#| light:dh        1    1086    1086    1.81 0.17828    
#| light:sm        1    9194    9194   15.35 9.2e-05 ***
#| light:sp       11  123386   11217   18.73 < 2e-16 ***
#| div:sp         11   39980    3635    6.07 6.8e-10 ***
#| light:div:sp   11   21638    1967    3.28 0.00017 ***
#| com            24   49942    2081    3.47 2.4e-08 ***
#| light:com      24   37143    1548    2.58 3.8e-05 ***
#| plot          188  380816    2026    3.38 < 2e-16 ***
#| sp:com          9   30246    3361    5.61 1.0e-07 ***
#| light:sp:com    9   12404    1378    2.30 0.01426 *  
#| Residuals    2777 1663438     599                    
#| 1024 observations deleted due to missingness

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
#|             nom          den  df  ddf         SS         MS      F     P   s
#| 1         block         plot   3  188    5922.24    1974.08   0.97 0.406    
#| 2         light    light:com   1   24   26102.78   26102.78  16.87 0.001 ***
#| 3            ed          com   1   24  738340.34  738340.34 354.82 0.001 ***
#| 4            dh          com   1   24  156311.51  156311.51  75.12 0.001 ***
#| 5            sm          com   1   24 1454090.24 1454090.24 698.78 0.001 ***
#| 6           div          com   1   24    5807.20    5807.20   2.79 0.108    
#| 7            sp       com:sp  11    9 3164058.89  287641.72  85.59 0.001 ***
#| 8      light:ed    light:com   1   24    9706.03    9706.03   6.27 0.020   *
#| 9      light:dh    light:com   1   24    1085.91    1085.91   0.70 0.411    
#| 10     light:sm    light:com   1   24    9193.76    9193.76   5.94 0.023   *
#| 11    light:div    light:com   1   24    6549.09    6549.09   4.23 0.051   .
#| 12     light:sp light:com:sp  11    9  123386.37   11216.94   8.14 0.002  **
#| 13       div:sp       com:sp  11    9   39979.78    3634.53   1.08 0.461    
#| 14 light:sp:div light:sp:com  11    9   21638.11    1967.10   1.43 0.302    
#| 15          com    light:com  24   24   49941.52    2080.90   1.34 0.237    
#| 16    light:com         plot  24  188   37142.88    1547.62   0.76 0.779    
#| 17         plot    Residuals 188 2777  380816.15    2025.62   3.38 0.001 ***
#| 18       com:sp light:com:sp   9    9   30246.10    3360.68   2.44 0.101    
#| 19 light:com:sp    Residuals   9 2777   12404.47    1378.27   2.30 0.015   *

## Inspect the effect of div after adjusting for dominant species:
m <- aov(height ~ block + light + ed + dh + sm + div,
           data = d14)
summary.lm(m)
#| Call:
#| aov(formula = height ~ block + light + ed + dh + sm + div, data = d14)
#| 
#| Residuals:
#|     Min      1Q  Median      3Q     Max 
#| -108.03  -28.19   -4.58   27.47  149.48 
#| 
#| Coefficients:
#|             Estimate Std. Error t value Pr(>|t|)    
#| (Intercept)   72.901      2.361   30.88   <2e-16 ***
#| blockB2       -0.105      2.156   -0.05   0.9611    
#| blockB3        1.194      2.134    0.56   0.5758    
#| blockB4        2.447      2.140    1.14   0.2529    
#| lights        -4.910      1.532   -3.20   0.0014 ** 
#| ed            65.158      2.198   29.64   <2e-16 ***
#| dh            41.000      2.353   17.42   <2e-16 ***
#| sm            64.537      2.272   28.40   <2e-16 ***
#| div          -11.801      0.993  -11.89   <2e-16 ***
#| 
#| Residual standard error: 42.5 on 3079 degrees of freedom
#|   (1024 observations deleted due to missingness)
#| Multiple R-squared:  0.301,	Adjusted R-squared:  0.299 
#| F-statistic:  166 on 8 and 3079 DF,  p-value: <2e-16
#| 
#| ## The effect of div (slope) is -11.801 

## Inspect the effect of div before adjusting for dominant species:
m <- aov(height ~ block + light + div,
           data = d14)
summary.lm(m)
#| Call:
#| aov(formula = height ~ block + light + ed + dh + sm + div, data = d14)
#| 
#| Residuals:
#|     Min      1Q  Median      3Q     Max 
#| -108.03  -28.19   -4.58   27.47  149.48 
#| 
#| Coefficients:
#|             Estimate Std. Error t value Pr(>|t|)    
#| (Intercept)   72.901      2.361   30.88   <2e-16 ***
#| blockB2       -0.105      2.156   -0.05   0.9611    
#| blockB3        1.194      2.134    0.56   0.5758    
#| blockB4        2.447      2.140    1.14   0.2529    
#| lights        -4.910      1.532   -3.20   0.0014 ** 
#| ed            65.158      2.198   29.64   <2e-16 ***
#| dh            41.000      2.353   17.42   <2e-16 ***
#| sm            64.537      2.272   28.40   <2e-16 ***
#| div          -11.801      0.993  -11.89   <2e-16 ***
#| 
#| Residual standard error: 42.5 on 3079 degrees of freedom
#|   (1024 observations deleted due to missingness)
#| Multiple R-squared:  0.301,	Adjusted R-squared:  0.299 
#| F-statistic:  166 on 8 and 3079 DF,  p-value: <2e-16
#| 
#| ## The effect of div (slope) is -11.801 

###
### End of script.
###
######################################################################
