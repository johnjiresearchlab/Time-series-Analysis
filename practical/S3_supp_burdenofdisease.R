## ################################################################################
## # S3_BURDEN OF DISEASE DUE TO ENVIRONMENTAL EXPOSURE
## # e-ASIA SUMMER SCHOOL 2023, TOKYO
## ################################################################################

## ################################################################################
## # PREPARATION

# LOAD THE PACKAGES
# install.packages(c("dlnm","splines","foreign","tsModel","lubridate"))
library(dlnm) ; library(splines) ; library(foreign) ; library(tsModel);
library(lubridate);library(MASS); library(abind);library(data.table)

# CHECK VERSION OF THE PACKAGE
if(packageVersion("dlnm")<"2.2.0")
  stop("update dlnm package to version >= 2.2.0")

# LOAD THE DATA
lndn <- read.csv("london.csv")

# INSPECT THE DATA
head(lndn, 3)

# SUMMARISE THE DATA
summary(lndn)

# CHECK THE STRUCTURE OF THE DATA
str(lndn)

# CREATE DATE
lndn$date <- as.Date(lndn$date, format="%Y/%m/%d")

# CREATE OTHER TIME VARIABLES
# lndn <- transform(lndn,
#                  year = year(lndn$date),
#                  month = month(lndn$date),
#                  day = day(lndn$date),           # DAY-OF-MONTH
#                  yday = yday(lndn$date),         # DAY-OF-YEAR
#                  dow = wday(lndn$date, label=T)  # DAY-OF-WEEK
# )

# LOAD THE FUNCTION attrdl
# (READ THE ACCOMPANYING PDF FOR DOCUMENTATION)
source("attrdl.R")
## ################################################################################
## OZONE AND MORTALITY

# PRODUCE THE CROSS-BASIS FOR OZONE (SCALING NOT NEEDED)
# A SIMPLE LINEAR TERM AND THE UNCONSTRAINED LAG STRUCTURE
cbo3unc <- crossbasis(lndn$ozone,lag=c(0,1),argvar=list(fun="lin"),
                      arglag=list(fun="integer"))

# PRODUCE THE CROSS-BASIS FOR TEMPERATURE
cb.temp <- crossbasis(lndn$tmean, lag=21,
                      argvar=list(fun="ns",knots=quantile(lndn$tmean,c(.10,.75,.90),na.rm=T)),
                      arglag=list(fun="ns",knots=logknots(21,3)))


# FIT THE REGRESSION MODEL
fit<-glm(all~cbo3unc+cb.temp+ns(date,df=7*length(unique(year)))+ dow,
           family=quasipoisson(), lndn, na.action="na.exclude")

# EXTRACT THE ESTIMATES
ind <- grep("cbo3unc", names(coef(fit)))
coefall <- coef(fit)[ind]
vcovall <- vcov(fit)[ind,ind]

# ATTRIBUTABLE NUMBERS AND FRACTION DUE TO O3 below 60𝜇𝑔/𝑚3
attrdl(lndn$ozone,cbo3unc,lndn$all,fit,type="an",dir="forw",cen=0,range=c(0,60))
attrdl(lndn$ozone,cbo3unc,lndn$all,fit,dir="forw",cen=0,range=c(0,60))*100


## ################################################################################
## TEMPERATURE AND MORTALITY

# KNOTS FOR EXPOSURE-RESPONSE FUNCTION
vk <- equalknots(lndn$tmean,fun="bs",df=4,degree=2)

# KNOTS FOR THE LAG-RESPONSE FUNCTION
maxlag <- 25
lk <- logknots(maxlag,3)

# CENTERING VALUE (AND PERCENTILE)
cen <- 28
# sum(lndn$tmean<cen,na.rm=T)/sum(!is.na(lndn$tmean))

# COMPUTE THE CROSS-BASIS FOR TEMPERATURE
cb <- crossbasis(lndn$tmean, lag=maxlag, argvar=list(fun="bs",degree=2,
                                                     knots=vk), arglag=list(knots=lk))

# RUN THE MODEL
model <- glm(all~cb+ns(date,7*length(unique(year)))+dow,family=quasipoisson(),lndn, na.action="na.exclude")

# ATTRIBUTABLE FRACTION COMPONENT DUE TO temperature above 28℃
attrdl(lndn$tmean,cb,lndn$all,model,cen=cen,range=c(cen,100))*100
attrdl(lndn$tmean,cb,lndn$all,model,type="an",cen=cen,range=c(cen,100))

# WITH EMPIRICAL CONFIDENCE INTERVALS
# (NB: eCI ARE DIFFERENT AS OBTAINED EMPIRACALLY FROM RANDOM SAMPLES)
# If sim=TRUE, the function computes samples of the attributable risk measures by simulating from
# the assumed normal distribution of the estimated coefficients (only implemented for total estimates).
# These samples can be used to defined empirical confidence intervals.
quantile(attrdl(lndn$tmean,cb,lndn$all,model,cen=cen,range=c(cen,100))*100,c(0.025,0.975))
quantile(attrdl(lndn$tmean,cb,lndn$all,model,type="an",cen=cen,range=c(cen,100)),c(0.025,0.975))


