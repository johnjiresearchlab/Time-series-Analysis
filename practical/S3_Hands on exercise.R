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

## ################################################################################
## # DEMO 1 PM10 AND MORTALITY

# THE MOVING AVERAGE OF PM10, TEMPERATURE, AND HUMIDITY
lndn$pm1001 <- runMean(lndn$pm10, lags=0:1)
lndn$tmean04 <- runMean(lndn$tmean, lags=0:4)
lndn$rh03 <- runMean(lndn$rh, lags=0:3)

# FIT THE REGRESSION MODEL
fitpm<-glm(all~pm1001+ns(tmean04,df=6)+ns(rh03,df=3)+ns(date,df=8*length(unique(year)))+ dow,
           family=quasipoisson(), lndn, na.action="na.exclude")

# EXTRACT THE ESTIMATES
ind <- grep("pm1001", names(coef(fitpm)))
coefall <- coef(fitpm)[ind]
vcovall <- vcov(fitpm)[ind,ind]

#Q1: What is your interpretation of the estimates (coef) above?

# DEFINE ONEBASIS FOR PM10 (I.E., X0, CONTERFACTURAL EXPOSURE LEVEL)
onepm10 <- onebasis(lndn$pm10, "lin")   # CONTERFACTUAL SCENARIO PM10 LEVEL (X0) IS 0
# onepm10<- onebasis(pmax(lndn$pm10-45,0), "lin") # ALTERNATIVE CONTERFACTUAL SCENARIO PM10 LEVEL (X0) IS WHO GUIDELINE 45µg/m3

# FORWARD MOVING AVERAGE OF DEATHS
y <- rowMeans(as.matrix(Lag(lndn$all, -1:0)))

# COMPUTE THE EXCESS DEATHS
anday <- (1-exp(-onepm10%*%coefall))*y
antot=sum(anday,na.rm=T)

# SIMULATE THE ASSUMED NORMAL DISTRIBUTION OF THE ESTIMATED COEFFICIENTS FOR PM10 (I.E.,"COEFALL")
set.seed(12345)
coefsim <- mvrnorm(1000, coefall, vcovall)

# SIMULATED DISTRIBUTION OF EXCESS DEATHS
ansim<- lapply(coefsim, function(b) {
  anday<- (1-exp(-onepm10%*%b))*y # Daily AN
  tot<-sum(anday, na.rm=T)         # Total AN for the study period
})

andata<- as.data.table(ansim)

# TOTAL AN WITH 95% EMPIRICAL CI
result_an<-c(est=antot,ci.low=quantile(andata,0.025,na.rm=T),ci.high=quantile(andata,0.975,na.rm=T))
# TOTAL AF WITH 95% EMPIRICAL CI
n<-sum(lndn$all,na.rm = T)
result_af<-c(est=antot/n,ci.low=quantile(andata,0.025,na.rm=T)/n,ci.high=quantile(andata,0.975,na.rm=T)/n)

# Q2: What is the attributable risk of PM10 to mortality in London between 2002 and 2006?

## ################################################################################
## # DEMO 2 TEMPERATURE AND MORTALITY

# DERIVE THE CROSS-BASIS FOR TEMPERATURE

# SPECIFICATION OF THE EXPOSURE FUNCTION
varfun = "bs"
vardegree = 2
varper <- c(10,75,90)

# SPECIFICATION OF THE LAG FUNCTION
lag <- 21
lagnk <- 3

# DEGREE OF FREEDOM FOR SEASONALITY
dfseas <- 8

# DEFINE THE CROSSBASIS
argvar <- list(fun=varfun,knots=quantile(lndn$tmean,varper/100,na.rm=T),
               degree=vardegree)
cb <- crossbasis(lndn$tmean,lag=lag,argvar=argvar,
                 arglag=list(knots=logknots(lag,lagnk)))

# RUN THE MODEL
model <- glm(all~cb+ns(date,7*length(unique(year)))+dow,
             family=quasipoisson(),lndn, na.action="na.exclude")

# PREDICT EFFECT SUMMARY FROM DLNM
cptmean <- crosspred(cb, model, cen=20)

# VISUALIZE THE TEMPERATURE-MORTALITY ASSOCIATION
# 3-D PLOT
d3<-plot(cptmean, xlab="Temperature", zlab="RR", phi=35, theta=205, ltheta=170,
           main="Exposure-lag-response risk surface")
lines(trans3d(x=25, y=seq(0, 21), z=cptmean$matRRfit["25",], pmat=d3),
      col=2, lwd=2)
lines(trans3d(x=cptmean$predvar, y=4, z=cptmean$matRRfit[,"lag4"], pmat=d3),
      col=3, lwd=2)

# PLOT CUMULATIVE ASSOCIATION
plot(cptmean, "overall", col=4, ylab="RR", xlab="Temperature", ylim=c(0.5,4),
     lwd=1.5, main="Overall cumulative exposure-response")

# Q3: What can you learn from the figures?

# LOAD THE FUNCTION attrdl
# (READ THE ACCOMPANYING PDF FOR DOCUMENTATION)
source("attrdl.R")

# SET THE COUNTERFACTUAL CONDITION (REFERENCE VALUE FOR EXPOSURE)
cen=20

# BACKWARD ATTRIBUTABLE RISK OF TEMPERATURE (AN AND AF)
attrdl(lndn$tmean,cb,lndn$all,model,type="an",cen=cen)
attrdl(lndn$tmean,cb,lndn$all,model,cen=cen)

# FORWARD ATTRIBUTABLE RISK OF TEMPERATURE (AN AND AF)
attrdl(lndn$tmean,cb,lndn$all,model,dir="forw",type="an",cen=cen)
attrdl(lndn$tmean,cb,lndn$all,model,dir="forw",cen=cen)

# Q4: What can you learn from the results?

# WITH EMPIRICAL CONFIDENCE INTERVALS
# (NB: eCI ARE DIFFERENT AS OBTAINED EMPIRACALLY FROM RANDOM SAMPLES)
# If sim=TRUE, the function computes samples of the attributable risk measures by simulating from
# the assumed normal distribution of the estimated coefficients (only implemented for total estimates).
# These samples can be used to defined empirical confidence intervals.
quantile(attrdl(lndn$tmean,cb,lndn$all,model,sim=T,nsim=1000,cen=cen),c(0.025,0.975))

# ATTRIBUTABLE FRACTION COMPONENT DUE TO HEAT
attrdl(lndn$tmean,cb,lndn$all,model,cen=cen,range=c(cen,100))*100

# Q5: Can you derive AN and AF due to heat with a reference level at 25℃? Is it the same with the results above?

# ATTRIBUTABLE FRACTION COMPONENT DUE TO MILD COLD
perc1<-quantile(lndn$tmean,0.01)
attrdl(lndn$tmean,cb,lndn$all,model,cen=cen,range=c(perc1,cen))*100

# DAILY ATTRIBUTABLE DEATHS DUE TO COLD IN SECOND MONTH, FORWARD & BACKWARD
attrdl(lndn$tmean,cb,lndn$all,model,tot=F,type="an",cen=cen,range=c(-100,cen))[31:60] # [31:60] RETURNS RESULTS IN SECOND MONTH
attrdl(lndn$tmean,cb,lndn$all,model,tot=F,type="an",dir="forw",cen=cen,
       range=c(-100,20))[31:60]



