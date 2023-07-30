################################################################################
##                                                                            ##
##  e-Asia Summer School: Time-series Regression for Public Health, 2023      ##
##  Session 2: Distributed lag nonlinear model                                ##
##  Last update: 03/07/2023                                                   ##
##                                                                            ##
################################################################################

# Remove all objects in current session
rm(list = ls())

# Exclude rows with missing in estimation but reinsert as NAs for prediction
options(na.action="na.exclude")

# Install packages
# install.packages("tsModel", "splines", "dlnm", "lubridate", "dplyr")

# Load useful functions
source("qAIC.R")

# Load libraries
library(tsModel); library(splines); library(dlnm)
library(lubridate); library(dplyr)

# Data input
data <- read.csv("london.csv")
names(data)
head(data)
str(data)

# Format date and day-of-week
data$date <- as.Date(data$date)
data$dow <- as.factor(data$dow)
str(data)

################################################################################
### Nonlinear exposure-response association
################################################################################

# Visual inspection
with(data, plot(tmean, all, ylim=c(100, 300),
                xlab=expression(paste("Temperature (",degree,"C)")),
                ylab="Num. of Deaths", pch=4, col="blue"))

# Create spline for season & long-term trend: natural cubic spline 6 df/year
numyears <- length(unique(data$year))
data$time <- seq(nrow(data))
nstime <- ns(data$time, df=6*numyears)

# Create basis variable for temperature spline: natural cubic spline 3 df
tmean.b <- onebasis(data$tmean, fun="ns", df=4)
attributes(tmean.b) # 3 equal-spaced internal knots at 25%, 50%, 75%

# Model with nonlinear temperature-mortality association
model1 <- glm(all ~ nstime + dow + tmean.b, data=data, family=quasipoisson)
pred.nl <- crosspred(tmean.b, model1, by=1)
plot(pred.nl, lwd=2, col="red",
     xlab=expression(paste("Temperature (",degree,"C)")),
     ylab="Relative Risk")

################################################################################
### Lag effects
################################################################################

##### Single-lag model #####
tlag <- data.frame(Lag(data$tmean, 0:27))
names.tlag <- c(sprintf("tlag%d", 0:27))
names(tlag) <- names.tlag
data <- cbind(data, tlag)
head(data,20)

# Single-lag model with linear temperature in summer (June-September)
# subset data to June - September
data.s <- subset(data, month>5 & month<10)
# create a day-of-season variable for seasonal adjustment
data.s$dos <- sequence(tapply(data.s$date, data.s$year, length))
nl <- length(names.tlag)

# matrix for storing results
mat1 <- matrix(numeric(nl*5), ncol=5, byrow=TRUE,
               dimnames=list(c(names.tlag),
                             c("logRR", "se", "RR", "CI.low", "CI.up")))

for(i in 1:nl){
  slag <- data.s[, names.tlag[i]]
  mod.sl <- glm(all ~ slag + dow + ns(dos, df=4):factor(year) + year,
                data=data.s, family=quasipoisson)
  # store estimates
  logb <- mod.sl$coefficients[2]
  se <- sqrt(diag(vcov(mod.sl)))[2]
  mat1[i,1] <- logb
  mat1[i,2] <- se
  mat1[i,3] <- exp(logb)
  mat1[i,4] <- exp(logb - 1.96*se)
  mat1[i,5] <- exp(logb + 1.96*se)
}
sldat <- data.frame(lagday=c(0:27), mat1)

# plotting the single lags
with(sldat, plot(RR~lagday, xlab="Lag day", ylab="RR",
                 main="Single-lag Model",
                 las=1, ylim=c(.98,1.04), type="n", xaxt="n"))
axis(1, at=c(0:27))
abline(h=1, lty=2, col="gray50")
with(sldat, segments(lagday, CI.low, lagday, CI.up))
with(sldat, points(lagday, RR, pch=19, col="red"))

##### Distributed lag model (DLM, unconstrained) #####
# using cross-basis
cbtemp <- with(data.s, crossbasis(tmean, lag=27,
                                  argvar=list(fun="lin"),
                                  arglag=list(fun="integer")), group=year)
summary(cbtemp) # check the basis functions
unconst <- glm(all ~ cbtemp + dow + ns(dos, df=4):factor(year) + year,
               family=quasipoisson(), data.s)

# temperature summary and percentiles in summer
summary(data.s$tmean)
sd(data.s$tmean)
round(with(data.s,
           quantile(tmean, probs=c(.01, .025, .1, .5, .9, .975, .99))), 1)

# prediction and plotting
pd.temp.hot <- crosspred(cbtemp, unconst, by=1, cen=10)
  # minimum summer tmean is about 10C, so reference is set at 10
plot(pd.temp.hot, var=24, # at about 97.5th percentile temperature
     type="p", ci="bars", pch=16, cex=1,
     col="red", las=1, xlab="Lag Day", ylab="RR",
     xlim=c(0,30), ylim=c(0.75,1.5),
     main="Unconstrained Distributed Lag Model")

# Distributed lag model (lags constrained using moving average)
cbtempma <- with(data.s, crossbasis(tmean, lag=27,
                            argvar=list(fun="lin"),
                            arglag=list(fun="strata", df=1)), group=year)

constma <- glm(all ~ cbtempma + dow + ns(dos, df=4):factor(year) + year,
               family=quasipoisson(), data.s)

pd.temp.hot <- crosspred(cbtempma, constma, by=1, cen=10)
  # minimum summer tmean is about 10C, so reference is set a 10
plot(pd.temp.hot, var=24, # at about 97.5th pct
     type="l", ci="area",
     col="red", las=1, xlab="Lag Day", ylab="RR",
     xlim=c(0,30), ylim=c(0.98,1.04),
     main="Lags Constrained using Moving Average")

##### DLM (lags constrained using natural cubic spline) #####
maxlag <- 27
lknots <- logknots(maxlag, nk=3) # set 3 knots at log scale

cbtempns <- with(data.s, crossbasis(tmean, lag=27,
                            argvar=list(fun="lin"),
                            arglag=list(fun="ns", knots=lknots)), group=year)

constns <- glm(all ~ cbtempns  + dow + ns(dos, df=4):factor(year) + year,
               family=quasipoisson(), data.s)

pd.temp.hot <- crosspred(cbtempns, constns, by=1, cen=10)

plot(pd.temp.hot, var=24, # at about 97.5th pct
     type="l", ci="area",
     col="red", las=1, xlab="Lag Day", ylab="RR",
     xlim=c(0,30), ylim=c(0.8,1.4),
     main="Lags Constrained using Natural Cubic Spline")

##### Distributed lag model using complete full-year data ####
cbtemp <- with(data, crossbasis(tmean, lag=27,
                                argvar=list(fun="lin"), # linear assumption
                                arglag=list(fun="ns", knots=lknots)))

model2 <- glm(all ~ cbtemp + nstime + dow, family=quasipoisson(), data)

# temperature summary and percentiles for whole year
summary(data$tmean)
sd(data$tmean)
round(with(data,
           quantile(tmean, prob=c(.01, .025, .1, .5, .9, .975, .99))), 1)
  # restrict the reference to between 10th & 90th percentile
  # between 4C - 19C

pd.temp <- crosspred(cbtemp, model2, by=1)
which.min(pd.temp$allRRfit) # although lowest risk at 28C, choose limit of 19C
pd.temp <- crosspred(cbtemp, model2, by=1, cen=19)

plot(pd.temp, var=22, # at about 97.5th pct
     type="l", ci="area",
     col="red", lwd=2, las=1, xlab="Lag Day", ylab="RR",
     xlim=c(0,30), ylim=c(0.95,1.1),
     main="Heat")

plot(pd.temp, var=2, # at about 2.5th pct
     type="l", ci="area",
     col="blue", lwd=2, las=1, xlab="Lag Day", ylab="RR",
     xlim=c(0,30), ylim=c(0.85,1.1),
     main="Cold")

# Check the exposure-response shape - linear
plot(crossreduce(cbtemp, model2, type="overall", cen=19),lwd=2,
     xlab=expression(paste("Temperature (",degree,"C)")), ylab="RR")

# combine plot for heat and cold lags
plot(pd.temp, var=22,
     type="l", ci="n",
     col="red", lwd=2, las=1, xlab="Lag Day", ylab="RR",
     xlim=c(0,30), ylim=c(0.85,1.1),
     main="Heat and Cold Lags")
lines(pd.temp, var=2, col="blue",lwd=2)

################################################################################
##### Distributed lag non-linear model (DLNM)
################################################################################

# Nonlinear exposure-response association using quadratic B-spline

#  temperature time series
with(data, plot(date, tmean,  xlab="Day",
                ylab=expression(paste("Temperature (",degree,"C)")),
                main="Daily Mean Temperature", pch=3, col=4))

#  daily number of deaths
with(data, plot(date, all, ylim=c(100, 300), xlab="Day", ylab="Deaths",
                main="Daily Number of Deaths", pch=4, col=2))

# create a heatwave indicator for heatwave in 2003 August 4-13
data$hw <- ifelse(data$date %within%
                    interval(as.Date("2003-08-04"), as.Date("2003-08-13")),1,0)
data$hw <- as.factor(data$hw)

# reduce number of lag days to 21
maxlag2 <- 21
lknots2 <- logknots(maxlag2, nk=3)

cbtemp <- with(data, crossbasis(tmean, lag=21,
                          argvar=list(fun="bs",degree=2, df=4),# non-linear
                          arglag=list(fun="ns", knots=lknots2)))
model2 <- glm(all ~ cbtemp + nstime + dow + hw, family=quasipoisson(), data)

# prediction and minimum mortality temperature (MMT, at which risk is minimum)
pd.temp <- crosspred(cbtemp, model2, by=1)
mmt <- as.numeric(names(which.min(pd.temp$allRRfit)))
pd.temp <- crosspred(cbtemp, model2, by=1, cen=mmt) # re-center at MMT
pct <- ecdf(data$tmean)
pct(mmt) # 85th percentile

# 3D surface plot
plot3d <- plot(pd.temp, shade=0.05, theta=235, phi=25, ltheta=-135, col="antiquewhite",
               xlab="Temperature", zlab="RR")
# lag-exposure curve at temperature 22C and 2C
lines(trans3d(x=22, y=0:21, z=pd.temp$matRRfit[as.character(22),], pmat=plot3d),
      lwd=3, col="red")
lines(trans3d(x=2, y=0:21, z=pd.temp$matRRfit[as.character(2),], pmat=plot3d),
      lwd=3, col="blue")
# exposure-response curve at lag 1
lines(trans3d(x=pd.temp$predvar, y=1, z=pd.temp$matRRfit[,"lag1"], pmat=plot3d),
      lwd=3, col="green3")

# lag-exposure
plot(pd.temp, var=22, # at about 97.5th pct
     type="l", ci="area",
     col="red", lwd=3, las=1, xlab="Lag Day", ylab="RR",
     xlim=c(0,21), ylim=c(0.95,1.2),
     main="At Temperature 22ºC (Heat)")

plot(pd.temp, var=2, # at about 12.5th pct
     type="l", ci="area",
     col="blue", lwd=3, las=1, xlab="Lag Day", ylab="RR",
     xlim=c(0,21), ylim=c(0.85,1.2),
     main="At Temperature 2ºC (Cold)")

# exposure-response
plot(pd.temp, lag=1,
     type="l", ci="area",
     col="green3", lwd=3, las=1, xlab="Temperature (ºC)", ylab="RR",
     ylim=c(0.85,1.2),
     main="At Lag 1")

# overall cumulative exposure-response
plot(crossreduce(cbtemp, model2, type="overall", cen=mmt),
     lwd=3, col=1, ci.arg=list(col="antiquewhite"),
     xlab=expression(paste("Temperature (",degree,"C)")), ylab="RR")
abline(v=c(2,mmt,22), lty=(c(2,1,2)), lwd=c(1.2,1,1.2))
text(x=c(2,mmt,22)+1.3, y=c(rep(1.88,3)), labels=c("2.5th","MMT","97.5th"))
box(lty=1)

# heat risk
with(pd.temp, cbind(allRRfit, allRRlow, allRRhigh))["22",]
# cold risk
with(pd.temp, cbind(allRRfit, allRRlow, allRRhigh))["2",]
