################################################################################
##  e-Asia Summer School: Time-series Regression for Public Health, 2023
##  Session 2: Distributed lag nonlinear model
##  Last update: 03/07/2023
################################################################################

# Remove all objects in current session
rm(list = ls())

# Install packages
# install.packages("tsModel", "splines", "dlnm", "Epi")

# Load libraries
library(tsModel); library(splines); library(dlnm); library(Epi)

# Load useful functions
source("qAIC.R")

# Load data
data <- read.csv("valencia.csv")
names(data)
head(data)
str(data)

# Format variable
data$date <- as.Date(data$date)
data$hol <- as.factor(data$hol)

# Generate time variables
data$time <- seq(nrow(data))
data$dow  <- as.factor(weekdays(data$date))

################################################################################
## Descriptive Analysis
################################################################################

# Visualize the time series
pdf("Figure 1.pdf", width=12, height=9)
par(mex=0.8,mfrow=c(2,1))
  # All-cause mortality
  plot(data$date, data$all, xaxt="n",
       type="l", col= "red",
       ylim=c(0,50),
       ylab="Num. of deaths", xlab="Date")
  axis(1, at =data$date[c(1,diff(data$year))==1], labels = 2000+1:7)
  # Mean temperature
  plot(data$date, data$tmean, xaxt="n",
       type="l", col= "blue",
       ylim=c(0,40),
       ylab="Temperature (ºC)", xlab="Date")
  axis(1, at =data$date[c(1,diff(data$year))==1], labels = 2000+1:7)
dev.off()
layout(1)

pdf("Figure 2.pdf", width=6, height=5)
plot(data$tmean, data$all, ylim=c(0, 50),
     xlab="Temperature (ºC)",
     ylab="Num. of Deaths", pch=4, col="red")
dev.off()

# temperature summary and percentiles
summary(data$tmean)
sd(data$tmean)
round(with(data,quantile(tmean,
                         probs=c(.01, .025, .1, .5, .9, .975, .99))), 1)

################################################################################
## Simple Time-series Regression Model
################################################################################

# Linear model
  # number of years in the data
  numyears <- length(unique(data$year))
  # Natural cubic splines for seasonal and long-term trend with 6df per year
  nstime <- ns(data$time, df=6*numyears)

  modlin <- glm(all ~ tmean + nstime + dow + hol,
              family=quasipoisson(), data=data, na.action="na.exclude")
  summary(modlin)
  ci.exp(modlin, subset="tmean")
  # What does the result mean?

# Moving average models
  # 7-day moving average of mean temperature
  data$tmean06 <- runMean(data$tmean, lags=0:6)

  modma3 <- glm(all ~ tmean06 + nstime + dow + hol,
                family=quasipoisson(), data=data, na.action="na.exclude")
  summary(modma3)
  ci.exp(modma3, subset="tmean06")

  # 2-day and 6-day moving averages of mean temperature
  data$tmean01 <- runMean(data$tmean, lags=0:1)
  data$tmean16 <- runMean(data$tmean, lags=1:6)

  modma3b <- glm(all ~ tmean01 + tmean16 + nstime + dow + hol,
                family=quasipoisson(), data=data, na.action="na.exclude")
  summary(modma3b)
  ci.exp(modma3b, subset="tmean01")
  ci.exp(modma3b, subset="tmean16")

  # Can you explain why the results change when two types of moving averages
  # are added?

################################################################################
## Distributed Lag Model (DLM)
################################################################################

# Distributed lag model (unconstrained)
  maxlag <- 21 # set number of lags, excluding current day lag 0
  cbtmean.unc <- crossbasis(data$tmean, lag=maxlag,
                            argvar=list(fun="lin"),
                            arglag=list(fun="integer"))
  summary(cbtmean.unc)
  dim(cbtmean.unc)

  mod.unc <- update(modma3, .~. -tmean06+cbtmean.unc)
  predtmean.unc <- crosspred(cbtmean.unc, mod.unc, by=1, cen=18)
  # when "lin" is used, the default reference value is 0
  # here we arbitrarily selected the median temperature 18C as the reference

  # for 1C increase above reference, var=19
  plot(predtmean.unc, var=19,
       type="p", ci="bars", pch=16, cex=1,
       col="red", las=1, xlab="Lag Day", ylab="RR",
       xlim=c(0,21), ylim=c(0.98,1.02),
       main="Unconstrained Distributed Lag Model")

  # Can you describe the patterns in the plot?

# Distributed lag model (constrained using natural cubic spline)
  cbtmean.spl <- crossbasis(data$tmean, lag=maxlag,
                            argvar=list(fun="lin"),
                            arglag=list(fun="ns", df=5)) # 3 knots
  mod.spl <- update(modma3, .~. -tmean06+cbtmean.spl)
  predtmean.spl <- crosspred(cbtmean.spl, mod.spl, by=1, cen=18)
  # for 1C increase above reference, var=19
  plot(predtmean.spl, var=19,
       type="p", ci="bar", pch=16, cex=1,
       col="red", las=1, xlab="Lag Day", ylab="RR",
       xlim=c(0,21), ylim=c(0.99,1.01),
       main="Lags Constrained using Moving Average")

  # Can you describe what happen after applying spline constraint to lags?

  # Obtain RR for 1C increase above reference=18C at lag 0
  predtmean.spl$matRRfit["19","lag0"]
  predtmean.spl$matRRlow["19","lag0"]; predtmean.spl$matRRhigh["19","lag0"]

################################################################################
## Distributed Lag Non-linear Model (DLNM)
################################################################################
# Distributed lag non-linear model
  # Number of lag days
  maxlag <- 21
  # number of years in the data
  numyears <- length(unique(data$year))
  # Natural cubic splines for seasonal and long-term trend with df per year
  nstime <- ns(data$time, df=6*numyears)

  # Crossbasis for temperature
  varknot <- equalknots(data$tmean, nk=2) # 2 equally spaced knots
  lknots <- logknots(maxlag, 3) # 3 knots at log scale

  cbtemp <- crossbasis(data$tmean, lag=maxlag,
                       argvar=list(fun="bs", degree=2, knots=varknot),
                       arglag=list(fun="ns", knots=lknots))
  summary(cbtemp)

  # quasi-Poisson regression model - main model
  mod <- glm(all ~ cbtemp + nstime + dow + hol, data=data,
             family=quasipoisson(), na.action="na.exclude")
  pdtemp <- crosspred(cbtemp, mod, by=1)

  # Get prediction centered at the MMT (re-center)
  mmt <- pdtemp$predvar[which.min(pdtemp$allRRfit)]
  pdtempcen <- crosspred(cbtemp, mod, cen=mmt, by=1)

  # Plot exposure-lag-response risk surface
  pdf("Figure 3.pdf", width=8, height=8)
    plot3d <- plot(pdtempcen, shade=0.05, theta=235, phi=30, ltheta=-135,
                   col="antiquewhite", xlab="Temperature", zlab="RR",
                   main="Exposure-lag-response risk surface")
    # lag-exposure curve at temperature 28C and 7C
    lines(trans3d(x=28, y=0:maxlag, z=pdtempcen$matRRfit["28",],
                  pmat=plot3d), lwd=2.5, col="red")
    lines(trans3d(x=7, y=0:maxlag, z=pdtempcen$matRRfit["7",],
                  pmat=plot3d), lwd=2.5, col="blue")
  dev.off()

  # Plot lag-response curve at specific temperature
  pdf("Figure 4.pdf", width=12, height=5)
  par(mfrow=c(1,2))
    plot(pdtempcen, var=28, type="l", ci="area",
         col="red", lwd=2, las=1, xlab="Lag", ylab="RR",
         xlim=c(0,maxlag), ylim=c(0.95,1.2), main="Lag-response at 28ºC")

    plot(pdtempcen, var=7, type="l", ci="area",
         col="blue", lwd=2, las=1, xlab="Lag", ylab="RR",
         xlim=c(0,maxlag), ylim=c(0.85,1.2), main="Lag-response at 7ºC")
  dev.off()

  # Plot overall cumulative exposure-response curve
  pdf("Figure 5.pdf", width=8, height=6)
    plot(pdtempcen, "overall", lwd=2, col=1, ci.arg=list(col="antiquewhite"),
         xlab="Temperature (ºC)", ylab="RR",
         main="Overall cumulative exposure-response")
    abline(v=c(7,mmt,28), lty=(c(2,1,2)), lwd=c(1.2,1,1.2))
    text(x=c(7,mmt,28)+1.3, y=c(rep(2.7,3)), labels=c("2.5th","MMT","97.5th"))
  dev.off()

  names(pdtempcen)

  # Extract the effect estimates
  # heat risk
    with(pdtempcen, cbind(allRRfit, allRRlow, allRRhigh))["28",]
  # cold risk
    with(pdtempcen, cbind(allRRfit, allRRlow, allRRhigh))["7",]

  # Model fit
  # qAIC
    qAIC(mod, type="dev")
  # Sum of absolute PACF of residuals
    res.mod <- residuals(mod, type="response")
    pacf.mod <- pacf(res.mod, na.action=na.omit)
    sum(abs(pacf.mod$acf))

  # Sensitivity to number of knots with a different placement
    pct <- quantile(data$tmean, prob=c(.10, .75, .90), na.rm=TRUE)
    varknot <- pct[c(1,2,3)] # knots at 10th, 75th, 90th percentile

    cbtemp2 <- crossbasis(data$tmean, lag=maxlag,
                         argvar=list(fun="bs", degree=2, knots=varknot),
                         arglag=list(fun="ns", knots=lknots) )
    summary(cbtemp2)

    # quasi-Poisson regression model
      mod2 <- glm(all ~ cbtemp2 + nstime + dow + hol, data=data,
                 family=quasipoisson(), na.action="na.exclude")
      pdtemp2 <- crosspred(cbtemp2, mod2, by=1)

      # Model fit
      # qAIC
        qAIC(mod2, type="dev")
      # Sum of absolute PACF of residuals
        res.mod2 <- residuals(mod2, type="response")
        pacf.mod2 <- pacf(res.mod2, na.action=na.omit)
        sum(abs(pacf.mod2$acf))

      # Get prediction centered at the MMT
        mmt2 <- pdtemp2$predvar[which.min(pdtemp2$allRRfit)]
        pdtempcen2 <- crosspred(cbtemp2, mod2, cen=mmt2, by=1)

      # Compare main model to alternative model
      plot(pdtempcen, "overall", lwd=3, col=1, ci.arg=list(col="antiquewhite"),
           xlab="Temperature (ºC)", ylab="RR",
           main="Overall cumulative exposure-response")
      lines(pdtempcen2, "overall", lwd=2.5, lty=2, col=4, ci="area",
            ci.arg=list(density=20, col=4, angle=75))

  # Compare the two curves. Can you identify the parts of curves that differ
  # between the two models? Why?

  # Re-center
    # Heat (99th vs 90th percentile) based on the second model
      # re-center at 90th percentile
      pdtempcen3 <- crosspred(cbtemp2, mod2, cen=26, by=1)
      with(pdtempcen3, cbind(allRRfit, allRRlow, allRRhigh))["29",]
    # Cold (1st vs 10th percentile)
      # re-center at 10th percentile
      pdtempcen4 <- crosspred(cbtemp2, mod2, cen=10, by=1)
      with(pdtempcen4, cbind(allRRfit, allRRlow, allRRhigh))["6",]

  # Compare these results to those obtained above in lines 216-218?
  # Can you explain the difference?

  # Reducing a DLNM to one dimension, use the main model
    crall <- crossreduce(cbtemp, mod, cen=mmt)
      plot(crall, xlab="Temperature", ylab="RR",
           main="Overall cumulative exposure-response")
    crlag <- crossreduce(cbtemp, mod, type="lag", value=5, cen=mmt)
      plot(crlag, xlab="Temperature", ylab="RR",
           main="Exposure-response at lag 5")
    crvar <- crossreduce(cbtemp, mod, type="var", value=28, cen=mmt)
      plot(crvar, xlab="Lag", ylab="RR", main="Lag-response at 28ºC")

# end
