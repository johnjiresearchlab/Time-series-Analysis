####################################################################################
####################################################################################
###  e-Asia Summer School. Time-series Regression for Public Health
###  Session 1: Introduction to the Time-series Design
###  Last update: 10/07/2023
####################################################################################
####################################################################################

# Remove all objects in workspace.
rm(list = ls())
options(na.action="na.exclude")

# Install packages.
#install.packages("tsModel","splines","dlnm","gnm")

# Load packages.
library("tsModel"); library("splines"); library("dlnm"); library("gnm")

# Load functions.
source("qAIC.R"); source("findmin.R")

####################################################################################
###  Load London dataset
####################################################################################

# Load and inspect the dataset.
data <- read.csv('london.csv')
names(data)
head(data)

# Formatting date.
data$date <- as.Date(data$date)

# Generate time variables.
# data <- transform(data,
#                  year = year(data$date),         # Year   
#                  month = month(data$date),       # Month
#                  day = day(data$date),           # Day-of-month
#                  yday = yday(data$date),         # Day-of-year
#                  dow = wday(data$date,label=T)   # Day-of-week
# )

####################################################################################
###  Time-series plots
####################################################################################

# Figure for Time series plots.
par(mex=0.8,mfrow=c(2,1))

# Deaths. 
plot(data$date, data$all, type="l", col= "blue", ylab="Num. of Deaths", xlab="Date")

# Temperature.
plot(data$date, data$tmean, type="l", col= "red", ylab="Temperature (ºC)", xlab="Date")

# Close figure.
layout(1)

####################################################################################
###  Time-series decomposition
####################################################################################

# Figure for decomposition plots.
par(mex=0.8,mfrow=c(2,3))

# Deaths. 
all.dc <- tsdecomp(data$all, c(1,2,15,1825))
plot(data$date, all.dc[,1], xlab="Date", ylab="Long-term trend", main="Deaths")
plot(data$date, all.dc[,2], xlab="Date", ylab="Seasonal")
plot(data$date, all.dc[,3], xlab="Date",  ylab="Residual")

# Temperature.
tmean.dc <- tsdecomp(data$tmean, c(1,2,15,1825))
plot(data$date, tmean.dc[,1], xlab="Date", ylab="Long-term trend", main="Temperature (ºC)")
plot(data$date, tmean.dc[,2], xlab="Date", ylab="Seasonal")
plot(data$date, tmean.dc[,3], xlab="Date", ylab="Residual")

# Close figure.
layout(1)

####################################################################################
###  Boxplots
####################################################################################

# Figure for boxplots.
par(mex=0.8,mfrow=c(2,3))

# Deaths. 
boxplot(all~year, data=data, ylab="Num. of Deaths", xlab="Year")
boxplot(all~month, data=data,ylab="Num. of Deaths", xlab="Year")
boxplot(all~dow, data=data, ylab="Num. of Deaths", xlab="Day-of-week")

# Temperature.
boxplot(tmean~year, data=data, ylab="Temperature (ºC)", xlab="Year")
boxplot(tmean~month, data=data, ylab="Temperature (ºC)", xlab="Month")
boxplot(tmean~dow, data=data, ylab="Temperature (ºC)", xlab="Day-of-week")

# Close figure.
layout(1)

####################################################################################
###  Modeling framework for time-series regression
####################################################################################

# Overdispersion and autocorrelation of the mortality time-series
#################################################################

# Check overdispersion parameter.
model0 <- glm(all ~ 1, data=data, family=quasipoisson)
summary(model0)

# Autocorrelation plot.
pacf(data$all, lag=30, ylim=c(0,1))

# Time-stratified models
########################

# Figure for time stratified models.
par(mex=0.8,mfrow=c(3,1))

# Model with year.
model1a <- glm(all ~ factor(year), data=data, family=quasipoisson)
summary(model1a)
pred1a <- predict(model1a, type="response")

plot(data$date, data$all, ylim=c(100,300), pch=19, cex=0.5, col=grey(0.6),
    ylab="Num. of all", xlab="Date")
lines(data$date, pred1a, lwd=5, col="blue")

# Model with month.
model1b <- glm(all ~ factor(month), data=data, family=quasipoisson)
summary(model1b)
pred1b <- predict(model1b, type="response")

plot(data$date, data$all, ylim=c(100,300), pch=19, cex=0.5, col=grey(0.6),
     ylab="Num. of all", xlab="Date")
lines(data$date, pred1b, lwd=5, col="blue")

# Model with year and month.
model1c <- glm(all ~ factor(year) + factor(month), data=data, family=quasipoisson)
summary(model1c)
pred1c <- predict(model1c, type="response")

plot(data$date, data$all, ylim=c(100,300), pch=19, cex=0.5, col=grey(0.6),
     ylab="Num. of all", xlab="Date")
lines(data$date, pred1c, lwd=5, col="blue")

# Close figure.
layout(1)

# Periodic functions
####################

# Figure for periodic functions.
par(mex=0.8,mfrow=c(3,1))

# Model with one sine-cosine pairs.
data$time <- seq(nrow(data))
fourier <- harmonic(data$time, nfreq=1, period=365.25)

model2a <- glm(all ~ fourier + time, data=data, family=quasipoisson)
summary(model2a)
pred2a <- predict(model2a, type="response")

plot(data$date, data$all, ylim=c(100,300), pch=19, cex=0.5, col=grey(0.6),
     ylab="Num. of all", xlab="Date")
lines(data$date, pred2a, lwd=5, col="blue")

# Model with two sine-cosine pairs.
fourier <- harmonic(data$time, nfreq=2, period=365.25)

model2b <- glm(all ~ fourier + time, data=data, family=quasipoisson)
summary(model2b)
pred2b <- predict(model2b, type="response")

plot(data$date, data$all, ylim=c(100,300), pch=19, cex=0.5, col=grey(0.6),
     ylab="Num. of all", xlab="Date")
lines(data$date, pred2b, lwd=5, col="blue")

# Model with four sine-cosine pairs.
fourier <- harmonic(data$time, nfreq=4, period=365.25)

model2c <- glm(all ~ fourier + time, data=data, family=quasipoisson)
summary(model2c)
pred2c <- predict(model2c, type="response")

plot(data$date, data$all, ylim=c(100,300), pch=19, cex=0.5, col=grey(0.6),
     ylab="Num. of all", xlab="Date")
lines(data$date, pred2c, lwd=5, col="blue")

# Close figure.
layout(1)

# Spline functions
##################

# Figure for splines.
par(mex=0.8,mfrow=c(3,1))

# Model with natural cubit splits with 1 df/year.
numyears <- length(unique(data$year))
spl <- ns(data$time, df=1*numyears)

model3a <- glm(all ~ spl , data=data, family=quasipoisson)
summary(model3a)
pred3a <- predict(model3a, type="response")

plot(data$date, data$all, ylim=c(100,300), pch=19, cex=0.5, col=grey(0.6),
     ylab="Num. of all", xlab="Date")
lines(data$date, pred3a, lwd=5, col="blue")

# Model with natural cubit splits with 6 df/year.
spl <- ns(data$time, df=6*numyears)

model3b <- glm(all ~ spl , data=data, family=quasipoisson)
summary(model3b)
pred3b <- predict(model3b, type="response")

plot(data$date, data$all, ylim=c(100,300), pch=19, cex=0.5, col=grey(0.6),
     ylab="Num. of all", xlab="Date")
lines(data$date, pred3b, lwd=5, col="blue")

# Model with natural cubit splits with 12 df/year.
spl <- ns(data$time, df=12*numyears)

model3c <- glm(all ~ spl , data=data, family=quasipoisson)
summary(model3c)
pred3c <- predict(model3c, type="response")

plot(data$date, data$all, ylim=c(100,300), pch=19, cex=0.5, col=grey(0.6),
     ylab="Num. of all", xlab="Date")
lines(data$date, pred3c, lwd=5, col="blue")

# Close figure.
layout(1)

####################################################################################
###  Comparing modeling strategies
####################################################################################

# Joint figure for Autocorrelation functions.
par(mex=0.8,mfrow=c(2,2))

# Unadjusted model.
qaic0 <- qAIC(model0, type="dev")
disp0 <- sqrt(summary(model0)$dispersion)
res0 <- residuals(model0, type="response")
pacf(res0, na.action=na.omit, ylim=c(0,1), 
     main=paste("Unadjusted | ", "qAIC=" , round(qaic0,1) , ", overdisp=" , round(disp0,2)))

# Time-stratified model.
qaic1 <- qAIC(model1b, type="dev")
disp1 <- sqrt(summary(model1c)$dispersion)
res1 <- residuals(model1c, type="response")
pacf(res1, na.action=na.omit, ylim=c(0,1), 
    main=paste("Time stratified | ", "qAIC=" , round(qaic1,1) , ", overdisp=" , round(disp1,2)))

# Periodic functions.
qaic2 <- qAIC(model2c, type="dev")
disp2 <- sqrt(summary(model2c)$dispersion)
res2 <- residuals(model2c, type="response")
pacf(res2, na.action=na.omit, ylim=c(0,1), 
    main=paste("Periodic functions | ", "qAIC=" , round(qaic2,1) , ", overdisp=" , round(disp2,2)))

# Spline functions.
qaic3 <- qAIC(model3b, type="dev")
disp3 <- sqrt(summary(model3b)$dispersion)
res3 <- residuals(model3b, type="response")
pacf(res3, na.action=na.omit, ylim=c(0,1), 
    main=paste("Splines | ", "qAIC=" , round(qaic3,1) , ", overdisp=" , round(disp3,2)))

# Close figure.
layout(1)

####################################################################################
###  Exposure-response
####################################################################################

# Linear association.
tmean.b <- onebasis(data$tmean, fun="lin")
model <- glm(all ~ spl + factor(dow) + tmean.b, data=data, family=quasipoisson)
summary(model)
pred.lin <- crosspred(tmean.b, model, by=1)
plot(pred.lin)

# Linear association centered at MMT.
pred.lin <- crosspred(tmean.b, model, by=1, cen=min(data$tmean))
plot(pred.lin)

# Non-linear association.
tmean.b <- onebasis(data$tmean, fun="ns", df=4)
model <- glm(all ~ spl + factor(dow) + tmean.b, data=data, family=quasipoisson)
summary(model)
pred.nl <- crosspred(tmean.b, model, by=1)
plot(pred.nl)

# Non-linear association centered at MMT.
mmt <- findmin(tmean.b, model)
pred.nl <- crosspred(tmean.b, model, by=1, cen=mmt)
plot(pred.nl)

# Close figure.
layout(1)

####################################################################################
###  Case-crossover
####################################################################################

# Generate time-stratified strata.
data$month <- as.factor(data$month)
data$year  <- as.factor(data$year)
data$dow   <- as.factor(data$dow)
data$stratum <- with(data, as.factor(year:month:dow))

# Fit fixed-effects conditional quasi-Poisson regression.
model <- gnm(all ~ tmean.b, data=data, family=quasipoisson, eliminate=stratum)
mmt <- findmin(tmean.b, model)
pred.cc <- crosspred(tmean.b, model, by=1, cen=mmt)
plot(pred.cc)

# Joint figure for Autocorrelation functions.
par(mex=0.8,mfrow=c(1,2))

# Time-series.
plot(pred.nl, main="Time-series")

# Case-crossover.
plot(pred.cc, main="Case-crossover")

# Close figure.
layout(1)

####################################################################################
####################################################################################
###  End of script
####################################################################################
####################################################################################
