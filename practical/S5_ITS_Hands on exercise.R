################################################################################
# e-Asia Summer School: Time-series Regression for Public Health, 2023
# Session 5: Interrupted time-series regression
# Last update: 09/07/2023
#
# The R codes were adapted from supplementary materials of the key reference paper:   
# Interrupted time series regression for the evaluation of public health interventions: a tutorial
# IJE 2017; J. Lopez Bernal, S. Cummins, A. Gasparrini
################################################################################

library(Epi)
library(tsModel)

# read data from csv file
data <- read.csv("sicily.csv")
head(data)


################################################################################
# Step 3: Descriptive analyses
################################################################################

# Compute the standardized rates
data$rate <- with(data, aces/stdpop*10^5)

### Scatter plot ###
# whether a linear model is appropriate
# whether there was a seasonal trend

plot(data$rate,type="n",ylim=c(0,300),xlab="Year", ylab="Std rate x 100,000", bty="l",xaxt="n")
rect(36,0,60,300,col=grey(0.9),border=F)
points(data$rate[data$smokban==0], cex=0.7) # the observed rate in the pre-intervention
axis(1, at=0:5*12, labels=F) # the x-axis (i.e. time units)
axis(1, at=0:4*12+6, tick=F, labels=2002:2006)
title("Sicily, 2002-2006")

### Summary statistics ###
summary(data)

# Simple comparisons of the outcome (ACE) before and after the smoking ban
summary(data$aces[data$smokban==0]) # 810.9
summary(data$aces[data$smokban==1]) # 857.5 

summary(data$rate[data$smokban==0]) # 213.3
summary(data$rate[data$smokban==1]) # 220.5 


################################################################################
# Step 4: Poisson regression model
################################################################################

# (A) Level change
# (B) Slope change
# (C) Level and slope change

### A Poisson regression model with the age-standardized population as an offset ###
# Used a Poisson model as the outcome is count data
# To model the count data directly (instead of the rate) and a Poisson distribution,
# incorporate the population (log transformed) as an offset variable 

# Naive model  
model.naive <- glm(aces ~ offset(log(stdpop)) + smokban, family=poisson, data)
ci.exp(model.naive)

# (A) Level change 
model.a <- glm(aces ~ offset(log(stdpop)) + time + smokban, family=poisson, data)
ci.exp(model.a)

# (B) Slope change
pre.itv <- tail( data$time[ data$smokban == 0 ], 1) # 36, the last month before the intervention
data$timeban <- ifelse(data$smokban == 0, 0, data$time - pre.itv)

model.b <- glm(aces ~ offset(log(stdpop)) + time + timeban, family=poisson, data)
ci.exp(model.b)

# (C) Level and slope change
model.c <- glm(aces ~ offset(log(stdpop)) + time + smokban + timeban, family=poisson, data)
ci.exp(model.c)

### Test which model fitted better ###
AIC(model.a, model.b, model.c)
anova(model.a, model.c, test="Chisq") # a likelihood ratio (LR) test 

# Print estimates 
round( ci.exp(model.naive), 3)
round( ci.exp(model.a), 3)
round( ci.exp(model.b), 4)
round( ci.exp(model.c), 4)

summary(model.a)
summary(model.b)
summary(model.c)

# --------------------------------------------------------------------------------
# Visualization
# --------------------------------------------------------------------------------

# (A) Level change 
pred.m.a <- predict(model.a, type="response") / data$stdpop * 10^5
pred.ctf <- predict(model.a, type="response", transform(data, smokban=0)) / data$stdpop * 10^5

plot(data$rate,type="n",ylim=c(0,300),xlab="Year", ylab="Std rate x 100,000", bty="l",xaxt="n")
rect(36,0,60,300,col=grey(0.9),border=F)
points(data$rate,cex=0.7)
axis(1,at=0:5*12,labels=F)
axis(1,at=0:4*12+6,tick=F,labels=2002:2006)
title("Sicily, 2002-2006")

lines(data$time, pred.m.a, col=2)
lines(data$time, pred.ctf, col=2, lty=2)


# (B) Slope change & (C) Level and slope change
pred.m.b <- predict(model.b, type="response") / data$stdpop * 10^5
pred.m.c <- predict(model.c, type="response") / data$stdpop * 10^5

plot(data$rate,type="n",ylim=c(0,300),xlab="Year", ylab="Std rate x 100,000", bty="l",xaxt="n")
rect(36,0,60,300,col=grey(0.9),border=F)
points(data$rate,cex=0.7)
axis(1,at=0:5*12,labels=F)
axis(1,at=0:4*12+6,tick=F,labels=2002:2006)
title("Sicily, 2002-2006")

lines(data$time, pred.m.b, col=3, lty=1, lwd=2)
lines(data$time, pred.m.c, col=4, lty=2, lwd=2)


# All three prediction lines  
plot(data$rate,type="n",ylim=c(0,300),xlab="Year", ylab="Std rate x 100,000", bty="l",xaxt="n")
rect(36,0,60,300,col=grey(0.9),border=F)
points(data$rate,cex=0.7)
axis(1,at=0:5*12,labels=F)
axis(1,at=0:4*12+6,tick=F,labels=2002:2006)
title("Sicily, 2002-2006")

lines(data$time, pred.m.a, col=2, lwd=1.5)
lines(data$time, pred.ctf, col=2, lty=2, lwd=1.5)
lines(data$time, pred.m.b, col=3, lty=1, lwd=1.5)
lines(data$time, pred.m.c, col=4, lty=2, lwd=1.5)
legend("bottomleft", c("(A) Level change", "(B) Slope change", "(C) Level and slope change"),
       lty=c(1,1,2), col=2:4, lwd=1.5, ncol=1, bty="n")


################################################################################
# Step 5: Methodological issues
################################################################################

### Adjusting for seasonality ###
# There are various ways of adjusting for seasonality 
# In this example, harmonic (Fourier) terms specifying the number of sine and cosine pairs, 
# including the length of the period (12 months)

model.a.adjs <- update(model.a, .~. + harmonic(month, nfreq=2, period=12))
summary(model.a.adjs)
# nfreq: the number of sine-cosine pairs representing different harmonics with period 12 months 
# the higher number, the higher harmonics (more fluctuations)

ci.exp(model.a.adjs, subset="smokban")

### Predict and plot of the season adjusted model ###
# "pred.adjs.de" is the de-seasonalised prediction based on the constant month of June (month=6)  
# the de-seasonalized predictions will help visualize the trends 
pred.adjs    <- predict(model.a.adjs, type="response") / data$stdpop * 10^5
pred.adjs.de <- predict(model.a.adjs, type="response", transform(data, month=6)) / data$stdpop * 10^5 

plot(data$rate,type="n",ylim=c(0,300),xlab="Year", ylab="Std rate x 100,000", bty="l",xaxt="n")
rect(36,0,60,300,col=grey(0.9),border=F)
points(data$rate,cex=0.7)
axis(1,at=0:5*12,labels=F)
axis(1,at=0:4*12+6,tick=F,labels=2002:2006)
title("(A) Level change model, adjusting for seasonality")

lines(data$time, pred.adjs, col=2)
lines(data$time, pred.adjs.de, col=3, lty=2, lwd=2)

legend("bottomleft", 
       c("Predictions w/ all 12 months", "Predictions w/ a constant month (de-seasonalized)"),
       lty=c(1,2), col=2:3, lwd=1.5, ncol=1, bty="n")

### Test which model fitted better ###
AIC(model.a, model.a.adjs)
anova(model.a, model.a.adjs, test="Chisq")


# --------------------------------------------------------
### Over-dispersion ###
model.a.adjs.q <- update(model.a.adjs, family=quasipoisson)
summary(model.a.adjs)$dispersion
summary(model.a.adjs.q)$dispersion

ci.exp(model.a.adjs.q, subset="smokban")

# Print estimates 
round( ci.exp(model.a,        subset="smokban"), 3)
round( ci.exp(model.a.adjs,   subset="smokban"), 3)
round( ci.exp(model.a.adjs.q, subset="smokban"), 3)

# end
