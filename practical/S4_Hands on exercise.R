## ################################################################################
## # S4_Analysing seasonality and the role of environmental exposure
## # e-ASIA SUMMER SCHOOL 2023, TOKYO
## ################################################################################

## ################################################################################
## # PREPARATION

# LOAD THE PACKAGES
# install.packages(c("dlnm","splines","foreign","tsModel","lubridate","MASS","abind","data.table","pbs))
library(mgcv);library(tsModel);library(dlnm);library(splines);
library(mondate);library(ggplot2); library(devtools);library(pbs)

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

# CREATE SEASONALITY INDICATOR (DAY OF YEAR)
years <- unique(as.factor(as.character(lndn$year)))
for(t in years){
  tempind <- lndn$year == t
  temp <- as.numeric( strftime( lndn$date[tempind] , format = "%j") )
  if( length(temp) == 365 ){ temp[60:365] <- 61:366 }
  else if( length(temp) == 364 ){ temp[60:364] <- 61:365 }

  lndn$seasonal[tempind] <- temp
}

# Q1: What are the other options to fit seasonality?

# DERIVE THE CYCLIC SPLINE FOR SEASONALITY
spline.season <- pbs(lndn$seasonal, df = 4  )

# RUN THE MODEL
# TEMPERATURE UNADJUSTED
fit1<- glm(all ~ spline.season + year+factor(dow), data=lndn, family=quasipoisson(),na.action="na.exclude")

# TEMPERATURE ADJUSTED
# DERIVE THE CROSS-BASIS FOR TEMPERATURE
cb.temp <- crossbasis(lndn$tmean, lag=21,
                      argvar=list(fun="ns",knots=quantile(lndn$tmean,c(.25,.50,.75),na.rm=T)),
                      arglag=list(fun="ns",knots=logknots(21,3)))
fit2<- glm(all ~ spline.season + year+factor(dow) +cb.temp,  data=lndn, family=quasipoisson(),na.action="na.exclude")

# COEFFICIENTS FROM CYCLIC SPLINE
ind <- grep("spline.season", names(coef(fit1)))
coef1<- coef(fit1)[ ind ]
vcov1<- vcov(fit1)[ ind,ind]
coef2<- coef(fit2)[ ind ]
vcov2<- vcov(fit2)[ ind,ind]

# COEFFICIENTS FOR EACH DAY OF YEAR (LOGRR, SE)
logrr1 <- matrix(NA, 1, 366 )
logrr_se1<- matrix(NA, 1, 366 )
logrr2 <- matrix(NA, 1, 366 )
logrr_se2 <- matrix(NA, 1, 366 )

Trange <- 1:366

bvar<- pbs( lndn$seasonal, df = 4  )

# LOGRR
l <- length(Trange)
for(j in 1:l){ logrr1[1,j] <-  t(bvar[j,]) %*% as.numeric( coef1)  }
for(j in 1:l){ logrr2[1,j] <-  t(bvar[j,]) %*% as.numeric( coef2)  }

# SE OF LOGRR
for(j in 1:l){ logrr_se1[1,j] <- sqrt( as.numeric( t(bvar[j,]) %*% vcov1 %*% (bvar[j,]) ) ) }
for(j in 1:l){ logrr_se2[1,j] <- sqrt( as.numeric( t(bvar[j,]) %*% vcov2%*% (bvar[j,]) ) ) }

# CENTERING LOGRR AT TROUGH
mincen<- apply(logrr1, 1, which.min )
mincen<- apply(logrr2, 1, which.min )

bvar_cen<- bvar[mincen,]
bvar_cen2<- bvar[mincen,]

# CENTERED LOGRR
l <- length(Trange)
for(j in 1:l){ logrr1[1,j] <- as.numeric( t(bvar[j,]-bvar_cen) %*% as.numeric( coef1) ) }
for(j in 1:l){ logrr2[1,j] <- as.numeric( t(bvar[j,]-bvar_cen2) %*% as.numeric( coef2) )}

# SD FOR CENTERED LOGRR
for(j in 1:l){ logrr_se1[1,j] <- sqrt( as.numeric( t(bvar[j,]-bvar_cen) %*% vcov1%*% (bvar[j,]-bvar_cen) ) ) }
for(j in 1:l){ logrr_se2[1,j] <- sqrt( as.numeric( t(bvar[j,]-bvar_cen2) %*% vcov2%*% (bvar[j,]-bvar_cen2) ) ) }

##############################################################################
# SUMMARY AND COMPARISON OF KEY FEATURES

# PLOT SEASONAL CURVE
# CREATE DAY OF YEAR FOR X AXIS
#lndn$data<-as.Date(lndn$date,"%Y-/%m-/%d")
DY <-as.Date( lndn$date[lndn$year == 2002])
dateY <- format(DY, format="%b %d")
M <- mondate("1-1-2000")
FDM <- as.Date( rev(M - 1:12) )
FDM <- c( FDM, as.Date( DY[length(DY)] ) )
dateM <- format(FDM, format="%b %d")
monIND <- which( as.character(dateY) %in% as.character(dateM) )

plot(Trange, exp(logrr1), type='l', xlab="Day of year", ylab="RR (95% CI)",  ylim=c(0.9,1.4), lwd=2, xaxt='n',col="black" )
axis(1, at=Trange[monIND], labels=FALSE)
text(x=Trange[monIND], y=par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]), labels=dateM, srt=45, adj=1, xpd=TRUE, cex=0.8)
polygon(c(Trange, rev(Trange)), c(exp(logrr1-1.96*logrr_se1), rev(exp(logrr1+1.96*logrr_se1))),
        col = adjustcolor("black", alpha.f = 0.10) ,border = NA)
# 95%CI FOR RR
lines(Trange, exp(logrr1-1.96*logrr_se1),lty=2,col="black" )
lines(Trange, exp(logrr1+1.96*logrr_se1),lty=2,col="black" )

lines(Trange, exp(logrr2),lwd=2,col="red")
polygon(c(Trange, rev(Trange)), c(exp(logrr2-1.96*logrr_se2), rev(exp(logrr2+1.96*logrr_se2))),
        col = adjustcolor("red", alpha.f = 0.10) ,border = NA)
# 95%CI FOR RR
lines(Trange, exp(logrr2-1.96*logrr_se2),lty=2,col="red")
lines(Trange, exp(logrr2+1.96*logrr_se2),lty=2,col="red")

abline(h=1, col="black")

# Q2: What can you learn from the figure?

# TROUGH (95%eCI)
source("sfindmin.R")

trough1_est <- findmin(spline.season,fit1)
trough1<-findmin(spline.season,fit1,sim=T)

trough2_est <- findmin(spline.season,fit2)
trough2<- findmin(spline.season,fit2,sim=T)

trough<-rbind(trough1,trough2)
rownames(trough)<-c("Unadjusted","Adjusted")

# PEAK (95%eCI)
source("findmax.R")

peak1_est <- findmax(spline.season,fit1)
peak1<-findmax(spline.season,fit1,sim=T)

peak2_est <- findmax(spline.season,fit2)
peak2<-findmax(spline.season,fit2,sim=T)

peak<-rbind(peak1,peak2)
rownames(peak)<-c("Unadjusted","Adjusted")

# TIMINGS
timings<-cbind(peak,trough)
print(timings)

# Q3: What is your interpretation of the results above?

# PEAK-TO-TROUGH (95%CI)
max1<-apply(logrr1, 1, which.max ) # identify the peak
logrr1_max<-apply(logrr1, 1, max,na.rm=TRUE) # identify the coefficient at the peak
se1<- logrr_se1[cbind(seq_along(max1), max1)] # identify the se at the peak
ptr1 <-cbind(ptr=exp(logrr1_max), ptr.low=exp(logrr1_max-1.96*se1),ptr.high=exp( logrr1_max+1.96*se1))

max2<-apply(logrr2, 1, which.max )
logrr2_max <-apply(logrr2, 1, max,na.rm=TRUE)
se2 <- logrr_se2 [cbind(seq_along(max2), max2)]
ptr2 <-cbind(ptr=exp(logrr2_max), ptr.low=exp( logrr2_max-1.96*se2),ptr.high=exp( logrr2_max+1.96*se2))

ptr<-rbind(Unadjusted=ptr1,Adjusted=ptr2)
print(ptr)

# DIFFERENCE BETWEEN PTR BEFORE AND AFTER ADJUSTMENT
# ABSOLUTE DIFFERENCE
delta<-logrr2_max-logrr1_max
delta_se<-sqrt(se1^2+se2^2)
ad_ptr<-c(change=delta,change.low=delta-1.96*delta_se,change.high=delta+1.96*delta_se)
print(ad_ptr)
#RELATIVE DIFFERENCE
rd_ptr<-exp(ad_ptr)
print(rd_ptr)

# Q4: What is your interpretation of the results above?

# ATTRIBUTABLE FRACTION
source("attrs.R")
af1_est <- attrs(lndn$seasonal,spline.season,lndn$all,lndn,fit1,type="af",tot=T)
simaf1 <- attrs (lndn$seasonal,spline.season,lndn$all,lndn,fit1,type="af",tot=T,sim = T)
af1.low <- quantile(simaf1,c(2.5)/100)
af1.up <- quantile(simaf1,c(97.5)/100)
af1<-c(af=af1_est,af1.low,af1.up)

af2_est <- attrs (lndn$seasonal,spline.season,lndn$all,lndn,fit2,type="af",tot=T)
simaf2 <- attrs (lndn$seasonal,spline.season,lndn$all,lndn,fit2,type="af",tot=T,sim = T)
af2.low <- quantile(simaf2,c(2.5)/100)
af2.up <- quantile(simaf2,c(97.5)/100)
af2<-c(af=af2_est,af2.low,af2.up)

af<-rbind(Unadjusted=af1,Adjusted=af2)
colnames(af)[2:3]<-c("af.low","af.up")

# Q5: What is your interpretation of the results above?

# DIFFERENCE BETWEEN AF BEFORE AND AFTER ADJUSTMENT
delta<-af2_est-af1_est
delta_se<-sqrt(((af1.up-af1.low)/1.96)^2+((af2.up-af2.low)/1.96)^2)
delta_af<-c(change=delta,change.low=delta-1.96*delta_se,change.high=delta+1.96*delta_se)

#SUMMARY AND COMPARISON OF KEY FEATURES (TABLE 2)
table<-rbind(cbind(timings,ptr,af),Changes=c(peak2_est-peak1_est,NA,NA,trough2_est-trough1_est,NA,NA,rd_ptr,delta_af))
table<-round(table,digits = 3)
print(table)

# Q6: What is your conclusion from all the results above?

