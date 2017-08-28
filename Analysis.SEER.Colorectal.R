

rm(list=ls())

#######
#######
#######
####### read in data from SEER
####### 
#######
#######
####### 

data001 <- data001[ !is.na( data001$Age.at.diagnosis) & !is.na(data001$Race.recode..White..Black..Other.) & !is.na(data001$Sex) & !is.na(data001$Survival.months) ,] 
dim(data001)

data001 <- data001[data001$Race %in% c("White", "Black"),]
 
Age  <- data001$Age.at.diagnosis
Race <- data001$Race.recode..White..Black..Other.
Race[Race ==" "] <- NA
Race[Race ==""] <- NA
table(Race)


##### Outcome 
Survival.months <- data001$Survival.months
Dead <-  data001$Vital.status.recode..study.cutoff.used.=="Dead"
mean( ( Survival.months < 3) & !Dead)

DeathWithin3 <- (Survival.months < 3) & Dead 
DeathWithin3[ is.na(Survival.months) ] <- NA
mean(DeathWithin3, na.rm=TRUE)
 
Age10 <- Age/10 - 40/10 
Age102 <- Age10^2
 
Sex <- data001$Sex
DxYr  <- data001$Year.of.diagnosis 

######################################################
####################################################
##################  Estimates  Black Race:  Logistic 
######################################################
######################################################

####### Odds ratios, confidence intervals, and p-values
x1 <- glm(DeathWithin3 ~ Age10 + Sex   , family=binomial, subset=(Race=="Black")) 
x2 <- glm(DeathWithin3 ~ Age10 + Sex + Age102 , family=binomial, subset=(Race=="Black")) 
x3 <- glm(DeathWithin3 ~ Age10 + Sex + Age10*Sex, family=binomial, subset=(Race=="Black")) 
x4 <- glm(DeathWithin3 ~ Age10 + Sex + Age102 +Age10*Sex, family=binomial, subset=(Race=="Black"))   

rnd <- 2
row1 <- paste(
round(exp(x1$coefficients),2), " (",
round(exp(confint(x1)[,1]),2),"-",
round(exp(confint(x1)[,2]),2),")",sep="") 
g1 <- paste(row1[-1],round(summary(x1)$coefficients[-1,4],3), sep=" & ")


row2 <- paste(
round(exp(x2$coefficients),2), " (",
round(exp(confint(x2)[,1]),2),"-",
round(exp(confint(x2)[,2]),2),")",sep="")
row2 
summary(x2)
g2 <- paste(row2[-1],round(summary(x2)$coefficients[-1,4],3), sep=" & ")


row3 <- paste(
round(exp(x3$coefficients),2), " (",
round(exp(confint(x3)[,1]),2),"-",
round(exp(confint(x3)[,2]),2),")",sep="")
row3 
summary(x3)
g3 <- paste(row3[-1],round(summary(x3)$coefficients[-1,4],3), sep=" & ")


row4 <- paste(
round(exp(x4$coefficients),2), " (",
round(exp(confint(x4)[,1]),2),"-",
round(exp(confint(x4)[,2]),2),")",sep="")
row4 
summary(x4)
g4 <- paste(row4[-1],round(summary(x4)$coefficients[-1,4],3), sep=" & ")


### Additional Metrics 
### Brier Score
round(mean((x1$fitted.values-x1$y)^2),4)
round(mean((x2$fitted.values-x1$y)^2),4)
round(mean((x3$fitted.values-x1$y)^2),4)
round(mean((x4$fitted.values-x1$y)^2),4)

library(clinfun)

### AUC 
round(roc.area.test( x1$fitted.values, x1$y)[1]$area,4)
round(roc.area.test( x2$fitted.values, x2$y)[1]$area,4)
round(roc.area.test( x3$fitted.values, x3$y)[1]$area,4)
round(roc.area.test( x4$fitted.values, x4$y)[1]$area,4)

### AIC
AIC(x1) 
AIC(x2) 
AIC(x3) 
AIC(x4) 

######################################################
####################################################
##################  Estimates  Black Race:  Complementary Log-Log 
######################################################
######################################################

####### Odds ratios, confidence intervals, and p-values
x1 <- glm(DeathWithin3 ~ Age10 + Sex   ,  subset=(Race=="Black"), family=binomial(link="cloglog") ) 

rnd <- 2

row1 <- paste(
round(exp(x1$coefficients),2), " (",
round(exp(confint(x1)[,1]),2),"-",
round(exp(confint(x1)[,2]),2),")",sep="")
  
g1 <- paste(row1[-1],round(summary(x1)$coefficients[-1,4],3), sep=" & ")

### Additional Metrics 
### Brier Score
round(mean((x1$fitted.values-x1$y)^2),4)

### AUC 
round(roc.area.test( x1$fitted.values, x1$y)[1]$area,4)

### AIC
AIC(x1) 


######################################################
####################################################
##################  Estimates  Black Race:  G-J Link
######################################################
######################################################


allmatrix <- model.matrix(~ Age10 +Sex ) 

allmatrixBlack <- allmatrix[Race=="Black", ] 

allmatrixUniqueBlack <- unique(allmatrixBlack)

DeathAmongBlack <- DeathWithin3[Race=="Black"]

n.val.Black  <- m.val.Black <- NULL
for(ii in 1:dim(allmatrixUniqueBlack )[1]){
	allunique <- matrix(rep(allmatrixUniqueBlack [ii, ], dim(allmatrixBlack )[1]), byrow=TRUE, nrow= dim(allmatrixBlack )[1])
	subset <- apply( allunique  == allmatrixBlack, 1, all)
	n.val.Black  <- c(n.val.Black  , sum(DeathAmongBlack [subset]))
	m.val.Black  <- c(m.val.Black , sum(!DeathAmongBlack [subset]))
	}

 
my.lambda <- seq(-0.5,0.5,by=0.01)

loglikvect <- NULL 

for(ii in 1:50){
	
	my.result =  get.exact.beta.est(n.val.Black, m.val.Black,allmatrixUniqueBlack, lambda=my.lambda[ii]) 
 
	if(any(is.na(my.result$beta.sd))) next
	loglikvect <- rbind(loglikvect ,c(my.lambda[ii], my.result$log.likeli))
	print(ii)
	}
 
plot(loglikvect[,1],loglikvect[,2])
	
selectedlambda <- loglikvect[which.max(loglikvect[,2]),1] 
selectedmod <- get.exact.beta.est(n.val.Black, m.val.Black,allmatrixUniqueBlack, lambda=selectedlambda  ) 
 
xb <- allmatrixBlack  %*%  selectedmod$beta.est

outcome <- DeathAmongBlack 
 
phat <-  ( 1+ (1 + selectedlambda*xb)^(-1/selectedlambda)   )^(-1)
 
# Brier Score 
round(mean(  (outcome - phat )^2 ),4)


## AUC
library(clinfun)
roc.area.test(phat , outcome )
 
## AIC 
AIC  <- 2*(length(selectedmod$beta.est)+1)- 2* selectedmod$log.likeli
round(AIC  )


######################################################
######################################################
##################  Estimates  White Race:  Logistic 
######################################################
######################################################

####### Odds ratios, confidence intervals, and p-values
x1 <- glm(DeathWithin3 ~ Age10 + Sex   , family=binomial, subset=(Race=="White")) 
x2 <- glm(DeathWithin3 ~ Age10 + Sex + Age102 , family=binomial, subset=(Race=="White")) 
x3 <- glm(DeathWithin3 ~ Age10 + Sex + Age10*Sex, family=binomial, subset=(Race=="White")) 
x4 <- glm(DeathWithin3 ~ Age10 + Sex + Age102 +Age10*Sex, family=binomial, subset=(Race=="White"))   

rnd <- 2
row1 <- paste(
round(exp(x1$coefficients),2), " (",
round(exp(confint(x1)[,1]),2),"-",
round(exp(confint(x1)[,2]),2),")",sep="") 
g1 <- paste(row1[-1],round(summary(x1)$coefficients[-1,4],3), sep=" & ")


row2 <- paste(
round(exp(x2$coefficients),2), " (",
round(exp(confint(x2)[,1]),2),"-",
round(exp(confint(x2)[,2]),2),")",sep="")
row2 
summary(x2)
g2 <- paste(row2[-1],round(summary(x2)$coefficients[-1,4],3), sep=" & ")


row3 <- paste(
round(exp(x3$coefficients),2), " (",
round(exp(confint(x3)[,1]),2),"-",
round(exp(confint(x3)[,2]),2),")",sep="")
row3 
summary(x3)
g3 <- paste(row3[-1],round(summary(x3)$coefficients[-1,4],3), sep=" & ")


row4 <- paste(
round(exp(x4$coefficients),2), " (",
round(exp(confint(x4)[,1]),2),"-",
round(exp(confint(x4)[,2]),2),")",sep="")
row4 
summary(x4)
g4 <- paste(row4[-1],round(summary(x4)$coefficients[-1,4],3), sep=" & ")


### Additional Metrics 
### Brier Score
round(mean((x1$fitted.values-x1$y)^2),4)
round(mean((x2$fitted.values-x1$y)^2),4)
round(mean((x3$fitted.values-x1$y)^2),4)
round(mean((x4$fitted.values-x1$y)^2),4)

library(clinfun)

### AUC 
round(roc.area.test( x1$fitted.values, x1$y)[1]$area,4)
round(roc.area.test( x2$fitted.values, x2$y)[1]$area,4)
round(roc.area.test( x3$fitted.values, x3$y)[1]$area,4)
round(roc.area.test( x4$fitted.values, x4$y)[1]$area,4)

### AIC
AIC(x1) 
AIC(x2) 
AIC(x3) 
AIC(x4) 

######################################################
#####################################################
##################  Estimates  White Race:  Complementary Log-Log 
######################################################
######################################################

####### Odds ratios, confidence intervals, and p-values
x1 <- glm(DeathWithin3 ~ Age10 + Sex   ,  subset=(Race=="White"), family=binomial(link="cloglog") ) 

rnd <- 2

row1 <- paste(
round(exp(x1$coefficients),2), " (",
round(exp(confint(x1)[,1]),2),"-",
round(exp(confint(x1)[,2]),2),")",sep="")
  
g1 <- paste(row1[-1],round(summary(x1)$coefficients[-1,4],3), sep=" & ")

### Additional Metrics 
### Brier Score
round(mean((x1$fitted.values-x1$y)^2),4)

### AUC 
round(roc.area.test( x1$fitted.values, x1$y)[1]$area,4)

### AIC
AIC(x1) 



######################################################
#####################################################
##################  Estimates  White Race:  G-J link 
######################################################
######################################################

 
allmatrix <- model.matrix(~ Age10 +Sex, contrasts=list(Sex="contr.helmert")) 
allmatrix <- model.matrix(~ Age10 +Sex ) 

allmatrixWhite  <-  allmatrix[Race=="White", ]

allmatrixUniqueWhite <- unique(allmatrixWhite )

DeathAmongWhite <- DeathWithin3[Race=="White"]

n.val.White  <- m.val.White <- NULL
for(ii in 1:dim(allmatrixUniqueWhite )[1]){
	allunique <- matrix(rep(allmatrixUniqueWhite [ii, ], dim(allmatrixWhite)[1]), byrow=TRUE, nrow= dim(allmatrixWhite)[1])
	subset <- apply( allunique  == allmatrixWhite, 1, all)
	n.val.White <- c(n.val.White, sum(DeathAmongWhite [subset]))
	m.val.White <- c(m.val.White , sum(!DeathAmongWhite [subset]))
	}

 
my.lambda <- seq(-1,1,by=0.01)

loglikvect <- NULL 

for(ii in 10:125){
	
	my.result =  get.exact.beta.est(n.val.White, m.val.White,allmatrixUniqueWhite, lambda=my.lambda[ii]) 
 
	if(any(is.na(my.result$beta.sd))) next
	loglikvect <- rbind(loglikvect ,c(my.lambda[ii], my.result$log.likeli))
	print(ii)
	}
 
 plot(loglikvect[,1],loglikvect[,2])
	
selectedlambda <- loglikvect[which.max(loglikvect[,2]),1] 
selectedmod <- get.exact.beta.est(n.val.White, m.val.White,allmatrixUniqueWhite,  lambda=selectedlambda  ) 
 
xb <- allmatrixWhite  %*%  selectedmod$beta.est

outcome <- DeathAmongWhite 
 
phat <-  ( 1+ (1 + selectedlambda*xb)^(-1/selectedlambda)   )^(-1)

# Brier Score  
round(mean((outcome - phat )^2 ),4)
#    0.1691

library(clinfun)
 roc.area.test(phat , outcome ) 
#  0.6892974

### extra parameter estimated: lambda
AIC  <- 2*(length(selectedmod$beta.est)+1)- 2* selectedmod$log.likeli
round(AIC)  
 


#################################################################
#################################################################
############## Figure 2, and Appendix Figure A2a
#################################################################
#################################################################
 
combineddata <- data.frame(DeathWithin3=1*DeathWithin3, Age10,Sex,Age102,Race)

x1W <- glm(DeathWithin3 ~ Age10 + Sex   , family=binomial,data= combineddata , subset=(Race=="White")) 
x2W <- glm(DeathWithin3 ~ Age10 + Sex + Age102 , family=binomial, data= combineddata , subset=(Race=="White")) 
x3W <- glm(DeathWithin3 ~ Age10 + Sex + Age10*Sex, family=binomial, data= combineddata , subset=(Race=="White")) 
x4W <- glm(DeathWithin3 ~ Age10 + Sex + Age102 +Age10*Sex, family=binomial, data= combineddata , subset=(Race=="White"))   

agenew <- seq(min(Age10), max(Age10), length.out = 10000)
sexMalenew <- rep("Male", length(agenew))
sexFemalenew <- rep("Female", length(agenew))

racenew <- rep("White", length(agenew))
outcomenew <- rep(1,  length(agenew))

datapictureWMale <- data.frame(DeathWithin3=outcomenew , Age10=agenew ,Sex=sexMalenew ,Age102=(agenew^2),Race=racenew )
datapictureWFemale <- data.frame(DeathWithin3=outcomenew , Age10=agenew ,Sex=sexFemalenew ,Age102=(agenew^2),Race=racenew )


x1B <- glm(DeathWithin3 ~ Age10 + Sex   , family=binomial,data= combineddata , subset=(Race=="Black")) 
x2B <- glm(DeathWithin3 ~ Age10 + Sex + Age102 , family=binomial, data= combineddata , subset=(Race=="Black")) 
x3B <- glm(DeathWithin3 ~ Age10 + Sex + Age10*Sex, family=binomial, data= combineddata , subset=(Race=="Black")) 
x4B <- glm(DeathWithin3 ~ Age10 + Sex + Age102 +Age10*Sex, family=binomial, data= combineddata , subset=(Race=="Black"))   


agenew <- seq(min(Age10), max(Age10), length.out = 10000)
sexMalenew <- rep("Male", length(agenew))
sexFemalenew <- rep("Female", length(agenew))

racenew <- rep("Black", length(agenew))
outcomenew <- rep(1,  length(agenew))

datapictureBMale <- data.frame(DeathWithin3=outcomenew , Age10=agenew ,Sex=sexMalenew ,Age102=(agenew^2),Race=racenew )
datapictureBFemale <- data.frame(DeathWithin3=outcomenew , Age10=agenew ,Sex=sexFemalenew ,Age102=(agenew^2),Race=racenew )


predx1WM <- predict(x1W, newdata=datapictureWMale ,type="response")
predx2WM <- predict(x2W, newdata=datapictureWMale ,type="response")
predx3WM <- predict(x3W, newdata=datapictureWMale ,type="response")
predx4WM <- predict(x4W, newdata=datapictureWMale ,type="response")

predx1WF <- predict(x1W, newdata=datapictureWFemale ,type="response")
predx2WF <- predict(x2W, newdata=datapictureWFemale ,type="response")
predx3WF <- predict(x3W, newdata=datapictureWFemale ,type="response")
predx4WF <- predict(x4W, newdata=datapictureWFemale ,type="response")

predx1BM <- predict(x1B, newdata=datapictureBMale ,type="response")
predx2BM <- predict(x2B, newdata=datapictureBMale ,type="response")
predx3BM <- predict(x3B, newdata=datapictureBMale ,type="response")
predx4BM <- predict(x4B, newdata=datapictureBMale ,type="response")

predx1BF <- predict(x1B, newdata=datapictureBFemale ,type="response")
predx2BF <- predict(x2B, newdata=datapictureBFemale ,type="response")
predx3BF <- predict(x3B, newdata=datapictureBFemale ,type="response")
predx4BF <- predict(x4B, newdata=datapictureBFemale ,type="response")



meanconfint <- function(RaceHH=Race, SexHH=Sex, AgeHH=Age, WindowHH=Window){
			sub1min <- Sex==SexHH & Race==RaceHH& Age10 >= AgeHH- WindowHH & Age10 <= AgeHH+ WindowHH
			mp1 <- mean(DeathWithin3[sub1min],na.rm=TRUE)
			mp1se <- sd(DeathWithin3[sub1min],na.rm=TRUE)/sqrt(length(DeathWithin3[sub1min & !is.na(DeathWithin3) ]))
			confintm1 <- c(mp1 -1.96*mp1se , mp1 +1.96*mp1se)
		as.vector((c(mp1,confintm1))) 
		}

 
 
timevect <- c(0.1, 0.5, 1.0,1.5,2,2.5,3,3.5,4,4.5,4.9)


pdf("ColorectalSEER.pdf", height=12, width=12)

par(mfrow=c(2,2))

plot( c(0,0)  , c(0,0),  xaxt="n",yaxt="n", ylab="",xlab="",typ="n",
ylim=c(0,1),xlim=c(0,5), lwd=2)
axis(1, c(0,1,2,3,4,5,6),c("0","10","20","30","40","50", "60"),cex.axis=1.5)
axis(2, c(0,0.25,0.50, 0.75, 1.0),cex.axis=1.5)
title(main=list("Colorectal Cancer \n White Males", cex=1.5), xlab=list("Age of Diagnosis from 40 yo", cex=1.5), ylab=list("Probability of Early Death", cex=1.5))

lines(agenew,predx1WM, lwd=2,col=2,lty=1,cex=2)
lines(agenew,predx2WM, lwd=2,col=3,lty=1,cex=2)
lines(agenew,predx3WM, lwd=2,col=4,lty=1,cex=2)
lines(agenew,predx4WM, lwd=2,col="salmon",lty=1,cex=2)

for(pointhh in timevect ){
points(pointhh , meanconfint(Race="White", Sex="Male", Age=pointhh , Window=0.1)[1] , col=1, cex=1,pch = 19)
lines(c(pointhh,pointhh) ,c(meanconfint(Race="White", Sex="Male", Age=pointhh , Window=0.1)[2:3]), col=1, cex=1)
}

legend(0,0.99,
c(
"Empirical",
"Model M1",
"Model M2",
"Model M3",
"Model M4"), 
col=c(1,2,3,4,"salmon"), lty=c(0,1,1,1,1),pch=c(19,NA,NA,NA, NA), lwd=2, bty="n", seg.len=3, cex=1.5)


plot( c(0,0)  , c(0,0),  xaxt="n",yaxt="n", ylab="",xlab="",typ="n",
ylim=c(0,1),xlim=c(0,5), lwd=2)
axis(1, c(0,1,2,3,4,5,6),c("0","10","20","30","40","50", "60"),cex.axis=1.5)
axis(2, c(0,0.25,0.50, 0.75, 1.0),cex.axis=1.5)
title(main=list("Colorectal Cancer \n White Females", cex=1.5), xlab=list("Age of Diagnosis from 40 yo", cex=1.5), ylab=list("Probability of Early Death", cex=1.5))

lines(agenew,predx1WF, lwd=2,col=2,lty=1,cex=2)
lines(agenew,predx2WF, lwd=2,col=3,lty=1,cex=2)
lines(agenew,predx3WF, lwd=2,col=4,lty=1,cex=2)
lines(agenew,predx4WF, lwd=2,col="salmon",lty=1,cex=2)

for(pointhh in timevect ){
points(pointhh , meanconfint(Race="White", Sex="Female", Age=pointhh , Window=0.1)[1] , col=1, cex=1,pch = 19)
lines(c(pointhh,pointhh) ,c(meanconfint(Race="White", Sex="Female", Age=pointhh , Window=0.1)[2:3]), col=1, cex=1)
}

legend(0,0.99,
c(
"Empirical",
"Model M1",
"Model M2",
"Model M3",
"Model M4"), 
col=c(1,2,3,4,"salmon"), lty=c(0,1,1,1,1),pch=c(19,NA,NA,NA, NA), lwd=2, bty="n", seg.len=3, cex=1.5)


plot( c(0,0)  , c(0,0),  xaxt="n",yaxt="n", ylab="",xlab="",typ="n",
ylim=c(0,1),xlim=c(0,5), lwd=2)
axis(1, c(0,1,2,3,4,5,6),c("0","10","20","30","40","50", "60"),cex.axis=1.5)
axis(2, c(0,0.25,0.50, 0.75, 1.0),cex.axis=1.5)
title(main=list("Colorectal Cancer \n Black Males", cex=1.5), xlab=list("Age of Diagnosis from 40 yo", cex=1.5), ylab=list("Probability of Early Death", cex=1.5))

lines(agenew,predx1BM, lwd=2,col=2,lty=1,cex=2)
lines(agenew,predx2BM, lwd=2,col=3,lty=1,cex=2)
lines(agenew,predx3BM, lwd=2,col=4,lty=1,cex=2)
lines(agenew,predx4BM, lwd=2,col="salmon",lty=1,cex=2)

for(pointhh in timevect ){
points(pointhh , meanconfint(Race="Black", Sex="Male", Age=pointhh , Window=0.1)[1] , col=1, cex=1,pch = 19)
lines(c(pointhh,pointhh) ,c(meanconfint(Race="Black", Sex="Male", Age=pointhh , Window=0.1)[2:3]), col=1, cex=1)
}

legend(0,0.99,
c(
"Empirical",
"Model M1",
"Model M2",
"Model M3",
"Model M4"), 
col=c(1,2,3,4,"salmon"), lty=c(0,1,1,1,1),pch=c(19,NA,NA,NA, NA), lwd=2, bty="n", seg.len=3, cex=1.5)


plot( c(0,0)  , c(0,0),  xaxt="n",yaxt="n", ylab="",xlab="",typ="n",
ylim=c(0,1),xlim=c(0,5), lwd=2)
axis(1, c(0,1,2,3,4,5,6),c("0","10","20","30","40","50", "60"),cex.axis=1.5)
axis(2, c(0,0.25,0.50, 0.75, 1.0),cex.axis=1.5)
title(main=list("Colorectal Cancer \n Black Females", cex=1.5), xlab=list("Age of Diagnosis from 40 yo", cex=1.5), ylab=list("Probability of Early Death", cex=1.5))

lines(agenew,predx1BF, lwd=2,col=2,lty=1,cex=2)
lines(agenew,predx2BF, lwd=2,col=3,lty=1,cex=2)
lines(agenew,predx3BF, lwd=2,col=4,lty=1,cex=2)
lines(agenew,predx4BF, lwd=2,col="salmon",lty=1,cex=2)

for(pointhh in timevect ){
points(pointhh , meanconfint(Race="Black", Sex="Female", Age=pointhh , Window=0.1)[1] , col=1, cex=1,pch = 19)
lines(c(pointhh,pointhh) ,c(meanconfint(Race="Black", Sex="Female", Age=pointhh , Window=0.1)[2:3]), col=1, cex=1)
}

legend(0,0.99,
c(
"Empirical",
"Model M1",
"Model M2",
"Model M3",
"Model M4"), 
col=c(1,2,3,4,"salmon"), lty=c(0,1,1,1,1),pch=c(19,NA,NA,NA, NA), lwd=2, bty="n", seg.len=3, cex=1.5)

dev.off()

 

subset1 <- Race=="White" & Sex == "Male"
subset2 <- Race=="White" & Sex == "Female"

outcomesyub1 <- 1*DeathWithin3[subset1]
agesub1 <- Age10[subset1]

outcomesyub2 <- 1*DeathWithin3[subset2]
agesub2 <- Age10[subset2]

pdf("ColorectalKMeans.pdf", height=6, width=12)

par(mfrow=c(1,2))

plot(agesub1 ,outcomesyub1 ,type="n",main="",xlab="",ylab="",xaxt="n")
axis(1, c(0:5),10*c(0:5),cex.axis=1.5) 
lines(ksmooth(agesub1 ,outcomesyub1,"normal" ,bandwidth = 2), lwd=2, col =1, cex=1.5)
lines(ksmooth(agesub2 ,outcomesyub2 ,"normal",bandwidth = 2 ), lwd=2, col =2, cex=1.5)
legend(0.1, 1, c("Males", "Females"), col=c(1,2),lty=1, bty="n", lwd=2, seg.len=3 , cex=1.5)
title(main=list("Colorectal White Cohort", cex=1.5), xlab=list("Age of Diagnosis from 40 yo", cex=1.5), ylab=list("Probability of Early Death", cex=1.5))


subset1 <- Race=="Black" & Sex == "Male"
subset2 <- Race=="Black" & Sex == "Female"

outcomesyub1 <- 1*DeathWithin3[subset1]
agesub1 <- Age10[subset1]

outcomesyub2 <- 1*DeathWithin3[subset2]
agesub2 <- Age10[subset2]

plot(agesub1 ,outcomesyub1 ,type="n",main="",xlab="",ylab="",xaxt="n")
axis(1, c(0:5),10*c(0:5),cex.axis=1.5)  
lines(ksmooth(agesub1 ,outcomesyub1,"normal" ,bandwidth =2),lwd=2,  col =1, cex=1.5)
lines(ksmooth(agesub2 ,outcomesyub2 ,"normal",bandwidth = 2),lwd=2,  col =2, cex=1.5)
legend(0.1, 1, c("Males", "Females"), col=c(1,2),lty=1, bty="n", lwd=2, seg.len=3 , cex=1.5)
title(main=list("Colorectal Black Cohort", cex=1.5), xlab=list("Age of Diagnosis from 40 yo", cex=1.5), ylab=list("Probability of Early Death", cex=1.5))

dev.off()
 
