



############################################################################################### Case-control Simulation
############################################################################################### N=300 
###############################################################################################
###############################################################################################

library(clinfun)

gompertz <- function(Xbeta, k,maxrisk) maxrisk*exp(-exp(-k*Xbeta))

beta1list <- c(0, 0.50, 1, 1.5, 2) 

resultsmatrix <- NULL

for(kk in 1:length(beta1list)){
 ytotal <- NULL
mselist1 <- mselist2 <- mselist3 <- mselist4 <- NULL 
auclist1 <- auclist2 <- auclist3 <- auclist4 <- NULL 
aiclist1 <- aiclist2 <- aiclist3 <- aiclist4  <- NULL
 p1quad <- p2quad <- p3quad <- p4quad <- NULL 
 p1inter <-  p2inter <-  p3inter <-  p4inter <-  NULL 

 
for(ii in 1:5000){

Nbig <- 10000
X1 <- rbinom(Nbig,1,p=0.5 )
X2  <- runif(Nbig,0,30)
Xmat <- cbind(rep(1,Nbig),X1 ,X2 )


beta <- c(-2,beta1list[kk],0.15)

betamat <- matrix(beta,ncol=1,nrow=3)
Xbeta  <- Xmat %*% betamat

ptrue <- gompertz(Xbeta,1,0.5)
 
Ybig <- rbinom(Nbig, 1, ptrue )

if(sum(Ybig) <= 300 | sum(1-Ybig) <= 300 ) print("increase sample size")
if(sum(Ybig) <= 300 | sum(1-Ybig) <= 300 ) stop

ndids <- sample((1:Nbig)[Ybig==0], 150, replace=FALSE)
dids <- sample((1:Nbig)[Ybig==1], 150, replace=FALSE)

selecth <- c(ndids, dids)

ytotal <- c(ytotal,sum(Ybig[selecth]))

X22 <- X2^2

xx1 <- glm(Ybig[selecth] ~ X1[selecth]+X2[selecth], family=binomial)
xx2 <-  glm(Ybig[selecth] ~ X1[selecth]+X2[selecth]+X22[selecth] , family=binomial)
xx3 <-  glm(Ybig[selecth] ~ X1[selecth]*X2[selecth], family=binomial)
xx4 <-  glm(Ybig[selecth] ~ X1[selecth]*X2[selecth]+ X22[selecth], family=binomial)

prev <- mean(Ybig)

xb1 <- xx1$linear.predictors + log(prev/(1-prev))
xb2 <- xx2$linear.predictors + log(prev/(1-prev))
xb3 <- xx3$linear.predictors + log(prev/(1-prev))
xb4 <- xx4$linear.predictors + log(prev/(1-prev))

phat1 <- exp(xb1)/(1+exp(xb1))
phat2 <- exp(xb2)/(1+exp(xb2))
phat3 <- exp(xb3)/(1+exp(xb3))
phat4 <- exp(xb4)/(1+exp(xb4))

mselist1 <- c(mselist1 ,mean((phat1-ptrue[selecth])^2))
mselist2 <- c(mselist2, mean((phat2-ptrue[selecth])^2))
mselist3 <- c(mselist3, mean((phat3-ptrue[selecth])^2))
mselist4 <- c(mselist4, mean((phat4-ptrue[selecth])^2))

auclist1  <- c(auclist1,roc.area.test(xx1$fitted.values, Ybig[selecth])[1]$area)
auclist2  <- c(auclist2,roc.area.test(xx2$fitted.values, Ybig[selecth])[1]$area)
auclist3  <- c(auclist3,roc.area.test(xx3$fitted.values, Ybig[selecth])[1]$area)
auclist4  <- c(auclist4,roc.area.test(xx4$fitted.values, Ybig[selecth])[1]$area)

aiclist1  <- c(aiclist1,AIC(xx1))
aiclist2  <- c(aiclist2,AIC(xx2))
aiclist3  <- c(aiclist3,AIC(xx3))
aiclist4  <- c(aiclist4,AIC(xx4))

p1quad <- 	c(p1quad,NA)
p1inter <- 	c(p1inter,NA) 
p2quad <- 	c(p2quad,summary(xx2)$coeff[4,4])
p2inter <- 	c(p2inter,NA) 
p3quad <- 	c(p3quad,NA)
p3inter <- 	c(p3inter,summary(xx3)$coeff[4,4]) 
p4quad <- 	c(p4quad,summary(xx4)$coeff[4,4])
p4inter <- 	c(p4inter,summary(xx4)$coeff[5,4]) 

}
 
 
resultsmatrix <- rbind(resultsmatrix, 
c(
mean(ytotal),
mean(mselist1)/mean(mselist1),mean(mselist1)/mean(mselist2),mean(mselist1)/mean(mselist3),mean(mselist1)/mean(mselist4),
mean(auclist1),mean(auclist2),mean(auclist3),mean(auclist4),
mean(aiclist1),mean(aiclist2),mean(aiclist3),mean(aiclist4),
mean(p1quad < 0.05),mean(p2quad < 0.05),mean(p3quad < 0.05),mean(p4quad < 0.05),
mean(p1inter < 0.05),mean(p2inter < 0.05),mean(p3inter < 0.05),mean(p4inter < 0.05)
))


print(kk)

}

lab1 <- paste("casecontrol.sim1.N=300.",Nbig,".csv",sep="")
 
write.csv(resultsmatrix,lab1)
 
