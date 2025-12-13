

### Simulate the acs type of retrospective study
### fit likelihood and check the estimators
### include 3-state and 5-state models
### test the robustness of models
### June, 2025


logl_acs <- function(theta,sdat){
  ss <- theta[1]  ### sensitivity
  bb <- theta[2]  ### hazard of exponential
  pp <-  ss*exp(-(sdat$dtime*bb))
  logl <-  sum(sdat$weight*ifelse(sdat$tpos==1,-sdat$tpos*log(pp), -(1-sdat$tpos)*log(1-pp)))
  return(logl)
}


####
#### simulate ACS type of study 
#### use three years data after blood draw
#### consider 100% detection among shedding cancers
#### 0% among non-shedding cancers
####


nsimulation <- 1000
npop <- 50000
out <- matrix(0,nsimulation,8)

for (i in 1:nsimulation){
  if (i%%10==0) cat(i,"..")
  JJ <- rexp(n=npop,rate=1)/0.01     ## baseline preclinical incidence
  tt <- rbinom(npop,size=1,prob=0.5) ## type 1 cancer shedding;  type II cancer not shedding
  ss1 <- rexp(n=npop,rate=1)         ## mean sojourn time about 1 years for type I cancer
  ss2 <- rexp(n=npop,rate=1)         ## mean sojourn time about 1 years for type II cancer
  
  cc <- ifelse(tt==1,JJ + ss1,JJ+ss2)
  
  JJ <- JJ[cc>=50]
  tt <- tt[cc>=50]
  cc <- cc[cc>=50]
  
  
  dd <- rep(0,length(cc))
  aa <- 1*(JJ<50 & cc>50) ## preclinical cancer or not at screening
  ## type of cancer
  
  dd[aa==1 & tt==1] <- rbinom(sum(aa[tt==1]),1,1)  ## preclinical sensitivity for type I cancer 
  dd[aa==1 & tt==0] <- rbinom(sum(aa[tt==0]),1,0)  ## preclinical sensitivity for type II cancer
  
  ## take the cancers in the three years after blood draw
  
  tpos <- dd[cc>50 & cc<53]
  dtime <- cc[cc>50 & cc<53]-50
  weight <- 1
  sdat <- data.frame(cbind(dtime,tpos,weight))
  names(sdat) <- c("dtime","tpos","weight")
  fit<-optim(c(0.5,0.5),logl_acs,method = "L-BFGS-B",hessian=T,lower = c(0.001, 0.001), upper = c(1, Inf),sdat=sdat)
  
  out[i,1:2] <- fit$par
  out[i,3:4] <- diag(solve(fit$hessian))
  out[i,5] <-  sum(dd)/sum(aa)
  out[i,6] <- mean(cc>50 & cc<51)
  out[i,7] <- mean(cc>51 & cc<52)
  out[i,8] <- mean(cc>52 & cc<53)
}  
apply(out[out[,4]<10,],2,mean)
sqrt(apply(out[out[,4]<10,],2,var))
apply(out,2,mean)
sqrt(apply(out,2,var))





nsimulation <- 1000
npop <- 50000
out <- matrix(0,nsimulation,8)

for (i in 1:nsimulation){
  if (i%%10==0) cat(i,"..")
  JJ <- rexp(n=npop,rate=1)/0.01     ## baseline preclinical incidence
  tt <- rbinom(npop,size=1,prob=0.5) ## type 1 cancer shedding;  type II cancer not shedding
  ss1 <- rexp(n=npop,rate=1)         ## mean sojourn time about 1 years for type I cancer
  ss2 <- rexp(n=npop,rate=0.5)       ## mean sojourn time about 2 years for type II cancer
  
  cc <- ifelse(tt==1,JJ + ss1,JJ+ss2)
  
  JJ <- JJ[cc>=50]
  tt <- tt[cc>=50]
  cc <- cc[cc>=50]
  
  
  dd <- rep(0,length(cc))
  aa <- 1*(JJ<50 & cc>50) ## preclinical cancer or not at screening
  ## type of cancer
  
  dd[aa==1 & tt==1] <- rbinom(sum(aa[tt==1]),1,1)  ## preclinical sensitivity for type I cancer 
  dd[aa==1 & tt==0] <- rbinom(sum(aa[tt==0]),1,0)  ## preclinical sensitivity for type II cancer
  
  ## take the cancers in the three years after blood draw
  
  tpos <- dd[cc>50 & cc<53]
  dtime <- cc[cc>50 & cc<53]-50
  weight <- 1
  sdat <- data.frame(cbind(dtime,tpos,weight))
  names(sdat) <- c("dtime","tpos","weight")
  fit<-optim(c(0.5,0.5),logl_acs,method = "L-BFGS-B",hessian=T,lower = c(0.001, 0.001), upper = c(1, Inf),sdat=sdat)
  
  out[i,1:2] <- fit$par
  out[i,3:4] <- diag(solve(fit$hessian))
  out[i,5] <-  sum(dd)/sum(aa)
  out[i,6] <- mean(cc>50 & cc<51)
  out[i,7] <- mean(cc>51 & cc<52)
  out[i,8] <- mean(cc>52 & cc<53)
}  
apply(out[out[,4]<10,],2,mean)





####
####
####
#### Simulate the five-state model in a retrospective study
#### 
####
####

rm(list=ls())
setwd("/Users/jdai/Documents/Sensitivity_estimation/")


sojourn_early_late_loglik <- function(theta,edat,ldat){
  
  ### four parameters of interest in the likelihood 
  ### edat: the data for early-stage cancers
  ### ldat: the data for late-stage cancers
  
  ss0 <- theta[1]  ### early-stage sensitivity
  bb0 <- theta[2]  ### early-stage hazard, sum of two hazards
  ss1 <- theta[3]  ### late-stage sensitivity
  bb1 <- theta[4]  ### late-stage hazard
  
  pp0 <-  ss0*exp(-(edat$dtime*bb0))
  logl <-  sum(ifelse(edat$tpos==1,-edat$tpos*log(pp0), -(1-edat$tpos)*log(1-pp0)))
  
  pp1 <-  ss0*(bb1/(bb0-bb1))*(exp(-(ldat$dtime*bb1))-exp(-ldat$dtime*bb0))
  pp2 <-  ss1*exp(-bb1*ldat$dtime)
  logl <-  logl + sum(ifelse(ldat$tpos==1,-ldat$tpos*log(pp1+pp2), -(1-ldat$tpos)*log(1-pp1-pp2)))
  return(logl)
}


nsimulation <- 1000
out <- matrix(0,nsimulation,4)
npop <- 50000
for (i in 1:nsimulation){
  cat(i,"..")
  JJ <- rexp(n=npop,rate=1)/0.01       ## baseline preclinical incidence to born into early stage
  ss1 <- rexp(n=npop,rate=1)           ## mean sojourn time about 1 years for early clinical cancer
  ss2 <- rexp(n=npop,rate=1)           ## mean sojourn time about 1 years for preclinical early going to late stage
  ss3 <- rexp(n=npop,rate=2)           ## mean sojourn time about 1/2 years for preclinical late going to clinical
  tt <- ifelse(ss1<ss2,1,0)              ## indicator whether early stage going into clinical or late stage
  ss <- ifelse(tt==1,ss1,ss2)
  
  cc <- ifelse(tt==1,JJ + ss1,JJ+ss2+ss3) # the observed cancer diagnosed time 
  
  ##
  ## restrict the cancer diagnosis between 50 and 53, blood test at age 50 and 3 year follow-up
  ##
  
  ss2 <- ss2[cc>50 & cc<53]
  JJ <- JJ[cc>=50 & cc<53]
  tt <- tt[cc>=50 & cc<53]
  cc <- cc[cc>=50 & cc<53]
  
  
  dd <- rep(0,length(cc))
  aa <- 1*(JJ<50 & cc>50)                    ## preclinical cancer or not at screening
  tt2 <- rep(0,length(tt))                   ## tt2 is the variable for cancer stage at time of bd, 0: no cancer; 1: early stage; 2: late stage
  tt2[tt==1 & aa==1] <- 1                    ## preclinical early-stage went on to be clinical early-stage 
  tt2[tt==0 & aa==1 & (JJ+ss2) >50 ] <- 1    ## preclinical early-stage went on to be clinical late-stage, stay on early-stage at blood draw
  tt2[tt==0 & aa==1 & (JJ+ss2) <50 ] <- 2    ## preclinical early-stage went on to be clinical late-stage, stay on late-stage at blood draw
  
  
  dd[aa==1 & tt2==1] <- rbinom(sum(aa[tt2==1]),1,0.3)    ## preclinical sensitivity for early stage cancer 
  dd[aa==1 & tt2==2] <- rbinom(sum(aa[tt2==2]),1,0.8)    ## preclinical sensitivity for late stage cancer
  
  tpos1 <- dd[tt==1]
  dtime1 <- cc[tt==1]-50
  
  tpos2 <- dd[tt==0]
  dtime2 <- cc[tt==0]-50
  
  cdat <- data.frame(cbind(tpos1,dtime1))
  names(cdat) <- c("tpos","dtime")
  
  sdat <- data.frame(cbind(tpos2,dtime2))
  names(sdat) <- c("tpos","dtime")
  
  fit<-optim(c(0.5,1,0.5,2),sojourn_early_late_loglik,method = "Nelder-Mead",edat=cdat,ldat=sdat)
  #fit<-optim(c(0.3,1,0.8,2),sojourn_early_late_loglik,method = "L-BFGS-B",hessian=F,lower = c(0.001,0.001,0.001,0.001), upper = c(1,10,1,10),edat=cdat,ldat=sdat)
  out[i,1] <- fit$par[1]
  out[i,2] <- fit$par[2]
  out[i,3] <- fit$par[3]
  out[i,4] <- fit$par[4]
}  
apply(out,2,mean)
#[1] 0.3147749 2.0828486 0.8049659 2.0418849
sqrt(apply(out,2,var))
#[1] 0.08433805 0.45886052 0.09459195 0.35786990


#####
##### now simulate a subgroup (type II) of low-shedding cancers, slow-growing, unable to detect at early stage, able to detect 70% at late stage
#####
#####


nsimulation <- 1000
out <- matrix(0,nsimulation,4)
npop <- 50000
for (i in 1:nsimulation){
  cat(i,"..")
  JJ <- rexp(n=npop,rate=1)/0.01                 ## baseline preclinical incidence to born into early stage
  type_cancer <- rbinom(n=npop,size=1,prob=0.3)  ## two types of cancers, type I 30%
  
  n_typeI  <- sum(type_cancer)
  n_typeII <- sum(1-type_cancer)
  
  ss1 <- rep(0,npop)
  ss2 <- rep(0,npop)
  ss3 <- rep(0,npop)
  
  ss1[type_cancer==1] <- rexp(n=n_typeI,rate=1)           ## mean sojourn time about 1 years for early clinical cancer
  ss2[type_cancer==1] <- rexp(n=n_typeI,rate=1)           ## mean sojourn time about 1 years for preclinical early going to late stage
  ss3[type_cancer==1] <- rexp(n=n_typeI,rate=2)           ## mean sojourn time about 1/2 years for preclinical late going to clinical
  
  ss1[type_cancer==0] <- rexp(n=n_typeII,rate=0.5)           ## mean sojourn time about 2 years for early clinical cancer
  ss2[type_cancer==0] <- rexp(n=n_typeII,rate=0.5)           ## mean sojourn time about 2 years for preclinical early going to late stage
  ss3[type_cancer==0] <- rexp(n=n_typeII,rate=2)           ## mean sojourn time about 0.66 years for preclinical late going to clinical
  
  
  tt <- ifelse(ss1<ss2,1,0)              ## indicator whether early stage going into clinical or late stage
  ss <- ifelse(tt==1,ss1,ss2)
  
  cc <- ifelse(tt==1,JJ + ss1,JJ+ss2+ss3) # the observed cancer diagnosed time 
  
  ##
  ## restrict the cancer diagnosis between 50 and 53, blood test at age 50 and 3 year follow-up
  ## the subset of cancers in 3 years as in ACS
  ##
  
  type_cancer_x <- type_cancer[cc>50 & cc<53]
  ss_x <- ss2[cc>50 & cc<53]
  JJ_x <- JJ[cc>=50 & cc<53]
  tt_x <- tt[cc>=50 & cc<53]
  cc_x <- cc[cc>=50 & cc<53]
  
  
  dd_x <- rep(0,length(cc_x))
  aa_x <- 1*(JJ_x<50 & cc_x>50)                    ## preclinical cancer or not at screening
  
  ### generate the test positive data at the time of blood draw, sensitivity_early=0.30, sensitivity_late=0.79
  ### tt2 denote whether the cancer is at early or late stage at the time of blood draw
  
  tt2 <- rep(0,length(tt_x))                        ## tt2 is the variable for cancer stage at time of bd, 0: no cancer; 1: early stage; 2: late stage
  tt2[tt_x==1 & aa_x==1] <- 1                       ## preclinical early-stage went on to be clinical early-stage 
  tt2[tt_x==0 & aa_x==1 & (JJ_x+ss_x) >50 ] <- 1    ## preclinical early-stage went on to be clinical late-stage, stay on early-stage at blood draw
  tt2[tt_x==0 & aa_x==1 & (JJ_x+ss_x) <50 ] <- 2    ## preclinical early-stage went on to be clinical late-stage, stay on late-stage at blood draw
  
  
  dd_x[aa_x==1 & tt_x==1 & type_cancer_x==1 ] <- 1    ## preclinical sensitivity for early stage cancer 
  dd_x[aa_x==1 & tt_x==0 & type_cancer_x==1 ] <- 1    ## preclinical sensitivity for late stage cancer
  dd_x[aa_x==1 & tt_x==0 & type_cancer_x==0 & tt2==2 ] <- rbinom(sum(aa_x==1 & tt_x==0 & type_cancer_x==0 & tt2==2),size=1,prob=0.7)    ## preclinical sensitivity for late stage cancer
  
  
  ### prepare the dataset to feed into likelihood
  
  tpos1 <- dd_x[tt_x==1]
  dtime1 <- cc_x[tt_x==1]-50
  
  tpos2 <- dd_x[tt_x==0]
  dtime2 <- cc_x[tt_x==0]-50
  
  cdat <- data.frame(cbind(tpos1,dtime1))
  names(cdat) <- c("tpos","dtime")
  
  sdat <- data.frame(cbind(tpos2,dtime2))
  names(sdat) <- c("tpos","dtime")
  
  fit<-optim(c(0.5,1,0.5,2),sojourn_early_late_loglik,method = "Nelder-Mead",edat=cdat,ldat=sdat)
  #fit<-optim(c(0.3,1,0.8,2),sojourn_early_late_loglik,method = "L-BFGS-B",hessian=F,lower = c(0.001,0.001,0.001,0.001), upper = c(1,10,1,10),edat=cdat,ldat=sdat)
  out[i,1] <- fit$par[1]
  out[i,2] <- fit$par[2]
  out[i,3] <- fit$par[3]
  out[i,4] <- fit$par[4]
}  
apply(out,2,mean)
#[1] 0.3066315 2.0322333 0.8019795 2.0147982
sqrt(apply(out,2,var))



