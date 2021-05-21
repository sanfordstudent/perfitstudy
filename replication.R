# Description ###################################################################################
# This file produces a replication of the simulations and results found in Sinharay (2017), with 
# the addition of four weighted PFSs and the inclusion of two asympototically-corrected versions
# of two PFSs from the original study, lz and eci4z.
#################################################################################################
set.seed(102856)
library(MESS)
library(tidyverse)
library(ltm)
library(irtoys)
library(PerFit)
library(ROCR)
library(metafor)


## from Magis (2012): most of the code to calculate lz, eci4z, and other stats in this family. 
## Response probabilities, first and second derivatives under the 3PL model (equation 36)
## th: ability value
## it: matrix of item parameters: one row per item, three columns
## (discrimination, difficulty, pseudo-guessing)
Pi<-function(th,it){
  res1<-res2<-res3<-NULL
  for (i in 1:nrow(it)) {
    res1[i]<-
      it[i,3]+(1-it[i,3])*exp(it[i,1]*(th-it[i,2]))/(1+exp(it[i,1]*(th-it[
        i,2])))
    res2[i]<-(1-it[i,3])*it[i,1]*exp(it[i,1]*(th-it[i,2]))/(1+exp(it[i,1]*(th-it[
      i,2])))^2
    res3[i]<-(1-it[i,3])*it[i,1]^2*exp(it[i,1]*(th-it[i,2]))*(1-exp(it[i,1]*(th-it[
      i,2])))/(1+exp(it[i,1]*(th-it[i,2])))^3
  }
  RES<-list(Pi=res1,dPi=res2,d2Pi=res3)
  return(RES)
}
## Functions ri and r0 (equations 24, 27 and 33)
## method: « ML » for maximum likelihood, « BM » for Bayesian modal (or MAp), « WL » for weighted likelihood
## mu, sigma: prior mean and standard deviation parameters of the normal distribution
ri<-function(th,it) Pi(th,it)$dPi/(Pi(th,it)$Pi*(1-Pi(th,it)$Pi))
r0<-function(method="ML",th,it,mu=0,sigma=1){
  res<-switch(method,
              ML=0,
              BM=(mu-th)/sigma^2,
              WL=sum(ri(th,it)*Pi(th,it)$d2Pi)/(2*sum(ri(th,it)*Pi(th,it)$dPi)))
  return(res)}
## Ability estimation (constrained to range [-4; 4]) (equation 14)
# thetaEst<-function(x,it,method="ML",mu=0,sigma=1){
#   f<-function(th){
#     r0(method=method,th,it=it,mu=mu,sigma=sigma)+sum((x-Pi(th,it)$Pi)*ri(th,it))
#   }
#   if (f(-4)<0 & f(4)<0) {res <- -4} else{
#     if (f(-4)>0 & f(4)>0) {res <- 4} else res<-uniroot(f,c(-4,4))$root
#   }
#   return(res)}

# Xia and Zheng, 2018: four different asym-weighted PFSs for detecting spuriously high scores
# SHa(1/2)
# responses: matrix of 0/1 item responses, columns are items and rows are persons
# Ability: vector of theta estimates
# IP: item parameters. disc, diff, guess.
sha1.2<- function(responses, theta, it){
  ## Weight function for SHa(1/2)
  wi<-function(th,it) (1-2* Pi(th,it)$Pi)/sqrt(Pi(th,it)$Pi)
  ## Function Wn(theta) (equation 9)
  ## x: response pattern (same length as nrow(it), zeros and ones as entries)
  Wn<-function(x,th,it) sum((x-Pi(th,it)$Pi)*wi(th,it))
  ## Function sig2n (equation 11)
  sig2n<-function(th,it)
    sum(wi(th,it)^2*Pi(th,it)$Pi*(1-Pi(th,it)$Pi))/nrow(it)
  ## Function cn (equation 15)
  cn<-function(th,it) sum(Pi(th,it)$dPi*wi(th,it))/sum(Pi(th,it)$dPi*ri(th,it))
  ## Function wi ‘tilde’ (equation 16)
  wiTilde<-function(th,it) wi(th,it)-cn(th,it)*ri(th,it)
  ## Function tau2n (equation 18)
  tau2n<-function(th,it)
    sum(wiTilde(th,it)^2*Pi(th,it)$Pi*(1-Pi(th,it)$Pi))/nrow(it)
  
  ## Indexes lz and lz*
  ## snijders: logical argument: FALSE returns lz, TRUE returns lz*
  ## so for SH's, we want it TRUE
  ## added th argument for pre-estimated theta
  SHa1.2<-function(x,it,th,method="ML",mu=0,sigma=1,snijders=TRUE){
    # th<-thetaEst(x=x,it=it,method=method,mu=mu,sigma=sigma)
    if (snijders==TRUE){
      EWn<-cn(th,it)*r0(method=method,th=th,it=it,mu=mu,sigma=sigma)
      VWn<-nrow(it)*tau2n(th,it)
    } else{
      EWn<-0
      VWn<-nrow(it)*sig2n(th,it)
    }
    res<-(Wn(x,th,it)-EWn)/sqrt(VWn)
    return(res)
  }
  
  stats <- vector()
  for(row in 1:nrow(responses)){
    th <- theta[row]
    x <- as.numeric(responses[row,][1,]) 
    stat <- SHa1.2(x, it, th)
    stats <- append(stats, stat)
  }
  
  return(stats)
}

# SHa(1)
# responses: matrix of 0/1 item responses, columns are items and rows are persons
# Ability: vector of theta estimates
# IP: item parameters. disc, diff, guess.
sha1<- function(responses, theta, it){
  ## Weight function for SHa(1)
  wi<-function(th,it) (1-2* Pi(th,it)$Pi)/(Pi(th,it)$Pi)
  ## Function Wn(theta) (equation 9)
  ## x: response pattern (same length as nrow(it), zeros and ones as entries)
  Wn<-function(x,th,it) sum((x-Pi(th,it)$Pi)*wi(th,it))
  ## Function sig2n (equation 11)
  sig2n<-function(th,it)
    sum(wi(th,it)^2*Pi(th,it)$Pi*(1-Pi(th,it)$Pi))/nrow(it)
  ## Function cn (equation 15)
  cn<-function(th,it) sum(Pi(th,it)$dPi*wi(th,it))/sum(Pi(th,it)$dPi*ri(th,it))
  ## Function wi ‘tilde’ (equation 16)
  wiTilde<-function(th,it) wi(th,it)-cn(th,it)*ri(th,it)
  ## Function tau2n (equation 18)
  tau2n<-function(th,it)
    sum(wiTilde(th,it)^2*Pi(th,it)$Pi*(1-Pi(th,it)$Pi))/nrow(it)
  
  ## Indexes lz and lz*
  ## snijders: logical argument: FALSE returns lz, TRUE returns lz*
  ## so for SH's, we want it TRUE
  ## added th argument for pre-estimated theta
  SHa1<-function(x,it,th,method="ML",mu=0,sigma=1,snijders=TRUE){
    # th<-thetaEst(x=x,it=it,method=method,mu=mu,sigma=sigma)
    if (snijders==TRUE){
      EWn<-cn(th,it)*r0(method=method,th=th,it=it,mu=mu,sigma=sigma)
      VWn<-nrow(it)*tau2n(th,it)
    } else{
      EWn<-0
      VWn<-nrow(it)*sig2n(th,it)
    }
    res<-(Wn(x,th,it)-EWn)/sqrt(VWn)
    return(res)
  }
  
  stats <- vector()
  for(row in 1:nrow(responses)){
    th <- theta[row]
    x <- as.numeric(responses[row,][1,]) 
    stat <- SHa1(x, it, th)
    stats <- append(stats, stat)
  }
  
  return(stats)
}

# SHb*(2)
# responses: matrix of 0/1 item responses, columns are items and rows are persons
# Ability: vector of theta estimates
# IP: item parameters. disc, diff, guess.
shb2<- function(responses, theta, it){
  ## Weight function for SHb(2)
  wi<-function(th,it) (1-2*Pi(th,it)$Pi)/exp(-(2*Pi(th,it)$Pi)^2)
  ## Function Wn(theta) (equation 9)
  ## x: response pattern (same length as nrow(it), zeros and ones as entries)
  Wn<-function(x,th,it) sum((x-Pi(th,it)$Pi)*wi(th,it))
  ## Function sig2n (equation 11)
  sig2n<-function(th,it)
    sum(wi(th,it)^2*Pi(th,it)$Pi*(1-Pi(th,it)$Pi))/nrow(it)
  ## Function cn (equation 15)
  cn<-function(th,it) sum(Pi(th,it)$dPi*wi(th,it))/sum(Pi(th,it)$dPi*ri(th,it))
  ## Function wi ‘tilde’ (equation 16)
  wiTilde<-function(th,it) wi(th,it)-cn(th,it)*ri(th,it)
  ## Function tau2n (equation 18)
  tau2n<-function(th,it)
    sum(wiTilde(th,it)^2*Pi(th,it)$Pi*(1-Pi(th,it)$Pi))/nrow(it)
  
  ## Indexes lz and lz*
  ## snijders: logical argument: FALSE returns lz, TRUE returns lz*
  ## so for SH's, we want it TRUE
  ## added th argument for pre-estimated theta
  SHb2<-function(x,it,th,method="ML",mu=0,sigma=1,snijders=TRUE){
    # th<-thetaEst(x=x,it=it,method=method,mu=mu,sigma=sigma)
    if (snijders==TRUE){
      EWn<-cn(th,it)*r0(method=method,th=th,it=it,mu=mu,sigma=sigma)
      VWn<-nrow(it)*tau2n(th,it)
    } else{
      EWn<-0
      VWn<-nrow(it)*sig2n(th,it)
    }
    res<-(Wn(x,th,it)-EWn)/sqrt(VWn)
    return(res)
  }
  
  stats <- vector()
  for(row in 1:nrow(responses)){
    th <- theta[row]
    x <- as.numeric(responses[row,][1,]) 
    stat <- SHb2(x, it, th)
    stats <- append(stats, stat)
  }
  
  return(stats)
}

# SHb*(3)
# responses: matrix of 0/1 item responses, columns are items and rows are persons
# Ability: vector of theta estimates
# IP: item parameters. disc, diff, guess.
shb3<- function(responses, theta, it){
  ## Weight function for SHb(3)
  wi<-function(th,it) (1-2*Pi(th,it)$Pi)/exp(-(3*Pi(th,it)$Pi)^2)
  ## Function Wn(theta) (equation 9)
  ## x: response pattern (same length as nrow(it), zeros and ones as entries)
  Wn<-function(x,th,it) sum((x-Pi(th,it)$Pi)*wi(th,it))
  ## Function sig2n (equation 11)
  sig2n<-function(th,it)
    sum(wi(th,it)^2*Pi(th,it)$Pi*(1-Pi(th,it)$Pi))/nrow(it)
  ## Function cn (equation 15)
  cn<-function(th,it) sum(Pi(th,it)$dPi*wi(th,it))/sum(Pi(th,it)$dPi*ri(th,it))
  ## Function wi ‘tilde’ (equation 16)
  wiTilde<-function(th,it) wi(th,it)-cn(th,it)*ri(th,it)
  ## Function tau2n (equation 18)
  tau2n<-function(th,it)
    sum(wiTilde(th,it)^2*Pi(th,it)$Pi*(1-Pi(th,it)$Pi))/nrow(it)
  
  ## Indexes lz and lz*
  ## snijders: logical argument: FALSE returns lz, TRUE returns lz*
  ## so for SH's, we want it TRUE
  ## added th argument for pre-estimated theta
  SHb3<-function(x,it,th,method="ML",mu=0,sigma=1,snijders=TRUE){
    # th<-thetaEst(x=x,it=it,method=method,mu=mu,sigma=sigma)
    if (snijders==TRUE){
      EWn<-cn(th,it)*r0(method=method,th=th,it=it,mu=mu,sigma=sigma)
      VWn<-nrow(it)*tau2n(th,it)
    } else{
      EWn<-0
      VWn<-nrow(it)*sig2n(th,it)
    }
    res<-(Wn(x,th,it)-EWn)/sqrt(VWn)
    return(res)
  }
  
  stats <- vector()
  for(row in 1:nrow(responses)){
    th <- theta[row]
    x <- as.numeric(responses[row,][1,]) 
    stat <- SHb3(x, it, th)
    stats <- append(stats, stat)
  }
  
  return(stats)
}

# Sinharay, 2016: asymptotically standard normal correction of eci4z
# SHb*(3)
# responses: matrix of 0/1 item responses, columns are items and rows are persons
# Ability: vector of theta estimates
# IP: item parameters. disc, diff, guess.
eci4zstar<- function(responses, theta, it){
  Ttheta <- function(th, it){sum(Pi(th,it)$Pi)/nrow(it)}
  ## Weight function for ECI4z. from sinharay, this is
  ## -[Pj(theta_i-hat)-T(theta_i-hat)]
  wi<-function(th,it) {-(Pi(th,it)$Pi - Ttheta(th,it))}
  ## Function Wn(theta) (equation 9)
  ## x: response pattern (same length as nrow(it), zeros and ones as entries)
  Wn<-function(x,th,it) sum((x-Pi(th,it)$Pi)*wi(th,it))
  ## Function sig2n (equation 11)
  sig2n<-function(th,it)
    sum(wi(th,it)^2*Pi(th,it)$Pi*(1-Pi(th,it)$Pi))/nrow(it)
  ## Function cn (equation 15)
  cn<-function(th,it) sum(Pi(th,it)$dPi*wi(th,it))/sum(Pi(th,it)$dPi*ri(th,it))
  ## Function wi ‘tilde’ (equation 16)
  wiTilde<-function(th,it) wi(th,it)-cn(th,it)*ri(th,it)
  ## Function tau2n (equation 18)
  tau2n<-function(th,it)
    sum(wiTilde(th,it)^2*Pi(th,it)$Pi*(1-Pi(th,it)$Pi))/nrow(it)
  
  ## Indexes lz and lz*
  ## snijders: logical argument: FALSE returns lz, TRUE returns lz*
  ## so for SH's, we want it TRUE
  ## added th argument for pre-estimated theta
  ECI4zstar<-function(x,it,th,method="ML",mu=0,sigma=1,snijders=TRUE){
    # th<-thetaEst(x=x,it=it,method=method,mu=mu,sigma=sigma)
    if (snijders==TRUE){
      EWn<-cn(th,it)*r0(method=method,th=th,it=it,mu=mu,sigma=sigma)
      VWn<-nrow(it)*tau2n(th,it)
    } else{
      EWn<-0
      VWn<-nrow(it)*sig2n(th,it)
    }
    res<-(Wn(x,th,it)-EWn)/sqrt(VWn)
    return(res)
  }
  
  stats <- vector()
  for(row in 1:nrow(responses)){
    th <- theta[row]
    x <- as.numeric(responses[row,][1,]) 
    stat <- ECI4zstar(x, it, th)
    stats <- append(stats, stat)
  }
  
  return(stats)
}

# Sinharay, 2017: Function to compute an ROC Area and standard error
ROC = function(good,bad,x) { 
  Areas=rep(0,ncol(good))
  SE=Areas  
  F=matrix(0,length(x),ncol(good))
  H=F
  for (j in 1:ncol(good)){
    goodstat=good[,j]
    badstat=bad[,j]
    for (i in 1:length(x)){
      # here, x[i] is the cutoff. it starts at -5 and ends at 5. all stats are standardized, and stats
      # where large values indicate aberrance are flipped
      F[i,j]=length(goodstat[goodstat<x[i]])/length(goodstat)
      H[i,j]=length(badstat[badstat<x[i]])/length(badstat)
    }
    #Use the R function `auc' to compute the ROC area
    a=auc(F[,j],H[,j],type='spline')
    Areas[j]=ifelse(a>1,0.999,a)
    a=Areas[j]
    ng=nrow(good)
    nb=nrow(bad)
    q1=a/(2-a)
    q2=2*a*a/(1+a)
    SE[j]=sqrt(( a*(1-a)+(nb-1)*(q1-a*a)+(ng-1)*(q2-a*a) )/(ng*nb))
  }
  return(rbind(Areas,SE,F,H))
}

# get everything in place for the analysis
ItemResponses <- list()
# from Sinharay, 2017
Output.ch <- NULL
SEs.ch <- NULL
x=seq(-5, 5, length.out=100)

# from me
Output.all <- NULL
SEs.all <- NULL
n <- 10000
lengths = c(17, 33, 65)

# cheating only, cheating plus... creative, guessing, careless, random
types = c("ch", "cr","g","ca","r")

amounts = c(0.05, 0.1, 0.25, 0.50)

# construct the matrix of test conditions
ds_conditions <- data.frame(
  length = c(
    rep(lengths[1], 20),
    rep(lengths[2], 20),
    rep(lengths[3], 20)
  ),
  type = c(
    rep(
      c(
        rep(types[1], 4),
        rep(types[2], 4),
        rep(types[3], 4),
        rep(types[4], 4),
        rep(types[5], 4)
      ),
      3
    )
  ),
  amount = c(
    rep(amounts, 15)
  )
)

trueitemdifficulties.17 = seq(-2, 2, length.out = 17)
trueitemdifficulties.33 = seq(-2, 2, length.out = 33)
trueitemdifficulties.65 = seq(-2, 2, length.out = 65)
statsK03=list() #PFSs of all examinees over 60 data sets

# now, construct datasets and calculate stats
# now, construct datasets and calculate stats
start <- Sys.time()
for(i in 1:60){
  condition <- ds_conditions[i,]
  length <- condition$length
  diffs <- switch(
    as.character(length),
    "17"=trueitemdifficulties.17,
    "33"=trueitemdifficulties.33,
    "65"=trueitemdifficulties.65
  )
  type = condition$type
  amt = 10000*condition$amount
  
  examinees <- data.frame(
    aberrant = rep(FALSE, n),
    type = rep(type, n)
  )
  
  examinees$theta <- 0
  # insert columns for item responses
  for(j in 1:length){
    examinees[[paste("i", j, sep="")]] <- -1
  }
  
  examinees$type = "normal"
  
  # next, designate a certain number of examinees as aberrant according to amt,
  # and simulate their aberrant responses
  if(type=="ch"){
    # select amt cheaters from the examinees with ability below -0.5
    examinees[1:amt,]$aberrant <- TRUE
    # assign thetas to cheaters from uniform -2:-0.5
    examinees[1:amt,]$theta <- runif(amt, -2, -0.5)
    # assign thetas to noncheaters from uniform -2:2
    examinees[amt+1:10000,]$theta <- runif(amt, -2, 2)
    # sort by ability
    examinees <- examinees[order(examinees$theta),]
    
    # select hardest items as "compromised"
    compromised <- which(diffs>=1.5)
    
    # now, simulate responses
    for(k in 1:length){
      item <- paste("i", k, sep="")
      examinees[[item]] <- apply(
        examinees,
        1,
        function(drow){
          ifelse(
            drow["aberrant"],
            ifelse(
              k %in% compromised,
              1,
              rbernoulli(
                1, 
                1/(1+exp(diffs[k]-as.numeric(drow["theta"])))
              ) %>% as.numeric()
            ),
            rbernoulli(
              1, 
              1/(1+exp(diffs[k]-as.numeric(drow["theta"])))
            ) %>% as.numeric()
          )
        }
      )
    }
  } else if(type=="cr") { # creative
    # select amt creative from the examinees with ability >= 0.5
    examinees[1:amt,]$aberrant <- TRUE
    # assign thetas to creative from uniform 0.5:2
    examinees[1:amt,]$theta <- runif(amt, 0.5, 2)
    # assign thetas to normal from uniform -2:2
    examinees[amt+1:10000,]$theta <- runif(amt, -2, 2)
    # sort by ability
    examinees <- examinees[order(examinees$theta),]
    
    # now, simulate  answering creatively/cheating
    # then simulate the normal respondents via the Rasch model
    for(k in 1:length){
      item <- paste("i", k, sep="")
      examinees[[item]] <- apply(
        examinees,
        1,
        function(drow){
          ifelse(
            drow["aberrant"],
            ifelse(
              diffs[k] <= -1.5,
              0,
              rbernoulli(
                1, 
                1/(1+exp(diffs[k]-as.numeric(drow["theta"])))
              ) %>% as.numeric()
            ),
            rbernoulli(
              1, 
              1/(1+exp(diffs[k]-as.numeric(drow["theta"])))
            ) %>% as.numeric()
          )
        }
      )
    }
  } else if(type=="g") { # guessing
    # select amt guessing from the examinees with ability below -0.5
    examinees[1:amt,]$aberrant <- TRUE
    # assign thetas to creative from uniform 0.5:2
    examinees[1:amt,]$theta <- runif(amt, -2, -0.5)
    # assign thetas to normal from uniform -2:2
    examinees[amt+1:10000,]$theta <- runif(amt, -2, 2)
    # sort by ability
    examinees <- examinees[order(examinees$theta),]
    
    # now, simulate guessing
    # then simulate the normal respondents via the Rasch model
    for(k in 1:length){
      item <- paste("i", k, sep="")
      examinees[[item]] <- apply(
        examinees,
        1,
        function(drow){
          ifelse(
            drow["aberrant"],
            ifelse(
              diffs[k] >= 0.5,
              rbernoulli(1, 0.25) %>% as.numeric(),
              rbernoulli(
                1, 
                1/(1+exp(diffs[k]-as.numeric(drow["theta"])))
              ) %>% as.numeric()
            ),
            rbernoulli(
              1, 
              1/(1+exp(diffs[k]-as.numeric(drow["theta"])))
            ) %>% as.numeric()
          )
        }
      )
    }
  } else if(type=="ca") { # careless
    # select amt/2 careless from the examinees with ability >= 0.5
    examinees[1:amt,]$aberrant <- TRUE
    # assign thetas to careless from uniform 0.5:2
    examinees[1:amt,]$theta <- runif(amt, 0.5, 2)
    # assign thetas to normal from uniform -2:2
    examinees[amt+1:10000,]$theta <- runif(amt, -2, 2)
    # sort by ability
    examinees <- examinees[order(examinees$theta),]
    
    # now, simulate  answering carelessly
    # then simulate the normal respondents via the Rasch model
    for(k in 1:length){
      item <- paste("i", k, sep="")
      examinees[[item]] <- apply(
        examinees,
        1,
        function(drow){
          ifelse(
            drow["aberrant"],
            ifelse(
              diffs[k] <= -0.5,
              rbernoulli(
                1, 
                0.5
              ) %>% as.numeric(),
              rbernoulli(
                1, 
                1/(1+exp(diffs[k]-as.numeric(drow["theta"])))
              ) %>% as.numeric()
            ),
            rbernoulli(
              1, 
              1/(1+exp(diffs[k]-as.numeric(drow["theta"])))
            ) %>% as.numeric()
          )
        }
      )
    }
  } else { # random
    # select amt random from the examinees with ability below -0.5
    examinees[1:amt,]$aberrant <- TRUE
    # assign thetas to creative from uniform
    examinees[1:amt,]$theta <- runif(amt, -2, -0.5)
    # assign thetas to normal from uniform -2:2
    examinees[amt+1:10000,]$theta <- runif(amt, -2, 2)
    # sort by ability
    examinees <- examinees[order(examinees$theta),]
    
    # now, simulate them answering at random
    # then simulate the normal respondents via the Rasch model
    for(k in 1:length){
      item <- paste("i", k, sep="")
      examinees[[item]] <- apply(
        examinees,
        1,
        function(drow){
          ifelse(
            drow["aberrant"],
            rbernoulli(1, 0.25) %>% as.numeric(),
            rbernoulli(
              1, 
              1/(1+exp(diffs[k]-as.numeric(drow["theta"])))
            ) %>% as.numeric()
          )
        }
      )
    }
  }
  
  # now analyze simulated responses
  # fit the Rasch model using ltm
  responses <- examinees[,4:(length+3)]
  scores.rasch = rasch(responses, constraint = cbind(length + 1, 1)) #fit the model
  itparm=coef(scores.rasch)#Estimated item parameters
  b=itparm[,1]#Difficulty parameters
  a=itparm[,2]#Common slope parameter
  raw=apply(responses,1,sum)
  
  responses <- responses[raw>0 & raw<length,]#Remove those with 0 or full score: HT undefined for them
  examinees <- examinees[raw>0 & raw<length,]
  theta=examinees$theta[raw>0 & raw<length]
  #Estimate the examinee abilities using package 'irtoys'
  itparms=cbind(a,b,rep(0,length))
  thetaest=mlebme(responses,itparms)[,1]
  #Compute lz*, HT, U3 using R package Perfit
  lzs=lzstar(responses,Ability=thetaest,IP=itparms)
  Hts = Ht(responses)
  U3s=U3(responses,Ability=thetaest,IP=itparms)#U3 is large for aberrant examinees
  # compute eci4z* 
  eci4zs <- eci4zstar(responses, thetaest, itparms)
  # and xia & zheng's asymmetrically-weighted PFSs
  sha1.2s <- sha1.2(responses, thetaest, itparms)
  sha1s <- sha1(responses, thetaest, itparms)
  shb2s <- shb2(responses, thetaest, itparms)
  shb3s <- shb3(responses, thetaest, itparms)
  
  # combine them all together; negative for PFSs where aberrance is indicated with a large value
  stats=cbind(
    lzs[[1]]$PFscores,
    Hts[[1]]$PFscores,
    -U3s[[1]]$PFscores,
    -eci4zs,
    -sha1.2s,
    -sha1s,
    -shb2s,
    -shb3s
  )
  
  statsK03[[paste(type, length, amt)]] <- stats
  stats=scale(stats,center=TRUE,scale=TRUE)#Standardize each PFS
  
  ## general aberrance detection rates
  good=stats[which(examinees$aberrant == FALSE),]#PFSs for the non-aberrant examinees
  bad=stats[which(examinees$aberrant == TRUE),]#PFSs for the aberrant examinees
  ROCOut=ROC(good,bad,x)#Compute ROC Areas using the function `ROC'
  Areas=ROCOut[1,]#ROC Areas of the PFSs
  F=ROCOut[3:(length(x)+2),]#False-alarm rates
  H=ROCOut[(length(x)+3):nrow(ROCOut),]#Hit rates
  SEs.all=rbind(SEs.all,ROCOut[2,])
  Output.all=rbind(Output.all,c(length,as.character(type),100*condition$amount,Areas))
  
  ## save item responses for later
  ItemResponses[[paste(type, length, amt)]] <- examinees
}

write.csv(Output.all, "/Users/sanfordstudent/Documents/Other Work/Person fit simulation project/rep_output_all.csv")

# mean AUROC over all conditions

out.all <- read.csv("/Users/sanfordstudent/Documents/Other Work/Person fit simulation project/rep_output_all.csv")[,2:12]
names(out.all) <- c("length", "type", "percent", "lz*", "ht", "u3", "eci4z*", "sha1/2", "sha1", "shb2", "shb3")

out.all.summary <- out.all %>% pivot_longer(
  names_to = "stat",
  values_to = "AUROC",
  cols = c(
    `lz*`, ht, u3, `eci4z*`, `sha1/2`, sha1, shb2, shb3
  )
) %>% group_by(stat) %>%
  summarise(
    mean = mean(AUROC),
    sd = sd(AUROC),
    min = min(AUROC),
    max = max(AUROC)
  )

out.all.summary.length <- out.all %>% pivot_longer(
  names_to = "stat",
  values_to = "AUROC",
  cols = c(
    `lz*`, ht, u3, `eci4z*`, `sha1/2`, sha1, shb2, shb3
  )
) %>% group_by(length, stat) %>%
  summarise(
    mean = mean(AUROC),
    sd = sd(AUROC),
    min = min(AUROC),
    max = max(AUROC)
  )

out.all.summary.type <- out.all %>% pivot_longer(
  names_to = "stat",
  values_to = "AUROC",
  cols = c(
    `lz*`, ht, u3, `eci4z*`, `sha1/2`, sha1, shb2, shb3
  )
) %>% group_by(type, stat) %>%
  summarise(
    mean = mean(AUROC),
    sd = sd(AUROC),
    min = min(AUROC),
    max = max(AUROC)
  )

out.all.summary.pct <- out.all %>% pivot_longer(
  names_to = "stat",
  values_to = "AUROC",
  cols = c(
    `lz*`, ht, u3, `eci4z*`, `sha1/2`, sha1, shb2, shb3
  )
) %>% group_by(percent, stat) %>%
  summarise(
    mean = mean(AUROC),
    sd = sd(AUROC),
    min = min(AUROC),
    max = max(AUROC)
  )

end<-Sys.time()
print(end - start)





## Dersimonian-Laird version: overall for all
out.all.overall.dsl <- NULL
Areas <- out.all
Areas <- Areas %>% pivot_longer(names_to = "stat", values_to = "AUC", cols = c("lz*", "ht", "u3", "eci4z*", "sha1/2", "sha1", "shb2", "shb3"))
se.df <- data.frame(SEs.all)
names(se.df) <- c("lz*", "ht", "u3", "eci4z*", "sha1/2", "sha1", "shb2", "shb3")
se.df <- se.df %>% pivot_longer(names_to = "st", values_to = "se", cols = c("lz*", "ht", "u3", "eci4z*", "sha1/2", "sha1", "shb2", "shb3"))

dsl.dat <- cbind(Areas, se.df)

#Use DerSimonian-Laird algorithm
for (stat in unique(dsl.dat$stat)){
  out.all.overall.dsl=c(
    out.all.overall.dsl,
    rma(
      yi = as.numeric(dsl.dat[which(dsl.dat$stat==stat),]$AUC),
      sei = as.numeric(dsl.dat[which(dsl.dat$stat==stat),]$se),
      # data = data.frame(yi = Areas[,i], sei = SEs.all[,i]),
      method="DL"
    )$beta
  )
}

out.all.overall.dsl <- data.frame(out.all.overall.dsl %>% t())

names(out.all.overall.dsl) <- c("lz*", "ht", "u3", "eci4z*", "sha1/2", "sha1", "shb2", "shb3")
out.all.overall.dsl %>% round(digits=2) %>% write.csv(
  "/Users/sanfordstudent/Documents/Other Work/Person fit simulation project/tables/rep/out_all_overall.csv"
)

## test length
out.all.length.dsl <- NULL
Areas <- out.all
Areas <- Areas %>% pivot_longer(names_to = "stat", values_to = "AUC", cols = c("lz*", "ht", "u3", "eci4z*", "sha1/2", "sha1", "shb2", "shb3"))
se.df <- data.frame(SEs.all)
names(se.df) <- c("lz*", "ht", "u3", "eci4z*", "sha1/2", "sha1", "shb2", "shb3")
se.df <- se.df %>% pivot_longer(names_to = "st", values_to = "se", cols = c("lz*", "ht", "u3", "eci4z*", "sha1/2", "sha1", "shb2", "shb3"))

dsl.dat <- cbind(Areas, se.df)

#Use DerSimonian-Laird algorithm
for(length in unique(dsl.dat$length)) {
  for (stat in unique(dsl.dat$stat)){
    auc <- rma(
      yi = as.numeric(dsl.dat[which(dsl.dat$stat==stat&dsl.dat$length==length),]$AUC),
      sei = as.numeric(dsl.dat[which(dsl.dat$stat==stat&dsl.dat$length==length),]$se),
      # data = data.frame(yi = Areas[,i], sei = SEs.all[,i]),
      method="DL"
    )$beta
    out.all.length.dsl <- rbind(out.all.length.dsl, c(length, stat, auc))
  }
}

out.all.length.dsl <- data.frame(out.all.length.dsl)

names(out.all.length.dsl) <- c("length", "stat", "auc")
out.all.length.dsl %>% mutate(auc = round(as.numeric(as.character(auc)), digits = 2)) %>%
  pivot_wider(names_from = stat, values_from = auc) %>% write.csv(
    "/Users/sanfordstudent/Documents/Other Work/Person fit simulation project/tables/rep/out_all_length.csv"
  )

## test type
out.all.type.dsl <- NULL
Areas <- out.all
Areas <- Areas %>% pivot_longer(names_to = "stat", values_to = "AUC", cols = c("lz*", "ht", "u3", "eci4z*", "sha1/2", "sha1", "shb2", "shb3"))
se.df <- data.frame(SEs.all)
names(se.df) <- c("lz*", "ht", "u3", "eci4z*", "sha1/2", "sha1", "shb2", "shb3")
se.df <- se.df %>% pivot_longer(names_to = "st", values_to = "se", cols = c("lz*", "ht", "u3", "eci4z*", "sha1/2", "sha1", "shb2", "shb3"))

dsl.dat <- cbind(Areas, se.df)

#Use DerSimonian-Laird algorithm
for(type in unique(dsl.dat$type)) {
  for (stat in unique(dsl.dat$stat)){
    auc <- rma(
      yi = as.numeric(dsl.dat[which(dsl.dat$stat==stat),]$AUC),
      sei = as.numeric(dsl.dat[which(dsl.dat$stat==stat),]$se),
      # data = data.frame(yi = Areas[,i], sei = SEs.all[,i]),
      method="DL"
    )$beta
    out.all.type.dsl <- rbind(out.all.type.dsl, c(type, stat, auc))
  }
}

out.all.type.dsl <- data.frame(out.all.type.dsl)

names(out.all.type.dsl) <- c("type", "stat", "auc")
out.all.type.dsl %>% mutate(auc = round(as.numeric(as.character(auc)), digits = 2)) %>%
  pivot_wider(names_from = stat, values_from = auc) %>% write.csv(
    "/Users/sanfordstudent/Documents/Other Work/Person fit simulation project/tables/rep/out_all_type.csv"
  )

## test pct
out.all.pct.dsl <- NULL
Areas <- out.all
Areas <- Areas %>% pivot_longer(names_to = "stat", values_to = "AUC", cols = c("lz*", "ht", "u3", "eci4z*", "sha1/2", "sha1", "shb2", "shb3"))
se.df <- data.frame(SEs.all)
names(se.df) <- c("lz*", "ht", "u3", "eci4z*", "sha1/2", "sha1", "shb2", "shb3")
se.df <- se.df %>% pivot_longer(names_to = "st", values_to = "se", cols = c("lz*", "ht", "u3", "eci4z*", "sha1/2", "sha1", "shb2", "shb3"))

dsl.dat <- cbind(Areas, se.df)

#Use DerSimonian-Laird algorithm
for(pct in unique(dsl.dat$percent)) {
  for (stat in unique(dsl.dat$stat)){
    auc <- rma(
      yi = as.numeric(dsl.dat[which(dsl.dat$stat==stat&dsl.dat$percent==pct),]$AUC),
      sei = as.numeric(dsl.dat[which(dsl.dat$stat==stat&dsl.dat$percent==pct),]$se),
      # data = data.frame(yi = Areas[,i], sei = SEs.all[,i]),
      method="DL"
    )$beta
    out.all.pct.dsl <- rbind(out.all.pct.dsl, c(pct, stat, auc))
  }
}

out.all.pct.dsl <- data.frame(out.all.pct.dsl)

names(out.all.pct.dsl) <- c("pct", "stat", "auc")
out.all.pct.dsl %>% mutate(auc = round(as.numeric(as.character(auc)), digits = 2)) %>%
  pivot_wider(names_from = stat, values_from = auc) %>% write.csv(
    "/Users/sanfordstudent/Documents/Other Work/Person fit simulation project/tables/rep/out_all_pct.csv"
  )













#### Further investigation: counts of true and false positives

thetastats <- list()
counts.flagged <- list()
for(r in names(statsK03)){
  # 1) join examinees' theta and response type to table of raw statistics
  st <- statsK03[[r]]
  th <- ItemResponses[[r]][,1:3]
  thetastats[[r]] <- cbind(th, st)
  names(thetastats[[r]]) <- c("aberrant", "type", "theta", "lz*", "ht", "u3", "eci4z*", "sha1/2", "sha1", "shb2", "shb3")
  condition <- strsplit(r, " ")
  counts.flagged[[r]] <- data.frame(
    type = condition[[1]][1],
    naber = condition[[1]][3],
    length = condition[[1]][2],
    sha1.2.cheating = length(which(thetastats[[r]]$type=="cheating" & thetastats[[r]]$`sha1/2` < -1.64)),
    sha1.2.not = length(which(thetastats[[r]]$type!="cheating" & thetastats[[r]]$`sha1/2` < -1.64)),
    
    sha1.cheating = length(which(thetastats[[r]]$type=="cheating" & thetastats[[r]]$sha1 < -1.64)),
    sha1.not = length(which(thetastats[[r]]$type!="cheating" & thetastats[[r]]$sha1 < -1.64)),
    
    shb2.cheating = length(which(thetastats[[r]]$type=="cheating" & thetastats[[r]]$shb2 < -1.64)),
    shb2.not = length(which(thetastats[[r]]$type!="cheating" & thetastats[[r]]$shb2 < -1.64)),
    
    shb3.cheating = length(which(thetastats[[r]]$type=="cheating" & thetastats[[r]]$shb3 < -1.64)),
    shb3.not = length(which(thetastats[[r]]$type!="cheating" & thetastats[[r]]$shb3 < -1.64)),
    
    lzstar.cheating = length(which(thetastats[[r]]$type=="cheating" & thetastats[[r]]$`lz*` < -1.64)),
    lzstar.not = length(which(thetastats[[r]]$type!="cheating" & thetastats[[r]]$`lz*` < -1.64)),
    
    ht.cheating = length(which(thetastats[[r]]$type=="cheating" & scale(thetastats[[r]]$ht) < -1.64)),
    ht.not = length(which(thetastats[[r]]$type!="cheating" & scale(thetastats[[r]]$ht) < -1.64)),
    
    eci.cheating = length(which(thetastats[[r]]$type=="cheating" & thetastats[[r]]$`eci4z*` < -1.64)),
    eci.not = length(which(thetastats[[r]]$type!="cheating" & thetastats[[r]]$`eci4z*` < -1.64)),
    
    u3.cheating = length(which(thetastats[[r]]$type=="cheating" & scale(thetastats[[r]]$u3) < -1.64)),
    u3.not = length(which(thetastats[[r]]$type!="cheating" & scale(thetastats[[r]]$u3) < -1.64))
  )
}

counts.concat <- do.call(what = "rbind", counts.flagged)

counts.bypct <- counts.concat %>%
  mutate(
    sha1.2 = sha1.2.cheating + sha1.2.not, 
    sha1  = sha1.cheating + sha1.not, 
    shb2 = shb2.cheating + shb2.not,
    shb3 = shb3.cheating + shb3.not,
    lzstar = lzstar.cheating + lzstar.not,
    ht = ht.cheating + ht.not,
    eci = eci.cheating + eci.not,
    u3 = u3.cheating + u3.not,
    sha1.2pct = sha1.2.cheating/sha1.2, 
    sha1pct = sha1.cheating/sha1, 
    shb2pct = shb2.cheating/shb2, 
    shb3pct = shb3.cheating/shb3, 
    lzstarpct = lzstar.cheating/lzstar,
    htpct = ht.cheating/ht,
    ecipct = eci.cheating/eci,
    u3pct = u3.cheating/u3
  ) %>% dplyr::select(
    type, naber, length, sha1.2, sha1, shb2, shb3, lzstar, ht, eci, u3,
    sha1.2pct, sha1pct, shb2pct, shb3pct, lzstarpct, htpct, ecipct, u3pct
  ) %>% 
  pivot_longer(
    names_to = "key",
    values_to = "value",
    cols = c(
      sha1.2, sha1, shb2, shb3, lzstar, ht, eci, u3,
      sha1.2pct, sha1pct, shb2pct, shb3pct, lzstarpct, htpct, ecipct, u3pct
    )
  ) %>% group_by(
    naber, key
  ) %>% summarise(
    mean = mean(value)
  ) %>% pivot_wider(
    names_from = naber, values_from = mean
  )


write.csv(
  counts.bypct, 
  "/Users/sanfordstudent/Documents/Other Work/Person fit simulation project/tables/rep/counts_by_pct.csv"
)




## replicate Sandip's approach to aggregating over test length
l17 <- NULL
l17.cheaters <- c()
l17.noncheaters <- c()
l33 <- NULL
l33.cheaters <- c()
l33.noncheaters <- c()
l65 <- NULL
l65.cheaters <- c()
l65.noncheaters <- c()

# create combined datasets
for(dataset in names(ItemResponses)){
  length <- (strsplit(dataset, " "))[[1]][2]
  type = (strsplit(dataset, " "))[[1]][1]
  dat <- drop_na(ItemResponses[[dataset]])
  
  if(length=="17") {
    # note cheaters
    if(type=="ch"){
      l17.cheaters <- append(l17.cheaters, (ifelse(is.null(nrow(l17)), 0, nrow(l17))+which(dat$aberrant==TRUE)))
      l17.noncheaters <- append(l17.noncheaters, (ifelse(is.null(nrow(l17)), 0, nrow(l17))+which(dat$aberrant==FALSE)))
    } else{
      l17.noncheaters <- append(l17.noncheaters, (ifelse(is.null(nrow(l17)), 0, nrow(l17))+1:nrow(dat)))
    }
    l17 <- rbind(l17, dat)
  } else if (length=="33") {
    # note cheaters
    if(type=="ch"){
      l33.cheaters <- append(l33.cheaters, (ifelse(is.null(nrow(l33)), 0, nrow(l33))+which(dat$aberrant==TRUE)))
      l33.noncheaters <- append(l33.noncheaters, (ifelse(is.null(nrow(l33)), 0, nrow(l33))+which(dat$aberrant==FALSE)))
    } else{
      l33.noncheaters <- append(l33.noncheaters, (ifelse(is.null(nrow(l33)), 0, nrow(l33))+1:nrow(dat)))
    }
    l33 <- rbind(l33, dat)
  } else if (length=="65") {
    # note cheaters
    if(type=="ch"){
      l65.cheaters <- append(l65.cheaters, (ifelse(is.null(nrow(l65)), 0, nrow(l65))+which(dat$aberrant==TRUE)))
      l65.noncheaters <- append(l65.noncheaters, (ifelse(is.null(nrow(l65)), 0, nrow(l65))+which(dat$aberrant==FALSE)))
    } else{
      l65.noncheaters <- append(l65.noncheaters, (ifelse(is.null(nrow(l65)), 0, nrow(l65))+1:nrow(dat)))
    }
    l65 <- rbind(l65, dat)
  }
}

# estimate difficulties
responses.17 <- l17[,4:20]
responses.33 <- l33[,4:36]
responses.65 <- l65[,4:68]

scores.rasch.17 = rasch(responses.17, constraint = cbind(18, 1)) #fit the model
scores.rasch.33 = rasch(responses.33, constraint = cbind(34, 1)) #fit the model
scores.rasch.65 = rasch(responses.65, constraint = cbind(66, 1)) #fit the model

# 17 items
itparm.17=coef(scores.rasch.17)#Estimated item parameters
b.17=itparm.17[,1]#Difficulty parameters
a.17=itparm.17[,2]#Common slope parameter

#Estimate the examinee abilities using package 'irtoys'
itparms.17=cbind(a.17,b.17,rep(0,17))
thetaest.17=mlebme(responses.17,itparms.17)[,1]

#Compute lz*, HT, U3 using R package Perfit
lzs=lzstar(responses.17,Ability=thetaest.17,IP=itparms.17)
Hts = Ht(responses.17)
U3s=U3(responses.17,Ability=thetaest.17,IP=itparms.17)#U3 is large for aberrant examinees
# compute eci4z* 
eci4zs <- eci4zstar(responses.17, thetaest.17, itparms.17)
# and xia & zheng's asymmetrically-weighted PFSs
sha1.2s <- sha1.2(responses.17, thetaest.17, itparms.17)
sha1s <- sha1(responses.17, thetaest.17, itparms.17)
shb2s <- shb2(responses.17, thetaest.17, itparms.17)
shb3s <- shb3(responses.17, thetaest.17, itparms.17)

# combine them all together; negative for PFSs where aberrance is indicated with a large value
stats.17=cbind(
  lzs[[1]]$PFscores,
  Hts[[1]]$PFscores,
  -U3s[[1]]$PFscores,
  -eci4zs,
  -sha1.2s,
  -sha1s,
  -shb2s,
  -shb3s
)

stats.17=scale(stats.17,center=TRUE,scale=TRUE)#Standardize each PFS

## general aberrance detection rates
good=stats.17[which(l17$aberrant == FALSE),]#PFSs for the non-aberrant examinees
bad=stats.17[which(l17$aberrant == TRUE),]#PFSs for the aberrant examinees
ROCOut=ROC(good,bad,x)#Compute ROC Areas using the function `ROC'
Areas=ROCOut[1,]#ROC Areas of the PFSs
F=ROCOut[3:19,]#False-alarm rates
H=ROCOut[20:nrow(ROCOut),]#Hit rates
SEs.all.17=ROCOut[2,]
Output.all.17=Areas

## cheating detection
good=stats.17[l17.noncheaters,]#PFSs for the non-aberrant examinees
bad=stats.17[l17.cheaters,]#PFSs for the aberrant examinees
ROCOut=ROC(good,bad,x)#Compute ROC Areas using the function `ROC'
Areas=ROCOut[1,]#ROC Areas of the PFSs
F=ROCOut[3:19,]#False-alarm rates
H=ROCOut[20:nrow(ROCOut),]#Hit rates
SEs.ch.17=ROCOut[2,]
Output.ch.17=Areas

# 33 items
itparm.33=coef(scores.rasch.33)#Estimated item parameters
b.33=itparm.33[,1]#Difficulty parameters
a.33=itparm.33[,2]#Common slope parameter

#Estimate the examinee abilities using package 'irtoys'
itparms.33=cbind(a.33,b.33,rep(0,33))
thetaest.33=mlebme(responses.33,itparms.33)[,1]

#Compute lz*, HT, U3 using R package Perfit
lzs=lzstar(responses.33,Ability=thetaest.33,IP=itparms.33)
Hts = Ht(responses.33)
U3s=U3(responses.33,Ability=thetaest.33,IP=itparms.33)#U3 is large for aberrant examinees
# compute eci4z* 
eci4zs <- eci4zstar(responses.33, thetaest.33, itparms.33)
# and xia & zheng's asymmetrically-weighted PFSs
sha1.2s <- sha1.2(responses.33, thetaest.33, itparms.33)
sha1s <- sha1(responses.33, thetaest.33, itparms.33)
shb2s <- shb2(responses.33, thetaest.33, itparms.33)
shb3s <- shb3(responses.33, thetaest.33, itparms.33)

# combine them all together; negative for PFSs where aberrance is indicated with a large value
stats.33=cbind(
  lzs[[1]]$PFscores,
  Hts[[1]]$PFscores,
  -U3s[[1]]$PFscores,
  -eci4zs,
  -sha1.2s,
  -sha1s,
  -shb2s,
  -shb3s
)

stats.33=scale(stats.33,center=TRUE,scale=TRUE)#Standardize each PFS

## general aberrance detection rates
good=stats.33[which(l33$aberrant == FALSE),]#PFSs for the non-aberrant examinees
bad=stats.33[which(l33$aberrant == TRUE),]#PFSs for the aberrant examinees
ROCOut=ROC(good,bad,x)#Compute ROC Areas using the function `ROC'
Areas=ROCOut[1,]#ROC Areas of the PFSs
F=ROCOut[3:19,]#False-alarm rates
H=ROCOut[20:nrow(ROCOut),]#Hit rates
SEs.all.33=ROCOut[2,]
Output.all.33=Areas

## cheating detection
good=stats.33[l33.noncheaters,]#PFSs for the non-aberrant examinees
bad=stats.33[l33.cheaters,]#PFSs for the aberrant examinees
ROCOut=ROC(good,bad,x)#Compute ROC Areas using the function `ROC'
Areas=ROCOut[1,]#ROC Areas of the PFSs
F=ROCOut[3:19,]#False-alarm rates
H=ROCOut[20:nrow(ROCOut),]#Hit rates
SEs.ch.33=ROCOut[2,]
Output.ch.33=Areas

# 65 items
itparm.65=coef(scores.rasch.65)#Estimated item parameters
b.65=itparm.65[,1]#Difficulty parameters
a.65=itparm.65[,2]#Common slope parameter

#Estimate the examinee abilities using package 'irtoys'
itparms.65=cbind(a.65,b.65,rep(0,65))
thetaest.65=mlebme(responses.65,itparms.65)[,1]

#Compute lz*, HT, U3 using R package Perfit
lzs=lzstar(responses.65,Ability=thetaest.65,IP=itparms.65)
Hts = Ht(responses.65)
U3s=U3(responses.65,Ability=thetaest.65,IP=itparms.65)#U3 is large for aberrant examinees
# compute eci4z* 
eci4zs <- eci4zstar(responses.65, thetaest.65, itparms.65)
# and xia & zheng's asymmetrically-weighted PFSs
sha1.2s <- sha1.2(responses.65, thetaest.65, itparms.65)
sha1s <- sha1(responses.65, thetaest.65, itparms.65)
shb2s <- shb2(responses.65, thetaest.65, itparms.65)
shb3s <- shb3(responses.65, thetaest.65, itparms.65)

# combine them all together; negative for PFSs where aberrance is indicated with a large value
stats.65=cbind(
  lzs[[1]]$PFscores,
  Hts[[1]]$PFscores,
  -U3s[[1]]$PFscores,
  -eci4zs,
  -sha1.2s,
  -sha1s,
  -shb2s,
  -shb3s
)

stats.65=scale(stats.65,center=TRUE,scale=TRUE)#Standardize each PFS

## general aberrance detection rates
good=stats.65[which(l65$aberrant == FALSE),]#PFSs for the non-aberrant examinees
bad=stats.65[which(l65$aberrant == TRUE),]#PFSs for the aberrant examinees
ROCOut=ROC(good,bad,x)#Compute ROC Areas using the function `ROC'
Areas=ROCOut[1,]#ROC Areas of the PFSs
F=ROCOut[3:19,]#False-alarm rates
H=ROCOut[20:nrow(ROCOut),]#Hit rates
SEs.all.65=ROCOut[2,]
Output.all.65=Areas

## cheating detection
good=stats.65[l65.noncheaters,]#PFSs for the non-aberrant examinees
bad=stats.65[l65.cheaters,]#PFSs for the aberrant examinees
ROCOut=ROC(good,bad,x)#Compute ROC Areas using the function `ROC'
Areas=ROCOut[1,]#ROC Areas of the PFSs
F=ROCOut[3:19,]#False-alarm rates
H=ROCOut[20:nrow(ROCOut),]#Hit rates
SEs.ch.65=ROCOut[2,]
Output.ch.65=Areas

# write out tables
all <- data.frame(rbind(Output.all.17, Output.all.33, Output.all.65))
names(all) <- c("lz*", "ht", "u3", "eci4z*", "sha1/2", "sha1", "shb2", "shb3")
write.csv(all, "/Users/sanfordstudent/Documents/Other Work/Person fit simulation project/tables/rep/overall_length_all.csv")


ch <- data.frame(rbind(Output.ch.17, Output.ch.33, Output.ch.65))
names(ch) <- c("lz*", "ht", "u3", "eci4z*", "sha1/2", "sha1", "shb2", "shb3")
write.csv(ch, "/Users/sanfordstudent/Documents/Other Work/Person fit simulation project/tables/rep/overall_length_cheating.csv")


## last thing: # flagged by each statistic and % cheating by test length
combined.17 <- stats.17 %>% data.frame()
names(combined.17) <- c("lz*", "ht", "u3", "eci4z*", "sha1/2", "sha1", "shb2", "shb3")
combined.17$cheating <- FALSE
combined.17[l17.cheaters,]$cheating <- TRUE
combined.17 <- combined.17 %>% pivot_longer(names_to = "stat", values_to = "value", cols = c(
  "lz*", "ht", "u3", "eci4z*", "sha1/2", "sha1", "shb2", "shb3"
))

combined.17.summary <- combined.17 %>% 
  mutate(
    flagged = as.numeric(value < -1.64)
  ) %>% group_by(stat) %>%
  summarise(
    flagged = sum(flagged),
    cheating = sum(flagged & cheating)/sum(flagged)
  )

combined.33 <- stats.33 %>% data.frame()
names(combined.33) <- c("lz*", "ht", "u3", "eci4z*", "sha1/2", "sha1", "shb2", "shb3")
combined.33$cheating <- FALSE
combined.33[l33.cheaters,]$cheating <- TRUE
combined.33 <- combined.33 %>% pivot_longer(names_to = "stat", values_to = "value", cols = c(
  "lz*", "ht", "u3", "eci4z*", "sha1/2", "sha1", "shb2", "shb3"
))

combined.33.summary <- combined.33 %>% 
  mutate(
    flagged = as.numeric(value < -1.64)
  ) %>% group_by(stat) %>%
  summarise(
    flagged = sum(flagged),
    cheating = sum(flagged & cheating)/sum(flagged)
  )

combined.65 <- stats.65 %>% data.frame()
names(combined.65) <- c("lz*", "ht", "u3", "eci4z*", "sha1/2", "sha1", "shb2", "shb3")
combined.65$cheating <- FALSE
combined.65[l65.cheaters,]$cheating <- TRUE
combined.65 <- combined.65 %>% pivot_longer(names_to = "stat", values_to = "value", cols = c(
  "lz*", "ht", "u3", "eci4z*", "sha1/2", "sha1", "shb2", "shb3"
))

combined.65.summary <- combined.65 %>% 
  mutate(
    flagged = as.numeric(value < -1.64)
  ) %>% group_by(stat) %>%
  summarise(
    flagged = sum(flagged),
    cheating = sum(flagged & cheating)/sum(flagged)
  )


combined.17.summary %>% write.csv("/Users/sanfordstudent/Documents/Other Work/Person fit simulation project/tables/rep/counts_17.csv")
combined.33.summary %>% write.csv("/Users/sanfordstudent/Documents/Other Work/Person fit simulation project/tables/rep/counts_33.csv")
combined.65.summary %>% write.csv("/Users/sanfordstudent/Documents/Other Work/Person fit simulation project/tables/rep/counts_65.csv")


