library(BayesRB)
library(Rcpp)

#The following libraries are used for constructing the dataset 
library(mgrp)
library(faraway)
source('mod_SimulateGenotypes.R')
library(methods)

sample.size <- 1000
# Define the odds ratios and risk allele frequencies for the risk SNPs
beta = c(rnorm(5,0,0.001),rnorm(2,0,0.01),rnorm(13,0,0.0001),rep(0,30),seq(0.1,1.2,by=0.1))
ngenes = length(beta)
true_beta = beta
or <- exp(beta)
fpool <-  c(0.297,0.636,0.857,0.202,0.743,0.83,0.097,0.764,0.787,0.512,0.485,0.212,0.312,0.478,0.728,0.614,0.637,0.458,0.443)
f <- sample(fpool,ngenes,replace = T)
prev=0.2
g3 <- sg(or=or,f=f,p=prev,nog=ngenes,n=sample.size,varyEffects=TRUE)
gdat <- data.frame(g3$g)
gdat <- -1*(gdat-3)
x <- as.matrix(gdat)
y <- as.vector(g3$disease)

#initialize beta and mu calculated by logistic regression:
mu=rep(NA,ngenes)
beta_lr=rep(NA,ngenes)
for(i in 1:ngenes){
  eq <- paste("y ~ ",paste0("X",i),collapse="")
  m3 <- glm(eq, family="binomial", data=gdat)
  mu[i] = m3$coef[1]
  if(is.na(m3$coef[2])){
    beta_lr[i] = 0.001
  }else{
    beta_lr[i] = m3$coef[2]
  }
  print(paste("finished",i,sep=""))
}
beta_init = c(mean(mu,na.rm=T),beta_lr)

mcmc_all = 41000
burn_in_num = 1000
thin_num= 50

result <- BayesRB(123, mcmc_all, burn_in_num, thin_num, x, y, beta_init)

