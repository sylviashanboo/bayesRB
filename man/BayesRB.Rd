\name{BayesRB}
\alias{BayesRB}
\alias{BRB}
\title{
	A Bayesian approach for SNP effect estimation and genetic risk prediction for dichotomous traits
}
\description{
BayesRB is a Bayesian approach for SNP effect estimation and genetic risk prediction for dichotomous traits. It is an extension of BayesR proposed by Moser et al. BayesR method performs well on the SNP effect estimation and genetic risk prediction, but it is designed to be appied to the data with quantitative traits. BayesRB allows the dichotomous outcomes. It inherited the characteristics of unbiasedness, accuracy, sparseness, robustness and powerfulness. 
}
\usage{
BayesRB(seed, MCMC_iteration, burn_iteration, thin_value, X, Y, beta_initial)
}
\arguments{
\tabular{ll}{
seed: \tab an integer. \cr
MCMC_iteration: \tab the number of MCMC interations. It should be an integer.\cr
burn_iteration: \tab the number of interations as warming up. It should be an integer.\cr
thin_value: \tab thin value. It should be an integer.\cr
X: \tab a n*p matrix with n samples and p variables.\cr
Y: \tab a vector with length of n, where n is the sample size. Y=1 indicates a case, otherwise a control.\cr
beta_initial: \tab a vector containing the intial value of the grant mean and the SNP effects. grant mean goes first. SNP order stays the same as the columns of X matrix.
}
}
\details{
BayesRB is used for SNP effect estimation and genetic risk prediction
}
\author{
Ying Shan <yis29@pitt.edu>, Daniel E. Weeks <weeks@pitt.edu>
}
\seealso{
This package is based on the R package \code{\link[Rcpp]{Rcpp}} and the c++ library \code{\link[gsl]{gsl}}

Moser, G., Lee, S. H., Hayes, B. J., Goddard, M. E., Wray, N. R. and Visscher, P. M. (2015). Simultaneous discovery, estimation and prediction analysis of complex traits using a bayesian mixture model. PLoS genetics 11, e1004969.
}
\examples{
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


}
