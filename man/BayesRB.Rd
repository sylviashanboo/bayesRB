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
MCMC_iteration: \tab the number of MCMC iterations. It should be an integer.\cr
burn_iteration: \tab the number of iterations as warming up. It should be an integer.\cr
thin_value: \tab thin value. It should be an integer.\cr
X: \tab a n rows x p columns matrix with n samples and p SNPs. For each SNP, the numbers of risk alleles are standardized as Norm(0,1). \cr
Y: \tab a vector with length of n, where n is the sample size. Y = 1 indicates a case, otherwise a control.\cr
beta_initial: \tab a vector with length p+1, containing the initial values of the grand mean and the SNP effects.The initial value of the grand mean goes first, followed by the initial values of the SNP effects. SNP order stays the same as the columns in the X matrix.\cr
}
}
\details{
BayesRB is used for SNP effect estimation and genetic risk prediction
}
\Values{
\tabular{ll}{
r_beta: \tab the estimated coefficients of all the SNPs in each MCMC loop after burn-in and thin (valid MCMC loop). The matrix contains I rows and P columns, where I is the number of valid MCMC loops, and P is the number of SNPs.\cr
r_bj: \tab the variance groups (1-4) to which the SNPs are estimated to belong in each valid MCMC loop. The matrix contains I rows and P columns, where I is the number of valid MCMC loops, and P is the number of SNPs.\cr
r_lambda:\tab the estimated lambda parameter for each sample in each valid MCMC loop.The matrix contains I rows and N columns, where I is the number of valid MCMC loops, and N is the number of samples.\cr
r_Z:\tab the estimated auxiliary variable Z for each sample in each valid MCMC loop.The matrix contains I rows and N columns, where I is the number of valid MCMC loops, and N is the number of samples.\cr
r_mu:\tab the estimated grand mean variable mu in each valid MCMC loop.\cr
r_sigma2:\tab the estimated genetic variance sigma2 in each valid MCMC loop.\cr
r_pi:\tab the estimated mixture proportions pi1-pi4 in each valid MCMC loop. The matrix contains I rows and 4 columns, where I is the number of valid MCMC loops.\cr
r_m:\tab the estimated SNPs being in the category 1-4 in each valid MCMC loop. The matrix contains I rows and 4 columns, where I is the number of valid MCMC loops.\cr
}
}
\author{
Ying Shan <sylvia.shanboo@gmail.com>, Daniel E. Weeks <weeks@pitt.edu>
}
\seealso{
This package is based on the R package \code{\link[Rcpp]{Rcpp}} and the c++ library \code{\link[gsl]{gsl}} \cr
	
Mollandin F, Rau A, Croiseau P. An evaluation of the predictive performance and mapping power of the BayesR model for genomic prediction. G3 (Bethesda). 2021 Oct 19;11(11):jkab225. doi: 10.1093/g3journal/jkab225. PMID: 34849780; PMCID: PMC8527474. \cr
	
Moser G, Lee SH, Hayes BJ, Goddard ME, Wray NR, Visscher PM. Simultaneous discovery, estimation and prediction analysis of complex traits using a bayesian mixture model. PLoS Genet. 2015 Apr 7;11(4):e1004969. doi: 10.1371/journal.pgen.1004969. PMID: 25849665; PMCID: PMC4388571. \cr
}
\examples{
library(BayesRB)

#The following libraries are used for constructing the dataset 
library(faraway)
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
#PRS:
est_Z = mean(result$r_mu) + x%*%apply(result$r_beta,2,mean)
#The probability of being a case:
prs_hat = apply(est_Z,1,ilogit)

}
