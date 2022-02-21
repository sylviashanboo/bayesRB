calcBaselineLOddsFromSample <-
function(or,f,nog,p,sampGenotypes) {

  # the population is assumed to be a matrix of 1,2,3 where 1 represents EE, 
  #     2 represents Ee, and 3 represents ee.  (E is the risk allele).
  #    So 3 minus sampGenotypes is the number of risk alleles for that gene.

  # this version does not consider all possible vectors m.
  # Instead, we use a population sample (with known prevalence)
  #    to "empirically" estimate P(genotype) ... ie P(genotype) = 
  #    either 0 or 1/samplesize
	
  # check consistency
  j = length(or)
  if((j != length(f)) | (j != length(nog))) {
    cat("parameter inconsistency in function calcBaselineLOddsFromSample")
    return
  } 
  
  # we don't use the whole sample for really big sample sizes...takes too long
  sampsize <- dim(sampGenotypes)[1]
  if(sampsize > 10000) {
    sampsize = 10000
  } 
    
  # now calculate logOR(m vs reference) for each possible m
  logORi = log(or)
  logORM = apply(sampGenotypes[1:sampsize,1:j],1,function(x) {
            sum=0 	  
            for(i in 1:j) sum = sum + (3-x[i] - 2*f[i]*nog[i])*logORi[i]
	    sum
	    })	    

  probM <- 1/sampsize
  
  expit <- function(x) exp(x) / (1+ exp(x)) 
  ff<-function(g){ 
    sum(expit(logORM + g)*probM) - p
  }
  g = uniroot(ff,c(100*log(p/(1-p)),0),tol=1e-20)$root

  list(or=or,f=f,p=p,bLO=g)  # g = log(baseline odds)
}
