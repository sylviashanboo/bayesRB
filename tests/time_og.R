setwd(".")
# renv::load(".")
library(BayesRB, verbose=T)

load("data/sample_data.rda")

mcmc_all = 41000
burn_in_num = 1000
thin_num= 50

# ptm <- proc.time()
# result <- BayesRB(123, mcmc_all, burn_in_num, thin_num, x, y, beta_init)
# time <- proc.time() - ptm
# print(time)

ptm2 <- proc.time()
result_og <- BayesRB_OG(123, mcmc_all, burn_in_num, thin_num, x, y, beta_init)
time2 <- proc.time() - ptm2
print(time2)

save(result, file="data/og_result.rda")