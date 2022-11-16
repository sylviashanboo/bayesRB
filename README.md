# Overview
BayesRB is a novel MCMC based polygenic genetic risk score algorithm for dichotomous traits. The R package under the same name is for SNP effect estimation and genetic risk prediction for dichotomous traits. BayesRB is an extension of BayesR proposed by Moser et al. BayesR method performs well on the SNP effect estimation and genetic risk prediction, but it is designed to be applied to the data with quantitative traits. BayesRB allows the dichotomous outcomes. It inherited the characteristics of unbiasedness, accuracy, sparseness, robustness and powerfulness. 

BayesRB R package is written by Rcpp.

# Installation

The package is written in R language. To install, proceed as follows:

1. Install the `devtools` and `faraway` packages by starting up R and issuing this command:

```
install.packages("devtools")
install.packages("faraway")
```

2. Load the `devtools` library to make its commands available:

```
library(devtools)
```

3. Install the GSL. Here is the [instruction](https://gist.github.com/TysonRayJones/af7bedcdb8dc59868c7966232b4da903).


4. Install the `BayesRB` R package from the github repository(the `main` branch) via this command:

```
install_github("sylviashanboo/bayesRB@main", build_vignettes = TRUE)
```

5. After the `BayesRB` R package has been installed, you can start to use the package:

```
library(BayesRB)
```

6. Run the following command for the details of how to use the package:
```
print(vignette("BayesRB-vignette",package="BayesRB"))
```

# Data Format

## The input data
* seed: the seed set for MCMC loop. It should be an integer.
* MCMC_iteration: the number of MCMC iterations. It should be an integer.
* burn_iteration: the number of interations as warming up. It should be an integer.
* thin_value: thin value. It should be an integer.
* X: the genotype file. It is an N rows x P columns matrix with N samples and P SNPs.
* Y: the phenotype file. It is a vector a length of N. Again N is the sample size. Y = 1 indicates a case, otherwise a control.
* beta_initial: a vector with a length of P+1.The initial value of the grand mean goes first, followed by the initial values of the P SNP effects. SNP order stays the same as the columns in the X matrix.

## The genotype file
The X matrix described above is the genotype file. It is an n rows x p columns matrix with N samples and P SNPs. For each SNP, the numbers of risk alleles are standardized as Norm(0,1).

| X1    | X2     | X3     | X4     | X5     | X6     |
|-------|--------|--------|--------|--------|--------|
| 0.553 | 1.338  | -0.944 | -0.391 | -1.257 | -0.463 |
| 0.553 | -0.072 | 0.766  | -0.391 | 0.129  | 1.935  |
| 0.553 | -1.482 | -0.944 | -0.391 | -1.257 | -0.463 |
| 0.553 | -0.072 | 0.766  | -1.845 | 0.129  | 1.935  |
| 0.553 | -0.072 | 0.766  | -1.845 | 0.129  | 1.935  |
| 0.553 | -0.072 | 0.766  | -0.391 | -1.257 | 1.935  |

## The phenotype file
The Y matrix described above is the phenotype file. It's a vector with N elements, which indicate the disease status (0,1) of the N samples. 

| V1 | V2 | V3 | V4 | V5 | V6 | V7 | V8 | V9 | V10 | V11 | V12 | V13 | V14 | V15 | V16 | V17 | V18 |
|----|----|----|----|----|----|----|----|----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
| 0  | 0  | 1  | 0  | 1  | 0  | 0  | 0  | 1  | 0   | 0   | 1   | 0   | 0   | 0   | 1   | 0   | 0   |

# How to Run the Program
The main functions in this package which implement the statistics described in Shan et al. (2022):

`BayesRB`

After installing the `BayesRB` R package, you can access BayesRB help page easily, which also contains example code.  
```
library(BayesRB)
print(vignette("BayesRB-vignette",package="BayesRB"))
```
OR
```
?BayesRB
```

# Explanation of the Results

## The output parameters
The results from BayesRB function are the estimated parameters in each MCMC loop after burn-in and thin (valid MCMC loop):

* r_beta: the estimated coefficients of all the SNPs in the valid MCMC loop. The matrix contains I rows and P columns, where I is the number of the valid MCMC loops, and P is the number of SNPs.
* r_bj: the variance groups (1-4) to which the SNPs are estimated to belong in each valid MCMC loop. The matrix contains I rows and P columns, where I is the number of the valid MCMC loops, and P is the number of SNPs.
* r_lambda: the estimated lambda parameter for each sample in each valid MCMC loop.The matrix contains I rows and N columns, where I is the number of the valid MCMC loops, and N is the number of samples.
* r_Z: the estimated auxiliary variable Z for each sample in each valid MCMC loop.The matrix contains I rows and N columns, where I is the number of the valid MCMC loops, and N is the number of samples.
* r_mu: the estimated grand mean variable mu in each valid MCMC loop.
* r_sigma2: the estimated genetic variance sigma2 in each valid MCMC loop.
* r_pi: the estimated mixture proportions pi1-pi4 in each valid MCMC loop. The matrix contains I rows and 4 columns, where I is the number of the valid MCMC loops.
* r_m: the estimated SNPs being in the category 1-4 in each valid MCMC loop. The matrix contains I rows and 4 columns, where I is the number of the valid MCMC loops.

## The calculation of PRS and probability of being a case
With the estimated parameters provided by BayesRB, PRS and probability of being a case can be easily calculated:

PRS:
```
data(result)
PRS = mean(result$r_mu) + x%*%apply(result$r_beta,2,mean)
hist(PRS)
```
The probability of being a case:
```
library(faraway) #for ilogit function
pro_hat = apply(PRS,1,ilogit)
hist(pro_hat)
```

# Citation
[Shan Y, Weeks D E. BayesRB: a markov chain Monte Carlo-based polygenic genetic risk score algorithm for dichotomous traits. bioRxiv, 2022, 2022.02.27.482193.](https://www.biorxiv.org/content/10.1101/2022.02.27.482193v1)


# Contact
Ying Shan, PhD

Clinical Research Academy, Peking University Shenzhen Hospital, Peking University, Shenzhen, China

BGI-Research, Shenzhen, China

Tel.: +86-755-83923333-6646

E-mail: sylvia.shanboo@gmail.com
