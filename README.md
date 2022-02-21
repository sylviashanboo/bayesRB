# BayesRB
BayesRB is a novel MCMC based polygenic genetic risk score algorithm for dichotomous traits. The R package under the same name is for SNP effect estimation and genetic risk prediction for dichotomous traits. BayesRB is an extension of BayesR proposed by Moser et al. BayesR method performs well on the SNP effect estimation and genetic risk prediction, but it is designed to be applied to the data with quantitative traits. BayesRB allows the dichotomous outcomes. It inherited the characteristics of unbiasedness, accuracy, sparseness, robustness and powerfulness. 

BayesRB R package is written by Rcpp.

## Installation
```
install.packages("devtools")
library(devtools)
devtools::install_github("sylviashanboo/bayesRB")
```
Or
```
R CMD INSTALL BayesRB_1.0.tar.gz
```

**Note**: You may be asked to download and install GNU Scientific Library (gsl) for C++ program. 

For MacOSX users:
Download gsl-latest.tar.gz from the [GSL ftp site](https://www.gnu.org/software/gsl/) or [this bayesRB GitHub site](https://github.com/sylviashanboo/bayesRB/tree/main/appendix%20files) and unzip it anywhere (e.g. /Downloads)

Open the unzipped gsl folder in Terminal (e.g. cd ~/Downloads/gsl-2.7.1)
```
sudo make clean
sudo chown -R $USER .
./configure && make
sudo make install
```

For other users: 

check [this website](https://gist.github.com/TysonRayJones/af7bedcdb8dc59868c7966232b4da903) for instructions.


## Citation
Shan, Y., 2016. Statistical methods for genetic risk confidence intervals, Bayesian disease risk prediction, and estimating mutation screening saturation (Doctoral dissertation, University of Pittsburgh).

## Contact
Ying Shan, PhD

Clinical Research Academy, Peking University Shenzhen Hospital, Peking University, Shenzhen, China

BGI-Research, Shenzhen, China

Tel.: +86-755-83923333-6646

E-mail: sylvia.shanboo@gmail.com
