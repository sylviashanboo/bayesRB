#include <R.h>
#include <Rmath.h>
#include <Rcpp.h>
//[[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "EigenRand/EigenRand"

//#include <RcppGSL.h> 
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sum.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <numeric>

#include <omp.h>

#include <cmath>
#include <algorithm>
using namespace std;
using namespace Rcpp;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXi;
using Eigen::ArrayXd;
using Eigen::ArrayXXd;
using Eigen::ArrayXXi;
using Eigen::Dynamic;
using Eigen::PermutationMatrix;

 int rightmost(double u, double lam)
 {
   double z=1;
   double x = exp(-0.5*lam);
   int j = 0;
   int OK;
   while(1)
   {
     j ++;
     z = z - (j+1)*(j+1)*pow(x,(j+1)*(j+1)-1);
     
     if (z > u){
       OK = 1;
       break;
     }
     j ++;
     z = z + (j+1)*(j+1)*pow(x,(j+1)*(j+1)-1);
     if (z < u){
       OK = 0;
       break;
     }
     
   }
   return (OK);
 }
 
 int leftmost(double u, double lam)
 {
   const double c_pi  = 3.141592653589793238463;
   double H = 0.5*log(2) + 2.5*log(c_pi) - 2.5*log(lam) - pow(c_pi,2)/(2*lam) + 0.5*lam;
   double lu = log(u);
   double z=1;
   double x = exp(-1*pow(c_pi,2)/(2*lam));
   double k = lam/pow(c_pi,2);
   
   int j =0;
   int OK;
   
   while(1){
     j++;
     z = z - k*pow(x,j*j-1);
     
     if(H+log(z)>lu){
       OK=1;
       break;
     }
     j++;
     z = z + (j+1)*(j+1)*pow(x,(j+1)*(j+1)-1);
     if(H+log(z)<lu){
       OK=0;
       break;
     }
     
   }
   return (OK);
 }  
 
 //[[Rcpp::export]]
 List BayesRB(int seed, int MCMC_inte, int burn_intee,int thinn, ArrayXXd X_unorm, ArrayXi Y, 
              ArrayXd beta_initial){
     Rcpp::Rcout << "hey we're starting" << endl;
     Eigen::Rand::P8_mt19937_64 rng {(std::size_t) seed};
     
     // gsl_rng * r;
     // r = gsl_rng_alloc(gsl_rng_mt19937);
     // gsl_rng_set(r, seed/3);
     gsl_rng_env_setup();
     gsl_rng **r = new gsl_rng *[omp_get_max_threads()];
     for (int b = 0; b < omp_get_max_threads(); b++) {
       r[b] = gsl_rng_alloc(gsl_rng_mt19937);
       gsl_rng_set(r[b], (b+1) * seed);
     }
     
     int N = X_unorm.rows(), P = X_unorm.cols();	
     // normalization?
     // ArrayXd std_dev = (X_unorm.rowwise() - X_unorm.colwise().mean()).square().rowwise().sum().sqrt() / (N - 1);
     // ArrayXXd X = (X_unorm.rowwise() - X_unorm.colwise().mean()).colwise() / std_dev;
     
     ArrayXXd X {N, P};
     for(int i = 0; i < P; i++){
       double std_dev2 = sqrt((X_unorm.col(i) - X_unorm.col(i).mean()).square().sum() / N);
       X.col(i) = (X_unorm.col(i) - X_unorm.col(i).mean()) / std_dev2;
     }
     
     ArrayXd lambda = ArrayXd::Zero(N);
     ArrayXd Z = ArrayXd::Zero(N);
     ArrayXd beta = beta_initial.segment(1, P);
     
     double mu = beta_initial(0);
     double sigma2 = gsl_rng_uniform_pos (r[0]) * 200;
     double sigma2_temp, u, logrrr;
     
     ArrayXd Ck {4};
     Ck << 0, 0.0001, 0.001, 0.01;
     
     ArrayXd pi {4};
     auto sumup = 1/Ck(1) + 1/Ck(2) + 1/Ck(3);
     
     pi(0) = 0.5;
     for (int i=1;i<=3;i++){
       pi(i) = (0.5 / Ck(i)) / sumup;
     }
     
     double theta = 1;
     
     ArrayXi m = ArrayXi::Zero(4);
     // MatrixXd Z_til = MatrixXd::Zero(N, P);
     // VectorXd data = VectorXd::LinSpaced(P, 0, P-1);
     ArrayXi bj {P};
     
     //output:
     
     int rest = MCMC_inte - burn_intee;
     int mod = rest % thinn;
     int store_num = (rest - mod)/ thinn;
     cout << store_num << endl;
     
     MatrixXd r_beta (store_num, P);
     Eigen::MatrixXi r_bj (store_num, P);
     MatrixXd r_lambda (store_num, N);
     MatrixXd r_Z (store_num, N);
     VectorXd r_mu (store_num);
     VectorXd r_sigma2 (store_num);
     MatrixXd r_pi (store_num, 4);
     Eigen::MatrixXi r_m (store_num, 4);
     
     double temp_a=0, temp_b=0;
     
     int num_of_store = 0;
     cout << "# of individuals:" << N <<"\n";
     cout << "# of SNPs:" << P << "\n";
     // cout << X << "\n";
     
     PermutationMatrix<Dynamic,Dynamic> permX(P);
     // ArrayXXd X_shuffle{N, P}, Z_til{N, P};
     ArrayXd xbeta_indi{N}, rand_mat{N}, h{P};
     // ArrayXd temp1{P}, temp2{P}; 
     
     // ArrayXXd temp1_add_sigma2k_inv {3, P}, beta_mu {4, P}, beta_var {4, P}, logL {4, P};
     // ArrayXXd log_pi_k {3, P}, log_temp1_mul_sigma2k {3, P}, beta_mu_mul_temp2 {3, P};
     // ArrayXXd threshold {4, P};
     
     //Start MCMC loop:
     for(int mcnum=1; mcnum < MCMC_inte; mcnum++){
       // cout << "mcmcnum: "<< mcnum <<"\n";
       
       //calculate sum_j(X_ij*beta[j])
       // xbeta_indi = (X.colwise() * beta.array()).colwise().sum();
       // #pragma omp parallel for shared ( xbeta_indi, loct )
       // for(int indi=0; indi<N; indi++){
       //   xbeta_indi(indi) = (X.row(N) * beta.transpose().array()).sum();
       //   loct(indi) = mu + xbeta_indi(indi);
       // }
       // xbeta_indi = (X.rowwise() * beta.transpose().array()).rowwise().sum();
       // cout << xbeta_indi << endl;
       // cout << "xbeta_indi \n" << xbeta_indi.transpose() << endl;
       // loct = mu + xbeta_indi;
       // cout << "loct \n" << loct.transpose() << endl;
       
       // rand_mat = Eigen::Rand::uniformReal<ArrayXd>(N, 1, rng, 1.0e-10, 1.);
       // rsam = 1 + (rand_mat/(1 - rand_mat)).log();
       // Z = rsam + loct;
       // cout << "rsam \n" << rsam.transpose() << endl;
       
       // rr = rsam.abs();
       // rr = (rr < 0.00001).select(0.00001, rr);
       // cout << "Z :\n" << Z.transpose() << endl;
       
       #pragma omp parallel for shared ( Y, Z, lambda, xbeta_indi )
       for(int indi=0; indi < N; indi++){
         double yy=0, uu=0, lambda_temp=0;
         int thread_num = omp_get_thread_num();
         
         xbeta_indi(indi) = (X.row(indi) * beta.transpose().array()).sum();
         double loct = mu + xbeta_indi(indi);
         double rsam = gsl_ran_logistic (r[thread_num],1);
         double z_indi = rsam + loct;
         
         if (Y(indi) == 1){
           while(z_indi < 0){
             rsam = gsl_ran_logistic(r[thread_num],1);
             z_indi = rsam + loct;
           }
         } else{
           while(z_indi >= 0){
             rsam = gsl_ran_logistic(r[thread_num],1);
             z_indi = rsam + loct;
           }
         } //finished Z
         Z(indi) = z_indi;
         
         double rr = abs(rsam);
         if(rr < 0.00001){
           rr = 0.00001;
         }
         
         int OK=0;
         while(OK==0){
           yy = gsl_ran_gaussian(r[thread_num],1);
           yy = yy*yy;
           yy = 1 + (yy-sqrt(yy*(4*rr+yy)))/(2*rr);
           while (yy==-1 || yy== 0){
             yy = gsl_ran_gaussian(r[thread_num],1);
             yy = yy*yy;
             yy = 1 + (yy-sqrt(yy*(4*rr+yy)))/(2*rr);
           }
           if (yy <= 0){
             yy = 0.00001;
           }
           uu = gsl_rng_uniform_pos (r[thread_num]);
           
           if (uu <= 1/(1+yy)){
             lambda_temp = rr/yy;
           }else{
             lambda_temp = rr*yy;
           }
           
           uu = gsl_rng_uniform_pos (r[thread_num]);
           if (lambda_temp > 4/3){
             
             OK = rightmost(uu, lambda_temp);
           }else{
             
             OK = leftmost(uu, lambda_temp);
           }
         }
         lambda[indi] = lambda_temp;
         //finish lambda
         ////////////////////////////////////////////////
       }//finish loop of individual;
       //finish step 1;
       // cout << "Y : \n" << Y.transpose() << endl;
       
       // cout << "lambda : \n" << lambda.transpose() << endl;
       temp_a = (Z.sum() - xbeta_indi.sum()) / N;
       temp_b = lambda.sum() / pow(N,2);
       mu = gsl_ran_gaussian(r[0], temp_b) + temp_a;

       permX.setIdentity();
       std::shuffle(permX.indices().data(), permX.indices().data()+permX.indices().size(), rng);
       /*
       X_shuffle = (X.matrix() * permX).array();
       beta = (beta.matrix().transpose() * permX).array();
       // cout << X_shuffle << endl; 
       Z_til = (Z - mu - xbeta_indi).replicate(1, P) + (X_shuffle.rowwise() * beta.transpose());
       
       // size of P, entire operation can be seen as a matmul
       // ArrayXd temp1 = (X_shuffle.pow(2).colwise() * (1/lambda)).colwise().sum(); 
       temp1 = lambda.inverse().transpose().matrix() * X_shuffle.pow(2).matrix();
       // size of P, convert to array for element wise op
       // the paper contains a power that was missing in the original code, ooof...
       temp2 = ((Z_til * X_shuffle).colwise() * lambda.inverse()).colwise().sum();
       // cout << "temp2 eigen: " << temp2(1) << endl;
       // cout << "temp2 og: " << (Z_til(all,1) * X_shuffle(all,1) / lambda).sum() << endl;
       // size of 3
       
       sigma2k = sigma2 * Ck.segment(1, 3);
       temp1_add_sigma2k_inv = (temp1.replicate(1, 3).transpose() + sigma2k.inverse().replicate(1, P)).inverse(); // size of (3, P)
       
       beta_mu.row(0) = VectorXd::Zero(P);
       beta_mu.middleRows(1, 3) = temp2.replicate(1, 3).transpose() * (temp1_add_sigma2k_inv);
       
       beta_var.row(0) = VectorXd::Zero(P);
       beta_var.middleRows(1, 3) = temp1_add_sigma2k_inv;

       log_pi_k = pi.segment(1, 3).log().replicate(1, P);
       log_temp1_mul_sigma2k = (temp1.replicate(1, 3).transpose() * sigma2k.replicate(1, P) + 1).log();
       beta_mu_mul_temp2 = beta_mu.middleRows(1, 3) * temp2.replicate(1, 3).transpose();
       
       logL.row(0) = VectorXd::Ones(P) * log(pi(0));
       logL.middleRows(1, 3) = log_pi_k - 0.5 * log_temp1_mul_sigma2k + 0.5 * beta_mu_mul_temp2;
       // cout << "logL :\n" << logL << endl;
       
       for (int k = 0; k < 4; k++){
         threshold.row(k) = (logL - logL.row(k).replicate(4, 1)).exp().colwise().sum().inverse();
       }

       // assign the bjs so we can use them for if-statements later
       bj = (h <= (threshold.row(0) + threshold.row(1) + threshold.row(2)).transpose()).select(ArrayXi::Constant(P, 3), 4);
       bj = (h <= (threshold.row(0) + threshold.row(1)).transpose()).select(2, bj);
       bj = (h <= threshold.row(0).transpose()).select(1, bj);
       // cout << "bj: \n" << bj.transpose() << endl;
       */
       // h = Eigen::Rand::uniformReal<ArrayXd>(P, 1, rng, 1.0e-20, 1);
       double vara = 0;
       int m0=0, m1=0, m2=0, m3=0;
       
       #pragma omp parallel for shared ( Ck, X, Z, mu, xbeta_indi, bj, beta, permX, pi ) reduction(+:m0,m1,m2,m3,vara)
       for (int snpji=0; snpji < P; snpji++){
            int snpj = permX.indices()[snpji], thread_num = omp_get_thread_num();
            ArrayXd Z_til = Z - mu - xbeta_indi + X.col(snpj) * beta(snpj);
            
            double temp1 = (X.col(snpj).square() / lambda).sum();
            double temp2 = (Z_til * X.col(snpj) / lambda).sum();
            double thresh[4], beta_mu[4], beta_var[4];
            ArrayXd logL = ArrayXd::Zero(4);
            
            for (int k = 0; k < 4; k++){
                if(k==0){
                    logL[k] = log(pi[k]);
                    beta_mu[k] = 0;
                    beta_var[k] = 0;
                }else{
                    double sigma2k = sigma2 * Ck[k];
                    beta_mu[k] = temp2/(temp1 + 1/sigma2k);
                    beta_var[k] = 1/(temp1 + 1/sigma2k);
                    logL[k] = log(pi[k]) - 0.5*log(temp1*sigma2k+1) + 0.5*(beta_mu[k]*temp2);
                }
            }
            
            for (int k = 0; k < 4; k++){
                thresh[k] = 1 / (logL - logL[k]).exp().sum();
            }
            
            double h_snpj = gsl_rng_uniform_pos(r[thread_num]);
            double beta_pre = beta[snpj];
            if (h_snpj <= thresh[0]){
                beta[snpj] = 0;
                m0+=1;
                bj[snpj]=1;
            }else if(h_snpj <= thresh[0] + thresh[1]){
                beta[snpj] = gsl_ran_gaussian(r[thread_num],beta_var[1])+beta_mu[1];
                m1+=1;
                vara += pow(beta[snpj],2)/(2*Ck[1]);
                bj[snpj]=2;
            }else if(h_snpj <= thresh[0] + thresh[1] + thresh[2]){
                beta[snpj] = gsl_ran_gaussian(r[thread_num],beta_var[2])+beta_mu[2];
                m2+=1;
                vara += pow(beta[snpj],2)/(2*Ck[2]);
                bj[snpj]=3;
            }else{
                beta[snpj] = gsl_ran_gaussian(r[thread_num],beta_var[3])+beta_mu[3];
                m3+=1;
                vara += pow(beta[snpj],2)/(2*Ck[3]);
                bj[snpj]=4;
            }
           //cout << "beta_j" << beta[snpj] <<"\n";
              // xbeta_indi = xbeta_indi - X.col(snpj)*beta_pre + X.col(snpj)*beta[snpj];
           // xbeta_indi -= X(all, snpj) * (beta_pre + beta(snpj));
          }
          m << m0, m1, m2, m3;
          
          //cout << "m[0,1,2,3]:"<<m[0]<<","<<m[1]<<","<<m[2]<<","<<m[3]<<"\n";
          if (m[1]==0 & m[2]==0 & m[3]==0){
              Rcpp::Rcout << "break, since all the SNPs are assigned to category 1\n";
              cout << "break, since all the SNPs are assigned to category 1\n";
            
              break;
          }
       
       ///////////////////////////////////////////////////////
       //finished Step 3 and 4: updated beta and bj.
       
       // MH sampling:
       double sigma2_temp = gsl_ran_gaussian(r[0],theta) + sigma2;
       while (sigma2_temp <= 0){
         sigma2_temp = gsl_ran_gaussian(r[0],theta) + sigma2;
       }
       
       double pnorm_sigma2 = R::pnorm(sigma2/sqrt(theta), 0, 1, 1, 0);
       double pnorm_sigma2_temp = R::pnorm(sigma2_temp/sqrt(theta), 0, 1, 1, 0);
       
       logrrr = -0.5*(m[1]+m[2]+m[3])*(-log(sigma2) + log(sigma2_temp)) 
         - vara*(-1/sigma2 + 1/sigma2_temp) + (log(pnorm_sigma2) - log(pnorm_sigma2_temp));
       
       // cout << "logrrr:" << logrrr <<"\n";
       
       u = gsl_rng_uniform_pos (r[0]);
       
       if(logrrr >= 0){
         sigma2 = sigma2_temp;
       }else if(exp(logrrr)>=u){
         sigma2 = sigma2_temp;
       }
       
       if (sigma2==0){
         Rcpp::Rcout << "break, since sigma_g^2 = 0\n";
         break;
       }
       // cout << "sigma2:" << sigma2 << "\n";
       ///////////////////////////////////////////////////////
       //finished step 5: sigma_g^2.
       
       ArrayXd m_add_one = m.cast<double>() + 1; 
       gsl_ran_dirichlet(r[0], 4, m_add_one.data(), pi.data());
       
       ///////////////////////////////////////////////////////
       //finished step 6: pi
       // Rcpp::Rcout << "finished mcmc loop"<< mcnum <<"\n";
       
       if (mcnum == burn_intee + thinn * num_of_store + thinn - 1){
         //r_beta:
         r_beta.row(num_of_store) = beta.transpose();
         //r_bj
         r_bj.row(num_of_store) = bj.transpose();
         //r_lambda 
         r_lambda.row(num_of_store) = lambda.transpose();
         // r_Z:
         r_Z.row(num_of_store) = Z.transpose();
         //r_mu:
         r_mu(num_of_store) = mu;
         //r_sigma2:
         r_sigma2(num_of_store)  = sigma2;
         //r_pi
         r_pi.row(num_of_store) = pi.transpose();
         //r_m:
         r_m.row(num_of_store) = m.transpose();
         
         num_of_store++;
       } 
     }
     return List::create(_["r_beta"] = Rcpp::wrap(r_beta), _["r_bj"] = Rcpp::wrap(r_bj), _["r_lambda"] = Rcpp::wrap(r_lambda),
                         _["r_Z"] = Rcpp::wrap(r_Z), _["r_mu"] = Rcpp::wrap(r_mu), _["r_sigma2"] = Rcpp::wrap(r_sigma2), 
                           _["r_pi"] = Rcpp::wrap(r_pi), _["r_m"] = Rcpp::wrap(r_m));
   }