#include <R.h>
#include <Rmath.h>
#include <Rcpp.h>
//[[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

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

int rightmost_og(double u, double lam)
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

int leftmost_og(double u, double lam)
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

List BayesRB_OG(int seed, int MCMC_inte, int burn_intee,int thinn, NumericMatrix X, NumericVector Y, 
	NumericVector beta_initial){
  cout << omp_get_thread_num() << endl;
	gsl_rng * r;
	r = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(r, seed);

	size_t N = X.nrow();
  	size_t P = X.ncol();

  	
  	
  	for (int i = 0; i < P; i ++){
  		X(_,i) = (X(_,i) - mean(X(_,i)))/sd(X(_,i));
  		
  	}	
  	
  	NumericVector lambda(N);
  	NumericVector Z(N);
  	NumericVector beta(P);
  	for (int i =1; i<=P; i++){
		beta[i-1] = beta_initial[i];
	}
  	double mu=beta_initial[0];
  	double sigma2=gsl_rng_uniform_pos (r)*200;
  	double sigma2_temp, u, logrrr;
  	
  	NumericVector Ck(4);
  	Ck[0] = 0;
  	Ck[1] = 0.0001;
  	Ck[2] = 0.001;
  	Ck[3] = 0.01;
  	
	//NumericVector pi(4);
	double* pi  = NULL;
	pi  = new double [4];
	double temp, sigma2k;
	double theta = 1;
	double sumup = 1/Ck[1] + 1/Ck[2] + 1/Ck[3];
	pi[0]=0.5;
	for (int i=1;i<=3;i++){
		temp=1/Ck[i];
		pi[i]=0.5*temp/sumup;
	}
	NumericVector m(4);
	NumericVector xbeta_indi(N);
	NumericMatrix Z_til(N, P);
	NumericVector data(P);
	for (int i=0;i<P;i++){
		data[i] = i;
	}
	NumericVector bj(P);

	//output:
	int rest = MCMC_inte - burn_intee;
	int mod = rest% thinn;
	int store_num = (rest - mod)/ thinn;
	MatrixXd r_beta{store_num,P};
	NumericMatrix r_bj(store_num,P);
	NumericMatrix r_lambda(store_num,N);
	NumericMatrix r_Z(store_num,N);
	NumericVector r_mu(store_num);
	NumericVector r_sigma2(store_num);
	NumericMatrix r_pi(store_num,4);
	NumericMatrix r_m(store_num,4);

	double xbeta, all_xbeta, loct, rr,yy,uu,lambda_temp,temp_a,temp_b,vara,  
	 temp1, temp2,h,rsam, beta_pre;
	int  snpj, num_pi, num_m;
	NumericVector thresh(4);
	NumericVector logL(4);
	NumericVector beta_mu(4);
	NumericVector beta_var(4);
	double* temp_d  = NULL;
	temp_d  = new double [4];
	
	int num_of_store = 0;
	cout << "# of individuals:" <<N<<"\n";
	cout << "# of SNPs:" << P << "\n";
	gsl_permutation * permut_p = gsl_permutation_alloc (P);
    gsl_permutation_init (permut_p);
	//Start MCMC loop:
	for(int mcnum=1; mcnum < MCMC_inte; mcnum++){
		// cout << "mcmcnum: "<<mcnum <<"\n";
		
    // #pragma omp parallel for
		for(int indi=0; indi < N; indi++){
			//calculate sum_j(X_ij*beta[j])
			xbeta = sum(X(indi,_)*beta);
			xbeta_indi[indi] = xbeta;//sum_j(Xij*beta_j)
			
			// cout << "mu: " << mu[mcnum-1]<<"\n";
			loct = mu + xbeta;
			
			rsam = gsl_ran_logistic (r,1);

			Z[indi] = rsam + loct;

			if (Y[indi]==1){
				while(Z[indi]<0){
					rsam = gsl_ran_logistic (r,1);
					Z[indi] = rsam + loct;
				}
			}
			else{
				while(Z[indi]>=0){
					rsam = gsl_ran_logistic (r,1);
					Z[indi] = rsam + loct;
				}
			}//finished Z
			// cout << "Z: "<<Z[indi]<<"\n";
			
			////////////////////////////////////////////////
			
			int OK=0;
			rr = abs(rsam);
			if(rr < 0.00001){
				rr = 0.00001;
			}
			while(OK==0){
				yy = gsl_ran_gaussian(r,1);
				yy = yy*yy;
				yy = 1 + (yy-sqrt(yy*(4*rr+yy)))/(2*rr);
				while (yy==-1 || yy== 0){
					yy = gsl_ran_gaussian(r,1);
					yy = yy*yy;
					yy = 1 + (yy-sqrt(yy*(4*rr+yy)))/(2*rr);
				}
				if (yy <= 0){
					yy = 0.00001;
				}
				uu = gsl_rng_uniform_pos (r);
				
				if (uu <= 1/(1+yy)){
					lambda_temp = rr/yy;
				}else{
					lambda_temp = rr*yy;
				}
				
				uu = gsl_rng_uniform_pos (r);
				if (lambda_temp > 4/3){
					
					OK = rightmost_og(uu, lambda_temp);
				}else{
					
					OK = leftmost_og(uu, lambda_temp);
				}
			}
			lambda[indi] = lambda_temp;
			//finish lambda;
			
			////////////////////////////////////////////////

		}//finish loop of individual;
		//finish step 1;
		// cout << xbeta_indi << endl;
		// cout << lambda << endl;
		all_xbeta = sum(xbeta_indi);//sum_i(sum_j(Xij*beta_j))

		temp_a = (sum(Z)-all_xbeta)/N;
		temp_b = sum(lambda)/pow(N,2);

		mu = gsl_ran_gaussian(r,temp_b)+temp_a;
		
		//finish mu;
		////////////////////////////////////////////////

		vara=0;
		
		for (num_m=0;num_m<4;num_m++){
			m[num_m] = 0;
		}
		//m = 0; // m is a vector;
		gsl_ran_shuffle (r, permut_p->data, P, sizeof(size_t));
		for (int snpji=0; snpji<P; snpji++){
			snpj = gsl_permutation_get(permut_p,snpji);
			Z_til(_,snpj) = Z - mu -  xbeta_indi + X(_,snpj)*beta[snpj];
			temp1 = sum(1/lambda * pow(X(_,snpj),2));
			temp2 = sum(Z_til(_,snpj) * X(_,snpj) /lambda); 

			for (int k = 0; k < 4; k++){
			
				if(k==0){
					logL[k] = log(pi[k]);
					beta_mu[k] = 0;
					beta_var[k] = 0;

				}else{
					sigma2k = sigma2*Ck[k];
					beta_mu[k] = temp2/(temp1 + 1/sigma2k);
					beta_var[k] = 1/(temp1 + 1/sigma2k);
					//logL[k]  = log(pi[k]) -0.5*log(sum(pow(X(_,snpj),2)/lambda)*Ck[k]*sigma2+1)
				//+0.5*(pow(sum(Z_til(_,snpj)*X(_,snpj)/lambda),2))/(sum(pow(X(_,snpj),2)/lambda)+1/(Ck[k]*sigma2));
					logL[k] = log(pi[k]) -0.5*log(temp1*sigma2k+1)+0.5*(beta_mu[k]*temp2);
				}

			}
			for (int k = 0; k < 4; k++){
				thresh[k] = 1/sum(exp(logL - logL[k]));
				
			}
			beta_pre = beta[snpj];
			h = gsl_rng_uniform_pos (r);

			if (h<=thresh[0]){
				beta[snpj] = 0;
				m[0]+=1;
				bj[snpj]=1;
				
			}else if(h<= thresh[0] + thresh[1]){
				beta[snpj] = gsl_ran_gaussian(r,beta_var[1])+beta_mu[1];
				m[1]+=1;
				vara += pow(beta[snpj],2)/(2*Ck[1]);
				bj[snpj]=2;
				
			}else if(h<= thresh[0] + thresh[1] + thresh[2]){
				beta[snpj] = gsl_ran_gaussian(r,beta_var[2])+beta_mu[2];
				m[2]+=1;
				vara += pow(beta[snpj],2)/(2*Ck[2]);
				bj[snpj]=3;

			}else{
				beta[snpj] = gsl_ran_gaussian(r,beta_var[3])+beta_mu[3];
				m[3]+=1;
				vara += pow(beta[snpj],2)/(2*Ck[3]);
				bj[snpj]=4;
			}
			//cout << "beta_j" << beta[snpj] <<"\n";
			xbeta_indi = xbeta_indi - X(_,snpj)*beta_pre + X(_,snpj)*beta[snpj];

			
		}
		
		// cout << "m[0,1,2,3]:"<<m[0]<<","<<m[1]<<","<<m[2]<<","<<m[3]<<"\n";
		if (m[1]==0 & m[2]==0 & m[3]==0){
			cout << "break, since all the SNPs are assigned to category 1\n";
			break;
		}
		///////////////////////////////////////////////////////
		//finished Step 3 and 4: updated beta and Ck.
	
		// MH sampling:
		//if (sigma2 < 0.000001){
		//	theta = 0.00000001;
		//}
		sigma2_temp = gsl_ran_gaussian(r,theta)+sigma2;
		while (sigma2_temp <= 0){
			sigma2_temp = gsl_ran_gaussian(r,theta)+sigma2;
		}
		// cout << "sigma2_temp:" << sigma2_temp << "\n";

		
		temp1 = R::pnorm(sigma2/sqrt(theta),0,1,1,0);
		temp2 = R::pnorm(sigma2_temp/sqrt(theta),0,1,1,0);
		

		logrrr = -0.5*(m[1]+m[2]+m[3])*(-log(sigma2) + log(sigma2_temp)) 
		- vara*(-1/sigma2 + 1/sigma2_temp) + (log(temp1) - log(temp2));
		// cout << "logrrr:" << logrrr <<"\n";
		
		u = gsl_rng_uniform_pos (r);

		if(logrrr >= 0){
			sigma2 = sigma2_temp;
		}else if(exp(logrrr)>=u){
			sigma2 = sigma2_temp;
		}

		if (sigma2==0){
			cout << "break, since sigma_g^2 = 0\n";
			break;
		}
		// cout << "sigma2:" << sigma2 << "\n";
		///////////////////////////////////////////////////////
		//finished step 5: sigma_g^2.
		
		for(int k=0;k<4;k++){
			temp_d[k] = m[k]+1;
		}
		
		gsl_ran_dirichlet (r,4,temp_d,pi);
		
		
		///////////////////////////////////////////////////////
		//finished step 6: pi
		//cout << "finished mcmc loop"<< mcnum <<"\n";
		
		if (mcnum == burn_intee + thinn * num_of_store + thinn - 1){
		  // cout << "mcmcnum: "<<mcnum <<"\n";
			//r_beta:
			r_beta.row(num_of_store) = Rcpp::as<Eigen::VectorXd>(beta).transpose();
			//r_bj
			r_bj(num_of_store,_) = bj;
			
			//r_lambda 
			r_lambda(num_of_store,_) = lambda;
			// r_Z:
			r_Z(num_of_store,_) = Z;
			//r_mu:
			r_mu[num_of_store] = mu;
			//r_sigma2:
			r_sigma2[num_of_store] = sigma2;
			//r_pi
			for(num_pi = 0; num_pi<4;num_pi++){
				r_pi(num_of_store,num_pi) = pi[num_pi];
			}
			//r_m:
			r_m(num_of_store,_) = m;
			
			num_of_store++;

		}
	}

	gsl_permutation_free (permut_p);
	//cout << "finished mcmc loop!\n";
	return List::create(_["r_beta"] = r_beta, _["r_bj"] = r_bj, _["r_lambda"] = r_lambda,
	 _["r_Z"] = r_Z, _["r_mu"] = r_mu, _["r_sigma2"] = r_sigma2, _["r_pi"] = r_pi,
	 _["r_m"] = r_m);


}	
