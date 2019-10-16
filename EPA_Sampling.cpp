#include <RcppArmadilloExtensions/sample.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
float log_pi_EPA_cpp(int n, float alpha, arma::mat lambda_permuted, arma::vec Z_permuted){
  arma::vec p_t(n);
  p_t.ones();
  arma::rowvec tmp(n);
  tmp = lambda_permuted.row(1);
  //Rcpp::Rcout << "Inside first function" << std::endl;
  
  
  for(int t=1; t<n; t++){
    arma::uvec ids = find(Z_permuted.subvec(0,t-1)==Z_permuted(t) );
    if(ids.n_elem>0){
      //Rcpp::Rcout << "Index: "<< t+1 << std::endl;
      float num = 0;
      float den = 0;
      arma::rowvec lambda_row(n);
      lambda_row = lambda_permuted.row(t);
      num = sum(lambda_row.elem(ids));
      //Rcpp::Rcout << "the numerator is: " << num << std::endl;
      den = sum(lambda_row.subvec(0,t-1) );
      //Rcpp::Rcout << "the denominator is: " << den << std::endl;
      p_t(t) = t/(alpha + t)*(num/den);
    }else{
      p_t(t) = alpha / (alpha + t);
    }
  }
  
  return sum(log(p_t));
}

// [[Rcpp::export]]
float L_cpp(int t_T, arma::rowvec psi_i, arma::rowvec beta_i, int c, float sigma2_psi_i, arma::mat X_i, arma::mat phi ){
  arma::vec Lik(t_T);
  for(int t=0; t<t_T; t++){
    //double mu = beta_i[t] + X_i[,t].t*phi[,c];
    arma::mat mu = beta_i[t] + X_i.col(t).t()*phi.col(c);
    Lik[t] = -0.5*log(2*3.14159265359*sigma2_psi_i) - 0.5*(1/sigma2_psi_i)*pow( (psi_i[t] - mu(0,0)), 2.0);
    
  }
    return sum(Lik);
}

// [[Rcpp::export]]
Rcpp::List EPA_step(arma::uvec permutation, int n, float alpha, arma::mat lambda_permuted, arma::vec z, int t_T, arma::mat psi, arma::mat beta, arma::vec sigma2_psi, Rcpp::List X, arma::mat phi){
  //Rcpp::Rcout << "Inside wrapper function" << std::endl;
  Rcpp::List Retx;
  
  
  //Sample from the EPA partition for each 
  for(int ii=0; ii<n; ii++){
    int i = permutation[ii]-1;
    //Rcpp::Rcout << "i = " << i << std::endl;
    arma::mat X_i = X[i];
    //Rcpp::Rcout << "X_i = " << X_i << std::endl;
    
    arma::rowvec psi_i = psi.row(i);
    //Rcpp::Rcout << "psi_i = " << psi_i << std::endl;
    
    
    arma::rowvec beta_i = beta.row(i);
    //Rcpp::Rcout << "beta_i = " << beta_i << std::endl;
    
    
    int Z_current = z[i];
    //Rcpp::Rcout << "Z[i] = " << Z_current << std::endl;
    
    arma::vec z_unique = unique(z);
    int nZ = z_unique.n_elem;
    arma::uvec Same_indices = find(z==Z_current);
    int Z_count = Same_indices.n_elem;
    
    //Rcpp::Rcout << "nZ = " << nZ << "; and Z_count = " << Z_count << std::endl;
    
    
    arma::vec l_pi(nZ); l_pi.zeros();
    arma::vec z_permuted(n); z_permuted = z.elem(permutation-1);
    
    arma::vec z_proposed = z_permuted;

    for(int c = 0; c<nZ; c++){
      z_proposed[ii] = c+1;
      l_pi[c] = log_pi_EPA_cpp(n, alpha, lambda_permuted, z_proposed) + L_cpp(t_T, psi_i, beta_i, c, sigma2_psi(i), X_i, phi);
      //l_pi[c] = log_pi_EPA_cpp(n, alpha, lambda_permuted, z_proposed);
      //l_pi[c] = L_cpp(t_T, psi_i, beta_i, c, sigma2_psi(i), X_i, phi);
    }
    
    //Now for adding a state:
    z_proposed[ii] = nZ+1;
    
    //generate the proposed phi from the prior
    NumericVector m_phi = NumericVector::create(-8.28, 0.061, 0.061);
    //Rcpp::Rcout << "mean vector: " << m_phi << std::endl;
    
    NumericVector sigma_phi = NumericVector::create(.6324555, .01414214, .01414214);
    //Rcpp::Rcout << "sd vector: " << sigma_phi << std::endl;
    
    NumericVector Z_norm = rnorm(3, 0, 1);
    arma::vec phi_prop(3); phi_prop.zeros();
    
    phi_prop(0) = m_phi[0] + sigma_phi[0]*Z_norm[0];
    phi_prop(1) = m_phi[1] + sigma_phi[1]*Z_norm[1];
    phi_prop(2) = m_phi[2] + sigma_phi[2]*Z_norm[2];
    
    //Rcpp::Rcout << "phi_1: " << phi_prop(0) << std::endl;
    //Rcpp::Rcout << "phi_2: " << phi_prop(1) << std::endl;
    //Rcpp::Rcout << "phi_3: " << phi_prop(2) << std::endl;
    
    arma::mat phi_tmp = phi;
    phi_tmp.insert_cols(nZ,phi_prop);
    //Rcpp::Rcout << "phi tmp" << phi_tmp << std::endl;
    
    arma::vec l_pi_p1(1); 
    l_pi_p1[0] = log_pi_EPA_cpp(n, alpha, lambda_permuted, z_proposed) + L_cpp(t_T, psi_i, beta_i, nZ, sigma2_psi(i), X_i, phi_tmp);
    //l_pi_p1[0] = log_pi_EPA_cpp(n, alpha, lambda_permuted, z_proposed);
    //l_pi_p1[0] = L_cpp(t_T, psi_i, beta_i, nZ, sigma2_psi(i), X_i, phi_tmp);
    l_pi.insert_rows(nZ,l_pi_p1);
    //Rcpp::Rcout << "l_pi = " << l_pi << std::endl;
    
    float s = l_pi.max();
    l_pi = l_pi - (s + log(sum(arma::exp(l_pi - s))));
    arma::vec pi_z = arma::exp(l_pi);
    //Rcpp::Rcout << "pi_z" << pi_z << std::endl;
    
    //Now sample the next value of z[i]
    NumericVector prob(nZ+1);
    NumericVector x(nZ+1);
    for(int c = 0; c <=nZ; c++){
      prob[c] = pi_z[c];
      x[c] = c+1; 
    }
    
    //Rcpp::Rcout << "x = " << x << std::endl;
    //Rcpp::Rcout << "Prob: " << prob << std::endl;
    
    int size = 1;
    NumericVector z_tmp = RcppArmadillo::sample(x,size, false, prob);
    //Rcpp::Rcout << "Current Z: " << z[i] << std::endl;
    z[i] = z_tmp[0];
    //Rcpp::Rcout << "Sampled Z: " << z[i] << std::endl;
    
    if(z[i] == (nZ+1)){
      phi = phi_tmp;
    }
    
    //Remove any unused states
    if(Z_count == 1 && z[i]!=Z_current){
      arma::uvec index_star = find(z>Z_current);
      z.elem(index_star) = z.elem(index_star) - 1;
      phi.shed_col(Z_current-1);
    }
    
    
  }
  
  return Rcpp::List::create(Rcpp::Named("z") = z,
                            Rcpp::Named("phi") = phi);
}

