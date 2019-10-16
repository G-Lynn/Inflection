#Data file
rm(list = ls())
library(BayesLogit)
library(truncnorm)
library(ggplot2)
library(rmutil) #needed for beta binomial density
library(parallel)
library(MASS)
library(invgamma)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("~/Box/Homeless/Code/EPA_Sampling.cpp")


#N_sims = 10
#steps = seq(from = -.3, to = .3, len = N_sims)
#steps = rep(0,N_sims)
#steps = c(.02)
#simulations 1:10 all have no increase to check reproducibility
#simulations 11:17 correspond to (-.3, -.2, -.1, 0, .1, .2, .3) linear increase from 2011-2017
#simulation 21 corresponds to single step of .02 at 2014
#load("~/Desktop/Research/Homeless/Data/Synthetic.RData")

main<-function(nsim,batch, var){
  N.MC = 10000
  B =  0
  Temp = .3
  
  initialization = batch + nsim #"synthetic" #"linear", "step", 1
  message(paste("initialization = ",initialization,sep=""))
  stem = "~/Box/Homeless/"
  load(paste(stem,"Data/ALL_CoC_Affordability_POV.RData",sep="") )
  #load(paste(stem,"Data/ALL_CoC_ZRI.RData",sep="") )
  source(paste(stem,"Code/MCMC_functions_DPM.R",sep="") )
  
  p_cov = dim(X[[1]])[1] # number of covariates
  begin_Year = 2011
  end_Year = 2017
  Year = begin_Year:end_Year
  t.T = length(Year)
  nPlaces = dim(Population)[1]
  
  
  #parameters for the population dynamic process
  m0 = 0
  C0 = .0001
  sigma2_eta = .0001
  
  #parameters for count dynamic proccess
  #lower_level = .75
  
  lower_level = .6*pct.unsheltered.2010 + .95*(1-pct.unsheltered.2010)
  delta = ( steps[nsim] / (t.T - 1) ) * pct.unsheltered.2010
  
  m = matrix(nrow = nPlaces, ncol = t.T)
  m[,1] = lower_level
  Step = F
  tau = 4
  for(t in 2:t.T){
    for(i in 1:nPlaces){
      m[i,t] = min(1,m[i,t-1] + delta[i])
      if(Step==T & t < tau){
        m[i,t] = m[i,t-1]
      }
      if(Step==T & t>=tau)
        m[i,t] = m[i,1] + (steps[nsim])*pct.unsheltered.2010[i]
    }
  }
  
  #m = c(.6,.6,.6,.7,.7,.7)
  #var = 5e-4 #converges with N.MC = 25000 and B = 15000.  inits 1-40
  #var = 5e-3 #provides much wider uncertainty interval.  inits 41 - 50
  
  a = m*((1-m)*m/var - 1)
  b = (var/m^2)*(a^2/m + a)
  var(rbeta(N.MC,a[15,t.T], b[15,t.T]))
  pi_count = list()
  for(i in 1:nPlaces){
    pi_count[[i]] = matrix(nrow = N.MC, ncol = t.T)
    for(t in 1:t.T) pi_count[[i]][,t] = (rbeta(N.MC,a[i,t], b[i,t]))
  }
  
  #hist(pi_count[[15]][,6])
  
  m.1 = apply(pi_count[[15]],2,mean)
  lw.1 = apply(pi_count[[15]],2,quantile,.005)
  up.1 = apply(pi_count[[15]],2,quantile,.995)
  df = data.frame(Year = begin_Year:end_Year,m.1,lw.1,up.1)
  
  
  H_t_supp = list()
  for(i in 1:nPlaces){
    H_t_supp[[i]] = list()
    for(t in 1:t.T){
      #H_t_supp[[i]][[t]] = Count[i,t]:Population[i,t]
      #lower_bound = quantile(1/rbeta(100000,a[i,t],b[i,t]),.0001)
      #upper_bound = quantile(1/rbeta(100000,a[i,t],b[i,t]),.9999)
      lower_bound = quantile(1/rbeta(100000,a[i,t],b[i,t]),.001)
      upper_bound = quantile(1/rbeta(100000,a[i,t],b[i,t]),.999)
      H_t_supp[[i]][[t]] = floor(lower_bound*Count[i,t]):floor(upper_bound*Count[i,t]) 
    }
  }
  
  beta_binom_density = list()
  for(i in 1:nPlaces){
    beta_binom_density[[i]] = list()
    for(t in 1:t.T){
      a_star = a[i,t]/(a[i,t] + b[i,t])
      b_star = (a[i,t] + b[i,t])
      beta_binom_density[[i]][[t]] = dbetabinom(Count[i,t],H_t_supp[[i]][[t]], m = a_star, s = b_star, log = TRUE)
    }
  }
  
  
  range(H_t_supp[[1]][[6]])
  #parameters for the homeless dynamic process
  homeless.2010 = rep(NA,nPlaces)
  for(i in 1:nPlaces){
    homeless.2010[i] = floor(  mean( count.2010[i]/pi_count[[i]][,1] ) )
  }
  
  
  a_psi = 3 #3 #101 #2#11 #3 #10.1
  b_psi =  .1 #.2  #10 #.1 #.01 #1
  sigma2_psi = rinvgamma(nPlaces, a_psi, b_psi) #initialize sigma2_psi
  mean(sigma2_psi)
  hist(sigma2_psi, breaks = 50)
  
  #sigma2_psi = rep(.1, nPlaces)
  
  p.hat = homeless.2010/pop.2010
  psi_2010 = log(p.hat/(1-p.hat))
  
  
  ZRI.2010 = rep(NA, nPlaces)
  Poverty.2010 = rep(NA,nPlaces)
  for(i in 1:nPlaces){
    ZRI.2010[i] = X[[i]][2,1]
    Poverty.2010[i] = X[[i]][3,1]
  }
  
  #for affordability and poverty
  base_measure = rnorm(nPlaces*10,-8.28,sqrt(.4)) + rnorm(nPlaces*10, .061, sqrt(.0002))*ZRI.2010+ rnorm(nPlaces*10, .061, sqrt(.0002))*Poverty.2010 + rnorm(nPlaces*10, 0, sqrt(rinvgamma(nPlaces*10, a_psi, b_psi)))
  
  
  
  #for ZRI only
  #base_measure = rnorm(nPlaces*10,-8.28, sqrt(.4)) + rnorm(nPlaces*10, 0.00157, sqrt(5e-8))*rep(ZRI.2010,10) + rnorm(nPlaces*10, 0, sqrt(rinvgamma(nPlaces*10, a_psi, b_psi)))
  
  
  
  hist(rep(psi_2010,times = 10), col = "black")
  hist(base_measure , col = "red", add = T)
  
  mean(psi_2010)
  mean(base_measure)
  
  var(psi_2010)
  var(base_measure)
  
  #parameters for phi:  this will need to change!  N(m_phi_bar, .25)
  #sigma2_phi = .055 #base measure for phi (same marginal prior for phi_i in AOAS)
  #m_phi = .94
  
  r.hat = sum(count.2010) / (.8*sum(pop.2010))
  log(r.hat / (1-r.hat))
  m_phi = c(-8.28, .061, .061 )
  #m_phi = c(-8.28, 0.00157)
  sigma2_phi = diag(p_cov) #.1*diag(p_cov)
  
  #ZRI only settings
  #sigma2_phi[1,1] = .4   
  #sigma2_phi[2,2] = 5e-8 
  
  #Affordability settings
  sigma2_phi[1,1] =  .4
  sigma2_phi[2,2] = .0002
  sigma2_phi[3,3] = .0002
  
  
  alpha_a = 1 #hyper parameters of gamma prior for alpha
  alpha_b = 1
  
  #parameters for nu
  sigma2_nu_i = .01
  m_nu_bar = 0
  sigma2_nu_bar = .005
  lambda.tilde = 2*pop.2010
  
  ###################################################################################
  #initialize
  #load(file = paste(stem, "Data/DPM/Posterior_Summaries_Initiliazation.RData",sep=""))
  #psi_star = matrix(nrow = nPlaces, ncol = t.T)
  
  
  load(file = "~/Box/Homeless/Data/DPM/IS_Initialization.RData")
  #eta = matrix(0,nrow = nPlaces, ncol = t.T) #dynamics of population process
  eta = Eta_Init
  
  p.hat = (Count/m)/Population
  #psi = log(p.hat/(1-p.hat)) #dynamics of homeless process
  psi = Psi_Init
  
  lambda = lambda.tilde*(1/(1+exp(-eta)))
  p = 1/(1+exp(-psi))
  
  ptm = proc.time()
  H_t = H_t_update = H_t_step(nPlaces, t.T, beta_binom_density, H_t_supp, Population, Count, psi)
  proc.time() - ptm
  H_t = floor(Count/.6)
  
  #initialize DP Mixture
  #J_init = 10 #start by initializing
  J_init = J_Init

  #z = sample(size = nPlaces, 1:J_init, replace = T)
  z = Z_Init
  nZ = length(table(z))
  #phi = t( mvrnorm(J_init,m_phi,sigma2_phi) )
  phi = t(Phi_Init)
  #alpha = rgamma(1,shape = alpha_a, rate = alpha_b)
  alpha = ALPHA_Init
  g = rbeta(1,alpha + 1, nPlaces) #auxiliary parameter for inferring alpha
  #initialize beta
  psi_star = matrix(nrow = nPlaces, ncol = t.T)
  for(t in 1:t.T) psi_star[,t] = psi[,t] - sapply(1:nPlaces, function(i) sum( X[[i]][,t] * phi[, z[i] ] )   ) 
  #beta = beta_step(t.T, nPlaces, psi_star, sigma2_psi)
  beta = Beta_Init
  
  #nu = rep(0,nPlaces)
  nu = NU_Init
  #N_star = sapply(1:t.T, function(i) rpois(nPlaces, lambda.tilde) ) 
  N_star = N_star_step(t.T, nPlaces, lambda.tilde, eta, Population)
  #Omega = sapply(1:t.T, function(t) rpg(nPlaces, as.numeric(N_star[,t]), as.numeric(eta[,t])) )
  Omega = Omega_step(t.T, nPlaces, N_star, eta)
  #Zeta = sapply(1:t.T, function(t) rpg(nPlaces, as.numeric(Population[,t]), as.numeric(psi[,t])) )
  Zeta = Zeta_step(t.T, nPlaces, Population, psi)
  nu_bar = nu_bar_step(nPlaces, m_nu_bar, sigma2_nu_bar, sigma2_nu_i, nu)
  
  sigma2_psi = SIGMA2_PSI_Init
  
  
  #MCMC specifics and storage
  HOMELESS = ETA = N_STAR = N_PRED = H_PRED = C_PRED = PSI = BETA = list()
  
  Z = SIGMA2_PSI = matrix(nrow = N.MC, ncol = nPlaces)
  ALPHA = NZ = rep(NA, N.MC)
  
  PHI = list()
  for(j in 1:nPlaces) PHI[[j]] = matrix(nrow = N.MC, ncol = p_cov)  #initialize to 387, allowing a cluster for each metro.  this is the max.    
  
  for(i in 1:nPlaces){
    HOMELESS[[i]] = ETA[[i]] = PSI[[i]] = N_STAR[[i]] = N_PRED[[i]] = C_PRED[[i]] = BETA[[i]] = matrix(nrow = N.MC, ncol = t.T)
  }
  
  NU = matrix(nrow = N.MC, ncol = nPlaces)
  NU_BAR = rep(NA,N.MC)
  PERMUTATION = matrix(nrow = N.MC, ncol = nPlaces)
  MH_RATIO = rep(NA,N.MC)
  
  #order Z and Phi before MCMC starts
  re_index = order(phi[2,], decreasing = FALSE)
  phi = phi[,re_index]
  z_ordered = rep(NA,nPlaces)
  for(i in 1:nZ){
    z_ordered[z==re_index[i]] = i
  }
  z = z_ordered
  
  ##EPA specifics
  Temp = .35
  lambda = matrix(nrow = nPlaces, ncol = nPlaces)
  
  for(i in 1:nPlaces){
    for(j in 1:nPlaces){
      lambda[i,j] = exp(-Temp*( sqrt( ( X[[i]][2,7] - X[[j]][2,7])^2 + ( X[[i]][3,7] - X[[j]][3,7])^2 ) ) )
    }
  }
  
  
  #Beginning of MCMC
  CPP = TRUE
  permutation = sample(1:nPlaces, size = nPlaces, replace = FALSE)
  lambda_permuted = lambda[permutation, permutation]
  Z_permuted = z[permutation]
  
  ptm = proc.time()
  for(m in -(B):N.MC){

    #Parameters estimated from Census population
    N_star = N_star_step(t.T, nPlaces, lambda.tilde,eta, Population)
    eta = eta_step(t.T, nPlaces, m0, C0, sigma2_eta, Population, N_star, nu, Omega)
    Omega = Omega_step(t.T, nPlaces, N_star, eta)
    nu = nu_step(nPlaces,C0,sigma2_eta,t.T,sigma2_nu_i, eta, nu_bar)
    nu_bar = nu_bar_step(nPlaces, m_nu_bar, sigma2_nu_bar, sigma2_nu_i, nu)
    
    #Noisy population and thinning for homeless population.
    Noisy_Pop = t( sapply(1:nPlaces, function(i) rpois(t.T,lambda.tilde[i]/(1+exp(-eta[i,])) ) ) )
    H_t = H_t_step(nPlaces, t.T, beta_binom_density, H_t_supp, Noisy_Pop, Count, psi)
    
    
    #update log odds of homelessness  
    psi = psi_step(t.T, nPlaces, sigma2_psi, H_t, Noisy_Pop, X, phi, z, Zeta, beta)
    Zeta = Zeta_step(t.T, nPlaces, Noisy_Pop, psi)
    
    #update permutation
    perm_K = 90
    perm_index = sample(1:nPlaces, perm_K, replace = F)
    perm_subset = permutation[perm_index]
    
    #now shuffle these.  
    perm_index_subset = sample(1:perm_K, perm_K, replace = F)
    perm_subset_shuffled = perm_subset[perm_index_subset]
    
    #propose new permutation
    permutation_proposed = permutation
    permutation_proposed[perm_index] = perm_subset_shuffled
    
    Z_perm_proposed = z[permutation_proposed]
    lambda_perm_proposed = lambda[permutation_proposed, permutation_proposed]
    
    #evaluate the Metropolis Ratio
    MH_ratio = min(1, exp(log_pi_EPA(nPlaces, alpha, lambda_perm_proposed, Z_perm_proposed) - log_pi_EPA(nPlaces, alpha, lambda_permuted, Z_permuted) ) )
    accept = sample(0:1, 1, prob = c(1-MH_ratio, MH_ratio))
    if(accept == 1){
      permutation = permutation_proposed
      lambda_permuted = lambda_perm_proposed
      Z_permuted = Z_perm_proposed
    }
  
    MH_RATIO[m] = MH_ratio
    
    
    ##Now sample from EPA partition
    
    if(CPP){
      retx = EPA_step(permutation, nPlaces, alpha, lambda_permuted, z, t.T, psi, beta, sigma2_psi, X, phi)
      #z = as.numeric( retx[[1]] )
      #phi = retx[[2]]
    }else{
      
      for(ii in 1:nPlaces){
        i = permutation[ii] 
        Z_current = z[i]
        Z_tab = table(z)
        Z_count = Z_tab[names(Z_tab)==Z_current]
      
        nZ = length(Z_tab)
        l.pi = rep(NA,nZ)
        Z_permuted = z[permutation]
        for(c in 1:nZ){
          n_ex_i_c = Z_tab[c] - as.numeric(c == z[i])  
          Z_proposed = Z_permuted
          Z_proposed[ii] = c
          l.pi[c] = log_pi_EPA(nPlaces, alpha, lambda_permuted, Z_proposed) + L(t.T, psi[i,], beta[i,], c, sigma2_psi[i], X[[i]], phi)
          #l.pi[c] = log_pi_EPA_cpp(nPlaces, alpha, lambda_permuted, Z_proposed)
          #l.pi[c] = L_cpp(t.T, psi[i,], beta[i,], c-1, sigma2_psi[i], X[[i]], phi)
        }
      
        #now for adding a state
        Z_proposed = Z_permuted
        Z_proposed[ii] = nZ+1
        phi_prop = mvrnorm(1,m_phi,sigma2_phi)
        l.pi_p1 =  log_pi_EPA(nPlaces, alpha, lambda_permuted, Z_proposed) + L(t.T, psi[i,], beta[i,], nZ+1, sigma2_psi[i], X[[i]], cbind(phi, phi_prop ) )
        #l.pi_p1 =  log_pi_EPA_cpp(nPlaces, alpha, lambda_permuted, Z_proposed)
        #l.pi_p1 =  L_cpp(t.T, psi[i,], beta[i,], nZ, sigma2_psi[i], X[[i]], cbind(phi, phi_prop ) )
        
        l.pi = c(l.pi, l.pi_p1) #create the probability of a new class.

        #sum log exponential
        s = max(l.pi)
        l.pi_z = l.pi - (s + log(sum( exp(l.pi - s) ) ) )
        pi_z = exp(l.pi_z)
        #print(pi_z)
        z[i] = sample(1, x = 1:(nZ+1), prob = pi_z)
        if(z[i]==(nZ+1)){phi = cbind(phi,phi_prop)}
      
        #now remove any unused states.  This needs to be fixed.  
        if( Z_count==1 & z[i]!=Z_current ){
         z[z>Z_current] = z[z>Z_current] - 1
         phi = phi[,-Z_current]
        }
      
      
     }#end of loop for i
    }
    
    
    Z_tab = table(z)
    nZ = length(Z_tab)
    
    #update alpha given nZ and a latent g.
    g = rbeta(1, alpha + 1, nPlaces)
    p.zeta = ( (alpha_a + nZ - 1)/(nPlaces*(alpha_b - log(g))) ) / (1 + (alpha_a + nZ - 1)/(nPlaces*(alpha_b - log(g) ) ) ) 
    zeta = sample(size = 1, x = 0:1, prob = c(1-p.zeta, p.zeta) )
    alpha = zeta*rgamma(1, shape = alpha_a + nZ, rate = alpha_b - log(g) ) + (1-zeta)*rgamma(1,shape = alpha_a + nZ - 1, rate = alpha_b - log(g))
    
    #now infer phi
    for(j in 1:nZ) phi[,j] = phi_c_update(t.T, c=j, z, sigma2_psi, sigma2_phi, m_phi, X, psi, beta)
    
    
    #re-order clusterings based on the slope. identifies cluster with slope of regression line 
    re_index = order(phi[2,], decreasing = FALSE)
    phi = phi[,re_index]
    colnames(phi) = 1:nZ
    z_ordered = rep(NA,nPlaces)
    for(i in 1:nZ){
      z_ordered[z==re_index[i]] = i
    }
    z = z_ordered
    
    
    #now update beta
    psi_star[,1] = psi[,1] - sapply(1:nPlaces, function(i) sum( X[[i]][,1] * phi[, z[i] ] )   ) 
    for(t in 2:t.T) psi_star[,t] = psi[,t] - sapply(1:nPlaces, function(i) sum( X[[i]][,t] * phi[, z[i] ] )   ) 
    beta = beta_step(t.T, nPlaces, psi_star, sigma2_psi)
    
    #now update sigma2_psi 
    sigma2_psi = sigma2_psi_step(t.T, nPlaces, a_psi, b_psi, psi, beta, phi, z, X)
    
    if(m>0){
      for(i in 1:nPlaces){
        HOMELESS[[i]][m,] = H_t[i,]
        BETA[[i]][m,] = beta[i,] 
        ETA[[i]][m,] = eta[i,]
        PSI[[i]][m,] = psi[i,]
        N_PRED[[i]][m,] = Noisy_Pop[i,]
        pi_m = rbeta(t.T,a[i,], b[i,])
        C_PRED[[i]][m,] = rbinom(t.T, H_t[i,], pi_m)
      }
      
      for(j in 1:nZ){
        PHI[[j]][m,] = phi[,j]
      }
      
      
      Z[m,] = z
      SIGMA2_PSI[m,] = sigma2_psi
      NU[m,] = nu
      NU_BAR[m] = nu_bar
      ALPHA[m] = alpha
      NZ[m] = nZ
      PERMUTATION[m,] = permutation
    }
    #message(paste("Clusters: ",nZ, sep =""))
    if(m%%1000==0){
      message(paste("Iteration: ",m,sep=""))
      message(paste("Clusters: ",nZ, sep =""))
    }
    
  }
  #End of MCMC
  proc.time() - ptm
  #save(file = paste(stem, "Data/DPM/Posterior_Summaries_Initiliazation.RData",sep=""), H_t, beta, eta, psi, Noisy_Pop, sigma2_psi, phi, z, nu, nu_bar, alpha, nZ, Omega, Zeta, N_star)
  save(file = paste(stem,"Data/DPM/Posterior_Summaries_", initialization, ".RData", sep=""), SIGMA2_PSI, NZ, ALPHA, HOMELESS, ETA, PSI, N_PRED, C_PRED, PHI, NU, Z, NU_BAR,pi_count, BETA)
}

#for now only run the truncated / non-truncated simulations.  
for(i in c(12)){
  
  if(i==12){ #final analysis 2018: EPA
    var = 1.5e-3
    N_sims = 4
    steps = rep(0,N_sims)
    batch = 120
    truncate = FALSE 
  }
  
  
  mclapply(1:N_sims, main, batch, var, mc.cores = 4)
}
