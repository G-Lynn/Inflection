#PG_Gibbs_Functions
N_star_step<-function(t.T, nPlaces, lambda.tilde,eta, Population){
  theta = 1/(1+exp(-eta))
  return( sapply(1:t.T, function(t) rpois(nPlaces,lambda.tilde*(1-theta[,t]) ) ) + Population )
}

eta_step<-function(t.T, nPlaces, m0, C0, sigma2_eta, Population, N_star, nu, Omega){
  m = m.tilde = a = matrix(nrow=nPlaces,ncol=t.T)
  C = C.tilde = R = matrix(nrow=nPlaces,ncol=t.T)
  kappa = Population - N_star/2
  
  for(t in 1:t.T){
    if(t == 1){
      R[,t] = C0 + sigma2_eta
      a[,t] = m0 + nu
      
      C[,t] = (Omega[,t] + 1/R[,t] )^(-1)
      m[,t] = C[,t]*(kappa[,t] + a[,t]/R[,t] )
    }else{
      R[,t] = C[,t-1] + sigma2_eta
      a[,t] = m[,t-1] + nu
      
      C[,t] = (Omega[,t] + 1/R[,t] )^(-1)
      m[,t] = C[,t]*(kappa[,t] + a[,t]/R[,t] )
    }
  }
  
  eta_ret = matrix(nrow=nPlaces,ncol=t.T)
  for(t in t.T:1){
    if(t == t.T){
      eta_ret[,t] = rnorm(nPlaces, m[,t], sqrt(C[,t]) )
    }else{
      C.tilde[,t] = (1/sigma2_eta + 1/C[,t])^(-1)
      m.tilde[,t] = C.tilde[,t]*(m[,t]/C[,t] + (eta_ret[,t+1]-nu)/sigma2_eta )
      eta_ret[,t] = rnorm(nPlaces,m.tilde[,t], sqrt(C.tilde[,t]) )
    }
  }
  
  return(eta_ret)
}

Omega_step<-function(t.T, nPlaces, N_star, eta){
  return( sapply(1:t.T, function(t) rpg(nPlaces, as.numeric(N_star[,t]), as.numeric(eta[,t]) ) ) )
}

nu_step<-function(nPlaces,C0,sigma2_eta,t.T,sigma2_nu_i, eta, nu_bar){
  Delta.eta = eta[,2:t.T] - eta[,1:(t.T-1)]
  sigma2_tilde_nu_i = (1/(C0 + sigma2_eta) + (t.T - 1)/sigma2_eta + 1/sigma2_nu_i)^(-1)
  m_tilde_nu_i = sigma2_tilde_nu_i*(eta[,1]/(C0 + sigma2_eta) + (1/sigma2_eta)*apply(Delta.eta,1,sum) + nu_bar/sigma2_nu_i)
  return( rnorm( nPlaces,mean = m_tilde_nu_i, sd = sqrt(sigma2_tilde_nu_i) ) )
}

nu_bar_step<-function(nPlaces, m_nu_bar, sigma2_nu_bar, sigma2_nu_i, nu){
  Sigma = (nPlaces/sigma2_nu_i + 1/sigma2_nu_bar)^(-1)
  M = Sigma*(sum(nu)/sigma2_nu_i + m_nu_bar/sigma2_nu_bar)
  return(rnorm(1,M,sqrt(Sigma)))
}

psi_step <-function(t.T, nPlaces, sigma2_psi, H_t, Population, X, phi, z, Zeta, beta){
  kappa = H_t - Population/2
  
  f =  matrix(nrow = nPlaces, ncol = t.T)
  Sigma =  matrix(nrow = nPlaces, ncol = t.T)
  psi_ret = matrix(nrow = nPlaces, ncol = t.T)
  
  for(i in 1:nPlaces){
  for(t in 1:t.T){
      Sigma[i,t] = (Zeta[i,t] + 1/sigma2_psi[i] )^(-1)
      f[i,t] = Sigma[i,t]*(kappa[i,t] + (beta[i,t] + sum( phi[,z[i] ] * X[[i]][,t] ) )/ sigma2_psi[i] )
      psi_ret[i,t] = rnorm(1,f[i,t], sqrt(Sigma[i,t]) )
    }
  }
  return(psi_ret)
}

Zeta_step <-function(t.T, nPlaces, Population, psi){
  return( sapply(1:t.T, function(t) rpg(nPlaces, as.numeric(Population[,t]), as.numeric(psi[,t]) ) ) )
}


H_t_step<-function(nPlaces, t.T, beta_binom_density, H_t_supp, Population, Count, psi){
  p = 1/(1+exp(-psi) )
  H_t = matrix(nrow = nPlaces, ncol = t.T)
  for(i in 1:nPlaces){
    for(t in 1:t.T){
      #a_star = a[i,t]/(a[i,t] + b[i,t])
      #b_star = (a[i,t] + b[i,t])
      #weight = dbetabinom(Count[i,t],H_t_supp[[i]][[t]], m = a_star, s = b_star) * dbinom(H_t_supp[[i]][[t]], Population[i,t], p[i,t])
      l.weight = beta_binom_density[[i]][[t]] + dbinom(H_t_supp[[i]][[t]], Population[i,t], p[i,t], log = TRUE)
      s = max(l.weight)
      l.prob = l.weight - (s + log(sum( exp(l.weight - s) ) ) )
      prob = exp(l.prob)
      H_t[i,t] = sample(x = H_t_supp[[i]][[t]], size = 1, prob = prob)
    }
  }
  return(H_t)
}

#Integral F(beta)dG_0 from Neal's algo 2
LG_integral <-function(t.T, sigma2_psi_i, m_phi, sigma2_phi, psi_i, beta_i, X_i){
  
  A = a = s = 0
  for(t in 1:t.T){
    A = A + X_i[,t]%*%t(X_i[,t]) / sigma2_psi_i
    a = a + X_i[,t]*(psi_i[t] - beta_i[t]) / sigma2_psi_i
    s = s + (psi_i[t] - beta_i[t])^2 / sigma2_psi_i    
  }
  
  A = A + solve(sigma2_phi)
  a = a + solve(sigma2_phi)%*%m_phi
  s = s + t(m_phi)%*%solve(sigma2_phi)%*%m_phi
  
  const = (2*pi*sigma2_psi_i)^(-(t.T)/2) * det(2*pi*sigma2_phi)^(-1/2) * det(2*pi*solve(A))^(1/2)
  #dens = as.numeric( const * exp(-0.5*s) * exp(0.5*t(a)%*%solve(A)%*%a) )
  log.dens = as.numeric( log(const) -.5*s + .5*t(a)%*%solve(A)%*%a)
  return(log.dens)
}

#evaluate the normal likelihood for psi_{i,1:T}
L <-function(t.T, psi_i, beta_i, c, sigma2_psi_i, X_i, phi ){
  Lik = rep(NA,t.T)
  for(t in 1:t.T){
    Lik[t] = dnorm(psi_i[t], mean = beta_i[t] + as.numeric( t(X_i[,t])%*%phi[,c]), sd = sqrt(sigma2_psi_i), log = TRUE  )
    #message(Lik[t])
  }
  
  return(sum(Lik))
}

#updated atoms from Neal's algo 2
phi_c_update <-function(t.T, c, C, sigma2_psi, sigma2_phi, m_phi, X, psi, beta){
  index = which(C == c)
  p = dim(X[[1]])[1]
  nC = length(index)
  A = 0
  a = rep(0,p)
  for(i in 1:nC){
    ii = index[i]
    B = b = 0
    for(t in 1:t.T){
      B = B + X[[ii]][,t]%*%t(X[[ii]][,t])/(sigma2_psi[ii])
      b = b + X[[ii]][,t]*(psi[ii,t] - beta[ii,t]) / (sigma2_psi[ii])
    }
    A = A + B 
    a = a + b
  }
  
  A = A + solve(sigma2_phi)
  a = a + solve(sigma2_phi)%*%m_phi
  
  A_star = solve(A)
  a_star = A_star%*%a
  phi_c = mvrnorm(1,a_star,A_star)
  return(phi_c)
}

#beta_step<-function(t.T, nPlaces, psi_star, sigma2_psi){
#  m0_beta = 0
#  C0_beta = .1
#  W = .005
  
#  m = C = R = A = matrix(nrow = nPlaces, ncol = t.T)
#  R[,1] = C0_beta + W
#  A[,1] = R[,1] / (R[,1] + sigma2_psi)
#  m[,1] = m0_beta + A[,1] * (psi_star[,1] - m0_beta)
#  C[,1] = R[,1] - R[,1]^2/(R[,1] + sigma2_psi)
  
#  for(t in 2:t.T){
#    R[,t] = C[,t-1] + W
#    A[,t] = R[,t] / (R[,t] + sigma2_psi)
#    m[,t] = m[,t-1] + A[,t] * (psi_star[,t] - m[,t-1])
#    C[,t] = R[,t] - R[,t]^2 / (R[,t] + sigma2_psi)
#  }
  
#  beta = matrix(nrow = nPlaces, ncol = t.T)
  
#  beta[,t.T] = rnorm(nPlaces, m[,t.T], sqrt(C[,t.T]) )
  
#  for(t in (t.T-1):1){
#    C.tilde = (1/C[,t] + 1/W)^(-1)
#    m.tilde = C.tilde*(m[,t]/C[,t] + beta[,t+1]/W)
#    beta[,t] = rnorm(nPlaces,m.tilde,sqrt(C.tilde))
#  }
  
#  return(beta)
#}

beta_step<-function(t.T, nPlaces, psi_star, sigma2_psi){
  m0_beta = matrix(c(0,0), ncol = 1)
  C0_beta = diag(c(.1,.000001) )
  W = diag(c(.005,.000001))
  W_inv = solve(W)
  G = matrix(c(1,1,0,1), nrow = 2, byrow = T)
  G_t_W_inv = t(G)%*%W_inv
  retx = matrix(nrow = nPlaces, ncol = t.T)
  F_t = matrix(c(1,0), nrow = 2)

  for(i in 1:nPlaces){
    m = a = matrix(nrow = 2, ncol = t.T)
    q = f = e = rep(NA,t.T)
    C = R = A = list()
    
    
    a[,1]  = G%*%m0_beta
    f[1] = t(F_t)%*%a[,1]
    
    R[[1]] = C0_beta + W
    q[1] = t(F_t)%*%R[[1]]%*%F_t + sigma2_psi[i]
    A[[1]] = R[[1]]%*%F_t/q[1]
    
    e[1] = psi_star[i,1] - f[1]
    m[,1] = a[,1] + A[[1]]%*%e[1]
    C[[1]] = R[[1]] - A[[1]]%*%q[1]%*%t(A[[1]])

    for(t in 2:t.T){
      a[,t]  = G%*%m[,t-1]
      f[t] = t(F_t)%*%a[,t]
      
      R[[t]] = C[[t-1]] + W
      q[t] = t(F_t)%*%R[[t]]%*%F_t + sigma2_psi[i]
      A[[t]] = R[[t]]%*%F_t / q[t]
      
      e[t] = psi_star[i,t] - f[t]
      m[,t] = a[,t] + A[[t]]%*%e[t]
      C[[t]] = R[[t]] - A[[t]]%*%q[t]%*%t(A[[t]])
    }

    beta = matrix(nrow = 2, ncol = t.T)

    beta[,t.T] = mvrnorm(1, m[,t.T], C[[t.T]]) 
    retx[i,t.T] = beta[1,t.T]
    for(t in (t.T-1):1){
      C_inv = solve(C[[t]])
      C.tilde = solve(C_inv + W_inv)
      m.tilde = C.tilde%*%(C_inv%*%m[,t] + G_t_W_inv%*%beta[,t+1])
      beta[,t] = mvrnorm(1,m.tilde,C.tilde)
      retx[i,t] = beta[1,t]
    }
  }
  return(retx)
}



sigma2_psi_step <-function(t.T, nPlaces, a_psi, b_psi, psi, beta, phi, z, X){
  sigma2_psi_star = rep(NA,nPlaces)
  for(i in 1:nPlaces){
    
    psi_i = psi[i,]
    a_psi_star = (t.T/2 + a_psi)
    b_psi_star = b_psi
    for(t in 1:t.T){
      b_psi_star = b_psi_star + .5*(psi_i[t] - beta[i,t] - sum( phi[,z[i] ] * X[[i]][,t] ) )^2 
    }
    sigma2_psi_star[i] = rinvgamma(1,a_psi_star, b_psi_star)
  }
  
  return(sigma2_psi_star)
  
}


log_pi_EPA <-function(n, alpha, lambda_permuted, Z_permuted){
  p_t = rep(1,n)
  for(t in 2:n){
    #check to see if Z is a new subset
    if(is.element(Z_permuted[t], Z_permuted[1: (t-1)] ) ){
      num = sum( lambda_permuted[t,1:(t-1)][ Z_permuted[1:t-1]==Z_permuted[t] ]  )
      den = sum(lambda_permuted[t, 1:(t-1)] )
      p_t[t] = (t-1)/(alpha + t-1)*(num/den)
    }else{
      p_t[t] = alpha / (alpha + t - 1)
    }
  }
  return( sum( log(p_t) ) )
}
