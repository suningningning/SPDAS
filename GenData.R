#------------------------------------------------------------------------
#Generating design matrix x,scenario(I) kind = 1
get.x1 = function(n, p, rho){
  mean = rep(0, p)
  V = matrix(0, p, p)
  for(i in 1:p){
    for(j in 1:p){
      V[i,j] = rho^abs(i - j)
    }
  }
  library(MASS)
  x = mvrnorm(n, mean, V)
  return(x)
}
#------------------------------------------------------------------------
#Generating design matrix x, scenario(II) kind = 2
get.x2 = function(n, p, rho){
  X1 = matrix(rnorm(n*p), nr = n, nc = p)
  x = matrix(0,nr = n,nc = p)
  x[, 1] = X1[, 1]
  x[, p] = X1[, p]
  for(j in 2:(p-1)){
    x[, j] = X1[, j] + rho*(X1[, j+1] + X1[, j-1])
  }
  return(x)
}
#------------------------------------------------------------------------
#Generating beta
get.beta = function(p, T0, R, seednum){
  #------------------------------------------------------------------------
  # Generate data for simulation
  #------------------------------------------------------------------------
  # INPUT:
  # 	p                dimension 
  # 	T0               support size of true beta
  # 	R                ratio of the largest non-zero absolute element to the smallest non-zero absolute element in beta
  # 	seednum          the seed number for repeatability and reproducibility
  #                                          
  # OUTPUT:
  # 	betature          the true value of beta 
  #   suppture          the true support of betature 
  #------------------------------------------------------------------------
  set.seed(seednum)
  b0 = (rbinom(T0, 1, 0.5)*2 - 1)*(R^runif(T0))
  supptrue = sample(p, T0)
  betatrue = rep(0, p)
  betatrue[supptrue] = b0
  return(list(betatrue, supptrue))
}
#------------------------------------------------------------------------
get.data = function(n, p, betatrue, sigma, rho, seednum, kind, cs.ita){
  #------------------------------------------------------------------------
  # Generate data for simulation
  #------------------------------------------------------------------------
  # INPUT:
  # 	n                sample size
  # 	p                dimension 
  #   betatrue         p-dim vector 
  # 	sigma            noise level for observation
  #	  rho              correlation parameter of X 
  # 	seednum          the seed number for repeatability and reproducibility
  #   kind    	       design matrix type: 1 = Classical Gaussian matrix (small p); 2 = Random Gaussian matrix (big p); 
  #   cs.ita           censor time ~ Uniform(0,cs.ita)
  # OUTPUT:
  # 	X                design matrix 
  # 	y                response vector 
  #   c0               censor time
  #------------------------------------------------------------------------
  # fix seed 
  set.seed(seednum)
  #create noise and observation 
  if(kind == 1){
    x = get.x1(n, p, rho)
    y = x%*%betatrue + rnorm(n, mean = 0, sd = sigma)
  }else{
    x = get.x2(n, p, rho)
    y = x%*%betatrue + rnorm(n, mean = 0, sd = sigma)
  }
  c0 = log(runif(n, 0, cs.ita))
  return(list(x, y, c0))
}
#------------------------------------------------------------------------

censor.rate = function(cs.ita,n,p,betatrue,sigma,rho,kind,sim){
  #------------------------------------------------------------------------
  #  simulation to define tuning parameter cs.ita
  #------------------------------------------------------------------------
  # INPUT:
  #   cs.ita           censor time ~ Uniform(0,cs.ita)
  # 	n                sample size
  # 	p                dimension 
  #   betatrue         p-dim vector 
  # 	sigma            noise level for observation
  #   rho              correlation parameter of X 
  #   kind    	       design matrix type: 1 = Classical Gaussian matrix (small p); 2 = Random Gaussian matrix (big p); 
  # OUTPUT:
  # 	c.r              average censor rate
  #------------------------------------------------------------------------
  c.r = 0
  for(i in 1:sim){
    data = get.data(n,p,betatrue,sigma,rho,i*1000,kind,cs.ita)
    c.r = c.r + mean(data[[2]] > data[[3]])
  }
  c.r = c.r/sim
  return(list(c.r))
}
#------------------------------------------------------------------------


get.weight = function(y, c0, n){
  #------------------------------------------------------------------------
  # Generate kaplan-Meier weight
  #------------------------------------------------------------------------
  # INPUT:
  # 	y                response vector 
  #   c0               censor time
  # 	n                sample size
  # OUTPUT:
  # 	y.order          observation with censor increasing order y_(1), y(2), y_(3)...
  # 	delta            censor indicator associate with y.order
  # 	omega            Kaplan-Meier estimator 
  # 	index            increasing order of observation with censor
  #------------------------------------------------------------------------
  # Written by:
  # 	Lican Kang (kanglican@whu.edu.cn)  
  # 	Ning Su (suning@whu.edu.cn)
  # This version: Feb 17, 2023   
  #------------------------------------------------------------------------
  omega = rep(0, n)
  cs.indicator = as.numeric(y <= c0)
  y = apply(cbind(y, c0), 1, min)
  index = order(y,decreasing = F) 
  y.order = y[index]
  delta = cs.indicator[index]
  cs.indicator.order = delta/((n + 1) - c(1:n))
  omega[1] = cs.indicator.order[1]
  for(i in 2:n){
    omega[i] = cs.indicator.order[i]*prod((1 - cs.indicator.order[c(1:i-1)]))
  }
  return(list(y.order, delta, omega, index))
}

xy.process = function(X, y, c0, n){
  #------------------------------------------------------------------------
  # Data process
  #------------------------------------------------------------------------
  # INPUT:
  #    x                covariates
  #    y                observation without censoring
  #    c0               censor time
  #    n                sample size
  # OUTPUT:
  # 	x.tilde          covariate weighted by KM estimator
  # 	y.tilde          observation weighted by KM estimator
  # 	D                pXp matrix
  #------------------------------------------------------------------------
  # Written by:
  # 	Lican Kang (kanglican@whu.edu.cn)  
  # 	Ning Su (suning@whu.edu.cn)
  # This version: Feb 17, 2023  
  #------------------------------------------------------------------------
  weight.result = get.weight(y, c0, n)
  DW = diag(sqrt(weight.result[[3]]))
  y = DW%*%weight.result[[1]]
  X = X[weight.result[[4]],]
  DX = DW%*%X
  X.length = sqrt(apply((DX)**2,2,sum))
  D = diag((sqrt(n)/X.length))
  X = DX%*%D
  return(list(x.tilde = X, y.tilde = y, D = D))
}

library(ncvreg)
library(glmnet)
library(MASS)

