###AFT PDAS###

PDAS = function(ita, d, lambda, K, x, y, n, p){
  #------------------------------------------------------------------------
  #  PDAS algorithm
  #------------------------------------------------------------------------
  # INPUT:
  #   ita            initial value p-dim vector ita for PDAS
  # 	d              initial value p-dim vector d for PDAS  
  #   lambda         threshold to define active set  
  # 	K              The maximum number of iterations for PDAS 
  #       x              design matrix (after data processing)
  # 	y              response vector (after data processing)
  # 	n              sample size 
  # 	p              dimension 
  # OUTPUT:
  # 	ita            estimator value p-dim vector ita for PDAS
  # 	ita            estimator value p-dim vector d for PDAS
  #   HBIC           high dimensional BIC see ref  Wang et al. (2013b)
  #   suppnum        support size of estimator ita
  #------------------------------------------------------------------------
  # Reference:
  #   Wang L, Kim Y, Li R (2013b) Calibrating nonconvex penalized regression in ultra-high dimension. The Annals of Statistics 41(5)
  #------------------------------------------------------------------------
  # Written by:
  # 	Lican Kang (kanglican@whu.edu.cn)  
  # 	Ning Su (suning@whu.edu.cn)
  # This version: Feb 17, 2023   
  #------------------------------------------------------------------------
  for(k in 1:K){
    A = which(abs(ita + d) > lambda)
    A0 = A
    I = which(abs(ita + d) <= lambda)
    ita[I] = 0
    d[A] = 0
    X_A = x[,A,drop = FALSE];
    X_I = x[,I,drop = FALSE];
    ita[A] = solve((t(X_A)%*%X_A),t(X_A)%*%y)
    d[I] = t(X_I)%*%(y - X_A%*%ita[A])/n
    A = which(abs(ita + d) > lambda)
    if(identical(A0, A) | (length(A) == 0)){
      M.lambda = length(A0)
      SSE.lambda = sum((y - x%*%ita)**2)/n
      HBIC = log(SSE.lambda) + log(p)/n*M.lambda*log(n)
      return(list(ita = ita, d = d, HBIC = HBIC, suppnum = M.lambda))
    }
  }
  M.lambda = length(A0)
  SSE.lambda = sum((y - x%*%ita)**2)/n
  HBIC = log(SSE.lambda) + log(p)/n*M.lambda*log(n)
  return(list(ita = ita, d = d, HBIC = HBIC, suppnum = M.lambda))
}


SPDAS = function(K, x.tilde, y.tilde, D, n, p, m, alpha){
  #------------------------------------------------------------------------
  #  SPDAS algorithm
  #------------------------------------------------------------------------
  # INPUT:
  # 	K              the maximum number of iterations for PDAS
  #   x.tilde        design matrix (after data processing)
  # 	y.tilde        response vector (after data processing)
  #   D              pXp diagonal matrix 
  # 	n              sample size 
  # 	p              dimension 
  # 	m              the value that controls how quickly lambda decreases  
  #   alpha          
  # OUTPUT:
  # 	beta           estimator value p-dim vector beta for SPDAS (define by HBIC criterion)
  #------------------------------------------------------------------------
  # Written by:
  # 	Lican Kang (kanglican@whu.edu.cn)  
  # 	Ning Su (suning@whu.edu.cn)
  # This version: Feb 17, 2023   
  #------------------------------------------------------------------------
  
  ita = rep(0, p)
  d = t(x.tilde)%*%(y.tilde)/n
  lambda0 = max(abs(d))**2/2
  lambda = sqrt(2*lambda0*(alpha**c(1:m)))
  HBIC = rep(100, m)
  ita.record = matrix(0, p, m)
  for (i in 1:m){
    pdas.sim = PDAS(ita, d, lambda[i], K, x.tilde, y.tilde, n, p)
    ita = pdas.sim[[1]]
    ita.record[, i] = ita
    d = pdas.sim[[2]]
    HBIC[i] = pdas.sim[[3]]
    if(pdas.sim[[4]] > floor(n/log(p))) break
    
  }
  HBIC.min.index = which(HBIC == min(HBIC))[[1]]
  beta = D%*%ita.record[,HBIC.min.index]
  return(list(beta))
}




