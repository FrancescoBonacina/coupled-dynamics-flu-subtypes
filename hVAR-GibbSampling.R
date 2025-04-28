################################################################################
# 
# Fit a Bayesian Hierarchical VAR process by Gibbs Sampler     
# using auxiliary variables to break variable dependence.
#
# Code based on: 
# Lu et al., Bayesian hierarchical vector autoregressive models for 
# patient-level predictive modeling, PLOS ONE (2018), 
# DOI 10.1371/journal.pone.0208082
#
# The main difference from the code of Lu et al., is that we consider a VAR
# process with the intercept term. That is, we consider:
#    y_t = nu + A_1*y_t-1 + ... + A_l * y_t-l + eps_t
# rather than:
#    y_t = A_1*y_t-1 + ... + A_l * y_t-l + eps_t
#
# Moreover, we changed a few minor details
#
################################################################################

library(MASS)
library(statmod)
library(mvnfast)

# ------------------------------------------------------------------------------
#    function to explore the posterior distribution
# ------------------------------------------------------------------------------
MSGibbs <- function(H,        # list of N matrices, each one of shape (B*l+1) x TT-l. Each matrix stores coordinates of past points of the series
                    HH,       # list of the N matrices HH^t, of shape ((B*l+1) x (TT-l)) x ((T-l) x (B*l+1)) = (B*l+1) x (B*l+1)
                    eta0,     # list of N matrices, each one containing B rows and TT-l columns. Columns contain coord. of points from l+1 to TT
                    w_MLE,    # initial value of matrix w
                    TT,       # a vector containing the length of each trajectory. (For us all lengths are equal)
                    B,        # R in the paper, is the length of vectors constituting the trajectories (for us B=2, trajectories live in a 2-dimensional space)
                    l=2,      # lag of the autoregressive model
                    N,        # nb. of trajectories
                    seed=100, 
                    numIter=2000,  # nb. of samples of each chain
                    thin=1,        # estimates will be computed by considering one sample each 'thin' samples of the chain.
                    warmup=0,      # initial part of the chain to be removed
                                                   #c Lambda has a Wishart distribution:
                                                   #c    Wishart distribution is a generalization to multiple dimensions of the gamma distribution.
                                                   #c    These distributions are of great importance in the estimation of covariance matrices in 
                                                   #c    multivariate statistics. In Bayesian statistics, the Wishart distribution is the conjugate prior of 
                                                   #c    the inverse covariance-matrix of a multivariate-normal random-vector.
                                                   #c A Whisart distrib. has 2 parameters. W(V,n), with V the scale matrix and n the degrees of freedom
                    K_inv = diag(1/(B-1), B),   # inverse of scale matrix of Lambda prior
                    nu = 1,                     # df of Lambda prior
                    mu1 = 1,                    # mean of lambda1sq prior
                    nu1 = 0.001,                # df of lambda1sq prior
                    mu2 = 1,                    # mean of lambda2 prior
                    nu2 = 0.01,                 # df of lambda2 prior
                    k1 = 0.0005,                # shape prior of omega_v
                    theta1 = 200,               # scale prior of omega_v
                    a = 100,                        #c prior for the covariance matrix of vector alpha: alpha ~ N(0,a*I_R²p)
                    initial=NULL
){   
  ### 1) store values
  #c create lists to store parameters as soon as they are sampled
  sim.v_n <- sim.w <- vector(mode="list", length=(numIter-warmup)/thin+1) 
  sim.Lambda <- sim.lambda1sq <- sim.lambda2 <- sim.tausq <- vector(mode="list", length=(numIter-warmup)/thin+1)
  sim.U_v <- sim.alpha <- vector(mode="list", length=(numIter-warmup)/thin+1)
  sim.v_n.star <- vector(mode="list", length=(numIter-warmup)/thin+1)
  sim.omega_v <- vector(mode="list", length=(numIter-warmup)/thin+1)
  
  ### 2) initialize
  set.seed(seed)
                                                                                  #c if no values are assigned for the initial condition, compute them:
  if(is.null(initial)){                                                                 #c Remember that the prior distrib. are:
    sim.w[[1]] <- w <- rnorm(B^2*l+B, 0, 1)                                             #c w|D ~ N(0,D⁻1)
    sim.v_n[[1]] <- v_n <- lapply(vector("list", N), function(x){rnorm(B^2*l+B, 0, 1)}) #c v_n|Omega_v ~ N(0, Omega_v⁻1) 
    sim.Lambda[[1]] <- Lambda <- Posdef(n=B, ev=runif(B, 0, 0.001))                     #c Lambda ~ Wishart(k,nu) --> to generate 
    sim.lambda1sq[[1]] <- lambda1sq <- runif(B^2*l+B, 0.1, 1000)                        #c lambda1 ~ Gamma(mu1, nu1), lambda1² ~ ??
    sim.lambda2[[1]] <- lambda2 <- runif(B^2*l+B, 0.1, 10)                              #c lambda2 ~ Gamma(mu2, nu2)
    sim.tausq[[1]] <- tausq <- runif(B^2*l+B, 0.1, 10)                                  #c 2*tau²|2*csi²,lambda1 ~ exp(lambda1²/(2*csi²))
    sim.U_v[[1]] <- U_v <- runif(B^2*l+B, 0.1, 1)                                       #c these are the theta_v in the paper, theta_v ~ Gamma(k,s)
    sim.alpha[[1]] <- alpha <- runif(B^2*l+B, 0.1, 1)                                   #c alpha ~ N(0, a*I_R²p)
    sim.v_n.star[[1]] <- lapply(1:N, function(n){alpha*v_n[[n]]})                       #c v*_n = alpha*v_n 
    sim.omega_v[[1]] <- U_v/alpha^2                                                     #c omega_v = theta_v / alpha²
  } else {                                                                        #c otherwise, take initial conditions passed in input
    sim.w[[1]] <- w <- initial$w
    sim.v_n[[1]] <- v_n <- initial$v_n
    sim.Lambda[[1]] <- Lambda <- initial$Lambda
    sim.lambda1sq[[1]] <- lambda1sq <- initial$lambda1sq
    sim.lambda2[[1]] <- lambda2 <- initial$lambda2
    sim.tausq[[1]] <- tausq <- initial$tausq
    sim.U_v[[1]] <- U_v <- initial$U_v
    sim.alpha[[1]] <- alpha <- initial$alpha          
    sim.v_n.star[[1]] <- v_n.star <- initial$v_n.star
    sim.omega_v[[1]] <- omega_v <- initial$omega_v
  }
  
  ### 3) exploration step
  for(i in 1:numIter){
    
    # recover G, xisq and prec
    temp <- Gamma_recover(tausq=tausq, lambda1sq=lambda1sq, Lambda=Lambda, B=B, l=l)    #c Compute values for gamma_i and csi². gamma_i values are the elements of the Q matrix (S2appendix, pag 2)   
    G <- temp$G
    xisq <- temp$xisq
    prec <- prec_get(HH=HH, Lambda=Lambda)                                              #c for each series, compute kronecker(HH^t, Lambda). This term is needed in several calculations of S2appendix, pag 1.
    # update parameters
    w <- w_update(prec=prec, w_MLE=w_MLE, v_n=v_n, alpha=alpha, tausq=tausq, lambda2, N=N)
    v_n <- v_n_update(prec=prec, w_MLE=w_MLE, w=w, alpha=alpha, U_v=U_v, N=N)
    Lambda <- Lambda_update(eta0=eta0, w=w, v_n=v_n, alpha=alpha, H=H, G=G, K_inv=K_inv, nu=nu, TT=TT, B=B, l=l, N=N)
    lambda1sq <- lambda1sq_update(tausq=tausq, xisq=xisq, B=B, l=l, mu1=mu1, nu1=nu1)          
    lambda2 <- lambda2_update(w=w, B=B, l=l, mu2=mu2, nu2=nu2)                         
    tausq <- tausq_update(lambda1sq=lambda1sq, w=w, xisq=xisq, B=B, l=l) 
    U_v <- U_v_update(v_n=v_n, k1=k1, theta1=theta1, N=N, B=B, l=l)
    alpha <- alpha_update(prec=prec, w_MLE=w_MLE, w=w, v_n=v_n, a=a, N=N, B=B, l=l)
    v_n.star <- lapply(1:N, function(n){alpha*v_n[[n]]})
    omega_v <- U_v/alpha^2

    # store thinned simulations, after warmup                         
    if((i>warmup) & (i %% thin == 0)){  
      sim.w[[(i-warmup)/thin+1]] <- w
      sim.v_n[[(i-warmup)/thin+1]] <- v_n
      sim.Lambda[[(i-warmup)/thin+1]] <- Lambda
      sim.lambda1sq[[(i-warmup)/thin+1]] <- lambda1sq
      sim.lambda2[[(i-warmup)/thin+1]] <- lambda2
      sim.tausq[[(i-warmup)/thin+1]] <- tausq
      sim.U_v[[(i-warmup)/thin+1]] <- U_v 
      sim.alpha[[(i-warmup)/thin+1]] <- alpha 
      sim.v_n.star[[(i-warmup)/thin+1]] <- v_n.star
      sim.omega_v[[(i-warmup)/thin+1]] <- omega_v
    }
    
    # report iter
    if (i%%5000==0){
      cat("i=", i, "\n", sep="")
    }
  }
  
  ### 4) return results
  return(list(sim.w=sim.w[-1],    # remove the first sample, that was the initial condition
              sim.v_n=sim.v_n[-1],
              sim.Lambda=sim.Lambda[-1], 
              sim.lambda1sq=sim.lambda1sq[-1], 
              sim.lambda2=sim.lambda2[-1], 
              sim.tausq=sim.tausq[-1], 
              sim.U_v=sim.U_v[-1],
              sim.alpha=sim.alpha[-1],
              sim.v_n.star=sim.v_n.star[-1],
              sim.omega_v=sim.omega_v[-1]))   
}

# ------------------------------------------------------------------------------
#    Functions to update parameters and other auxiliary functions
# ------------------------------------------------------------------------------
Posdef <- function (n, ev = runif(n, 0, 10)) {
  # generate n-by-n positive definite matrix
  # used for Lambda ~ Wishart(k,nu), equation (18) in the paper
  # (when I call the function, n=B, that is the length of vectors constituting the trajectories)
  Z <- matrix(ncol=n, rnorm(n^2))   
  decomp <- qr(Z)                   #c QR decomposition: Z=QR, with Q=orthogonal, R=upper triangular matrix.
  Q <- qr.Q(decomp)                 
  R <- qr.R(decomp)                 
  d <- diag(R)                      #c vector containing the n elements of the eigenvalues of R 
  ph <- d / abs(d)                  #c (d1,d2,d3) / (|d1|, |d2|, |d3|) = (+-1, +-1, +-1)
  O <- Q %*% diag(ph)               
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}


Gamma_recover <- function(tausq, lambda1sq, Lambda, B, l){
  # (0.1) recover values of gamma_i e csi²_i (S2appendix, pag 2)
  
  M <- kronecker(diag(1, B*l+1), chol2inv(chol(Lambda))) #c M = I_B*p+1 x Lambda⁻1, M is a squared matrix of size (B*l+1)B = B²l+B
  k <- B^2*l+B
  
  # for xisq                                       
  xisq <- c(rep(NA, B-1), M[B,B])                        #c to store the B values of csi²_i (define a vector of NAs and append the element M_B,B)
  for(j in (B-1):1) xisq[j] <- M[j,j]- M[j,(j+1):B] %*% solve(M[(j+1):B,(j+1):B]) %*% M[(j+1):B,j]   #c calculation at the top of pag. 5 in the paper
  xisq <- rep(xisq, B*l+1)
  
  # for m_j
  m.all.possible <- vector("list", B-1)                  #c it's not clear what those m_j are... btw, they appear in the computation of csi². They are B-1 values (the Bth value is M_B,B, I think)
  for(j in (B-1):1) m.all.possible[[j]] <- M[j,(j+1):B] %*% solve(M[(j+1):B,(j+1):B])
  
  # for gamma
  G <- c(rep(NA, k-1), sqrt(tausq[k]*lambda1sq[k]))      #c B²p values, calculation nb. 7, at pag 2 of S2appendix --> computation of gamma_i in matrix Q
  for(j in (k-1):1){
    k1 <- (k-j)%%B        
    k2 <- floor((k-j)/B)
    if(k1>0){
      m <- m.all.possible[[B-k1]]
      G[j] <- sqrt(tausq[j]*lambda1sq[j]) - sum(m*G[(j+1):(j+k1)])
    } else {
      G[j] <- sqrt(tausq[j]*lambda1sq[j])
    }
  }
  return(list(G=G, xisq=xisq))                           #c return gamma_i and csi²_i values
}


prec_get <- function(HH, Lambda){lapply(HH, function(x){kronecker(x, Lambda)})}
   #c list of length N. For each series, compute kronecker(HH^t, lambda) --> this terms is needed in several calculations of S2appendix, pag 1.


w_update <- function(prec, w_MLE, v_n, alpha, tausq, lambda2, N){
  # (1) update w --> Computation of g(w). Point 1), pag 1, S2appendix.
  #c w is the same for all trajectories
  
  D_vec <- 2*tausq/(2*lambda2*tausq+1)                             #c diagonal elements of D⁻1. D⁻1 is the covariance matrix of w. w|D ~ N(0, D⁻1)
  mu_n <- lapply(1:N, function(n){w_MLE[[n]] - alpha*v_n[[n]] })   #c Estimation of w. To recap:
                                                                   #c VAR process for each trajectory (n) is: Y_n = w_n Z_n + noise_n = (w + v_n)Z_n + noise_n. The initial condition of w is its ML estimation.
                                                                   #c for computational reasons we also assumed: w_n = w + alpha*v_n 
                                                                   #c This implies that w = w_n - alpha*v_n. Therefore, for each serie I get a different estimate for w, that I call mu_n.
  
  Ssq_w <- chol2inv(chol(Reduce("+", prec) + diag(1/D_vec)))       #c S²_w, calculations of point 1), pag 1, S2appendix. S²_w = (sum_n=1..N {HH_n x Lambda} + D)⁻1
                                                                   #c prec = vector containing all matrices kronecker(HH_n, Lambda)
  M_w <- Ssq_w %*% apply(mapply("%*%", prec, mu_n), 1, sum)        #c M_w = S²_w * [sum_n=1..N{kronecker(HH_n, Lambda*mu_n)}]
  w <- as.vector(rmvn(n=1, mu=M_w, sigma=chol(Ssq_w), isChol=TRUE,  ncores = 8)) 
                                                                   #c Create n=1 matrices of Gaussian r.v., with mean defined by M_w and covariance by S²_w.
  return(w)
}


v_n_update <- function(prec, w_MLE, w, alpha, U_v, N){  
  # (2) update v_n. Computation of g(v_n), conditionally on the other params. Point 2), S2appendix
  # It is a MVN, defined by different hyperparameters for each series
  
  v_n <- lapply(1:N, function(n){                                       #c for each series:
    Ssq <- chol2inv(chol(t(t(prec[[n]]*alpha)*alpha) + diag(U_v)))      #c covariance matrix S²_v_n = [diag(alpha)*kronecker(HH,Lambda)*diag(alpha) + Theta_v]⁻1
    mu <- w_MLE[[n]]-w                                                  
    M <- Ssq %*% (alpha*(prec[[n]]%*%mu))                               #c expected value of the distrib. M = S²*[diag(alpha)*kronecker(HH,Lambda)*(w_n-w)]⁻1
    as.vector(rmvn(n=1, mu=M, sigma=chol(Ssq), isChol=TRUE, ncores=8))  #c generate a random matrix, by sampling from MVN(M,S²).
  })
  return(v_n)
}


Lambda_update <- function(eta0, w, v_n, alpha, H, G, K_inv, nu, TT, B, l, N){
  # (4) update Lambda. Calculation of g(Lambda). Point 7), S2appendix. 
  # g(Lambda) = Wishart(S_Lambda, nu_Lambda)
  
  S <- lapply(1:N, function(n){                                         #c for each series, compute S_n (Point 7, S2appendix)
    W <- matrix(w + alpha*v_n[[n]], B)                                  #c W_n = vec⁻1(w + alpha*v_n), Fill in a matrix of size (R x Rp) with the R²p coefs contained in (w + alpha*v_n)
    tcrossprod(eta0[[n]] - W%*%H[[n]])                                  #c S_n = (Y-WH)(Y-WH)^t  
  })
  S_Lambda <- chol2inv(chol(Reduce('+', S) + K_inv + 2*tcrossprod(matrix(G, nrow=B))))  #c S⁻1_Lambda = [sum_n{S_n} + K⁻1 + 2QQ^t]⁻1
  Sm <- min(eigen(S_Lambda)$values)                                     #c smallest eigenval of S_Lambda
  if(Sm<=0) S_Lambda <- S_Lambda+diag(abs(Sm)+0.1, B)                   #c if this eigenval is negative, then modify S_Lambda such that all eigenvals are positive (S_Lambda has to be positive defined) 
  Lambda <- rWishart(1, df=(TT-l)*N+2*(B*l+1)+nu, Sigma=S_Lambda)[,,1]  #c generate Lambda ~ Wishart(nu_Lambda=sum_n{T_n - p}+nu+2Rp, df=sum(TT-l)+2*B*l+nu)
  while(sum(eigen(Lambda)$value<=0) | (!isSymmetric(Lambda))){          #c sample again, until all eigenvals of Lambda are positives and Lambda is symmetric                           
    Lambda <- rWishart(1, df=(TT-l)*N+2*(B*l+1)+nu, Sigma=S_Lambda)[,,1]
  }
  return(Lambda)
}


lambda1sq_update <- function(tausq, xisq, B, l, mu1, nu1){
  # (5) update lambda1sq for j=1,...,B^2*p. Computation of g(lambda1²) = Gamma(mu_lambda1, nu_lambda1). Point 5), S2appendix.
  
  if(length(mu1)==1) mu1 <- rep(mu1, B^2*l+B)   #c We need R²p coefs mu1, and R²p coefs nu1
  if(length(nu1)==1) nu1 <- rep(nu1, B^2*l+B)
  nu_lambda1sq <- nu1+2
  mu_lambda1sq <- nu_lambda1sq*xisq*mu1/(2*tausq*mu1+nu1*xisq)
  lambda1sq <- abs(rgamma(n=rep(1, B^2*l+B), shape=nu_lambda1sq/2, scale=2*mu_lambda1sq/nu_lambda1sq))
  return(lambda1sq)
}


lambda2_update <- function(w, B, l, mu2, nu2){
  # (6) update lambda2 for j=1,...,B^2*p. Computation of g(lambda2) = Gamma(mu_lambda2, nu_lambda2). Point 6), S2appendix
  
  if(length(mu2)==1) mu2 <- rep(mu2, B^2*l+B)   #c We need R²p coefs mu2, and R²p coefs nu2
  if(length(nu2)==1) nu2 <- rep(nu2, B^2*l+B)
  nu_lambda2 <- nu2+2
  mu_lambda2 <- nu_lambda2*mu2/(w^2*mu2+nu2)
  lambda2 <- abs(rgamma(n=rep(1, B^2*l+B), shape=nu_lambda2/2, scale=2*mu_lambda2/nu_lambda2))
  return(lambda2)
}


tausq_update <- function(lambda1sq, w, xisq, B, l){
  # (7) update tausq for j=1,...,B^2*p. Computation of g(tau²). Point 8), S2appendix.
  # First, compute g(1/2/taus²) = InverseGaussian(M,S). Then, invert to obtain values for tau².
  # (N.B. The InverseGaussian distrib. is related with the Gaussin, but it's not its inverse)
  
  M <- sqrt(lambda1sq/(w^2)/xisq)
  S <- lambda1sq/xisq
  temp <- abs(rinvgauss(n=rep(1, B^2*l+B), mean=M, shape=S))
  tausq <- 1/2/temp
  return(tausq)
}


U_v_update <- function(v_n, k1, theta1, N, B, l){
  # (8) update U_v. Computation of g(theta_v) = Gamma(k_v, s_v). Point 4), S2 appendix
  # They are R²p coefs
  
  d <- Reduce("+", lapply(v_n, function(x){x^2}))
  kv <- rep(N/2+k1, B^2*l+B)
  thetav <- 1/(d/2+1/theta1)                                    #c theta1 is named s in calculations of S2appendix
  U_v <- abs(rgamma(n=rep(1, B^2*l+B), shape=kv, scale=thetav))  
  return(U_v)
}


alpha_update <- function(prec, w_MLE, w, v_n, a, N, B, l){
  # (10) update alpha. Calculation of g(alpha) = N(M_alpha, S²_alpha). Point 3), S2appendix
  
  Ssq_alpha <- chol2inv(chol(Reduce("+",lapply(1:N, function(n){
    t(t(prec[[n]]*v_n[[n]])*v_n[[n]])
    })) + diag(a, B^2*l+B)))
  M_alpha <- Ssq_alpha %*%  Reduce("+", lapply(1:N, function(n){(prec[[n]]%*%(w_MLE[[n]]-w))*v_n[[n]]}))
  alpha <- as.vector(rmvn(n=1, mu=M_alpha, sigma=chol(Ssq_alpha), isChol=TRUE, ncores=8))
  return(alpha)
}
