################################################################################
# Anaylize results of the Gibb sampling
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
################################################################################

library(rstan)

# ------------------------------------------------------------------------------
#    Inference of coefficients of the HVAR model
# ------------------------------------------------------------------------------
w_inference <- function(chains, B, l, n.sim, warmup=0, thin=1, filename, plot=FALSE, d=3){
  #c function to infer the estimation for w, given the samples got from the MC exploration
  
  n.chain <- length(chains)                                                   
  sel <- seq(warmup+1, n.sim, by=thin)                                          #c sequence: (warmup+1, warmup+1+thin, ..+2*thin, ..., nsim+1) = (252, ..., 501)
  ll <- length(sel)                                                             #c nb. of samples that we will consider for each chain
  sw <- array(NA, dim=c(ll, n.chain, B^2*l+B))                                  #c empty matrix of size (l x nchain x R²p) = (nb samples x nb chains x nb. of elements in w)
  for(i in 1:n.chain) sw[,i,] <- do.call("rbind", chains[[i]]$sim.w[sel])       #c for each chain, look at samples of w and take the ones with index 'sel'
  mon.w <- rstan::monitor(sw, warmup=warmup, print=FALSE, probs=c(0.025, 0.975, 0.05, 0.95))    #c function to print summary of a MC exploration:
                                                                                #c   sims=sw, A 3D array (iterations * chains * parameters) of MCMC simulations
                                                                                #c   warmup, nb of warmup iterations to be excluded when computing the summaries. 
                                                                                #c           The default is half of the total number of iterations. If sims doesn't 
                                                                                #c           include the warmup iterations then warmup should be set to zero.
                                                                                #c   probs, quantiles of interest
  w.mode <-  apply(sw, 3, function(x){find.mode(x, d=d)})                       #c compute mode of samples. (Apply the function along the third dimension of matrix sw)
  w.median <- mon.w$Q50                                               
  if(plot){
    pdf(filename, width=14)   # trace & hist
    par(mfrow=c(1, 2))
    for(j in 1:(B^2*l+B)) {
      matplot(sw[,,j], type="l", lty=1, main=paste("w", j, sep=""), ylab="")
      hist(sw[,,j], main=paste("w", j, sep=""), breaks=30, ylab="", xlab="")
      #abline(v= w.median[j], col=2)
      #axis(side=1, at= w.median[j], labels=round( w.median[j], d))
      abline(v=w.mode[j], col=4)
      axis(side=3, at=w.mode[j], labels=w.mode[j])
    }
    dev.off()
  }
  return(list(sw=sw, mon.w=mon.w, w.mode=w.mode, w.median=w.median))            #c return: list of samples for w, 
                                                                                #c         table with summary of statistics for sw,
                                                                                #c         mode and median
}


v_inference <- function(chains, B, l, N, n.sim, warmup=0, thin=1, filename, plot=FALSE, d=3){
  #c function to infer estimates for v_n, analogously to w_inference() 

  n.chain <- length(chains)
  sel <- seq(warmup+1, n.sim, by=thin)
  ll <- length(sel)
  sv <- array(NA, dim=c(ll, n.chain, (B^2*l+B)*N))
  for(i in 1:n.chain){sv[,i,] <-  do.call("rbind", lapply(chains[[i]]$sim.v_n.star[sel], unlist))}
  mon.v <- rstan::monitor(sv, warmup=warmup, print=F, probs=c(0.025, 0.975, 0.05, 0.95))
  v.mode <-  apply(sv, 3, function(x){find.mode(x, d=d)}) 
  v.median <- mon.v$Q50
  if(plot){
    pdf(filename, width=14)
    par(mfrow=c(1, 2))
    for(j in 1:((B^2*l+B)*N)) {
      matplot(sv[,,j], type="l", lty=1, main=paste("v", j, sep=""), ylab="")
      hist(sv[,,j], main=paste("v", j, sep=""), breaks=40, ylab="", xlab="")
      #abline(v=v.median[j], col=2)
      #axis(side=1, at=v.median[j], labels=round(v.median[j], d))
      abline(v=v.mode[j], col=4)
      axis(side=3, at=v.mode[j], labels=v.mode[j])
    }
    dev.off()
  }
  return(list(sv=sv, mon.v=mon.v, v.mode=v.mode, v.median=v.median))
}


u_inference <- function(chains, B, l, J, n.sim, warmup=0, thin=1, filename, plot=FALSE, d=3){
  n.chain <- length(chains)
  sel <- seq(warmup+1, n.sim, by=thin)
  ll <- length(sel)
  su <- array(NA, dim=c(ll, n.chain, (B^2*l+B)*N*J))
  for(i in 1:n.chain){su[,i,] <-  do.call("rbind", lapply(chains[[i]]$sim.u_nj.star[sel], unlist))}   # changed u_j to u_nj
  mon.u <- rstan::monitor(su, warmup=warmup, print=F)
  u.mode <-  apply(su, 3, function(x){find.mode(x, d=d)}) 
  u.median <- mon.u$Q50
  if(plot){
    pdf(filename, width=14)
    par(mfrow=c(1, 2))
    for(j in 1:((B^2*l+B)*J)){
      matplot(su[,,j], type="l", lty=1, main=paste("u", j, sep=""),ylab="")
      hist(su[,,j], main=paste("u", j, sep=""), breaks=40, ylab="", xlab="")
      #abline(v=u.median[j], col=2)
      #axis(side=1, at=u.median[j], labels=round(u.median[j], d))
      abline(v=u.mode[j], col=4)
      axis(side=3, at=u.mode[j], labels=u.mode[j])
    }
    dev.off()
  }
  return(list(su=su, mon.u=mon.u, u.mode=u.mode, u.median=u.median))
}


L_inference <- function(chains, B, l, n.sim, warmup=0, thin=1, filename, plot=FALSE, d=6){
  n.chain <- length(chains)
  sel <- seq(warmup+1, n.sim, by=thin)
  ll <- length(sel)
  sL <- array(NA, dim=c(ll, n.chain, B^2))
  for(i in 1:n.chain){sL[,i,] <-  do.call("rbind", lapply(chains[[i]]$sim.Lambda[sel], c))}
  mon.L <- rstan::monitor(sL, warmup=warmup, print=F)
  L.mode <-  apply(sL, 3, function(x){find.mode(x, d=d)}) 
  L.median <- mon.L$Q50
  if(plot){
    pdf(filename, width=14)
    par(mfrow=c(1, 2))
    for(j in 1:(B^2)){
      matplot(sL[,,j], type="l", lty=1, main=paste("lambda", j, sep=""),ylab="")
      hist(sL[,,j], main=paste("lambda", j, sep=""), breaks=40, ylab="", xlab="")
      #abline(v=L.median[j], col=2)
      #axis(side=1, at=L.median[j], labels=round(L.median[j], d))
      abline(v=L.mode[j], col=4)
      axis(side=3, at=L.mode[j], labels=L.mode[j])
    }
    dev.off()
  }
  return(list(sL=sL, mon.L=mon.L, L.mode=L.mode, L.median=L.median))
}


omega_v_inference <- function(chains, B, l, n.sim, warmup=0, thin=1, filename, plot=FALSE, d=3, log=T){
  n.chain <- length(chains)
  sel <- seq(warmup+1, n.sim, by=thin)
  ll <- length(sel)
  somega_v <- array(NA, dim=c(ll, n.chain, B^2*l+B))
  for(i in 1:n.chain){somega_v[,i,] <-  do.call("rbind", chains[[i]]$sim.omega_v[sel])}
  mon.omega_v <- rstan::monitor(somega_v, warmup=warmup, print=F)
  # omega_v.mode <-  apply(somega_v, 3, function(x){find.mode(x, d=-1)}) 
  omega_v.mode <-  apply(somega_v, 3, function(x){find.mode(x, d=d, log=log)}) 
  omega_v.median <- mon.omega_v$Q50 ###[,6]
  if(plot){
    pdf(filename, width=14)
    par(mfrow=c(1, 2))
    for(j in 1:(B^2*l+B)){
      matplot(somega_v[,,j], type="l", lty=1, log="y", main=paste("omega_v", j, sep=""), ylab="")
      hist(log(somega_v[,,j]), main=paste("log(omega_v", j, ")", sep=""), breaks=40, ylab="", xlab="")
      #abline(v=log(omega_v.median[j]), col=2)
      #axis(side=1, at=log(omega_v.median[j]), labels=round(log(omega_v.median[j]), 2))
      abline(v=log(omega_v.mode[j]), col=4)
      axis(side=3, at=log(omega_v.mode[j]), labels=round(log(omega_v.mode[j]), 2))
    }
    dev.off()
  }
  return(list(somega_v=somega_v, 
              mon.omega_v=mon.omega_v, 
              omega_v.mode=omega_v.mode, 
              omega_v.median=omega_v.median))
}


lambda1sq_inference <- function(chains, B, l,  n.sim, warmup=0, thin=1, filename, plot=FALSE){
  n.chain <- length(chains)
  sel <- seq(warmup+1, n.sim, by=thin)
  ll <- length(sel)
  slambda1sq <- array(NA, dim=c(ll, n.chain, B^2*l+B))
  for(i in 1:n.chain){slambda1sq[,i,] <-  do.call("rbind", chains[[i]]$sim.lambda1sq[sel])}
  mon.lambda1sq <- rstan::monitor(slambda1sq, warmup=warmup, print=F)
  lambda1sq.median <- mon.lambda1sq$Q50
  if(plot){
    pdf(filename, width=14)
    par(mfrow=c(1, 2))
    for(j in 1:(B^2*l+B)){
      matplot(slambda1sq[,,j], type="l", lty=1, log="y", main=paste("lambda1sq", j, sep=""), ylab="")
      hist(log(slambda1sq[,,j]), main=paste("log(lambda1sq", j, ")", sep=""), breaks=40, ylab="", xlab="")
      abline(v=log(lambda1sq.median[j]), col=2)
      axis(side=1, at=log(lambda1sq.median[j]), labels=round(log(lambda1sq.median[j]), 2))
    }
    dev.off()
  }
  return(list(slambda1sq=slambda1sq, mon.lambda1sq=mon.lambda1sq, lambda1sq.median=lambda1sq.median))
}


lambda2_inference <- function(chains, B, l,  n.sim, warmup=0, thin=1, filename, plot=FALSE){
  n.chain <- length(chains)
  sel <- seq(warmup+1, n.sim, by=thin)
  ll <- length(sel)
  slambda2<- array(NA, dim=c(ll, n.chain, B^2*l+B))
  for(i in 1:n.chain){slambda2[,i,] <-  do.call("rbind", chains[[i]]$sim.lambda2[sel])}
  mon.lambda2<- rstan::monitor(slambda2, warmup=warmup, print=F)
  lambda2.median <- mon.lambda2$Q50
  if(plot){
    pdf(filename, width=14)
    par(mfrow=c(1, 2))
    for(j in 1:(B^2*l+B)){
      matplot(slambda2[,,j], type="l", lty=1, log="y", main=paste("lambda2_", j, sep=""), ylab="")
      hist(log(slambda2[,,j]), main=paste("log(lambda2_", j, ")", sep=""), breaks=40, ylab="", xlab="")
      abline(v=log(lambda2.median[j]), col=2)
      axis(side=1, at=log(lambda2.median[j]), labels=round(log(lambda2.median[j]), 2))
    }
    dev.off()
  }
  return(list(slambda2=slambda2, mon.lambda2=mon.lambda2, lambda2.median=lambda2.median))
}


tausq_inference <- function(chains, B, l,  n.sim, warmup=0, thin=1, filename, plot=FALSE){
  n.chain <- length(chains)
  sel <- seq(warmup+1, n.sim, by=thin)
  ll <- length(sel)
  stausq <- array(NA, dim=c(ll, n.chain, B^2*l+B))
  for(i in 1:n.chain){stausq[,i,] <-  do.call("rbind", chains[[i]]$sim.tausq[sel])}
  mon.tausq <- rstan::monitor(stausq, warmup=warmup, print=F)
  tausq.median <- mon.tausq$Q50
  if(plot){
    pdf(filename, width=14)
    par(mfrow=c(1, 2))
    for(j in 1:(B^2*l+B)){
      matplot(stausq[,,j], type="l", lty=1, log="y", main=paste("tausq", j, sep=""), ylab="")
      hist(log(stausq[,,j]), main=paste("log(tausq", j, ")", sep=""), breaks=40, ylab="", xlab="")
      abline(v=log(tausq.median[j]), col=2)
      axis(side=1, at=log(tausq.median[j]), labels=round(log(tausq.median[j]), 2))
    }
    dev.off()
  }
  return(list(stausq=stausq, mon.tausq=mon.tausq, tausq.median=tausq.median))
}


# ------------------------------------------------------------------------------
#    Other functions
# ------------------------------------------------------------------------------
find.mode <- function(x, d=3, log=F){
  # compute the mode
  
  if(log){
    hh <- density(log(x))
    y <- exp(hh$x[which.max(hh$y)])
  } else {
    hh <- density(x)            #c estimate the density by using a Gaussian kernel (by default).
    y <- hh$x[which.max(hh$y)]  #c hh$x = the n coordinates of the points where the density is estimated.
                                #c hh$y = the estimated density values. These will be non-negative, but can be zero
                                #c therefore, y will be the max
  }
  return(y)
}

dss <- function(y,mu,Sigma){
  # Dawid-Sebastiani score to evaluate multivariate probabilistic forecasting.
  #
  # Problem: we have to predict a d-dimensional point y. Our prediction is given
  # by a forecast distribution F, with mean mu and covariance Sigma. The D-S
  # score provide a negative-oriented evaluation of the forecast, meaning
  # that smaller scores indicate better forecasts.
  #
  # Ref:
  # Dawid, A. P., and P. Sebastiani, 1999: Coherent dispersion criteria for 
  # optimal experimental design. Ann. Stat.. doi:10.1214/aos/1018031101.
  
  d <- length(y)
  Sig_inv <- solve(Sigma) # chol2inv(chol(Sigma)) --> faster method but more instable
  logdet_Sig <- determinant(Sigma, logarithm=TRUE)$modulus[1]
  ds_score <- logdet_Sig + matrix(y-mu,1,d) %*% Sig_inv %*% matrix(y-mu,d,1)
  return(as.numeric(ds_score))
}
