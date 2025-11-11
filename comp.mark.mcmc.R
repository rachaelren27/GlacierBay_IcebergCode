comp.mark.mcmc <- function(s.obs, n, u, X.full, full.win.idx, win.idx, ds, 
                           a.zeta, b.zeta, zeta.tune,
                           mu.beta, s2.beta, beta.tune, 
                           mu.alpha, s2.alpha, mu.gamma, s2.gamma, q, r,
                           n.mcmc){
  n <- length(u)
  p <- ncol(X.full)
  
  X.obs <- X.full[win.idx,]
  X.obs.plus <- cbind(rep(1, nrow(X.obs)), X.obs)
  
  alpha.save <- matrix(0, n.mcmc, p+1)
  gamma.save <- rep(0, n.mcmc)
  s2.u.save <- rep(0, n.mcmc)
  beta.0.save <- rep(0, n.mcmc)
  beta.save <- matrix(0, n.mcmc, p)
  
  # initialize values
  alpha <- rep(0, p+1)
  gamma <- 0
  s2.u <- 1
  zeta <- 1
  beta <- rep(0, p)
  
  X.full.beta <- X.full%*%beta
  X.obs.beta <- X.full.beta[win.idx,]
  X.full.win.beta <- X.full.beta[full.win.idx,]
  X.obs.plus.alpha <- X.obs.plus%*%alpha
  lambda <- exp(X.obs%*%beta)
  X.obs.plus.2 <- t(X.obs.plus)%*%X.obs.plus
  
  accept.zeta <- 0
  accept.beta <- 0
  
  # sub-routine
  spp.loglik <- function(zeta, X.obs.beta, X.full.win.beta, ds){
    llam <- log(zeta) + X.obs.beta
    lam.int <- sum(zeta*exp(log(ds) + X.full.win.beta))
    sum(llam) - lam.int
  }
  
  # MCMC algorithm
  for(k in 1:n.mcmc){
    if(k %% 1000 == 0){cat(k, " ")}
    # propose and update zeta (exp(beta_0))
    zeta.star <- rlnorm(1, log(zeta), zeta.tune)
    
    mh.1 <- spp.loglik(zeta.star, X.obs.beta, X.full.win.beta, ds) + 
            dgamma(zeta.star, a.zeta, b.zeta, log = TRUE) +
            dlnorm(zeta, log(zeta), zeta.tune, log = TRUE)
    
    mh.2 <- spp.loglik(zeta, X.obs.beta, X.full.win.beta, ds) + 
            dgamma(zeta, a.zeta, b.zeta, log = TRUE) +
            dlnorm(zeta.star, log(zeta.star), zeta.tune, log = TRUE)
    
    if(exp(mh.1 - mh.2) > runif(1)){
      zeta <- zeta.star
      accept.zeta <- accept.zeta + 1
    }
    
    # propose and update beta
    beta.star <- t(mvnfast::rmvn(1, beta, beta.tune*diag(p)))
    lambda.star <- exp(X.obs%*%beta.star)
    X.full.beta.star <- X.full%*%beta.star
    X.obs.beta.star <- X.full.beta.star[win.idx,]
    X.full.win.beta.star <- X.full.beta.star[full.win.idx,]
    
    mh.1 <- spp.loglik(zeta, X.obs.beta.star, X.full.win.beta.star, ds) +
            sum(dnorm(u, X.obs.plus.alpha + gamma*lambda.star, sqrt(s2.u), log = TRUE)) +
            sum(mvnfast::dmvn(t(beta.star), mu.beta, s2.beta*diag(p), log=TRUE))
    
    mh.2 <- spp.loglik(zeta, X.obs.beta, X.full.win.beta, ds) +
            sum(dnorm(u, X.obs.plus.alpha + gamma*lambda, sqrt(s2.u), log = TRUE)) +
            sum(mvnfast::dmvn(t(beta), mu.beta, s2.beta*diag(p), log=TRUE))
    
    if(exp(mh.1 - mh.2) > runif(1)){
      beta <- beta.star
      lambda <- exp(X.obs%*%beta)
      X.full.beta <- X.full%*%beta
      X.obs.beta <- X.full.beta[win.idx,]
      X.full.win.beta <- X.full.beta[full.win.idx,]
      accept.beta <- accept.beta + 1
    }
    
    # update gamma
    a.gamma <- sum(lambda^2)/s2.u + 1/s2.gamma
    b.gamma <- (t(u - X.obs.plus.alpha)%*%lambda)/s2.u + mu.gamma/s2.gamma
    gamma <- rnorm(1, b.gamma/a.gamma, sqrt(1/a.gamma))
    
    # update alpha
    A.alpha <- X.obs.plus.2/s2.u + diag(p+1)/s2.alpha
    b.alpha <- t(X.obs.plus)%*%(u - gamma*lambda)/s2.u + mu.alpha/s2.alpha
    A.alpha.inv <- solve(A.alpha)
    alpha <- t(rmvn(1, A.alpha.inv%*%b.alpha, A.alpha.inv))
    X.obs.plus.alpha <- X.obs.plus%*%alpha
    
    # update s2.u
    q.tilde <- (n/2) + q
    r.tilde <- 1/(sum((u - (X.obs.plus.alpha + gamma*lambda))^2)/2 + (1/r))
    s2.u <- 1/rgamma(1, q.tilde, scale = r.tilde)
    
    alpha.save[k,] <- alpha
    gamma.save[k] <- gamma
    beta.0.save[k] <- log(zeta)
    beta.save[k,] <- t(beta)
    s2.u.save[k] <- s2.u
  };cat("\n")
  
  print(paste("Zeta acceptance ratio:", accept.zeta/n.mcmc))
  print(paste("Beta acceptance ratio:", accept.beta/n.mcmc))
  
  return(list(alpha.save = alpha.save, gamma.save = gamma.save, beta.0.save = beta.0.save,
              beta.save = beta.save, s2.u.save = s2.u.save))
}