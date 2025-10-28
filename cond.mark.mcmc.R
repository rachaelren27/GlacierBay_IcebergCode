cond.mark.mcmc <- function(u, X, beta.prop, n.mcmc, mu.alpha, Sig.alpha,
                           mu.gamma, s2.gamma, q, r){
  n <- length(u)
  p <- ncol(X)
  
  X.plus <- cbind(rep(1, nrow(X)), X)
  
  alpha.save <- matrix(0, n.mcmc, p+1)
  gamma.save <- rep(0, n.mcmc)
  s2.u.save <- rep(0, n.mcmc)
  beta.save <- matrix(0, n.mcmc, p)
  
  # initialize values
  alpha <- rep(0, p+1)
  gamma <- 0
  s2.u <- 1
  beta <- beta.prop[, n.mcmc]
  
  X.plus.alpha <- X.plus%*%alpha
  X.beta <- X%*%beta
  s2.I.inv <- (1/s2.u)*diag(n)
  
  Sig.alpha.inv <- solve(Sig.alpha)
  
  accept <- 0
  
  for(k in 1:n.mcmc){
    if(k %% 1000 == 0){cat(k, " ")}
    # propose and update beta
    idx.star <- sample(1:n.mcmc, 1)
    beta.star <- c(beta.prop[,idx.star])
    
    mh.1 <- sum(dnorm(X.plus.alpha + gamma*X%*%beta.star, s2.u, log = TRUE))
    mh.2 <- sum(dnorm(X.plus.alpha + gamma*X%*%beta, s2.u, log = TRUE))
    
    if(exp(mh.1 - mh.2) > runif(1)){
      beta <- beta.star
      X.beta <- X%*%beta
      accept <- accept + 1
    }
    
    # update gamma
    a.gamma <- t(X.beta)%*%s2.I.inv%*%X.beta + (1/s2.gamma)
    b.gamma <- t(u - X.plus.alpha)%*%s2.I.inv%*%X.beta + (mu.gamma/s2.gamma)
    gamma <- rnorm(1, b.gamma/a.gamma, 1/a.gamma)
    
    # update alpha
    A.alpha <- t(X.plus)%*%s2.I.inv%*%X.plus + Sig.alpha.inv
    b.alpha <- t(t(u - gamma*X.beta)%*%s2.I.inv%*%X.plus + mu.alpha%*%Sig.alpha.inv)
    A.alpha.inv <- solve(A.alpha)
    alpha <- t(rmvn(1, A.alpha.inv%*%b.alpha, A.alpha.inv))
    X.plus.alpha <- X.plus%*%alpha
    
    # update s2.u
    q.tilde <- (n/2) + q
    r.tilde <- 1/(sum((u - (X.plus.alpha + gamma*X.beta))^2)/2 + (1/r))
    s2.u <- 1/rgamma(1, q.tilde, scale = r.tilde)
    s2.I.inv <- (1/s2.u)*diag(n)
    
    alpha.save[k,] <- alpha
    gamma.save[k] <- gamma
    beta.save[k,] <- t(beta)
    s2.u.save[k] <- s2.u
  };cat("\n")
  
  print(paste("Beta acceptance ratio:", accept/n.mcmc))
  
  return(list(alpha.save = alpha.save, gamma.save = gamma.save, beta.save = beta.save, s2.u.save = s2.u.save))
}