cond.mark.mcmc <- function(u, X, beta.prop, n.mcmc, mu.alpha, s2.alpha,
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
  beta <- beta.prop[, ncol(beta.prop)]
  
  X.plus.alpha <- X.plus%*%alpha
  lambda <- exp(X%*%beta)
  X.plus.2 <- t(X.plus)%*%X.plus
  
  accept <- 0
  
  for(k in 1:n.mcmc){
    if(k %% 1000 == 0){cat(k, " ")}
    # propose and update beta
    idx.star <- sample(1:n.mcmc, 1)
    beta.star <- c(beta.prop[,idx.star])
    
    mh.1 <- sum(dnorm(u, X.plus.alpha + gamma*exp(X%*%beta.star), sqrt(s2.u), log = TRUE))
    mh.2 <- sum(dnorm(u,X.plus.alpha + gamma*lambda, sqrt(s2.u), log = TRUE))
    
    if(exp(mh.1 - mh.2) > runif(1)){
      beta <- beta.star
      lambda <- exp(X%*%beta)
      accept <- accept + 1
    }
    
    # update gamma
    a.gamma <- sum(lambda^2)/s2.u + 1/s2.gamma
    b.gamma <- (t(u - X.plus.alpha)%*%lambda)/s2.u + mu.gamma/s2.gamma
    gamma <- rnorm(1, b.gamma/a.gamma, sqrt(1/a.gamma))
    
    # update alpha
    A.alpha <- X.plus.2/s2.u + diag(p+1)/s2.alpha
    b.alpha <- t(X.plus)%*%(u - gamma*lambda)/s2.u + mu.alpha/s2.alpha
    A.alpha.inv <- solve(A.alpha)
    alpha <- t(rmvn(1, A.alpha.inv%*%b.alpha, A.alpha.inv))
    X.plus.alpha <- X.plus%*%alpha
    
    # update s2.u
    q.tilde <- (n/2) + q
    r.tilde <- 1/(sum((u - (X.plus.alpha + gamma*lambda))^2)/2 + (1/r))
    s2.u <- 1/rgamma(1, q.tilde, scale = r.tilde)
    
    alpha.save[k,] <- alpha
    gamma.save[k] <- gamma
    beta.save[k,] <- t(beta)
    s2.u.save[k] <- s2.u
  };cat("\n")
  
  print(paste("Beta acceptance ratio:", accept/n.mcmc))
  
  return(list(alpha.save = alpha.save, gamma.save = gamma.save, beta.save = beta.save, s2.u.save = s2.u.save))
}