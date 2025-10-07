cond.mark.mcmc <- function(u, X.obs, n.mcmc, mu.alpha, Sig.alpha, q, r){
  n <- length(u)
  p <- ncol(X.obs)
  
  alpha.save <- matrix(0, n.mcmc, p)
  s2.save <- rep(0, n.mcmc)
  
  alpha <- rep(0, p)
  s2 <- 1
  
  for(k in 1:n.mcmc){
    A <- t(X.obs)%*%((1/s2)*diag(n))%*%X.obs + solve(Sig.alpha)
    b <- t(t(u)%*%((1/s2)*diag(n))%*%X.obs + mu.alpha%*%solve(Sig.alpha))
    A.inv <- solve(A)
    alpha <- t(rmvn(1, A.inv%*%b, A.inv))
    
    X.alpha <- X.obs%*%alpha
    q.tilde <- (n/2) + q
    r.tilde <- 1/(sum((u - X.alpha)^2)/2 + (1/r))
    s2 <- rinvgamma(1, q.tilde, scale = r.tilde)
    
    alpha.save[k,] <- alpha
    s2.save[k] <- s2
  }
  
  return(list(alpha.save = alpha.save, s2.save = s2.save))
}