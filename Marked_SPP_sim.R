setwd("/Users/rlr3795/Desktop/GlacierBay_Project")

library(tidyverse)
library(raster)
library(mvnfast)
library(viridis)
library(sf)
library(spatstat)
library(ggforce)
library(here)
library(tictoc)
library(parallel)
library(foreach)
library(doParallel)
library(coda)

set.seed(1234)

# --- Construct domain --------------------------------------------------------
domain.range <- c(0,1.05)
survey.win <- owin(xrange = domain.range, yrange = domain.range)

# define the coordinates for window squares
win.length <- 0.2
gap <- 0.05
domain.length <- domain.range[2] - domain.range[1]
win.coords <- expand.grid(x = seq(gap, domain.length - win.length - gap, by = win.length + gap), 
                      y = seq(gap, domain.length - win.length - gap, by = win.length + gap))

# create individual squares
wins <- lapply(1:nrow(win.coords), function(i) {
  x0 <- win.coords$x[i]
  y0 <- win.coords$y[i]
  owin(xrange = c(x0, x0 + win.length), yrange = c(y0, y0 + win.length))
})

# combine squares into single window
combined.window <- do.call(union.owin, wins)

# calculate areas
tot.area <- domain.length^2
tot.win.area <- area.owin(combined.window)
n.win <- 16
win.area <- tot.win.area/n.win
tot.nonwin.area <- tot.area - tot.win.area

# plot the window
domain <- owin(xrange = c(0, domain.length), yrange = c(0, domain.length))
plot(domain)
plot(combined.window, add = TRUE)

# construct quadrature grid
n.grid <- 500

x.full <- seq(domain.range[1], domain.range[2], length.out = n.grid)
y.full <- seq(domain.range[1], domain.range[2], length.out = n.grid)
s.full <- expand.grid(x = x.full, y = y.full)

ds <- (x.full[2] - x.full[1])^2

# --- Simulate GP covariate ----------------------------------------------------
GP_2D_sim <- function(n.grid, domain.range, mu = 0,
                                sig2 = 0.5, phi = 0.2, nu = 1e-5) {

  s <- seq(domain.range[1], domain.range[2], length.out = n.grid)
  D1 <- as.matrix(dist(s))
  K <- exp(-(D1^2) / (phi^2) / 2)
  K <- K + nu * diag(n.grid)
  R <- chol(K)
  L <- t(R)
  
  Z1 <- matrix(rnorm(n.grid * n.grid), n.grid, n.grid)
  Z2 <- rnorm(n.grid * n.grid)
  
  part1_mat <- L %*% Z1 %*% t(L)
  part1_vec <- as.vector(part1_mat)
  
  y.sim <- mu + sqrt(sig2) * part1_vec + sqrt(nu) * Z2
  
  locs <- expand.grid(x = s, y = s)
  data.frame(x = locs[,1], y = locs[,2], z = y.sim)
}

out.df <- GP_2D_sim(n.grid, domain.range, mu = 0,
                    sig2 = 0.5, phi = 0.15, nu = 1e-5)

ggplot(out.df, aes(x = x, y = y, fill = z)) +
  geom_tile() +
  scale_fill_viridis() +
  coord_fixed() +
  theme_minimal()


# --- Simulate from windowed IPP -----------------------------------------------
beta.0 <- 7
beta <- c(0.5, 0.5, 1)

X.full <- scale(cbind(s.full[,1], s.full[,2], out.df$z))

lam.full <- exp(beta.0 + X.full%*%beta)
lam.max <- max(lam.full)

full.df <- as.data.frame(cbind(x = s.full[,1], y = s.full[,2], z = lam.full))
full.raster <- rasterFromXYZ(full.df)

M <- rpois(1, lam.max) 
x.superpop <- runif(M, domain.range[1], domain.range[2])
y.superpop <- runif(M, domain.range[1], domain.range[2])
s.superpop <- cbind(x.superpop, y.superpop)
superpop.idx <- cellFromXY(full.raster, cbind(x = x.superpop, y = y.superpop))
X.superpop <- X.full[superpop.idx,]
lam.superpop <- lam.full[superpop.idx]

obs.idx <- rbinom(M, 1, lam.superpop/lam.max) == 1 # superpop -> observed
s.obs <- s.superpop[obs.idx,] # total observed points 
X.obs <- X.superpop[obs.idx,] 
lam.obs <- lam.superpop[obs.idx]
N <- nrow(s.obs) 

# get windowed data
obs.win <- inside.owin(s.obs[,1], s.obs[,2], combined.window)
obs.win.idx <- (1:N)[obs.win] # full observed -> observed windowed
n <- length(obs.win.idx)

full.win <- inside.owin(s.full[,1], s.full[,2], combined.window)
full.win.idx <- (1:n.grid^2)[full.win] # full grid -> full windowed

s.win <- s.obs[obs.win.idx,]
X.win <- X.obs[obs.win.idx,]

X.win.full <- X.full[full.win.idx,]
X.nonwin.full <- X.full[-full.win.idx,]

domain.sf <- st_as_sfc(as.polygonal(domain)) %>% st_sf()
combined.window.sf <- st_as_sfc(as.polygonal(combined.window)) %>% st_sf()
ggplot() +
  geom_sf(data = domain.sf) +
  geom_sf(data = combined.window.sf) +
  geom_point(aes(x = s.win[,1], y = 1.05 - s.win[,2]), col = "red", size = 0.1) +
  # geom_point(aes(x = s.obs[-obs.win.idx,1], y = s.obs[-obs.win.idx,2]), size = 0.5) +
  theme(axis.title = element_blank())


# # --- Fit SPP w/ complete likelihood -------------------------------------------
# n.mcmc <- 100000
# source(here("GlacierBay_Seal", "GlacierBay_SealCode", "spp_win_2D", "spp.comp.mcmc.R"))
# tic()
# out.comp.full <- spp.comp.mcmc(s.win, X.win, X.win.full, ds, n.mcmc, 0.1, 0.01)
# toc() # 385 sec
# 
# # discard burn-in
# n.burn <- 0.1*n.mcmc
# beta.0.save <- out.comp.full$beta.0.save[-(1:n.burn)]
# beta.save <- out.comp.full$beta.save[,-(1:n.burn)]
# 
# # trace plot
# layout(matrix(1:2,2,1))
# plot(beta.0.save,type="l")
# abline(h=beta.0,col=rgb(0,1,0,.8),lty=2)
# matplot(t(beta.save),lty=1,type="l")
# abline(h=beta,col=rgb(0,1,0,.8),lty=2)


# --- Berman Turner Device -----------------------------------------------------
n.bg <- 100000
bg.pts <- rpoint(n.bg, win = combined.window)

# prepare covariates for background sample 
bg.mat <- cbind(bg.pts$x, bg.pts$y)

bg.idx <- cellFromXY(full.raster, bg.mat)
X.bern <- rbind(X.win, X.full[bg.idx,])
y.bern <- rep(0, n + n.bg)
y.bern[1:n] <- 1
bern.df <- data.frame(y = y.bern, x1 = X.bern[,1], x2 = X.bern[,2], x3 = X.bern[,3])


# --- Fit SPP w/ GLM-E ---------------------------------------------------------
n.mcmc <- 100000
n.burn <- 0.1*n.mcmc

## stage 1
out.glm <- glm(y ~ x1 + x2 + x3, family = binomial(link="logit"), data = bern.df)
beta.glm <- coef(out.glm)[-1]
vcov.glm <- vcov(out.glm)[-1,-1]
# sample from glm estimated density
beta.save <- t(mvnfast::rmvn(n.mcmc, mu = beta.glm, sigma = vcov.glm))

## stage 2
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

tic()
lam.int.save <- foreach(k = 1:ncol(beta.save), .combine = c) %dopar% {
  lam.int <- sum(exp(log(ds) + X.win.full%*%beta.save[,k]))
  return(lam.int)
}

X.beta.sum.save <- foreach(k = 1:ncol(beta.save)) %dopar% {
  X.beta.sum <- sum(X.win%*%beta.save[,k])
  return(X.beta.sum)
}
toc() # 33.3 sec

stopCluster(cl) 

out.glm2 <- list(beta.save = beta.save, mu.beta = beta.glm, sigma.beta = vcov.glm,
                 n.mcmc = n.mcmc, n = n, ds = ds, X.full = X.win.full,
                 X.beta.sum.save = X.beta.sum.save, lam.int.save = lam.int.save)

## stage 3
source(here("GlacierBay_Seal", "GlacierBay_SealCode", "spp.stg3.mcmc.nb.R"))
tic()
out.glm3 <- spp.stg3.mcmc.nb(out.glm2)
toc() # 3.1 sec

beta.save <- out.glm3$beta.save
beta.0.save <- out.glm3$beta.0.save

layout(matrix(1:2,2,1))
plot(beta.0.save,type="l")
abline(h=beta.0,col=rgb(0,1,0,.8),lty=2)
matplot(t(beta.save),lty=1,type="l")
abline(h=beta,col=rgb(0,1,0,.8),lty=2)


# --- Generate marks -----------------------------------------------------------
alpha.0 <- -7
alpha <- c(-1, 0.1, -1)
gamma <- -0.1
xi.full <- alpha.0 + X.full%*%alpha + gamma*exp(X.full%*%beta)
xi.full.rast <- rasterFromXYZ(cbind(x = s.full[,1], y = s.full[,2],
                                          z = xi.full))

win.idx <- cellFromXY(full.raster, s.win) # full grid -> observed windowed
s2.u <- 0.5^2
u.win <- rnorm(n, xi.full[win.idx], sqrt(s2.u))
# u.nonwin <- rnorm(N - n, alpha.0 + X.obs.marks[-obs.win.idx,]%*%alpha, 0.75)

ggplot() +
  geom_sf(data = domain.sf) +
  geom_sf(data = combined.window.sf) +
  geom_circle(aes(x0 = s.win[,1], y0 = 1.05 - s.win[,2], r = exp(u.win)), color = "red", fill = NA) +
  # geom_circle(aes(x0 = s.obs[-obs.win.idx,1], y0 = 1.05 - s.obs[-obs.win.idx,2],
  #                 r = exp(u.nonwin)), color = "grey", fill = NA) +
  theme(axis.title = element_blank())


# --- Fit cond. mark model -----------------------------------------------------
source(here("GlacierBay_Iceberg", "GlacierBay_IcebergCode", "cond.mark.mcmc.R"))
p <- ncol(X.win)
mu.alpha <- rep(0, p + 1)
s2.alpha <- 10000
mu.gamma <- 0
s2.gamma <- 10000
q <- 0.0000001
r <- 10000000
tic()
out.mark.cond <- cond.mark.mcmc(u.win, X.win, out.glm3$beta.0.save, out.glm3$beta.save,
                                n.mcmc, mu.alpha, s2.alpha, mu.gamma, s2.gamma, q, r)
toc() # 12.6 sec

alpha.save <- out.mark.cond$alpha.save[-(1:n.burn),]
gamma.save <- out.mark.cond$gamma.save[-(1:n.burn)]
beta.0.save <- out.mark.cond$beta.0.save[-(1:n.burn)]
beta.save <- out.mark.cond$beta.save[-(1:n.burn),]
s2.u.save <- out.mark.cond$s2.u.save[-(1:n.burn)]

# effective sample size
effectiveSize(beta.0.save)
effectiveSize(beta.save)
effectiveSize(gamma.save)
effectiveSize(alpha.save)
effectiveSize(s2.u.save)

par(mfrow = c(2,5))
plot(alpha.save[,1], type = "l", ylab = "alpha_0")
abline(h = alpha.0, col = "green", lty = 2)

for(i in 1:p){
  plot(alpha.save[,i+1], type = "l", ylab = paste("alpha_", i))
  abline(h = alpha[i], col = "green", lty = 2)
}

plot(gamma.save, type = "l", ylab = "gamma")
abline(h = gamma, col = "green", lty = 2)

plot(beta.0.save, type = "l", ylab = "beta_0")
abline(h = beta.0, col = "green", lty = 2)

for(i in 1:p){
  plot(beta.save[,i], type = "l", ylab = paste("beta_", i))
  abline(h = beta[i], col = "green", lty = 2)
}

plot(s2.u.save, type = "l", ylab = "s2")
abline(h = s2.u, col = "green", lty = 2)


# --- Single stage MCMC --------------------------------------------------------
source(here("GlacierBay_Iceberg", "GlacierBay_IcebergCode", "comp.mark.mcmc.R"))
a.zeta <- 0.0000000001
b.zeta <- 0.0000000001
mu.beta <- rep(0, p)
s2.beta <- 10000
zeta.tune <- 0.01
beta.tune <- 0.00005
tic()
out.comp.mark <- comp.mark.mcmc(s.win, n, u.win, X.full, full.win.idx, win.idx, ds, 
                                a.zeta, b.zeta, zeta.tune,
                                mu.beta, s2.beta, beta.tune, 
                                mu.alpha, s2.alpha, mu.gamma, s2.gamma, q, r,
                                n.mcmc)
toc() # 559.2 sec (n.grid = 500)

# trace plots
n.burn <- 16000
alpha.save <- out.comp.mark$alpha.save[-(1:n.burn),]
gamma.save <- out.comp.mark$gamma.save[-(1:n.burn)]
beta.0.save <- out.comp.mark$beta.0.save[-(1:n.burn)]
beta.save <- out.comp.mark$beta.save[-(1:n.burn),]
s2.u.save <- out.comp.mark$s2.u.save[-(1:n.burn)]

# effective sample size
effectiveSize(beta.0.save)
effectiveSize(beta.save)
effectiveSize(gamma.save)
effectiveSize(alpha.save)
effectiveSize(s2.u.save)

par(mfrow = c(2,5))
plot(alpha.save[,1], type = "l", ylab = "alpha_0")
abline(h = alpha.0, col = "green", lty = 2)

for(i in 1:p){
  plot(alpha.save[,i+1], type = "l", ylab = paste("alpha_", i))
  abline(h = alpha[i], col = "green", lty = 2)
}

plot(gamma.save, type = "l", ylab = "gamma")
abline(h = gamma, col = "green", lty = 2)

plot(beta.0.save, type = "l", ylab = "beta_0")
abline(h = beta.0, col = "green", lty = 2)

for(i in 1:p){
  plot(beta.save[,i], type = "l", ylab = paste("beta_", i))
  abline(h = beta[i], col = "green", lty = 2)
}

plot(s2.u.save, type = "l", ylab = "s2")
abline(h = s2.u, col = "green", lty = 2)


# --- Compare Posterior Inference ----------------------------------------------
par(mfrow = c(2,5))
hist(out.comp.mark$beta.0.save[-(1:n.burn)], prob=TRUE, breaks=60,main="", 
     xlab=bquote(beta[0]), ylim = c(0,10))
lines(density(out.glm3$beta.0.save[-(1:n.burn)], n=1000, adj=2), col="green",
      lwd=2)
abline(v = beta.0, lty = "dashed", lwd = 2, col = "red")

for(i in 1:p){
hist(out.comp.mark$beta.save[-(1:n.burn),i], prob=TRUE, breaks=60,main="", 
     xlab=bquote(beta[i]), ylim = c(0,20))
lines(density(out.mark.cond$beta.save[-(1:n.burn),i], n=1000, adj=2), col="green",
      lwd=2)
abline(v = beta[i], lty = "dashed", lwd = 2, col = "red")
}

hist(out.comp.mark$gamma.save[-(1:n.burn)], prob=TRUE, breaks=60,main="", 
     xlab=bquote(gamma), ylim = c(0,50))
lines(density(out.mark.cond$gamma.save[-(1:n.burn)], n=1000, adj=2), col="green",
      lwd=2)
abline(v = gamma, lty = "dashed", lwd = 2, col = "red")

hist(out.comp.mark$alpha.save[-(1:n.burn),1], prob=TRUE, breaks=60,main="", 
     xlab=bquote(alpha[1]), ylim = c(0,20))
lines(density(out.mark.cond$alpha.save[-(1:n.burn),1], n=1000, adj=2), col="green",
      lwd=2)
abline(v = alpha.0, lty = "dashed", lwd = 2, col = "red")

for(i in 2:(p+1)){
  hist(out.comp.mark$alpha.save[-(1:n.burn),i], prob=TRUE, breaks=60,main="", 
       xlab=bquote(alpha[i]), ylim = c(0,20))
  lines(density(out.mark.cond$alpha.save[-(1:n.burn),i], n=1000, adj=2), col="green",
        lwd=2)
  abline(v = alpha[i-1], lty = "dashed", lwd = 2, col = "red")
}

hist(out.comp.mark$s2.u.save[-(1:n.burn)], prob=TRUE, breaks=60,main="", 
     xlab=bquote(sigma^2), ylim = c(0,50))
lines(density(out.mark.cond$s2.u.save[-(1:n.burn)], n=1000, adj=2), col="green",
      lwd=2)
abline(v = s2.u, lty = "dashed", lwd = 2, col = "red")


# --- N posterior predictive ---------------------------------------------------
N.comp.save <- rep(0, n.mcmc - n.burn)

beta.0.save <- out.comp.full$beta.0.save[-(1:n.burn)]
beta.save <- out.comp.full$beta.save[,-(1:n.burn)]

for(k in 1:(n.mcmc - n.burn)){
  if(k%%1000 == 0){cat(k," ")}
  beta.0.tmp <- beta.0.save[k]
  beta.tmp <- beta.save[,k]
  lam.nonwin.int <- sum(exp(log(ds) + beta.0.tmp + X.nonwin.full%*%beta.tmp))
  N.comp.save[k] <- n + rpois(1, lam.nonwin.int)
};cat("\n")

hist(N.comp.save)
abline(v = N, col = "green", lty = 2, lwd = 2)


# --- Mark Distribution Model Checking -----------------------------------------
# simulate circle centers outside footprints
sim_points <- function(lam, full.coord, win.idx, survey.win, footprint.win, nonwin = TRUE){
  if(nonwin){
    lam <- lam[-win.idx]
    coord <- full.coord[-win.idx,]
  } else{
    lam <- lam[win.idx]
    coord <- full.coord[win.idx,]
  }
  lam.max <- max(lam)
  M <- rpois(1, area.owin(survey.win)*lam.max)
  superpop.full <- rpoint(M, win = survey.win)
  
  if(nonwin){
    is.superpop.nonwin <- !inside.owin(superpop.full$x, superpop.full$y, footprint.win)
    superpop <- cbind(x = superpop.full$x, superpop.full$y)[which(is.superpop.nonwin == TRUE),]
    
  } else{
    superpop <- cbind(superpop.full$x, superpop.full$y)
  }
  
  lam.df <- data.frame(x = coord[,1], y = coord[,2], z = lam)
  lam.rast <- rasterFromXYZ(lam.df)
  superpop.idx <- cellFromXY(lam.rast, superpop)
  lam.superpop <- values(lam.rast)[superpop.idx]
  lam.superpop.mat <- na.omit(cbind(superpop, lam.superpop))
  M <- nrow(lam.superpop.mat)
  
  obs.idx <- rbinom(M, 1, lam.superpop.mat[,3]/lam.max)==1
  s.obs <- lam.superpop.mat[obs.idx,1:2] # total observed points 
  lam.obs <- lam.superpop.mat[obs.idx,3]
  
  return(list(s.obs, lam.obs))
} 

sim.point.list <- list()

for(k in 1:1000){
  lam.temp <- exp(beta.0.save[k + 8000] + X.full%*%beta.save[, k + 8000])
  sim.point <- sim_points(lam.temp, s.full, win.idx, survey.win, combined.window)
  sim.point.list[[k]] <- sim.point[[1]]
}

# simulate radii for each set of simulated points
sim.mark.list <- list()

# test_fun <- function(){
  for(k in 1:1000){
    sim.points <- sim.point.list[[k]]
    sim.idx <- cellFromXY(full.raster, sim.points)
    
    for(q in 1:1000){
      xi.temp <- alpha.save[1, q + 8000] + X.full.marks%*%alpha.save[-1, q + 8000]
      u.nonwin <- rnorm(length(sim.idx), xi.temp[sim.idx], sqrt(s2.save[q + 8000]))
      
      sim.mark.list[[(k-1)*1000 + q]] <- u.nonwin
    }
    
    print(k)
  }
# }

plot(density(sim.mark.list[[1]]), col = alpha(rgb(0,0,0), 0.2), ylim = c(0,0.8))
for(i in sample(1:length(sim.mark.list), 1000)){
  lines(density(sim.mark.list[[i]]), col = alpha(rgb(0,0,0), 0.2))
}
lines(density(u.nonwin), col = "red")
