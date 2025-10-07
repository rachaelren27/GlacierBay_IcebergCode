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
library(invgamma)

set.seed(677)

# --- Construct domain --------------------------------------------------------
domain.range <- c(0,1.05)

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
n.grid <- 100

x.full <- seq(domain.range[1], domain.range[2], length.out = n.grid)
y.full <- seq(domain.range[1], domain.range[2], length.out = n.grid)
s.full <- expand.grid(x = x.full, y = y.full)

ds <- (x.full[2] - x.full[1])^2

# --- Simulate GP covariate ----------------------------------------------------
GP_2D_sim <- function(n.grid, domain.range, mu, sig2, phi, nu){
  s1 <- seq(domain.range[1], domain.range[2], length.out = n.grid)
  s2 <- seq(domain.range[1], domain.range[2], length.out = n.grid)
  locs <- expand.grid(s1, s2)
  
  L <- n.grid^2
  
  D.full <- as.matrix(dist(locs))
  Sig <- nu*diag(L) + sig2*exp(-(D.full^2)/(phi^2)/2)
  
  y.sim <- c(rmvn(1, rep(mu, L), Sig))
  
  return(y.sim)
}

mu <- 0
sig2 <- 0.5
phi <- 0.2
nu <- 0.0001

s1 <- seq(domain.range[1], domain.range[2], length.out = n.grid)
s2 <- seq(domain.range[1], domain.range[2], length.out = n.grid)
locs <- expand.grid(s1, s2)

y.sim <- GP_2D_sim(n.grid, domain.range, mu, sig2, phi, nu)
out.df <- data.frame(cbind(s1 = locs[,1], s2 = locs[,2], y = y.sim))

ggplot() +
  geom_tile(data = out.df, aes(x = s1, y = s2, fill = y)) +
  scale_fill_viridis() +
  labs(fill = "y", x = "", y = "") +
  coord_fixed() +
  theme_minimal()

# --- Simulate from windowed IPP -----------------------------------------------
beta.0 <- 6
beta <- c(-0.5, -0.5, 1)

X.full <- scale(cbind(locs[,1], locs[,2], y.sim))

lam.full <- exp(beta.0 + X.full%*%beta)
lam.max <- max(lam.full)

full.df <- as.data.frame(cbind(x = locs[,1], y = locs[,2], z = lam.full))
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
obs.win.idx <- (1:N)[obs.win] # observed -> windowed
n <- length(obs.win.idx)

full.win <- inside.owin(s.full[,1], s.full[,2], combined.window)
full.win.idx <- (1:n.grid^2)[full.win]

s.win <- s.obs[obs.win.idx,]
X.win <- X.obs[obs.win.idx,]

X.win.full <- X.full[full.win.idx,]
X.nowin.full <- X.full[-full.win.idx,]

domain.sf <- st_as_sfc(as.polygonal(domain)) %>% st_sf()
combined.window.sf <- st_as_sfc(as.polygonal(combined.window)) %>% st_sf()
ggplot() +
  geom_sf(data = domain.sf) +
  geom_sf(data = combined.window.sf) +
  geom_point(aes(x = s.win[,1], y = 1.05 - s.win[,2]), col = "red", size = 0.1) +
  # geom_point(aes(x = s.obs[-obs.win.idx,1], y = s.obs[-obs.win.idx,2]), size = 0.5) +
  theme(axis.title = element_blank())

# --- Fit SPP w/ complete likelihood -------------------------------------------
n.mcmc <- 10000
source(here("GlacierBay_Seal", "GlacierBay_SealCode", "spp_win_2D", "spp.comp.mcmc.R"))
tic()
out.comp.full=spp.comp.mcmc(s.win, X.win, X.win.full, ds, n.mcmc, 0.1, 0.01)
toc() # 385 sec

# discard burn-in
n.burn <- 0.1*n.mcmc
beta.0.save <- out.comp.full$beta.0.save[-(1:n.burn)]
beta.save <- out.comp.full$beta.save[,-(1:n.burn)]

# trace plot
layout(matrix(1:2,2,1))
plot(beta.0.save,type="l")
abline(h=beta.0,col=rgb(0,1,0,.8),lty=2)
matplot(t(beta.save),lty=1,type="l")
abline(h=beta,col=rgb(0,1,0,.8),lty=2)

# --- Generate marks -----------------------------------------------------------
alpha.0 <- -6
alpha <- c(-0.1, 0.1, -0.1, -0.3)
X.full.marks <- scale(cbind(X.full, lam.full))
win.idx <- cellFromXY(full.raster, s.win) # full -> windowed
X.win.marks <- X.full.marks[win.idx,]

u.win <- rnorm(n, alpha.0 + X.win.marks%*%alpha, 0.75)

ggplot() +
  geom_sf(data = domain.sf) +
  geom_sf(data = combined.window.sf) +
  geom_circle(aes(x0 = s.win[,1], y0 = s.win[,2], r = exp(u.win)), color = "red", fill = NA) +
  # geom_point(aes(x = s.obs[-obs.win.idx,1], y = s.obs[-obs.win.idx,2]), size = 0.5) +
  theme(axis.title = element_blank())


# --- Fit cond. mark model -----------------------------------------------------
source(here("GlacierBay_Iceberg", "GlacierBay_IcebergCode", "cond.mark.mcmc.R"))
p.marks <- ncol(X.full.marks)
mu.alpha <- rep(0, p.marks + 1)
Sig.alpha <- 10*diag(p.marks + 1)
q <- 0.00001
r <- 1000000
out.mark.cond <- cond.mark.mcmc(u.win, cbind(rep(1, nrow(X.win.marks)), X.win.marks),
                                n.mcmc, mu.alpha, Sig.alpha, q, r)

plot(out.mark.cond$alpha.save[-(1:n.burn),4], type = "l",)
plot(out.mark.cond$s2.save, type = "l")
