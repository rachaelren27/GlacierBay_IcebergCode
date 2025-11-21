library(tidyverse)
library(viridis)
library(raster)
library(spatstat)
library(sf)
library(ggforce)

set.seed(677)

# --- Construct domain ---------------------------------------------------------
domain.range <- c(0,1)
domain.win <- owin(xrange = domain.range, yrange = domain.range)
domain.sf <- st_as_sfc(as.polygonal(domain.win)) %>% st_sf()

# construct grid
n.grid <- 100

x.full <- seq(domain.range[1], domain.range[2], length.out = n.grid)
y.full <- seq(domain.range[1], domain.range[2], length.out = n.grid)
s.full <- expand.grid(x = x.full, y = y.full)

ds <- (x.full[2] - x.full[1])^2


# --- Define variables and covariates ------------------------------------------
# simulate GP covariate
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

gp.df <- GP_2D_sim(n.grid, domain.range, mu = 0,
                   sig2 = 0.5, phi = 0.15, nu = 1e-5)
gp.rast <- rasterFromXYZ(gp.df)

ggplot(gp.df, aes(x = x, y = y, fill = z)) +
  geom_tile() +
  scale_fill_viridis() +
  coord_fixed() +
  theme_minimal()

# construct design matrix
X.full <- scale(gp.df)
p <- ncol(X.full)


# --- Simulate from HPP --------------------------------------------------------
lambda <- 1000
n.y <- rpois(1, lambda)
s.obs <- cbind(x = runif(n.y), y = runif(n.y))

obs.idx <- cellFromXY(gp.rast, s.obs)
X.obs <- X.full[obs.idx,]

plot(x = s.obs[,1], y = s.obs[,2], pch = 19)


# --- Simulate radii -----------------------------------------------------------
beta.0 <- -5
beta <- c(-0.5, 0.5, 1)

mu.u <- c(beta.0 + X.obs%*%beta)
s2.u <- 0.5

u <- rnorm(length(mu.u), mu.u, s2.u)
r <- exp(u)

ggplot() +
  geom_sf(data = domain.sf) +
  geom_circle(aes(x0 = s.obs[,1], y0 = 1 - s.obs[,2], r = r), color = "red", fill = NA) +
  theme(axis.title = element_blank())

marked.s.obs <- cbind(s.obs, r)

# --- Thin PP (Matern I) -------------------------------------------------------
# remove overlapping points
s.obs.dist <- as.matrix(dist(s.obs))

r.sum <- outer(r, r, "+")

overlap.mat <- s.obs.dist < r.sum & s.obs.dist > 0

overlap <- apply(overlap.mat, 1, any)

# additional thinning
p.0 <- 0.9
thin <- rbinom(n.y, 1, 1 - p.0)
keep <- overlap + thin == 0

marked.s.obs.thinned <- marked.s.obs[keep,]
r.thinned <- r[keep]

# plot thinned process
ggplot() +
  geom_sf(data = domain.sf) +
  geom_circle(aes(x0 = marked.s.obs.thinned[,1], y0 = 1 - marked.s.obs.thinned[,2],
                  r = r.thinned), color = "red", fill = NA) +
  theme(axis.title = element_blank())


# --- Thin PP (Matern II) ------------------------------------------------------
# remove overlapping points
s.obs.dist <- as.matrix(dist(s.obs))
s.obs.dist[lower.tri(s.obs.dist)] <- 0

r.sum <- outer(r, r, "+")
r.sum[lower.tri(r.sum)] <- 0

overlap.mat <- s.obs.dist < r.sum & s.obs.dist > 0

no.overlap <- !apply(overlap.mat, 1, any)

# additional thinning
p.0 <- 0.9
no.thin <- rbinom(n.y, 1, p.0)
keep <- no.overlap*no.thin == 1

marked.s.obs.thinned <- marked.s.obs[keep,]
r.thinned <- r[keep]

# plot thinned process
ggplot() +
  geom_sf(data = domain.sf) +
  geom_circle(aes(x0 = marked.s.obs.thinned[,1], y0 = 1 - marked.s.obs.thinned[,2],
                  r = r.thinned), color = "red", fill = NA) +
  theme(axis.title = element_blank())
