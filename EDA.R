library(terra)
library(here)
library(car)
library(sf)
library(tidyverse)
library(spatstat)
library(mvnfast)
library(viridis)
library(raster)
library(igraph)
library(concaveman)
library(lpSolve)
library(sp)
library(ggforce)

# --- Simulate jittery ovals ---------------------------------------------------
# create study domain
x.domain <- c(0,1)
y.domain <- c(0,1)

draw_oval <- function(x.domain, y.domain,
                      a.mu, a.sigma, b.mu, b.sigma,
                      npoints = 20, eps.sigma = 0.01,
                      border = "black", col = NA){
  
  x.center <- runif(1, x.domain[1], x.domain[2])
  y.center <- runif(1, y.domain[1], y.domain[2])
  
  theta <- seq(0, 2*pi, length.out = npoints)
  a <- rnorm(1, a.mu, a.sigma)
  b <- rnorm(1, b.mu, b.sigma)
  x <- a * cos(theta)
  y <- b * sin(theta)

  angle <- runif(1, 0, 359)
  ang <- angle * pi / 180
  x.rot <- x*cos(ang) - y*sin(ang)
  y.rot <- x*sin(ang) + y*cos(ang)

  x.final <- x.center + x.rot
  y.final <- y.center + y.rot
  
  dx <- x.final - x.center
  dy <- y.final - y.center
  
  dist <- sqrt(dx^2 + dy^2)
  dx.unit <- dx / dist
  dy.unit <- dy / dist

  eps.vec <- rnorm(npoints, mean = 0, sd = eps.sigma)
  
  x.final <- x.final + eps.vec * dx.unit
  y.final <- y.final + eps.vec * dy.unit
  
  polygon(x.final, y.final, border = border, col = col)
}

# plot
n.oval <- 10
plot(x = x.centers, y = y.centers, xlim = x.domain, ylim = y.domain, col = NA)
for(i in 1:n.oval){
  draw_oval(x.domain = c(0,1), y.domain = c(0,1),
            a.mu = 0.1, a.sigma = 0.02, b.mu = 0.1, b.sigma = 0.03)
}


# --- Read in iceberg data -----------------------------------------------------
ice <- rast("ice_filled.tif")

ice[ice != 255] <- NA

icebergs <- patches(ice, directions=8)

icebergs <- as.polygons(icebergs, dissolve=TRUE)

areas <- expanse(icebergs, unit="m")

iceberg <- icebergs[which(areas == max(areas)),]

centers <- centroids(icebergs)
center.coords <- crds(centers)

ggplot(data = iceberg.dat, aes(x = x, y = y, color = log(areas))) + 
  geom_point() + 
  scale_color_gradient()

# --- Fit statistical ellipse --------------------------------------------------
iceberg.coord <- st_coordinates(st_as_sf(iceberg, as_points=FALSE,
                                             merge=FALSE))[,1:2]

ell <- dataEllipse(x = iceberg.coord[,1], y = iceberg.coord[,2], levels = 0.8)

centroid <- st_coordinates(st_centroid(st_as_sf(iceberg)))

plot(iceberg.coord, asp=1)
lines(ell, col="blue", lwd=2)
points(centroid[1], centroid[2], col = "red")


# --- Marked Point Process -----------------------------------------------------
footprint <- ext(ice)
footprint.win <- owin(xrange = c(footprint[1], footprint[2]), 
                      yrange = c(footprint[3], footprint[4]))
iceberg.mpp <- ppp(x = center.coords[,1], y = center.coords[,2], 
                   window = footprint.win, marks = areas)

# variogram
plot(markvario(iceberg.mpp, "isotropic"))

# checking dependence between marks and locations
D.obs <- as.matrix(dist(center.coords))
bins <- seq(10, 320, by = 10)

calc_exp_var <- function(bins, areas, D.obs){
  results <- data.frame(
    radius = bins,
    mean = NA,
    sd = NA
  )
  
  for(i in 1:length(bins)){
    r <- bins[i]
    
    has.neighbor <- apply(D.obs <= r & D.obs > 0, 1, any)
    
    select.areas <- areas[has.neighbor]
    
    if(length(select.areas > 0)){
      results$mean[i] <- mean(select.areas)
      results$sd[i] <- sd(select.areas)
    }
  }
  return(results)
}

exp.var.out <- calc_exp_var(bins, areas, D.obs)

plot(x = exp.var.out[,1], y = exp.var.out[,2], type = 'l')
plot(x = exp.var.out[,1], y = exp.var.out[,3], type = 'l')


# --- Simulate Overlapping Circles ---------------------------------------------
## Simulate IPP points
# simulate GP covariate
GP_2D_sim <- function(n, mu, sig2, phi, nu){
  s1 <- seq(0, 1, length.out = n)
  s2 <- seq(0, 1, length.out = n)
  locs <- expand.grid(s1, s2)
  n.full <- nrow(locs)
  
  D.full <- as.matrix(dist(locs))
  Sig <- nu*diag(n.full) + sig2*exp(-(D.full^2)/(phi^2)/2)
  
  y.sim <- c(rmvn(1, rep(mu, n.full), Sig))
  
  return(y.sim)
}

n <- 50

mu <- 0
sig2 <- 1
phi <- 0.05
nu <- 0.001

y.sim <- GP_2D_sim(n, mu, sig2, phi, nu)
s1 <- seq(0, 1, length.out = n)
s2 <- seq(0, 1, length.out = n)
locs <- expand.grid(s1, s2)
out.df <- data.frame(cbind(s1 = locs[,1], s2 = locs[,2], y = y.sim))

ggplot() +
  geom_tile(data = out.df, aes(x = s1, y = s2, fill = y)) +
  scale_fill_viridis() +
  labs(fill = "y", x = "", y = "") +
  coord_fixed() +
  theme_minimal()

# simulate from IPP
beta.0 <- 4
beta <- 2

X.full <- y.sim

lam.full <- exp(beta.0 + beta*X.full)
lam.max <- max(lam.full)

full.df <- as.data.frame(cbind(x = s1, y = s2, z = lam.full))
full.raster <- rasterFromXYZ(full.df)

M <- rpois(1, lam.max) 
x.superpop <- runif(M, x.domain[1], x.domain[2])
y.superpop <- runif(M, y.domain[1], y.domain[2])
s.superpop <- cbind(x.superpop, y.superpop)
X.superpop <- y.sim[cellFromXY(full.raster, cbind(x = x.superpop, y = y.superpop))]
lam.superpop <- exp(beta.0 + beta*X.superpop)

obs.idx <- rbinom(M, 1,lam.superpop/lam.max) == 1
s.obs <- s.superpop[obs.idx,] # total observed points 
X.obs <- X.superpop[obs.idx] 
lam.obs <- lam.superpop[obs.idx]
N <- nrow(s.obs) 

# plot with circles
plot(s.obs[,1], s.obs[,2])
r <- 0.01
symbols(s.obs[,1], s.obs[,2], circles = rep(r, N), inches = FALSE,
        add = TRUE, bg = "black")

## Join overlapping circles
# building adjacency matrix
adj <- matrix(FALSE, N, N)
for(i in 1:N){
  for(j in 1:N){
    if(i != j){
      d <- sqrt((s.obs[i,1] - s.obs[j,1])^2 + (s.obs[i,2] - s.obs[j,2])^2)
      if(d <= 2*r) adj[i,j] <- TRUE
    }
  }
}

## join by adjacency matrix
g <- graph_from_adjacency_matrix(adj, mode = "undirected")
group <- components(g)$membership
plot(s.obs[,1], s.obs[,2], col = group)

groups <- split(1:vcount(g), clusters(g)$membership)

## draw concave hulls
theta <- seq(0, 2*pi, length.out = 100)

plot(NA, xlim = c(0,1), ylim = c(0,1))
for(g in groups){
  circle_pts <- do.call(rbind, lapply(g, function(i) {
    cbind(s.obs[i,1] + r * cos(theta), s.obs[i,2] + r * sin(theta))
  }))
  hull <- concaveman(circle_pts, concavity = 2, length_threshold = 0)
  polygon(hull, border = "darkblue", lwd = 2, col = adjustcolor("lightgreen", 0.3))
}


# --- Fit Observed Data --------------------------------------------------------
## circle covering using hexagonal grid
n.classes <- 6
size.classes <- c(14, 50, 100, 500, 1000, 5000, Inf)

find_centers <- function(n.classes, size.classes, icebergs){
  centers <- list()
  icebergs.sf <- st_as_sf(icebergs)
  
  for(i in 1:n.classes){
    icebergs.sub <- icebergs.sf[which(areas >= size.classes[i] & areas < size.classes[i + 1]),]
    
    r <- sqrt(size.classes[i]/(2*pi))
    cellsize <- sqrt(3) * r
    
    centers.sub <- list()
    for(j in 1:nrow(icebergs.sub)){
      iceberg_sp <- as(icebergs.sub[j,], "Spatial")
      
      # Hexagonal sampling inside bounding box
      hex_pts <- spsample(iceberg_sp, type="hexagonal", cellsize = cellsize, iter = 100)
      
      # Convert to sf
      centers_sf <- st_as_sf(hex_pts)
      
      # 2) Keep only centers whose circle intersects polygon (within r)
      keep <- lengths(st_is_within_distance(centers_sf, st_as_sf(icebergs.sub[j,]), dist = r)) > 0
      centers.sub[[j]] <- centers_sf[keep, ]
    }
    
    centers[[i]] <- do.call(rbind, lapply(seq_along(centers.sub), function(j) {
      cbind(centers.sub[[j]], poly.id = j)}))
  }
  
  centers.df <- do.call(rbind, lapply(seq_along(centers), function(i) {
    cbind(centers[[i]], class.id = i)}))
  
  return(centers.df)
}

centers.out <- find_centers(n.classes, size.classes, icebergs)


# 3) Plot
radii <- sqrt(size.classes[1:n.classes]/(2*pi))
  
ggplot() +
  geom_sf(data = st_as_sf(icebergs), fill = "grey85", color = "black") +
  # geom_sf(data = centers_sf, color = "blue", size = 2) +
  # geom_point(data = st_coordinates(centers_sf) |> as.data.frame(),
  #            aes(X, Y), color = "blue") +
  # coord_equal() +
  # Draw circles around centers
  geom_sf(data = st_buffer(centers.out[[1]], dist = radii[1]), color = "red", fill = NA) + 
  geom_sf(data = st_buffer(centers.out[[2]], dist = radii[2]), color = "blue", fill = NA) +
  geom_sf(data = st_buffer(centers.out[[3]], dist = radii[3]), color = "green", fill = NA) +
  geom_sf(data = st_buffer(centers.out[[4]], dist = radii[4]), color = "purple", fill = NA)
