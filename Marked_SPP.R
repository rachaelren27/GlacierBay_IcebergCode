setwd("/Users/rlr3795/Desktop/GlacierBay_Project")

library(tidyverse)
library(terra)
library(sf)
library(here)
library(spatstat)
library(parallel)
library(doParallel)
library(foreach)
library(raster)
library(tictoc)
library(mvnfast)

set.seed(1234)
# --- Read in data -------------------------------------------------------------
year <- "2007"
date <- "20070618"

## footprints
path <- here("NPS_data", paste0("HARBORSEAL_", year), "footprints")
footprint <- st_read(dsn = path, layer =  paste0("JHI_", date, "_footprint"))
footprint <- st_transform(footprint$geometry, 
                          CRS("+proj=longlat +datum=WGS84"))

## survey polygon
survey.poly <- st_read(dsn = here("cropped_survey_bounds"),
                       layer = "cropped_survey_poly")
survey.poly <- st_transform(survey.poly$geometry, 
                            CRS("+proj=longlat +datum=WGS84"))
survey.poly.mat <- survey.poly[[1]][[1]][[1]]

## icebergs
ice.tiff.path <- here("GlacierBay_Iceberg", "seal_ice", "data", date, "created",
                      "filled_tiffs")
ice.tiff.files <- list.files(ice.tiff.path, pattern = "\\.tif$", full.names = TRUE)

# ice.rast <- rast(ice.tiff.files[1])
# ice.rast[ice.rast != 255] <- NA
# icebergs <- patches(ice.rast, directions = 8)
# 
# ice.poly <- as.polygons(icebergs, dissolve = TRUE)
# ice.poly <- ice.poly[!is.na(ice.poly[[1]]), ]

tic()
ncores <- max(1, parallel::detectCores() - 2)
cl <- makeCluster(ncores)
registerDoParallel(cl)

# parallel loop
ice.list <- foreach(f = ice.tiff.files, .packages = c("terra","sf")) %dopar% {
  ice.rast <- rast(f)
  ice.rast[ice.rast != 255] <- NA
  icebergs <- patches(ice.rast, directions = 8)

  ice.poly <- as.polygons(icebergs, dissolve = TRUE)
  
  areas <- expanse(ice.poly, unit="m")
  
  centers <- centroids(ice.poly)
  center.coords <- crds(centers) 
  
  list(areas = areas, centers = center.coords)
}

stopCluster(cl)
toc()

ice.list <- list()
tic()
for(i in 1:10){
ice.rast <- rast(ice.tiff.files[i])
ice.rast[ice.rast != 255] <- NA
icebergs <- patches(ice.rast, directions = 8)

ice.poly <- as.polygons(icebergs, dissolve = TRUE)

areas <- expanse(ice.poly, unit="m")
ice.idx <- which(areas > 1.6)
areas <- areas[ice.idx]
ice.poly <- ice.poly[ice.idx]

centers <- centroids(ice.poly)
center.coords <- crds(centers) 

ice.list[[i]] <- list(areas = areas, centers = center.coords)
rm(ice.poly)
rm(icebergs)
rm(ice.rast)
gc()
print(i)
}
toc()

areas <- unlist(lapply(ice.list, `[[`, "areas"))
centers <- do.call(rbind, lapply(ice.list, `[[`, "centers"))


# --- Read in covariates -------------------------------------------------------
bath.rast <- raster(here("covariates", "bathymetry.tiff"))
bath.rast <- raster::crop(bath.rast, extent(st_bbox(footprint[[3]])))
bath.rast <- raster::mask(bath.rast, as(footprint[[3]], 'Spatial'))

glac.dist.rast <- raster(here("covariates", "glacier_dist.tiff"))
glac.dist.rast <- raster::crop(glac.dist.rast, extent(st_bbox(footprint[[3]])))
glac.dist.rast <- raster::mask(glac.dist.rast, as(footprint[[3]], 'Spatial'))

ds <- res(bath.rast)[1]*res(bath.rast)[2]

plot(bath.rast)
plot(footprint[[3]], add = TRUE)
points(center.coords)


# --- Construct X matrices -----------------------------------------------------
cell.idx <- which(!is.na(values(bath.rast)))
s.full <- coordinates(bath.rast)[cell.idx,]

X.full <- scale(cbind(values(bath.rast)[cell.idx], values(glac.dist.rast)[cell.idx]))
p <- ncol(X.full)

ice.center.idx.all <- cellFromXY(bath.rast, center.coords.all)
center.coords.bath <- cbind(center.coords.all, values(bath.rast)[ice.center.idx.all])
na.row.idx <- which(is.na(center.coords.bath[,3]))
ice.center.idx <- ice.center.idx.all[-na.row.idx]

center.coords <- center.coords.bath[-na.row.idx, -3]

radii <- sqrt(areas[-na.row.idx]/(2*pi))

n <- length(ice.center.idx)

X.obs <- scale(cbind(values(bath.rast)[ice.center.idx], 
                     values(glac.dist.rast)[ice.center.idx]))


# --- Fit SPP w/ complete likelihood -------------------------------------------
n.mcmc <- 10000
source(here("GlacierBay_Seal", "GlacierBay_SealCode", "spp_win_2D", "spp.comp.mcmc.R"))
theta.tune <- 0.1
beta.tune <- 0.1
tic()
out.comp.full <- spp.comp.mcmc(center.coords, X.obs, X.full, ds, n.mcmc, theta.tune, beta.tune)
toc() # 453.14 sec elapsed (~7.5 min)

# discard burn-in
n.burn <- 0.1*n.mcmc
beta.save.full.lik <- out.comp.full$beta.save[,-(1:n.burn)]
beta.0.save.full.lik <- out.comp.full$beta.0.save[-(1:n.burn)]

# trace plots
layout(matrix(1:2,2,1))
plot(beta.0.save.full.lik,type="l")
matplot(t(beta.save.full.lik),lty=1,type="l")


# --- Fit conditional mark likelihood ------------------------------------------
source(here("GlacierBay_Iceberg", "GlacierBay_IcebergCode", "cond.mark.mcmc.R"))
mu.alpha <- rep(0, p+1)
Sig.alpha <- 100*diag(p+1)
q <- 0.0001
r <- 0.0001
out.mark.cond <- cond.mark.mcmc(log(radii), cbind(rep(1, nrow(X.obs)), X.obs), n.mcmc,
                                mu.alpha, Sig.alpha, q, r)

plot(out.mark.cond$alpha.save[,3], type = "l")
plot(out.mark.cond$s2.save, type = "l")
