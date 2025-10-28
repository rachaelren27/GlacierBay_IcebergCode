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
library(ggforce)
library(geosphere)

set.seed(1234)
# --- Read in data -------------------------------------------------------------
year <- "2007"
date <- "20070618"

## survey polygon
survey.poly <- st_read(dsn = here("cropped_survey_bounds"),
                       layer = "cropped_survey_poly")
survey.poly <- st_transform(survey.poly$geometry, 
                            CRS("+proj=longlat +datum=WGS84"))
survey.poly.mat <- survey.poly[[1]][[1]][[1]]

## footprints
path <- here("NPS_data", paste0("HARBORSEAL_", year), "footprints")
footprint <- st_read(dsn = path, layer =  paste0("JHI_", date, "_footprint"))
footprint <- st_transform(footprint$geometry, 
                          CRS("+proj=longlat +datum=WGS84"))
footprint <- st_intersection(footprint, survey.poly)
# 
# temp.bounds.x <- c(-137.1144, -137.102)
# temp.bounds.y <- c(58.844, 58.853)
# 
# temp.bounds <- owin(xrange = temp.bounds.x, yrange = temp.bounds.y)

## icebergs
ice.tiff.path <- here("GlacierBay_Iceberg", "seal_ice", "data", date, "created",
                      "filled_tiffs")
ice.tiff.files <- list.files(ice.tiff.path, pattern = "\\.tif$", full.names = TRUE)

ice.rast <- rast(ice.tiff.files[1])
ice.rast[ice.rast != 255] <- NA
icebergs <- patches(ice.rast, directions = 8)

ice.poly <- as.polygons(icebergs, dissolve = TRUE)
ice.poly <- ice.poly[!is.na(ice.poly[[1]]), ]

ice.ext <- ext(ice.rast)

plot(footprint[3])
plot(ice.poly, add = TRUE)


ggplot() +
  geom_sf(data = survey.poly) + 
  geom_sf(data = footprint) + 
  xlim(temp.bounds.x) + 
  ylim(temp.bounds.y)
 
ice.list <- list()
tic()
count <- 1
for(i in 1:length(ice.tiff.files)){
  ice.rast <- rast(ice.tiff.files[i])
  
  # check if footprint lies within bounds
  ice.ext <- ext(ice.rast)
  ice.win <- as.owin(c(ice.ext[1], ice.ext[2], ice.ext[3], ice.ext[4]))
  union <- intersect.owin(temp.bounds, ice.win)
  is.nonempty <- !is.empty(union)
    
  if(is.nonempty){
    ice.rast[ice.rast != 255] <- NA
    icebergs <- patches(ice.rast, directions = 8)
    
    ice.poly <- as.polygons(icebergs, dissolve = TRUE)
    
    areas <- expanse(ice.poly, unit="m")
    ice.idx <- which(areas > 1.6)
    areas <- areas[ice.idx]
    ice.poly <- ice.poly[ice.idx]
    
    centers <- centroids(ice.poly)
    center.coords <- crds(centers) 
    
    ice.list[[count]] <- list(areas = areas, centers = center.coords)
    
    rm(ice.poly)
    rm(icebergs)
    rm(ice.rast)
    gc()
    print(i)
    count <- count + 1
  }
}
toc()

areas <- unlist(lapply(ice.list, `[[`, "areas"))
centers <- do.call(rbind, lapply(ice.list, `[[`, "centers"))

rad.m <- sqrt(areas/(2*pi))
rad.decdeg <- rad.m/(111320*cos(centers[,2]*pi/180))
  
ggplot() +
  geom_circle(aes(x0 = centers[,1], y0 = centers[,2], r = rad.decdeg), color = "red", fill = NA) +
  theme(axis.title = element_blank())


# --- Read in covariates -------------------------------------------------------
bath.rast <- raster(here("covariates", "bathymetry.tiff"))
temp.rast <- raster(extent(c(min(centers[,1]), max(centers[,1]), 
                             min(centers[,2]), max(centers[,2]))), crs = crs(bath.rast))
bath.rast <- raster::crop(bath.rast, temp.rast)
# bath.rast <- raster::mask(bath.rast, temp.rast)

# ## write glacier distance raster
# # calculate distance from southern boundary (glacier)
# ggplot() +
#    geom_sf(data = survey.poly) +
#    geom_point(aes(x = -137.1311, y = 58.84288), color = "red") # westmost point
# 
# survey.poly.df <- as.data.frame(survey.poly.mat)
# glacier.poly <- survey.poly.df %>% filter((V2 <= 58.84288) & (V1 <= -137.1035))
# 
# ggplot() +
#   geom_sf(data = survey.poly) +
#   geom_line(data = glacier.poly, aes(x = V1, y = V2), col = "blue")
# 
# # calculate glacier distance
# # seal.glac.dist <- dist2Line(seal.mat, glacier.poly) # in meters
# 
# bath.survey.idx <- which(!is.na(values(bath.rast)))
# full.coord <- xyFromCell(bath.rast, bath.survey.idx)
# 
# full.glac.dist <- dist2Line(s.full, glacier.poly) # takes a while
# 
# glac.dist.df <- data.frame(x = s.full[,1], y = s.full[,2],
#                            z = full.glac.dist[,1])
# glac.dist.rast <- rasterFromXYZ(glac.dist.df)
# writeRaster(glac.dist.rast, filename = "glacier_dist.tiff", format = "GTiff")

glac.dist.rast <- raster(here("covariates", "glacier_dist.tiff"))
glac.dist.rast <- raster::crop(glac.dist.rast, 
                               extent(c(min(centers[,1]), max(centers[,1]), 
                                        min(centers[,2]), max(centers[,2]))),
                                        crs = crs(bath.rast))
# glac.dist.rast <- raster::mask(glac.dist.rast, as(footprint[[3]], 'Spatial'))

ds <- res(bath.rast)[1]*res(bath.rast)[2]


# --- Construct X matrices -----------------------------------------------------
s.full <- coordinates(bath.rast)[cell.idx,]
X.full <- scale(na.omit(cbind(values(bath.rast), values(glac.dist.rast))))

p <- ncol(X.full)

ice.center.idx <- cellFromXY(bath.rast, centers)
centers.cov <- cbind(centers, values(bath.rast)[ice.center.idx],
                     values(glac.dist.rast)[ice.center.idx])
na.idx <- which(apply(is.na(centers.cov), 1, any))

if(length(na.idx) > 0){
  obs.idx <- ice.center.idx[-na.idx]
  centers <- centers.cov[-na.idx, 1:2]
}else{
  obs.idx <- ice.center.idx
}

X.obs <- scale(cbind(values(bath.rast)[obs.idx], 
                     values(glac.dist.rast)[obs.idx]))


radii <- sqrt(areas/(2*pi))

n <- length(obs.idx)

plot(bath.rast)
points(x = centers[,1], y = centers[,2])


# --- Fit SPP w/ complete likelihood -------------------------------------------
n.mcmc <- 10000
source(here("GlacierBay_Seal", "GlacierBay_SealCode", "spp_win_2D", "spp.comp.mcmc.R"))
theta.tune <- 0.01
beta.tune <- 0.005
tic()
out.comp.full <- spp.comp.mcmc(centers, X.obs, X.full, ds, n.mcmc, theta.tune, beta.tune)
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
