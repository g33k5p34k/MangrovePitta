library(PleistoDist)
library(basemaps)
library(ggplot2)
library(sf)
library(terra)
library(qgisprocess)
library(ENMeval)
library(raster)
library(ggpubr)
library(dplyr)
library(magick)
library(rlist)
library(nngeo)
library(dismo)
library(sfheaders)
library(spatialEco)
library(spThin)
library(tidyterra)
library(doParallel)
library(foreach)
library(smoothr)
library(windfetch)
library(purrr)
library(exactextractr)

###################################################

## Run PleistoDist and account for landmass subsidence
setwd("E:/Pittas/nichemodelling/GEBCO_10_Jul_2024_46870ad1c37c/")

#generate interval file
getintervals_time(time = 1000, intervals = 1000, outdir = "E:/Pittas/nichemodelling")

#bathymetry <- terra::rast("gebco_2023_SEA_clipped_reprojected_downsampled.asc") %>%
#terra::aggregate(fact=10)
#terra::writeRaster(bathymetry,"gebco_2023_SEA_clipped_reprojected_downsampled_10km_1.asc",overwrite=TRUE,filetype="AAIGrid")
PleistoDist::makemaps(inputraster="gebco_2023_SEA_clipped_reprojected_downsampled.asc",
         epsg = 6933,
         intervalfile = "E:/Pittas/nichemodelling/intervals.csv",
         offset = -0.25,
         outdir = "E:/Pittas/test/")

#read intervalfile, bathymetry file, and clipping file
intervalfile <- read.csv("E:/Pittas/nichemodelling/intervals.csv", header = TRUE)
basemap <- terra::rast("gebco_2023_SEA_clipped_reprojected_downsampled.asc")
#sea_clip <- terra::vect("E:/Pittas/nichemodelling/SEA_clipped.shp")
#basecrop <- terra::crop(basemap,sea_clip)

#prepare a shapefile of the present day land extent with some allowance for intertidal
presentday <- terra::clamp(basemap, lower = -3, value = FALSE)
#presentday_shp <- terra::as.polygons(presentday)

#for interval0, generate slope map
writeRaster(presentday,filename="E:/Pittas/nichemodelling/raster_topo_new/topo_interval0.tif",gdal=c("COMPRESS=DEFLATE"),overwrite=TRUE)
presentslope <- terrain(basemap, v = "slope", neighbors = 8, unit = "degrees")
writeRaster(presentslope, filename = "E:/Pittas/nichemodelling/raster_slope/slope_interval0.tif",gdal=c("COMPRESS=DEFLATE"),overwrite=TRUE)

#generate slope maps for all intervals
for (i in 1:1000) {
  sealvl <- dplyr::filter(intervalfile,Interval==i)
  #zero basemap to interval sea level and set any pixels < -10 m to NA
  basenow = basemap - sealvl$MeanDepth - (-0.25 * sum(intervalfile$TimeInterval[1:i]))
  slopenow <- terrain(basenow, v = "slope", neighbors = 8, unit = "degrees")
  writeRaster(slopenow,filename = paste0("E:/Pittas/nichemodelling/raster_slope/slope_interval",i,".tif"),gdal=c("COMPRESS=DEFLATE"),overwrite=TRUE)
  gc()
}


#generate raster_topo maps for all intervals
for (i in 1:1000) {
  sealvl <- dplyr::filter(intervalfile,Interval==i)
  #zero basemap to interval sea level and set any pixels < -3 m to NA
  basenow = terra::clamp(basemap- (-0.25 * sum(intervalfile$TimeInterval[1:i])) - sealvl$MeanDepth, lower = -3, value = FALSE)
  writeRaster(basenow,filename = paste0("E:/Pittas/nichemodelling/raster_topo_new/topo_interval",i,".tif"),gdal=c("COMPRESS=DEFLATE"),overwrite=TRUE)
  gc()
}

for (i in 601:1000) {
  landshape <- sf::st_read(paste0("E:/Pittas/nichemodelling/shapefile/interval",i,".shp"))
  landvect_noholes <- smoothr::fill_holes(landshape,threshold=units::set_units(100000000000, km^2))
  landvect_noholes <- sf::st_union(landvect_noholes)
  landvect_noholes <- smoothr::fill_holes(landvect_noholes,threshold=units::set_units(10, km^2))
  sf::st_write(landvect_noholes,paste0("E:/Pittas/nichemodelling/shapefile_clean/interval",i,".shp"),overwrite=TRUE)
}

## function for generating point grids, limited to a 10km buffer from the shore

#read baseline grid file
grid <- sf::st_read("E:/Pittas/nichemodelling/basegrid.shp")
#set CRS
epsg <- "epsg:6933"
makegrids <- function(interval) {
  #read shapefile of land extent
  # landshape <- sf::st_read(paste0("E:/Pittas/nichemodelling/shapefile_clean/interval",interval,".shp")) %>%
  #   dplyr::select(-FID)

  #create shapefile of land extent, allowing for intertidal zone
  raster_topo_new <- terra::rast(paste0("E:/Pittas/nichemodelling/raster_topo_new/topo_interval",interval,".tif"))
  #flatten out topo raster
  m <- c(-4,10000,1)
  convmat <- matrix(m,ncol=3)
  raster_topo_flat <- terra::classify(raster_topo_new,convmat)
  shapefile_with_intertidal <- terra::as.polygons(raster_topo_flat)

  #read shapefile of land extent
  #landshape <- sf::st_read(paste0("E:/Pittas/nichemodelling/shapefile/interval",interval,".shp"))
  #close all internal holes so distance is only calculated from shoreline
  landvect_noholes <- smoothr::fill_holes(sf::st_as_sf(shapefile_with_intertidal),threshold=units::set_units(100000000000, km^2))
  landvect_noholes <- sf::st_union(landvect_noholes)
  landvect_noholes <- smoothr::fill_holes(landvect_noholes,threshold=units::set_units(10, km^2))
  sf::st_write(landvect_noholes,paste0("E:/Pittas/nichemodelling/shapefile_with_intertidal/interval",interval,".shp"),overwrite=TRUE)

  landvect_noholes <- sf::st_read(paste0("E:/Pittas/nichemodelling/shapefile_with_intertidal/interval",interval,".shp"))
  landvect_noholes <- smoothr::fill_holes(landvect_noholes,threshold=units::set_units(100000000000, km^2))
  #draw linestring of map boundary
  boundbox <- as.numeric(sf::st_bbox(landvect_noholes))
  boundary <- sf::st_sfc(sf::st_linestring(c(sf::st_point(c(boundbox[1],boundbox[2])),sf::st_point(c(boundbox[3],boundbox[2])),sf::st_point(c(boundbox[3],boundbox[4])),sf::st_point(c(boundbox[1],boundbox[4])),sf::st_point(c(boundbox[1],boundbox[2])))), crs=epsg)

  landline <- sf::st_cast(landvect_noholes,"MULTILINESTRING")

  #use the boundary line to clip the landline vector
  landline <- sf::st_difference(landline,boundary)
  landsimp <- sf::st_simplify(landline)
  sf::st_write(landsimp,paste0("E:/Pittas/nichemodelling/landline/landline_interval",interval,".shp"))

  landbuff <- qgisprocess::qgis_run_algorithm(
    'gdal:buffervectors',
    INPUT=paste0("E:/Pittas/nichemodelling/landline/landline_interval",interval,".shp"),
    DISTANCE=30000,
    DISSOLVE=TRUE
  )

  landbuff <- sf::st_as_sf(landbuff)

  landmask <- qgisprocess::qgis_run_algorithm(
    'native:intersection',
    INPUT=landbuff,
    OVERLAY=shapefile_with_intertidal
  )

  landmask <- sf::st_as_sf(landmask)
  sf::st_write(landmask,paste0("E:/Pittas/nichemodelling/shoremask/shoremask_interval",interval,".shp"))

  #load masking shapefile
  # maskingshape <- sf::st_read("E:/Pittas/nichemodelling/mangrove_maskinglayer_final_final.shp")
  # maskingshape %>%
  #   dplyr::select(-fid) %>%
  #   dplyr::select(-FID_2)

  # #mask landshape file
  # landmasked <- qgisprocess::qgis_run_algorithm(
  #   'native:difference',
  #   INPUT=landshape,
  #   OVERLAY=maskingshape
  # )

  #landmasked <- st_as_sf(landmasked)

  #subsample grid file based on land extent
  landgrid <- qgisprocess::qgis_run_algorithm(
    'native:extractbylocation',
    INPUT=grid,
    INTERSECT=landmask
  )

  #convert qgis output into sf format
  landgrid <- sf::st_as_sf(landgrid)

  landgrid <- dplyr::mutate(landgrid,exposure=NA)
  sf::st_write(landgrid,paste0("E:/Pittas/nichemodelling/landgrids/landgrid_interval",interval,".shp"))
  gc()
}

for (i in 601:800) {
  makegrids(i)
}

#convert landline to points, sampled at one point per 1km
for (i in 601:800) {
  landline <- sf::st_read(paste0("E:/Pittas/nichemodelling/landline/landline_interval",i,".shp"))
  landline <- sf::st_geometry(sf::st_cast(sf::st_union(landline),"LINESTRING"))
  shorepoints <- sf::st_line_sample(landline,density=0.001)
  shorepoints <- sf::st_cast(sf::st_union(shorepoints),"POINT")
  sf::st_write(shorepoints,paste0("E:/Pittas/nichemodelling/shorepoints/shorepoints_interval",i,".shp"),overwrite=TRUE)
}

#original shore exposure calculation function I wrote. PLEASE DON'T USE THIS
# exposure <- function(interval) {
#
#   #set up parallelisation
#   cores = parallel::detectCores()
#   cl <- parallel::makeCluster(cores[1]-4)
#   doParallel::registerDoParallel(cl)
#
#   template <- terra::rast(paste0("E:/Pittas/nichemodelling/raster_flat/interval",interval,".tif"))
#   landgrid <- sf::st_read(paste0("E:/Pittas/nichemodelling/landgrids/landgrid_interval",interval,".shp"))
#   landline <- sf::st_read(paste0("E:/Pittas/nichemodelling/landline/landline_interval",interval,".shp"))
#   landline <- smoothr::smooth(landline,method = "ksmooth", smoothness = 2)
#   boundbox <- as.numeric(st_bbox(landline))
#   boundary <- st_sfc(st_linestring(c(st_point(c(boundbox[1],boundbox[2])),st_point(c(boundbox[3],boundbox[2])),st_point(c(boundbox[3],boundbox[4])),st_point(c(boundbox[1],boundbox[4])),st_point(c(boundbox[1],boundbox[2])))), crs=epsg)
#
#   #set CRS
#   epsg <- "epsg:6933"
#
#   #iterate through points and calculate exposure index
#   exposurevect <- foreach::foreach (i=1:20000, .combine = 'rbind', .errorhandling = "pass") %dopar% {
#     #nrow(landgrid)
#     #exposurevect <- for (i in 1:nrow(landgrid)) {
#     #identify point
#     keypoint <- landgrid$geometry[i]
#
#     #extract coordinates of selected point
#     x_coord <- sf::st_coordinates(keypoint)[1]
#     y_coord <- sf::st_coordinates(keypoint)[2]
#
#     #identify nearest point on coastline
#     #nearestlines <- sf::st_nearest_points(keypoint,landline)
#     nearestline_index <- sf::st_nearest_feature(keypoint,landline)
#     nearestline <- landline$geometry[nearestline_index]
#     shoreintersectline <- sf::st_nearest_points(keypoint,nearestline)
#     #find index of shortest line
#     #shortest_line <- which.min(sf::st_length(nearestlines))
#     #extract coordinates of nearest shore point, which we will use to calculate the extended line
#     secondpoint <- sf::st_cast(shoreintersectline,"POINT")[2]
#
#
#     linegradient <- (y_coord - sf::st_coordinates(secondpoint)[2])/(x_coord - sf::st_coordinates(secondpoint)[1])
#     if ((linegradient < 0 & linegradient > -Inf & sf::st_coordinates(secondpoint)[1] < x_coord) | (linegradient > 0 & linegradient < Inf & linegradient & sf::st_coordinates(secondpoint)[1] > x_coord)) { #line is moving towards NE or NW...
#       longline <- sf::st_sfc(sf::st_linestring(rbind(sf::st_coordinates(keypoint),sf::st_point(c((boundbox[4]-(y_coord-linegradient*x_coord))/linegradient,boundbox[4])))),crs=epsg) #extend existing line towards ymax
#     } else if ((linegradient < 0 & linegradient > -Inf & sf::st_coordinates(secondpoint)[1] > x_coord) | (linegradient > 0 & linegradient < Inf & linegradient & sf::st_coordinates(secondpoint)[1] < x_coord)) { #if line is moving towards SE or SW
#       longline <- sf::st_sfc(sf::st_linestring(rbind(sf::st_coordinates(keypoint),sf::st_point(c((boundbox[2]-(y_coord-linegradient*x_coord))/linegradient,boundbox[2])))),crs=epsg)
#     } else if (abs(linegradient) == Inf & y_coord < sf::st_coordinates(secondpoint)[2]) {
#       longline <- sf::st_sfc(sf::st_linestring(rbind(sf::st_coordinates(keypoint),sf::st_point(c(x_coord,boundbox[4])))), crs=epsg)
#     } else if (abs(linegradient) == Inf & y_coord > sf::st_coordinates(secondpoint)[2]) {
#       longline <- sf::st_sfc(sf::st_linestring(rbind(sf::st_coordinates(keypoint),sf::st_point(c(x_coord,boundbox[2])))), crs=epsg)
#     } else if (linegradient == 0 & x_coord < sf::st_coordinates(secondpoint)[1]) {
#       longline <- sf::st_sfc(sf::st_linestring(rbind(sf::st_coordinates(keypoint),sf::st_point(c(boundbox[3],y_coord)))), crs=epsg)
#     } else if (linegradient == 0 & x_coord > sf::st_coordinates(secondpoint)[1]) {
#       longline <- sf::st_sfc(sf::st_linestring(rbind(sf::st_coordinates(keypoint),sf::st_point(c(boundbox[1],y_coord)))), crs=epsg)
#     }
#
#     #calculate the distance between longline and 2nd nearest shore
#     longline_intersect <- sort(sf::st_distance(keypoint,sf::st_cast(sf::st_intersection(longline,landline),"POINT")))
#     if (length(longline_intersect) == 1) {
#       d2shore <- as.numeric(longline_intersect[1])
#       d2nextshore <- as.numeric(sf::st_distance(keypoint,sf::st_cast(sf::st_intersection(longline,boundary),"POINT"))[1])
#     } else {
#       d2shore <- as.numeric(longline_intersect[1])
#       d2nextshore <- as.numeric(longline_intersect[2])
#     }
#
#     #calculate weighted exposure score
#     ((d2nextshore-d2shore)/d2shore)
#   }
#   landgrid <- landgrid %>%
#     dplyr::mutate(exposure = as.numeric(exposurevect))
#   landgridrast <- terra::rasterize(terra::vect(landgrid),template,"exposure")
#   terra::writeRaster(landgridrast,paste0("E:/Pittas/nichemodelling/raster_exposure/exposure_interval",interval,".tif"),filetype="GTiff",overwrite=TRUE)
#   parallel::stopCluster(cl)
# }
#
# for (i in 1:1000) {
#   exposure(i)
# }


# #function to calculate windfetch along the shoreline
# calc_shorefetch <- function(interval) {
#   library(doParallel)
#   library(foreach)
#   library(windfetch)
#
#   #set up parallelisation
#   cores = parallel::detectCores()
#   cl <- parallel::makeCluster(cores[1]-4)
#   doParallel::registerDoParallel(cl)
#
#   #load shapefile of shore point samples
#   shorepoints <- sf::st_read(paste0("E:/Pittas/nichemodelling/shorepoints/shorepoints_interval",interval,".shp"))
#   #load shapefile of land polygons
#   landshape <- sf::st_read(paste0("E:/Pittas/nichemodelling/shapefile/interval",interval,".shp"))
#   #landshape_smooth <- smoothr::smooth(landshape,method = "ksmooth", smoothness = 2)
#   #start parallelised for loop that calculates windfetch for each point
#   fetchvector <- foreach::foreach(i=1:nrow(shorepoints), .combine = 'rbind') %dopar% {
#     #load individual shorepoint
#     keypoint <- shorepoints$geometry[i]
#     #calculate windfetch for the keypoint
#     fetch <- tryCatch(
#       {
#         fetch_estimate <- windfetch_fixed(landshape,keypoint)
#         #calculate mean windfetch from 3 directions
#         as.numeric(sum(round(fetch_estimate@.Data[[5]]))/(length(which(as.numeric(fetch_estimate@.Data[[5]]) != 0))))
#       },
#       error = function(fetch_error) {
#         return(NA)
#       }
#     )
#   }
#   shorepoints <- dplyr::mutate(shorepoints,fetchvector)
#   #write points shapefile with windfetch values
#   sf::st_write(shorepoints,paste0("E:/Pittas/nichemodelling/shorefetch/shorefetch_interval",interval,".shp"))
# }


# ## Estimating shorefetch using windfetch and furrr parallelisation
#
# library(furrr)
# library(future)
# interval = 0
# #load shapefile of shore point samples
# shorepoints <- sf::st_read(paste0("E:/Pittas/nichemodelling/shorepoints/shorepoints_interval",interval,".shp"))
# #load shapefile of land polygons
# landshape <- sf::st_read(paste0("E:/Pittas/nichemodelling/shapefile/interval",interval,".shp"))
# calc_shorefetch <- function(shorepoint) {
#   library(sf)
#   library(windfetch)
#   sf::st_crs(landshape) <- "epsg:6933"
#   shorepoint <- sf::st_sfc(shorepoint)
#   sf::st_crs(shorepoint) <- "epsg:6933"
#   windfetch::windfetch(landshape,shorepoint)
#   sum(summary(fetch_estimate)$avg_fetch)/3
#   #shorepoints <- sf::st_sfc(shorepoint)
#   # tryCatch(
#   #   {
#   #     fetch_estimate <- windfetch::windfetch(landshape,shorepoint)
#   #     sum(summary(fetch_estimate)$avg_fetch)/3
#   #   },
#   #   error = function(fetch_error) {
#   #     return(NA)
#   #   }
#   # )
# }
# plan(multicore)
# fetchvector <- furrr::future_map(sf::st_geometry(shorepoints),purrr::safely(calc_shorefetch))
# write.csv(fetchvector,)

#after multiple false starts, this seems to be the best approach to calculating shore energy



#calculate fetch raster for sea pixels using the fetchr package
calc_fetchraster <- function(interval) {
  library(terra)
  library(fetchr)
  land_raster <- terra::rast(paste0("/carc/scratch/projects/andersen2016005/pitta/nichemodelling/raster_flat/interval",interval,".tif"))
  fetchraster <- fetchr::get_fetch(r=land_raster, max_dist=300000, in_parallel=TRUE, verbose=TRUE, func = "mean")
  terra::writeRaster(fetchraster,paste0("/carc/scratch/projects/andersen2016005/pitta/nichemodelling/shorefetch/fetchraster_interval",interval,".tif"))
}

for (m in 0:600) {
  calc_fetchraster(m)
}

#sample fetch values along shoreline
calc_shorefetch <- function(interval) {
  library(sf)
  library(terra)
  library(exactextractr)
  shorepoints <- sf::st_read(paste0("E:/Pittas/nichemodelling/shorepoints/shorepoints_interval",interval,".shp"))
  fetchraster <- terra::rast(paste0("E:/Pittas/nichemodelling/shorefetch/fetchraster_interval",interval,".tif"))
  #buffer shorepoints by 500m
  shorepoints_buff <- sf::st_buffer(shorepoints,500)
  #extract the mean fetch value for each point along the shoreline
  shorefetch <- exactextractr::exact_extract(fetchraster,shorepoints_buff,fun="mean")
  #append fetch values to shorepoints vector
  shorepoints$fetch <- shorefetch
  #write file
  sf::st_write(shorepoints,paste0("E:/Pittas/nichemodelling/shorefetch/shorefetch_interval",interval,".shp"))
}


for (i in 601:800) {
  calc_shorefetch(i)
}

#and now for the final step, calculating the on-land fetch values
calc_landfetch <- function(interval) {
  library(sf)
  library(terra)
  template <- terra::rast(paste0("E:/Pittas/nichemodelling/raster_flat/interval",interval,".tif"))
  shorefetch <- sf::st_read(paste0("E:/Pittas/nichemodelling/shorefetch/shorefetch_interval",interval,".shp"))
  landpoints <- sf::st_read(paste0("E:/Pittas/nichemodelling/landgrids/landgrid_interval",interval,".shp"))
  nearestindex <- sf::st_nearest_feature(landpoints,shorefetch)
  dist <- as.numeric(sf::st_distance(landpoints,shorefetch[nearestindex,], by_element=TRUE))
  landpoints$FID <- c(1:nrow(landpoints))
  exposure <- function(i) {
    shorefetch$fetch[nearestindex[i]]/(dist[i]/1000)
  }
  landpoints <- landpoints %>% dplyr::mutate(exposure = exposure(landpoints$FID))
  landgridrast <- terra::rasterize(terra::vect(landpoints),template,"exposure")
  terra::writeRaster(landgridrast,paste0("E:/Pittas/nichemodelling/raster_exposure/exposure_interval",interval,".tif"),filetype="GTiff",overwrite=TRUE)
  gc()
}

for (i in 601:800) {
  calc_landfetch(i)
}

#load world mangrove shapefile
world_mangroves <- sf::st_read("E:/Pittas/nichemodelling/gmw_v3_f1996_t2020_vec.shp")
#reproject shapefile
world_mangroves <- sf::st_transform(world_mangroves, crs="epsg:6933")
#crop mangrove polygon to area of interest
mangroves_sea <- sf::st_crop(world_mangroves,landshape)
#write file
sf::st_write(mangroves_sea,"E:/Pittas/nichemodelling/mangroves_SEA.shp")

#generate a raster of present-day mangroves so we can simulate effects of contemporary mangrove loss
mangroves_sea <- terra::vect("E:/Pittas/nichemodelling/mangroves_SEA.shp")

#niche modelling script




#now this subroutine generates maps of distance from coastline

#for interval0
#presentvect <- sf::st_read("E:/Pittas/nichemodelling/shapefile/interval0.shp")
#close all internal holes so distance is only calculated from shoreline
#presentvect_noholes <- smoothr::fill_holes(presentvect,threshold=units::set_units(100000, km^2))
presentvect <- sf::st_read("E:/Pittas/nichemodelling/shapefile_with_intertidal/interval0.shp")
#sf::st_write(presentvect_noholes,"E:/Pittas/nichemodelling/shapefile_clean/interval0.shp")
#rasterise the nohole vector
template <- terra::rast(terra::vect(presentvect),res=1000)

mangroves_today_rast <- terra::rasterize(mangroves_sea,template)
terra::writeRaster(mangroves_today_rast,"E:/Pittas/nichemodelling/mangrove_predictions/predicted_intervaltoday.tif")

landraster <- terra::rasterize(vect(presentvect),template)
m <- c(NA,1,1,NA)
convmat <- matrix(m,ncol=2,byrow=TRUE)
landraster_inverse <- terra::classify(landraster,convmat)
shoredist <- terra::distance(landraster_inverse)
terra::writeRaster(shoredist,"E:/Pittas/nichemodelling/raster_shoredist/shoredist_interval0.tif",overwrite=TRUE)

#now let's do this for each timeslice
calc_shoredist <- function(interval) {
  landvect <- sf::st_read(paste0("E:/Pittas/nichemodelling/shapefile_clean/interval",interval,".shp"))
  #rasterise the nohole vector
  template <- terra::rast(terra::vect(landvect),res=1000)
  landraster <- terra::rasterize(vect(landvect),template)
  m <- c(NA,1,1,NA)
  convmat <- matrix(m,ncol=2,byrow=TRUE)
  landraster_inverse <- terra::classify(landraster,convmat)
  shoredist <- terra::distance(landraster_inverse)
  terra::writeRaster(shoredist,paste0("E:/Pittas/nichemodelling/raster_shoredist/shoredist_interval",interval,".tif"),overwrite=TRUE)
  gc()
}

for (i in 601:800) {
  calc_shoredist(i)
}

#generate curvature map
#present_curvature <- spatialEco::curvature(presentday,type="profile")
#writeRaster(present_curvature,"E:/Pittas/nichemodelling/raster_curvature/curvature_interval0.tif")

#generate interior buffer for present-day shoreline
presentvect_buff <- terra::buffer(presentvect_noholes_conv,width=-50000)
terra::writeVector(presentvect_buff,"E:/Pittas/nichemodelling/presentvect_buff.shp",overwrite=TRUE)

#then use QGIS to clean up the buffered shapefile to remove extraneous artefacts (fix polygons, multiparts to singleparts, then manually delete artefacts)
#Also use QGIS to filter out areas higher than 100m and add to the masking layer

#from QGIS, export final masking polygon as a Shapefile for use in MaxEnt
presentday <- terra::rast("E:/Pittas/nichemodelling/raster_topo/interval0.tif")
maskfile <- terra::vect("E:/Pittas/nichemodelling/shoremask/shoremask_interval0.shp")
masked_raster <- terra::mask(presentday,maskfile)

#convert spatRaster to rasterlayer
masked <- raster::raster(masked_raster)

#generate background points
bg <- dismo::randomPoints(masked,10000)
write.csv(bg,"E:/Pittas/nichemodelling/background_points.csv",row.names = FALSE)
bg <- read.csv("E:/Pittas/nichemodelling/background_points.csv")


#load predictor variable rasters
pred_topo <- terra::rast("E:/Pittas/nichemodelling/raster_topo_new/topo_interval0.tif")
pred_shoredist <- terra::rast("E:/Pittas/nichemodelling/raster_shoredist/shoredist_interval0.tif")
pred_slope <- terra::rast("E:/Pittas/nichemodelling/raster_slope/slope_interval0.tif")
pred_exposure <- terra::rast("E:/Pittas/nichemodelling/raster_exposure/exposure_interval0.tif")
#pred_curv <- raster::raster("E:/Pittas/nichemodelling/raster_curvature/curvature_interval0.tif")

# #mask predictor variable rasters
# pred_topo_mask <- terra::mask(pred_topo,maskfile,inverse=TRUE)
# pred_shoredist_mask <- terra::mask(pred_shoredist,maskfile,inverse=TRUE)
# pred_slope_mask <- terra::mask(pred_slope,maskfile,inverse=TRUE)
# pred_exposure_mask <- terra::mask(pred_exposure,maskfile,inverse=TRUE)
# #convert masked predictor variables to raster format
# pred_topo_mask <- raster::raster(pred_topo)
# pred_shoredist_mask <- raster::raster(pred_shoredist)
# pred_slope_mask <- raster::raster(pred_slope)
# pred_exposure_mask <- raster::raster(pred_exposure)

#convert unmasked predictor variables to raster format
pred_topo <- raster::raster(pred_topo)
pred_shoredist <- raster::raster(pred_shoredist)
pred_slope <- raster::raster(pred_slope)
pred_exposure <- raster::raster(pred_exposure)

#rename predictor raster layers
names(pred_topo) <- c("elevation")
names(pred_shoredist) <- c("shoredist")
names(pred_slope) <- c("slope")
names(pred_exposure) <- c("exposure")
# names(pred_topo_mask) <- c("elevation")
# names(pred_shoredist_mask) <- c("shoredist")
# names(pred_slope_mask) <- c("slope")
# names(pred_exposure_mask) <- c("exposure")
#names(pred_curv) <- c("curvature")

#stack predictor variables
#pred_stack_bg <- raster::stack(pred_topo_mask,pred_shoredist_mask,pred_slope_mask,pred_exposure_mask)

#stack unmasked predictor variables
pred_stack <- raster::stack(pred_topo,pred_shoredist,pred_slope,pred_exposure)

#load occurrence points
mangroves <- read.csv("E:/Pittas/nichemodelling/globalmangrovewatch_mangrove_points.csv",header=TRUE,sep=",")
#spThin::thin(loc.data=mangroves,lat.col = "y", long.col = "x",spec="spec",thin.par = 1, reps = 5, out.dir = "E:/Pittas/nichemodelling/")

#load thinned occurrence points
#mangroves_thin <- read.csv("E:/Pittas/nichemodelling/thinned_data_thin1.csv",header=T,sep=",")

#exclude first column
mangroves_thin <- dplyr::select(mangroves,x,y)
#mangroves_nothin <- dplyr::select(mangroves,longitude,latitude)
#mangroves_nothin <- mangroves_nothin[!duplicated(mangroves_nothin),]
colnames(bg) <- colnames(mangroves_thin)

#data partitioning
block <- ENMeval::get.block(mangroves_nothin,bg,orientation="lat_lon")
evalplot.grps(pts = mangroves_nothin, pts.grp = block$occs.grp, envs = pred_stack_bg)

#run Maxent
tune.args = list(fc= c("L","LQ","Q","P","LQP","LP"),rm=1:5)
mangrove_model <- ENMeval::ENMevaluate(occs = mangroves_thin,
                                       envs = pred_stack,
                                       bg = bg,
                                       algorithm = "maxnet",
                                       partitions = 'block',
                                       tune.args = tune.args)
eval.results(mangrove_model)
write.csv(eval.results(mangrove_model),"E:/Pittas/nichemodelling/mp_nichemodel/eval_results.csv")
eval.results.partitions(mangrove_model)
eval.predictions(mangrove_model)
eval.models(mangrove_model)[["fc.LQ_rm.1"]]$betas
#model evaluation

#LQP RM1 and LP RM1 both have the lowest AICc (146988.0, all other delta AICcs > 2).
#auc.val.avg = 0.8727914, 0.8727777
#or.10p.avg = 0.11359485, 0.11429269

#run null model
mod.null <- ENMeval::ENMnulls(mangrove_model,mod.settings = list(fc="LQ",rm=1),no.iter=100)
null.results(mod.null) %>% head()
null.emp.results(mod.null)
write.csv(null.emp.results(mod.null),"E:/Pittas/nichemodelling/mp_nichemodel/null_results.csv")
evalplot.nulls(mod.null, stats = c("or.10p", "auc.val"), plot.type = "histogram")

mangrove_results <- eval.results(mangrove_model)
#opt.aicc <- mangrove_results %>% filter(delta.AICc == 0)
opt.aicc_LQ1 <- mangrove_results %>% dplyr::filter(tune.args == "fc.LQ_rm.1")

opt.aicc_LQ1

opt.seq <- mangrove_results %>%
  dplyr::filter(or.10p.avg == min(or.10p.avg)) %>%
  dplyr::filter(auc.val.avg == max(auc.val.avg))

opt.seq

best_mangrove_aicc <- eval.models(mangrove_model)[[opt.aicc_LQ1[1,]$tune.args]]
#best_mangrove_aicc_alternative <- eval.models(mangrove_model)[[opt.aicc_LP1[1,]$tune.args]]
best_mangrove_aicc$betas
#best_mangrove_aicc_alternative$betas
plot(best_mangrove_aicc, type="cloglog")
#plot(best_mangrove_aicc_alternative, type="cloglog")

bestmod = which(mangrove_model@results$AICc==min(mangrove_model@results$AICc))

#plot model outputs
dev.off()
predicted_mangroves <- eval.predictions(mangrove_model)[[opt.aicc_LQ1[1,]$tune.args]]
plot(predicted_mangroves)
terra::writeRaster(predicted_mangroves,filename="E:/Pittas/nichemodelling/predicted_mangroves.tif",overwrite=TRUE)

#write enmevaluate model to file
rmm <- eval.rmm(mangrove_model)
rangeModelMetadata::rmmToCSV(rmm,"E:/Pittas/nichemodelling/mangrove_models.csv")

#calculate model threshold
pr <- dismo::predict(pred_stack,mangrove_model@models[[bestmod[1]]],type="cloglog")
est.loc <- raster::extract(pr,mangrove_model@occs[1:2])
est.bg <- raster::extract(pr,mangrove_model@bg[1:2])
ev <- dismo::evaluate(est.loc,est.bg)
thr <- dismo::threshold(ev)
thr
pr_thr <- pr > thr$sensitivity
plot(pr_thr)
pr_thr_filt <- raster::clamp(pr_thr,lower=0.5,useValues=FALSE)
terra::writeRaster(pr_thr,filename="E:/Pittas/nichemodelling/mangrove_predictions_withland/predicted_interval0.tif",overwrite=TRUE)
pr_thr_filt <- terra::rast(pr_thr_filt)
pr_thr_filt_clean <- terra::sieve(pr_thr_filt,threshold=4,directions=8)
pr_thr_filt_clean <- terra::classify(pr_thr_filt_clean,rcl = matrix(c(0,NA),ncol=2))
terra::writeRaster(pr_thr_filt_clean,filename="E:/Pittas/nichemodelling/mangrove_predictions/predicted_interval0.tif",overwrite=TRUE)

#test historical rasters

for (i in 601:800) {
  pred_topo_interval <- raster::raster(paste0("E:/Pittas/nichemodelling/raster_topo_new/topo_interval",i,".tif"))
  pred_slope_interval <- raster::raster(paste0("E:/Pittas/nichemodelling/raster_slope/slope_interval",i,".tif"))
  pred_shoredist_interval <- raster::raster(paste0("E:/Pittas/nichemodelling/raster_shoredist/shoredist_interval",i,".tif"))
  pred_exposure_interval <- raster::raster(paste0("E:/Pittas/nichemodelling/raster_exposure/exposure_interval",i,".tif"))
  names(pred_topo_interval) <- c("elevation")
  names(pred_slope_interval) <- c("slope")
  names(pred_shoredist_interval) <- c("shoredist")
  names(pred_exposure_interval) <- c("exposure")

  pred_stack_interval <- raster::stack(pred_topo_interval,pred_slope_interval,pred_shoredist_interval,pred_exposure_interval)
  pr_interval <- dismo::predict(pred_stack_interval,mangrove_model@models[[bestmod[1]]],type="cloglog")
  pr_interval_threshold <- pr_interval > 0.5006234
  terra::writeRaster(pr_interval_threshold,paste0("E:/Pittas/nichemodelling/mangrove_predictions_withland/predicted_interval",i,".tif"),overwrite=TRUE)
  pr_interval_threshold_filt <- raster::clamp(pr_interval_threshold,lower=0.5,useValues=FALSE)
  pr_interval_threshold_filt <- terra::rast(pr_interval_threshold_filt)
  pr_interval_threshold_filt <- terra::sieve(pr_interval_threshold_filt,threshold=4,directions=8)
  pr_interval_threshold_filt <- terra::classify(pr_interval_threshold_filt,rcl = matrix(c(0,NA),ncol=2))
  terra::writeRaster(pr_interval_threshold_filt,paste0("E:/Pittas/nichemodelling/mangrove_predictions/predicted_interval",i,".tif"),overwrite=TRUE)
  gc()
}

andaman <- terra::vect("E:/Pittas/nichemodelling/andaman_mangroves.shp")
scs <- terra::vect("E:/Pittas/nichemodelling/scs_mangroves.shp")
metrics <- data.frame()

interval_mangrove <- terra::rast(paste0("E:/Pittas/nichemodelling/mangrove_predictions/predicted_intervaltoday.tif"))
total_area <- terra::expanse(interval_mangrove,unit="m",transform=FALSE)
andaman_area <- terra::expanse(terra::crop(interval_mangrove,andaman,mask=TRUE),unit="m",transform=FALSE)
scs_area <- expanse(terra::crop(interval_mangrove,scs,mask=TRUE),unit="m",transform=FALSE)

#calculate area and overlap metrics
for (i in 599:799) {
  print(paste0("Evaluating interval ",i,"..."))
  interval_mangrove <- terra::rast(paste0("E:/Pittas/nichemodelling/mangrove_predictions/predicted_interval",i,".tif"))
  total_area <- terra::expanse(interval_mangrove,unit="m",transform=FALSE)
  andaman_area <- terra::expanse(terra::crop(interval_mangrove,andaman,mask=TRUE),unit="m",transform=FALSE)
  scs_area <- expanse(terra::crop(interval_mangrove,scs,mask=TRUE),unit="m",transform=FALSE)
  earlier_interval <- terra::rast(paste0("E:/Pittas/nichemodelling/mangrove_predictions/predicted_interval",i+1,".tif"))
  andaman_earlier <- terra::crop(earlier_interval,andaman,mask=TRUE)
  overlap_andaman <- terra::clamp(terra::crop(interval_mangrove,andaman,mask=TRUE) + andaman_earlier,lower = 1.5,value=FALSE)
  overlap_andaman_area <- terra::expanse(overlap_andaman,unit="m",transform=FALSE)
  andaman_percent_diff <- (overlap_andaman_area$area / andaman_area$area * 100)
  scs_earlier <- terra::crop(earlier_interval,scs,mask=TRUE)
  overlap_scs <- terra::clamp(terra::crop(interval_mangrove,scs,mask=TRUE) + scs_earlier, lower = 1.5, value = FALSE)
  overlap_scs_area <- terra::expanse(overlap_scs,transform=FALSE)
  scs_percent_diff <- (overlap_scs_area$area / scs_area$area * 100)
  metrics <- rbind(metrics,c(i,total_area$area,andaman_area$area,scs_area$area,overlap_andaman_area$area,andaman_percent_diff,overlap_scs_area$area,scs_percent_diff))
  gc()
}

write.csv(metrics,"E:/Pittas/nichemodelling/mangrove_metrics_599-799.csv")

#time to plot the maps and graphs
ext <- st_bbox(basemap)
cartolight <- basemap_ggplot(ext,map_service = "carto", map_type = "light")
sealvl <- read.csv("E:/Pittas/ddRAD/8_Biogeography/sealvl.csv")
mangrove_metrics <- read.csv("E:/Pittas/nichemodelling/mangrove_metrics.csv",header = TRUE)

for (i in 0:600) {
  interval_shp <- terra::vect(paste0("E:/Pittas/nichemodelling/shapefile/interval",i,".shp"))
  #interval_shp <- terra::crop(interval_shp,sea_clip)
  interval_shp <- terra::project(interval_shp,"epsg:3857")
  interval_mangroves <- terra::rast(paste0("E:/Pittas/nichemodelling/mangrove_predictions/predicted_interval",i,".tif"))
  #andaman_mangroves <- terra::crop(interval_mangroves,andaman,mask=TRUE)
  #scs_mangroves <- terra::crop(interval_mangroves,scs,mask=TRUE)
  interval_mangroves <- terra::project(interval_mangroves,"epsg:3857")
  #andaman_mangroves <- terra::project(andaman_mangroves,"epsg:3857")
  #scs_mangroves <- terra::project(scs_mangroves,"epsg:3857")
  interval_mangroves <- as.factor(interval_mangroves)
  interval_upptime <- filter(intervalfile,Interval==i)$UpperTimeBound
  interval_lowtime <- filter(intervalfile,Interval==i)$LowerTimeBound
  map <- cartolight+
    geom_spatvector(data=interval_shp,fill="green4",alpha=0.4,colour=NA)+
    geom_spatraster(data=interval_mangroves)+
    theme_void()+
    annotate("text",label=paste0("Time: ",interval_lowtime,"-",interval_upptime," kya\nMean Sea Level: ",filter(intervalfile,Interval==i)$MeanDepth,"m"),x=9600000,y=-0.8E06,hjust=0)
  sealvlgraph <- ggplot(sealvl,aes(x=Time, y=Sealevel_Corrected))+
    geom_line(colour="blue4",linewidth=1.2)+xlim(0,600)+
    xlab("Time (kya)")+
    ylab("Sea level (m)")+
    ggtitle("Sea level from 0 to 600 kya (Bintanja and van de Wal, 2008)")+
    annotate("rect",xmin=interval_lowtime,xmax=interval_upptime,ymin=-Inf,ymax=Inf,fill="lightblue3",alpha=0.9)
  mangrovegraph <- ggplot(mangrove_metrics,aes(x=Time,y=Total.Area))+
    geom_line(colour = "black",linewidth=1.2)+
    geom_line(aes(x=Time,y=SCS.Area),colour="royalblue3",linewidth=1.2,alpha=0.7)+
    geom_line(aes(x=Time,y=Andaman.Area),colour="orange2",linewidth=1.2,alpha=0.7)+
    annotate("rect",xmin=interval_lowtime,xmax=interval_upptime,ymin=-Inf,ymax=Inf,fill="lightblue3",alpha=0.9)+
    xlab("Time (kya)")+
    ylab(expression(paste("Mangrove Area (km "^"2",")")))+
    ggtitle("Mangrove area from 0 to 600 kya")
  percentgraph <- ggplot(mangrove_metrics)+
    geom_line(aes(x=Time,y=Andaman.Percent.Similarity,colour="andaman"),linewidth=1.2,alpha=0.7)+
    geom_line(aes(x=Time,y=SCS.Percent.Similarity,colour="scs"),linewidth=1.2,alpha=0.7)+
    xlab("Time (kya)")+
    ylab("Overlap (%)")+
    annotate("rect",xmin=interval_lowtime,xmax=interval_upptime,ymin=-Inf,ymax=Inf,fill="lightblue3",alpha=0.9)+
    ggtitle("Percentage mangrove area overlap with previous interval")+
    scale_colour_manual(values = c(andaman = "orange2", scs = "royalblue3"), labels = c(andaman = "Andaman Sea", scs = "South China Sea"), limits = c("andaman","scs"))+
    theme(legend.position="bottom",legend.title=element_blank())
  metric_graphs <- ggarrange(sealvlgraph,mangrovegraph,percentgraph,nrow=3,align="hv",common.legend = TRUE)
  ggarrange(map,metric_graphs,ncol=2,widths = c(1.5,1))
  ggsave(paste0("E:/Pittas/nichemodelling/mangrove_maps/interval",i,".png"),width = 20, height = 10, units = "in",bg="white")
  gc()
}

pr <- predict(pred_stack_interval2,mangrove_model@models[[bestmod]],type="cloglog")
pr_df = as.data.frame(pr, xy=T)

ggplot() +
  geom_raster(data = pr_df, aes(x = x, y = y, fill = layer)) +
  coord_quickmap() +
  theme_bw() +
  scale_fill_gradientn(colours=viridis::viridis(99),na.value = "black")

terra::writeRaster(pr,"E:/Pittas/nichemodelling/predicted_mangroves_interval2.tif")
#presentshore <- clamp(basemap,upper=0,value=FALSE)
#presentshore_patches <- patches(presentshore)
#presentshore_patches <- terra::rast("E:/Pittas/nichemodelling/presentshore_patches.tif")
#presentshore_clean <- mask(presentshore,presentshore_patches,inverse=TRUE,maskvalues=2)
#presentshoredist <- terra::distance(presentshore_clean, haversine=TRUE)

#baseshp <- st_read("E:/Pittas/nichemodelling/shapefile/interval0.shp")
#baseshp_union <- st_union(baseshp)
#baseshp_clean <- st_make_valid(baseshp_union)
#baseshp_clean <- sfheaders::sf_remove_holes(baseshp_clean)
#st_write(baseshp, "E:/Pittas/nichemodelling/test.shp")




###################################################

##old script for plotting animated sea level change map/chart

setwd("E:/Pittas/ddRAD/8_Biogeography")
getintervals_time(time = 600,intervals = 600,outdir="E:/Pittas/ddRAD/8_Biogeography/")
makemaps(inputraster="C:/Users/david/Dropbox/ResearchProjects/PhD_Research/pleistodist_testenv/SEA/Sundaland.asc",epsg = 4326,intervalfile = "intervals.csv",outdir = "E:/Pittas/ddRAD/8_Biogeography/")

setwd("E:/Pittas/ddRAD/8_Biogeography/batch1/shapefile/")

sealvl <- read.csv("E:/Pittas/ddRAD/8_Biogeography/sealvl.csv")
ggplot(sealvl,aes(x=Time, y=Sealevel_Corrected),)+
  geom_point(colour="blue4")+xlim(400,600)+
  xlab("Time (kya)")+
  ylab("Sea level relative to present day (m)")+
  ggtitle("Sea level change from 400 to 600 kya")+
  geom_vline(xintercept=i+400, linetype="dashed", size= 1.5)

sundaland <- rast("C:/Users/david/Dropbox/ResearchProjects/PhD_Research/pleistodist_testenv/SEA/Sundaland.asc")
crs(sundaland) <- "epsg:4326"
i=0
interval <- st_read(paste0("interval",i,".shp"))
ext <- st_bbox(interval)
cartolight <- basemap_ggplot(ext,map_service = "carto", map_type = "light")
intervalfile <- read.csv("E:/Pittas/ddRAD/8_Biogeography/intervals.csv")


for (i in 0:85) {
  interval <- st_read(paste0("interval",i,".shp"))
  interval_proj <- st_transform(interval,crs="epsg:3857")
  time <- i + 400
  map <- cartolight+
    geom_sf(data=interval_proj, fill="green4", alpha = 0.4,colour=NA)+
    theme_void()+
    theme(
      panel.border=element_blank(),
      panel.background=element_blank(),
      axis.text = element_blank(),
      plot.background = element_blank(),
      panel.grid=element_blank(),
      axis.ticks = element_blank())+
    annotate("text", label =paste0("Time: ",i+400," kya\nSea Level: ",filter(intervalfile, Interval == i+400)$MeanDepth,"m"),x=10300000,y=-1050000,hjust=0)
  graph <- ggplot(sealvl,aes(x=Time, y=Sealevel_Corrected),)+
    geom_point(colour="blue4")+xlim(400,600)+
    xlab("Time (kya)")+
    ylab("Sea level (m)")+
    ggtitle("Sea level change from 400 to 600 kya")+
    geom_vline(xintercept=i+400, linetype="dashed", size= 1.5)
  ggarrange(graph,map,nrow=2, heights = c(1,2.5))
  ggsave(paste0("interval",i+400,".png"),width = 10, height = 10, units = "in")
}

setwd("E:/Pittas/ddRAD/8_Biogeography/batch2/shapefile/")

for (i in 0:115) {
  interval <- st_read(paste0("interval",i,".shp"))
  interval_proj <- st_transform(interval,crs="epsg:3857")
  time <- i + 486
  map <- cartolight+
    geom_sf(data=interval_proj, fill="green4", alpha = 0.4,colour=NA)+
    theme_void()+
    theme(
      panel.border=element_blank(),
      panel.background=element_blank(),
      axis.text = element_blank(),
      plot.background = element_blank(),
      panel.grid=element_blank(),
      axis.ticks = element_blank())+
    annotate("text", label =paste0("Time: ",i+486," kya\nSea Level: ",filter(intervalfile, Interval == i+486)$MeanDepth,"m"),x=10300000,y=-1050000,hjust=0)
  graph <- ggplot(sealvl,aes(x=Time, y=Sealevel_Corrected),)+
    geom_point(colour="blue4")+xlim(400,600)+
    xlab("Time (kya)")+
    ylab("Sea level (m)")+
    ggtitle("Sea level change from 400 to 600 kya")+
    geom_vline(xintercept=i+486, linetype="dashed", size= 1.5)
  ggarrange(graph,map,nrow=2, heights = c(1,2.5))
  ggsave(paste0("interval",i+486,".png"),width = 10, height = 10, units = "in")
}

setwd("E:/Pittas/ddRAD/8_Biogeography/images/")

imgs <- list.files("E:/Pittas/ddRAD/8_Biogeography/images/", full.names = TRUE)
imgs <- list.reverse(imgs)
img_list <- lapply(imgs, image_read)
img_joined <- image_join(img_list)
img_animated <- image_animate(img_joined, fps = 10)

img_animated

image_write(image = img_animated,
            path = "500kya_sealvl.gif")




