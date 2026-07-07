## script for plotting Barreto et al 2023 paleoclimate data for Sundaland
options(java.parameters="-Xmx8g")
library(pastclim)
library(sf)
library(ggplot2)
library(dplyr)
library(gstat)
library(ENMeval)
library(dismo)
library(basemaps)
library(tidyterra)
library(ggnewscale)
library(biscale)
library(cowplot)
library(fitMaxnet)
library(SDMtune)

pastclim::set_data_path("E:/Pittas/nichemodelling/Paleoclimate_Barreto2023/")
landvect_noholes <- sf::st_transform(sf::st_read("E:/Pittas/nichemodelling/shapefile_with_intertidal/interval0.shp"),4326)
boundbox <-sf::st_bbox(landvect_noholes)
boundary <- sf::st_polygon(sf::st_sfc(sf::st_linestring(c(sf::st_point(c(boundbox[1],boundbox[2])),sf::st_point(c(boundbox[3],boundbox[2])),sf::st_point(c(boundbox[3],boundbox[4])),sf::st_point(c(boundbox[1],boundbox[4])),sf::st_point(c(boundbox[1],boundbox[2]))))))
boundary <- vect(boundary)
plot(boundary)

#SEA_timeseries <- region_series(bio_variables=c("bio01","bio12","bio13","bio14"),dataset="Barreto2023", time_bp = list(min=-800000,max=0),crop=boundary)
download_dataset("Krapp2021")
#SEA_timeseries_krapp <- region_series(bio_variables=c("temperature_05","temperature_06","temperature_07","temperature_08","precip-itation_05","precipitation_06","precipitation_07","precipitation_08"),dataset="Krapp2021", time_bp = list(min=-800000,max=0),crop=boundary)
#extract bioclim variables to model distribution of tropical monsoon forest (BWPitta habitat)
#bio01 - mean annual temperature
#bio05 - max temp of warmest month
#bio06 - min temp of coldest month
#bio12 - annual precip
#bio13 - precip of wettest month
#bio14 - precip of driest month
#bio15 - precip seasonality
#bio16 - precip of wettest quarter
#bio17 - precip of driest quarter
SEA_timeseries_krapp <- region_series(bio_variables=c("bio01","bio05","bio06","bio12","bio13","bio14","bio15","bio16","bio17"),dataset="Krapp2021", time_bp = list(min=-800000,max=0),crop=boundary)
plot(SEA_timeseries_krapp$bio01)

krapp_biome <- region_slice(bio_variables="biome",dataset="Krapp2021",time_bp=list(min=-800000,max=0),crop=boundary)
plot(krapp_biome)

elevationmap <- terra::rast("E:/Pittas/nichemodelling/raster_topo_new/topo_interval799.tif")

#generate paleoclimate rasters for pitta core breeding season (June to July) with interpolation for missing areas
for (i in 1:800) {
  meantemp <- (SEA_timeseries_krapp$temperature_05[[i]] + SEA_timeseries_krapp$temperature_06[[i]] + SEA_timeseries_krapp$temperature_07[[i]] + SEA_timeseries_krapp$temperature_08[[i]])/4
  meantemp <- terra::focal(meantemp,NAonly=T,fun=mean,na.rm=T,w=5,na.policy="only")
  meantemp <- terra::project(meantemp,"epsg:6933")
  meantemp <- terra::resample(meantemp,elevationmap,method="bilinear")
  terra::writeRaster(meantemp,paste0("E:/Pittas/nichemodelling/krapp_paleoclim_meantemp/meantemp_interval",800-i,".tif"),gdal=c("COMPRESS=DEFLATE"),overwrite=TRUE)
  maxtemp <- max(SEA_timeseries_krapp$temperature_05[[i]],SEA_timeseries_krapp$temperature_06[[i]],SEA_timeseries_krapp$temperature_07[[i]],SEA_timeseries_krapp$temperature_08[[i]])
  maxtemp <- terra::focal(maxtemp,NAonly=T,fun=mean,na.rm=T,w=5,na.policy="only")
  maxtemp <- terra::project(maxtemp,"epsg:6933")
  maxtemp <- terra::resample(maxtemp,elevationmap,method="bilinear")
  terra::writeRaster(maxtemp,paste0("E:/Pittas/nichemodelling/krapp_paleoclim_maxtemp/maxtemp_interval",800-i,".tif"),gdal=c("COMPRESS=DEFLATE"),overwrite=TRUE)
  mintemp <- min(SEA_timeseries_krapp$temperature_05[[i]],SEA_timeseries_krapp$temperature_06[[i]],SEA_timeseries_krapp$temperature_07[[i]],SEA_timeseries_krapp$temperature_08[[i]])
  mintemp <- terra::focal(mintemp,NAonly=T,fun=mean,na.rm=T,w=5,na.policy="only")
  mintemp <- terra::project(mintemp,"epsg:6933")
  mintemp <- terra::resample(mintemp,elevationmap,method="bilinear")
  terra::writeRaster(mintemp,paste0("E:/Pittas/nichemodelling/krapp_paleoclim_mintemp/mintemp_interval",800-i,".tif"),gdal=c("COMPRESS=DEFLATE"),overwrite=TRUE)
  gc()
  meanprecip <- (SEA_timeseries_krapp$precipitation_05[[i]] + SEA_timeseries_krapp$precipitation_06[[i]] + SEA_timeseries_krapp$precipitation_07[[i]] + SEA_timeseries_krapp$precipitation_08[[i]])/4
  meanprecip <- terra::focal(meanprecip,NAonly=T,fun=mean,na.rm=T,w=5,na.policy="only")
  meanprecip <- terra::project(meanprecip,"epsg:6933")
  meanprecip <- terra::resample(meanprecip,elevationmap,method="bilinear")
  terra::writeRaster(meanprecip,paste0("E:/Pittas/nichemodelling/krapp_paleoclim_meanprecip/meanprecip_interval",800-i,".tif"),gdal=c("COMPRESS=DEFLATE"),overwrite=TRUE)
  maxprecip <- max(SEA_timeseries_krapp$precipitation_05[[i]],SEA_timeseries_krapp$precipitation_06[[i]],SEA_timeseries_krapp$precipitation_07[[i]],SEA_timeseries_krapp$precipitation_08[[i]])
  maxprecip <- terra::focal(maxprecip,NAonly=T,fun=mean,na.rm=T,w=5,na.policy="only")
  maxprecip <- terra::project(maxprecip,"epsg:6933")
  maxprecip <- terra::resample(maxprecip,elevationmap,method="bilinear")
  terra::writeRaster(maxprecip,paste0("E:/Pittas/nichemodelling/krapp_paleoclim_maxprecip/maxprecip_interval",800-i,".tif"),gdal=c("COMPRESS=DEFLATE"),overwrite=TRUE)
  minprecip <- min(SEA_timeseries_krapp$precipitation_05[[i]],SEA_timeseries_krapp$precipitation_06[[i]],SEA_timeseries_krapp$precipitation_07[[i]],SEA_timeseries_krapp$precipitation_08[[i]])
  minprecip <- terra::focal(minprecip,NAonly=T,fun=mean,na.rm=T,w=5,na.policy="only")
  minprecip <- terra::project(minprecip,"epsg:6933")
  minprecip <- terra::resample(minprecip,elevationmap,method="bilinear")
  terra::writeRaster(minprecip,paste0("E:/Pittas/nichemodelling/krapp_paleoclim_minprecip/minprecip_interval",800-i,".tif"),gdal=c("COMPRESS=DEFLATE"),overwrite=TRUE)
  gc()
}


#generate paleoclimate rasters for blue-winged pitta breeding habitat (tropical monsoon forest) with interpolation for missing areas

# for (i in 1:800) {
#   meantemp <- (SEA_timeseries_krapp$temperature_05[[i]] + SEA_timeseries_krapp$temperature_06[[i]] + SEA_timeseries_krapp$temperature_07[[i]] + SEA_timeseries_krapp$temperature_08[[i]])/4
#   meantemp <- terra::focal(meantemp,NAonly=T,fun=mean,na.rm=T,w=5,na.policy="only")
#   meantemp <- terra::project(meantemp,"epsg:6933")
#   meantemp <- terra::resample(meantemp,elevationmap,method="bilinear")
#   terra::writeRaster(meantemp,paste0("E:/Pittas/nichemodelling/krapp_paleoclim_meantemp/meantemp_interval",800-i,".tif"),gdal=c("COMPRESS=DEFLATE"),overwrite=TRUE)
#   maxtemp <- max(SEA_timeseries_krapp$temperature_05[[i]],SEA_timeseries_krapp$temperature_06[[i]],SEA_timeseries_krapp$temperature_07[[i]],SEA_timeseries_krapp$temperature_08[[i]])
#   maxtemp <- terra::focal(maxtemp,NAonly=T,fun=mean,na.rm=T,w=5,na.policy="only")
#   maxtemp <- terra::project(maxtemp,"epsg:6933")
#   maxtemp <- terra::resample(maxtemp,elevationmap,method="bilinear")
#   terra::writeRaster(maxtemp,paste0("E:/Pittas/nichemodelling/krapp_paleoclim_maxtemp/maxtemp_interval",800-i,".tif"),gdal=c("COMPRESS=DEFLATE"),overwrite=TRUE)
#   mintemp <- min(SEA_timeseries_krapp$temperature_05[[i]],SEA_timeseries_krapp$temperature_06[[i]],SEA_timeseries_krapp$temperature_07[[i]],SEA_timeseries_krapp$temperature_08[[i]])
#   mintemp <- terra::focal(mintemp,NAonly=T,fun=mean,na.rm=T,w=5,na.policy="only")
#   mintemp <- terra::project(mintemp,"epsg:6933")
#   mintemp <- terra::resample(mintemp,elevationmap,method="bilinear")
#   terra::writeRaster(mintemp,paste0("E:/Pittas/nichemodelling/krapp_paleoclim_mintemp/mintemp_interval",800-i,".tif"),gdal=c("COMPRESS=DEFLATE"),overwrite=TRUE)
#   gc()
#   meanprecip <- (SEA_timeseries_krapp$precipitation_05[[i]] + SEA_timeseries_krapp$precipitation_06[[i]] + SEA_timeseries_krapp$precipitation_07[[i]] + SEA_timeseries_krapp$precipitation_08[[i]])/4
#   meanprecip <- terra::focal(meanprecip,NAonly=T,fun=mean,na.rm=T,w=5,na.policy="only")
#   meanprecip <- terra::project(meanprecip,"epsg:6933")
#   meanprecip <- terra::resample(meanprecip,elevationmap,method="bilinear")
#   terra::writeRaster(meanprecip,paste0("E:/Pittas/nichemodelling/krapp_paleoclim_meanprecip/meanprecip_interval",800-i,".tif"),gdal=c("COMPRESS=DEFLATE"),overwrite=TRUE)
#   maxprecip <- max(SEA_timeseries_krapp$precipitation_05[[i]],SEA_timeseries_krapp$precipitation_06[[i]],SEA_timeseries_krapp$precipitation_07[[i]],SEA_timeseries_krapp$precipitation_08[[i]])
#   maxprecip <- terra::focal(maxprecip,NAonly=T,fun=mean,na.rm=T,w=5,na.policy="only")
#   maxprecip <- terra::project(maxprecip,"epsg:6933")
#   maxprecip <- terra::resample(maxprecip,elevationmap,method="bilinear")
#   terra::writeRaster(maxprecip,paste0("E:/Pittas/nichemodelling/krapp_paleoclim_maxprecip/maxprecip_interval",800-i,".tif"),gdal=c("COMPRESS=DEFLATE"),overwrite=TRUE)
#   minprecip <- min(SEA_timeseries_krapp$precipitation_05[[i]],SEA_timeseries_krapp$precipitation_06[[i]],SEA_timeseries_krapp$precipitation_07[[i]],SEA_timeseries_krapp$precipitation_08[[i]])
#   minprecip <- terra::focal(minprecip,NAonly=T,fun=mean,na.rm=T,w=5,na.policy="only")
#   minprecip <- terra::project(minprecip,"epsg:6933")
#   minprecip <- terra::resample(minprecip,elevationmap,method="bilinear")
#   terra::writeRaster(minprecip,paste0("E:/Pittas/nichemodelling/krapp_paleoclim_minprecip/minprecip_interval",800-i,".tif"),gdal=c("COMPRESS=DEFLATE"),overwrite=TRUE)
#   gc()
# }

for (i in 1:800) {
  bio01 <- terra::focal(SEA_timeseries_krapp$bio01[[i]],NAonly=T,fun=mean,na.rm=T,w=5,na.policy="only")
  bio01 <- terra::project(bio01,"epsg:6933")
  bio01 <- terra::resample(bio01,elevationmap,method="bilinear")
  terra::writeRaster(bio01,paste0("E:/Pittas/nichemodelling/krapp_paleoclim_bio01/bio01_interval",800-i,".tif"),gdal=c("COMPRESS=DEFLATE"),overwrite=TRUE)
  bio05 <- terra::focal(SEA_timeseries_krapp$bio05[[i]],NAonly=T,fun=mean,na.rm=T,w=5,na.policy="only")
  bio05 <- terra::project(bio05,"epsg:6933")
  bio05 <- terra::resample(bio05,elevationmap,method="bilinear")
  terra::writeRaster(bio05,paste0("E:/Pittas/nichemodelling/krapp_paleoclim_bio05/bio05_interval",800-i,".tif"),gdal=c("COMPRESS=DEFLATE"),overwrite=TRUE)
  bio06 <- terra::focal(SEA_timeseries_krapp$bio06[[i]],NAonly=T,fun=mean,na.rm=T,w=5,na.policy="only")
  bio06 <- terra::project(bio06,"epsg:6933")
  bio06 <- terra::resample(bio06,elevationmap,method="bilinear")
  terra::writeRaster(bio06,paste0("E:/Pittas/nichemodelling/krapp_paleoclim_bio06/bio06_interval",800-i,".tif"),gdal=c("COMPRESS=DEFLATE"),overwrite=TRUE)
  gc()
  bio12 <- terra::focal(SEA_timeseries_krapp$bio12[[i]],NAonly=T,fun=mean,na.rm=T,w=5,na.policy="only")
  bio12 <- terra::project(bio12,"epsg:6933")
  bio12 <- terra::resample(bio12,elevationmap,method="bilinear")
  terra::writeRaster(bio12,paste0("E:/Pittas/nichemodelling/krapp_paleoclim_bio12/bio12_interval",800-i,".tif"),gdal=c("COMPRESS=DEFLATE"),overwrite=TRUE)
  bio13 <- terra::focal(SEA_timeseries_krapp$bio13[[i]],NAonly=T,fun=mean,na.rm=T,w=5,na.policy="only")
  bio13 <- terra::project(bio13,"epsg:6933")
  bio13 <- terra::resample(bio13,elevationmap,method="bilinear")
  terra::writeRaster(bio13,paste0("E:/Pittas/nichemodelling/krapp_paleoclim_bio13/bio13_interval",800-i,".tif"),gdal=c("COMPRESS=DEFLATE"),overwrite=TRUE)
  bio14 <- terra::focal(SEA_timeseries_krapp$bio14[[i]],NAonly=T,fun=mean,na.rm=T,w=5,na.policy="only")
  bio14 <- terra::project(bio14,"epsg:6933")
  bio14 <- terra::resample(bio14,elevationmap,method="bilinear")
  terra::writeRaster(bio15,paste0("E:/Pittas/nichemodelling/krapp_paleoclim_bio14/bio14_interval",800-i,".tif"),gdal=c("COMPRESS=DEFLATE"),overwrite=TRUE)
  bio17 <- terra::focal(SEA_timeseries_krapp$bio15[[i]],NAonly=T,fun=mean,na.rm=T,w=5,na.policy="only")
  bio17 <- terra::project(bio15,"epsg:6933")
  bio17 <- terra::resample(bio15,elevationmap,method="bilinear")
  terra::writeRaster(bio15,paste0("E:/Pittas/nichemodelling/krapp_paleoclim_bio15/bio15_interval",800-i,".tif"),gdal=c("COMPRESS=DEFLATE"),overwrite=TRUE)
}

for (i in 1:800) {
  bio16 <- terra::focal(SEA_timeseries_krapp$bio16[[i]],NAonly=T,fun=mean,na.rm=T,w=5,na.policy="only")
  bio16 <- terra::project(bio16,"epsg:6933")
  bio16 <- terra::resample(bio16,elevationmap,method="bilinear")
  terra::writeRaster(bio16,paste0("E:/Pittas/nichemodelling/krapp_paleoclim_bio16/bio16_interval",800-i,".tif"),gdal=c("COMPRESS=DEFLATE"),overwrite=TRUE)
  bio17 <- terra::focal(SEA_timeseries_krapp$bio17[[i]],NAonly=T,fun=mean,na.rm=T,w=5,na.policy="only")
  bio17 <- terra::project(bio17,"epsg:6933")
  bio17 <- terra::resample(bio17,elevationmap,method="bilinear")
  terra::writeRaster(bio17,paste0("E:/Pittas/nichemodelling/krapp_paleoclim_bio17/bio17_interval",800-i,".tif"),gdal=c("COMPRESS=DEFLATE"),overwrite=TRUE)
  gc()
}

for (i in 1:800) {
  writeRaster(krapp_biome$biome[[i]],paste0("E:/Pittas/nichemodelling/krapp_paleoclim_biomes/biomes_interval",800-i,".tif"),gdal=c("COMPRESS=DEFLATE"),overwrite=TRUE,datatype="INT1U")
}

##Niche modelling script

#from QGIS, export final masking polygon as a Shapefile for use in MaxEnt
presentday <- terra::rast("E:/Pittas/nichemodelling/raster_topo_new/topo_interval0.tif")
maskfile <- terra::vect("E:/Pittas/nichemodelling/BWPitta_BreedingRange.shp")
masked_raster <- terra::mask(presentday,maskfile)

#convert spatRaster to rasterlayer
masked <- raster::raster(masked_raster)

#generate background points
bg <- dismo::randomPoints(masked,10000)
write.csv(bg,"E:/Pittas/nichemodelling/bwp_background_points.csv",row.names = FALSE)
bg <- read.csv("E:/Pittas/nichemodelling/bwp_background_points.csv")

terra::crs(boundary) <- "epsg:4326"
# worldclim_temp_may <- terra::project(terra::crop(terra::rast("E:/Pittas/nichemodelling/worldclim/wc2.1_30s_tavg/wc2.1_30s_tavg_05.tif"),boundary),"epsg:6933")
# worldclim_temp_may <- terra::resample(worldclim_temp_may,elevationmap,method="bilinear")
# worldclim_temp_jun <- terra::project(terra::crop(terra::rast("E:/Pittas/nichemodelling/worldclim/wc2.1_30s_tavg/wc2.1_30s_tavg_06.tif"),boundary),"epsg:6933")
# worldclim_temp_jun <- terra::resample(worldclim_temp_jun,elevationmap,method="bilinear")
# worldclim_temp_jul <- terra::project(terra::crop(terra::rast("E:/Pittas/nichemodelling/worldclim/wc2.1_30s_tavg/wc2.1_30s_tavg_07.tif"),boundary),"epsg:6933")
# worldclim_temp_jul <- terra::resample(worldclim_temp_jul,elevationmap,method="bilinear")
# worldclim_temp_aug <- terra::project(terra::crop(terra::rast("E:/Pittas/nichemodelling/worldclim/wc2.1_30s_tavg/wc2.1_30s_tavg_08.tif"),boundary),"epsg:6933")
# worldclim_temp_aug <- terra::resample(worldclim_temp_aug,elevationmap,method="bilinear")
# worldclim_precip_may <- terra::project(terra::crop(terra::rast("E:/Pittas/nichemodelling/worldclim/wc2.1_30s_prec/wc2.1_30s_prec_05.tif"),boundary),"epsg:6933")
# worldclim_precip_may <- terra::resample(worldclim_precip_may,elevationmap,method="bilinear")
# worldclim_precip_jun <- terra::project(terra::crop(terra::rast("E:/Pittas/nichemodelling/worldclim/wc2.1_30s_prec/wc2.1_30s_prec_06.tif"),boundary),"epsg:6933")
# worldclim_precip_jun <- terra::resample(worldclim_precip_jun,elevationmap,method="bilinear")
# worldclim_precip_jul <- terra::project(terra::crop(terra::rast("E:/Pittas/nichemodelling/worldclim/wc2.1_30s_prec/wc2.1_30s_prec_07.tif"),boundary),"epsg:6933")
# worldclim_precip_jul <- terra::resample(worldclim_precip_jul,elevationmap,method="bilinear")
# worldclim_precip_aug <- terra::project(terra::crop(terra::rast("E:/Pittas/nichemodelling/worldclim/wc2.1_30s_prec/wc2.1_30s_prec_08.tif"),boundary),"epsg:6933")
# worldclim_precip_aug <- terra::resample(worldclim_precip_aug,elevationmap,method="bilinear")


#load predictor variable rasters
# pred_meantemp <- terra::rast("E:/Pittas/nichemodelling/krapp_paleoclim_meantemp/meantemp_interval0.tif")
# pred_temprange <- terra::rast("E:/Pittas/nichemodelling/krapp_paleoclim_maxtemp/maxtemp_interval0.tif") - terra::rast("E:/Pittas/nichemodelling/krapp_paleoclim_mintemp/mintemp_interval0.tif")
# pred_meanprecip <- terra::rast("E:/Pittas/nichemodelling/krapp_paleoclim_meanprecip/meanprecip_interval0.tif")
# pred_preciprange <- terra::rast("E:/Pittas/nichemodelling/krapp_paleoclim_maxprecip/maxprecip_interval0.tif") - terra::rast("E:/Pittas/nichemodelling/krapp_paleoclim_minprecip/minprecip_interval0.tif")
# pred_meantemp <- mean(worldclim_temp_may,worldclim_temp_jun,worldclim_temp_jul,worldclim_temp_aug)
# pred_temprange <- max(worldclim_temp_may,worldclim_temp_jun,worldclim_temp_jul,worldclim_temp_aug) - min(worldclim_temp_may,worldclim_temp_jun,worldclim_temp_jul,worldclim_temp_aug)
# pred_meanprecip <- mean(worldclim_precip_may,worldclim_precip_jun,worldclim_precip_jul,worldclim_precip_aug)
# pred_preciprange <- max(worldclim_precip_may,worldclim_precip_jun,worldclim_precip_jul,worldclim_precip_aug) - min(worldclim_precip_may,worldclim_precip_jun,worldclim_precip_jul,worldclim_precip_aug)
pred_bio01 <- terra::rast("E:/Pittas/nichemodelling/krapp_paleoclim_bio01/bio01_interval0.tif")
#pred_bio05 <- terra::rast("E:/Pittas/nichemodelling/krapp_paleoclim_bio05/bio05_interval0.tif")
#pred_bio06 <- terra::rast("E:/Pittas/nichemodelling/krapp_paleoclim_bio06/bio06_interval0.tif")
pred_bio12 <- terra::rast("E:/Pittas/nichemodelling/krapp_paleoclim_bio12/bio12_interval0.tif")
pred_bio13 <- terra::rast("E:/Pittas/nichemodelling/krapp_paleoclim_bio13/bio13_interval0.tif")
pred_bio14 <- terra::rast("E:/Pittas/nichemodelling/krapp_paleoclim_bio14/bio14_interval0.tif")
pred_bio15 <- terra::rast("E:/Pittas/nichemodelling/krapp_paleoclim_bio15/bio15_interval0.tif")
pred_elevation <- terra::rast("E:/Pittas/nichemodelling/raster_topo_new/topo_interval0.tif")
pred_bio16 <- terra::rast("E:/Pittas/nichemodelling/krapp_paleoclim_bio16/bio16_interval0.tif")
pred_bio17 <- terra::rast("E:/Pittas/nichemodelling/krapp_paleoclim_bio17/bio17_interval0.tif")
#pred_elevation <- terra::resample(pred_elevation,downsamp,method="bilinear")
pred_slope <- terra::rast("E:/Pittas/nichemodelling/raster_slope/slope_interval0.tif")
#pred_slope <- terra::resample(pred_slope,downsamp,method="bilinear")
pred_stack <- c(pred_bio01,pred_bio12,pred_bio13,pred_bio14,pred_bio15,pred_bio16,pred_bio17,pred_slope,pred_elevation)
names(pred_stack) <- c("bio01","bio12","bio13","bio14","bio15","bio16","bio17","slope","elevation")
#check for correlation between rasters
#totalpred <- c(pred_elevation,pred_bio01,pred_bio16,pred_bio17,pred_bio12,pred_bio15,pred_slope)
pred_corr <- terra::layerCor(pred_stack,fun="cor")
#convert unmasked predictor variables to raster format
#pred_elevation <- raster::raster(pred_elevation)
# pred_bio01 <- raster::raster(pred_bio01)
# pred_bio05 <- raster::raster(pred_bio05)
# pred_bio06 <- raster::raster(pred_bio06)
# pred_bio12 <- raster::raster(pred_bio12)
# pred_bio13 <- raster::raster(pred_bio13)
# pred_bio14 <- raster::raster(pred_bio14)
# pred_bio17 <- raster::raster(pred_bio17)
# pred_slope <- raster::raster(pred_slope)

# names(pred_elevation) <- c("Elevation")
# names(pred_bio01) <- c("Mean_AnnualTemp")
# names(pred_bio05) <- c("maxtemp_warmestmonth")
# names(pred_bio06) <- c("mintemp_coldestmonth")
# names(pred_bio12) <- c("Annual_Precip")
# names(pred_bio13) <- c("Precip_wettestmonth")
# names(pred_bio14) <- c("Precip_driestmonth")
# names(pred_bio15) <- c("Precip_Seasonality")
# names(pred_bio16) <- c("Precip_WettestQuarter")
# names(pred_bio17) <- c("Precip_DirestQuarter")
# names(pred_slope) <- c("Slope")

maskfile <- sf::st_read("E:/Pittas/nichemodelling/BWPitta_BreedingRange.shp")
#exclude meantemp since it's strongly correlated with elevation
pred_stack_masked <- terra::mask(pred_stack,maskfile)
#pred_stack_masked <- terra::rast(pred_stack_masked)

bwp_occ <- read.csv("E:/Pittas/nichemodelling/gbif_bwp_jun_jul/gbif_bwp_jun_jul.csv",header=T,sep=",")
bwp_occ <- bwp_occ %>% dplyr::select(x,y)

#run Maxent
tune.args = list(fc= c("L","LQ","Q","P","LQP","LP"),rm=1:5)
bwp_model <- ENMeval::ENMevaluate(occs = bwp_occ,
                                       envs = pred_stack,
                                       bg = bg,
                                       algorithm = "maxnet",
                                       partitions = 'block',
                                       tune.args = tune.args)
gc()

eval.results(bwp_model)
write.csv(eval.results(bwp_model),"E:/Pittas/nichemodelling/bwp_nichemodel/eval_results_newmodel.csv")
#best fitting model appears to be fc.LQP_rm.1, and results look non-ideal but not a deal-breaker
#auc.val.avg = 0.7480139
#cbi.val.avg = 0.89850 
#or.10p.avg = 0.1985882 
eval.results.partitions(bwp_model)
eval.predictions(bwp_model)
eval.models(bwp_model)[["fc.LQP_rm.1"]]$betas

bwp.mod.null <- ENMeval::ENMnulls(bwp_model,mod.settings = list(fc="LQP",rm=1),no.iter=100)
null.results(bwp.mod.null) %>% head()
null.emp.results(bwp.mod.null)
write.csv(null.emp.results(bwp.mod.null),"E:/Pittas/nichemodelling/bwp_nichemodel/null_results_newmodel.csv")
evalplot.nulls(bwp.mod.null, stats = c("or.10p", "auc.val"), plot.type = "histogram")

#empirical results don't seem significantly different from the null model

bwp_results <- eval.results(bwp_model)
#opt.aicc <- mangrove_results %>% filter(delta.AICc == 0)
opt.aicc_LQP1 <- bwp_results %>% dplyr::filter(tune.args == "fc.LQP_rm.1")

opt.aicc_LQP1

opt.seq <- mangrove_results %>%
  dplyr::filter(or.10p.avg == min(or.10p.avg)) %>%
  dplyr::filter(auc.val.avg == max(auc.val.avg))

opt.seq

best_bwp_aicc <- eval.models(bwp_model)[[as.character(opt.aicc_LQP1[1,]$tune.args)]]
#best_mangrove_aicc_alternative <- eval.models(mangrove_model)[[opt.aicc_LP1[1,]$tune.args]]
best_bwp_aicc$betas
#best_mangrove_aicc_alternative$betas
plot(best_bwp_aicc, type="cloglog")
#plot(best_mangrove_aicc_alternative, type="cloglog")

bestmod = which(bwp_model@results$AICc==min(bwp_model@results$AICc))

#plot model outputs
dev.off()
predicted_bwp <- eval.predictions(bwp_model)[[as.character(opt.aicc_LQP1[1,]$tune.args)]]
plot(predicted_bwp)
#terra::writeRaster(predicted_bwp,filename="E:/Pittas/nichemodelling/predicted_bwp/predicted_interval0.tif",overwrite=TRUE)
terra::writeRaster(predicted_bwp,filename="E:/Pittas/nichemodelling/predicted_bwp_new/predicted_interval0.tif",overwrite=TRUE)

#calculate model threshold
pr <- dismo::predict(pred_stack,bwp_model@models[[bestmod[1]]],type="cloglog",na.rm=T)
est.loc <- terra::extract(pr,bwp_model@occs[1:2])
est.bg <- terra::extract(pr,bwp_model@bg[1:2])
ev <- dismo::evaluate(p=est.loc[,2],a=est.bg[,2])
thr <- dismo::threshold(ev)
thr
#use the spec_sens threshold of 0.476742
pr_thr <- pr > 0.7952567
plot(pr_thr)
pr_thr_filt <- raster::clamp(pr_thr,lower=0.5,useValues=FALSE)
#terra::writeRaster(pr_thr,filename="E:/Pittas/nichemodelling/bwp_breedingrange_withland/predicted_interval0.tif",overwrite=TRUE)
pr_thr_filt <- terra::rast(pr_thr_filt)
pr_thr_filt_clean <- terra::sieve(pr_thr_filt,threshold=4,directions=8)
pr_thr_filt_clean <- terra::classify(pr_thr_filt_clean,rcl = matrix(c(0,NA),ncol=2))
terra::writeRaster(pr_thr_filt_clean,filename="E:/Pittas/nichemodelling/bwp_predicted_breedingrange/predicted_interval0.tif",overwrite=TRUE)

for (i in 0:799) {
  pred_elevation <- raster::raster(paste0("E:/Pittas/nichemodelling/raster_topo_new/topo_interval",800-i,".tif"))
  pred_meantemp <- raster::raster(paste0("E:/Pittas/nichemodelling/krapp_paleoclim_meantemp/meantemp_interval",800-i,".tif"))
  pred_meanprecip <- raster::raster(paste0("E:/Pittas/nichemodelling/krapp_paleoclim_meanprecip/meanprecip_interval",800-i,".tif"))
  pred_temprange <- raster::raster(paste0("E:/Pittas/nichemodelling/krapp_paleoclim_maxtemp/maxtemp_interval",800-i,".tif")) - raster::raster(paste0("E:/Pittas/nichemodelling/krapp_paleoclim_mintemp/mintemp_interval",800-i,".tif"))
  pred_preciprange <- raster::raster(paste0("E:/Pittas/nichemodelling/krapp_paleoclim_maxprecip/maxprecip_interval",800-i,".tif")) - raster::raster(paste0("E:/Pittas/nichemodelling/krapp_paleoclim_minprecip/minprecip_interval",800-i,".tif"))
  pred_slope <- raster::raster(paste0("E:/Pittas/nichemodelling/raster_slope/slope_interval",800-i,".tif"))
  names(pred_elevation) <- c("elevation")
  names(pred_meantemp) <- c("meantemp")
  names(pred_meanprecip) <- c("meanprecip")
  names(pred_temprange) <- c("temprange")
  names(pred_preciprange) <- c("preciprange")
  names(pred_slope) <- c("slope")
  
  pred_stack_interval <- raster::stack(pred_elevation,pred_meantemp,pred_meanprecip,pred_temprange,pred_preciprange,pred_slope)
  pr_interval <- predict(pred_stack_interval,bwp_model@models[[bestmod[1]]],type="cloglog")
  pr_interval_threshold <- pr_interval > 0.3773045
  #terra::writeRaster(pr_interval_threshold,paste0("E:/Pittas/nichemodelling/mangrove_predictions_withland/predicted_interval",i,".tif"),overwrite=TRUE)
  pr_interval_threshold_filt <- raster::clamp(pr_interval_threshold,lower=0.5,useValues=FALSE)
  pr_interval_threshold_filt <- terra::rast(pr_interval_threshold_filt)
  pr_interval_threshold_filt <- terra::sieve(pr_interval_threshold_filt,threshold=4,directions=8)
  pr_interval_threshold_filt <- terra::classify(pr_interval_threshold_filt,rcl = matrix(c(0,NA),ncol=2))
  terra::writeRaster(pr_interval_threshold_filt,paste0("E:/Pittas/nichemodelling/bwp_predicted_breedingrange/predicted_interval",800-i,".tif"),overwrite=TRUE)
  terra::writeRaster(pr_interval,paste0("E:/Pittas/nichemodelling/predicted_bwp/predicted_interval",800-i,".tif"),overwrite=TRUE)
  gc()
}

#time to plot the maps and graphs
cartolight <- basemap_ggplot(boundbox,map_service = "carto", map_type = "light")
sealvl <- read.csv("E:/Pittas/ddRAD/8_Biogeography/sealvl.csv")
intervalfile <- read.csv("E:/Pittas/nichemodelling/intervals.csv", header = TRUE)

#make mean temp, mean precip, mean temp range, and mean precip range file
climatemetrics <- data.frame()
sundaland_shapefile <- terra::vect("E:/Pittas/nichemodelling/Sundaland_hotspot.shp")
sundaland_shapefile <- terra::project(sundaland_shapefile,"epsg:6933")
indoburma_shapefile <- terra::vect("E:/Pittas/nichemodelling/IndoBurma_hotspot.shp")
indoburma_shapefile <- terra::project(indoburma_shapefile,"epsg:6933")
total_shapefile <- terra::project(terra::vect("E:/Pittas/nichemodelling/Indoburma_Sundaland_hotspots.shp"),"epsg:6933")

# for (i in 0:600) {
#   meantemp_sundaland <- as.numeric(global(terra::mask(terra::rast(paste0("E:/Pittas/nichemodelling/krapp_paleoclim_meantemp/meantemp_interval",i,".tif")),sundaland_shapefile),"mean",na.rm=TRUE))
#   meantemp_indoburma <- as.numeric(global(terra::mask(terra::rast(paste0("E:/Pittas/nichemodelling/krapp_paleoclim_meantemp/meantemp_interval",i,".tif")),indoburma_shapefile),"mean",na.rm=TRUE))
#   maxtemp_sundaland <- as.numeric(global(terra::mask(terra::rast(paste0("E:/Pittas/nichemodelling/krapp_paleoclim_maxtemp/maxtemp_interval",i,".tif")),sundaland_shapefile),"max",na.rm=TRUE))
#   maxtemp_indoburma <- as.numeric(global(terra::mask(terra::rast(paste0("E:/Pittas/nichemodelling/krapp_paleoclim_maxtemp/maxtemp_interval",i,".tif")),indoburma_shapefile),"max",na.rm=TRUE))
#   mintemp_sundaland <- as.numeric(global(terra::mask(terra::rast(paste0("E:/Pittas/nichemodelling/krapp_paleoclim_mintemp/mintemp_interval",i,".tif")),sundaland_shapefile),"min",na.rm=TRUE))
#   mintemp_indoburma <- as.numeric(global(terra::mask(terra::rast(paste0("E:/Pittas/nichemodelling/krapp_paleoclim_mintemp/mintemp_interval",i,".tif")),indoburma_shapefile),"min",na.rm=TRUE))
#   meanprecip_val <- as.numeric(global(terra::rast(paste0("E:/Pittas/nichemodelling/krapp_paleoclim_meanprecip/meanprecip_interval",i,".tif")),"mean",na.rm=TRUE))
#   maxprecip_val <- as.numeric(global(terra::rast(paste0("E:/Pittas/nichemodelling/krapp_paleoclim_maxprecip/maxprecip_interval",i,".tif")),"max",na.rm=TRUE))
#   minprecip_val <- as.numeric(global(terra::rast(paste0("E:/Pittas/nichemodelling/krapp_paleoclim_minprecip/minprecip_interval",i,".tif")),"min",na.rm=TRUE))
#   climatemetrics <- rbind(climatemetrics,c(i,meantemp_val,maxtemp_val,mintemp_val,meanprecip_val,maxprecip_val,minprecip_val))
# }

# write.csv(climatemetrics,"E:/Pittas/nichemodelling/climate_metrics.csv")
# climatemetrics <- read.csv("E:/Pittas/nichemodelling/climate_metrics.csv",header = TRUE, sep=",")

interval_bwp <- terra::mask(terra::rast("E:/Pittas/nichemodelling/worldclim/wc2.1_30s_tavg/"),total_shapefile)
interval_bwp <- terra::project(interval_bwp,"epsg:3857")
interval_bwp <- as.factor(interval_bwp)
interval_bwp <- sf::st_union(sf::st_as_sf(terra::as.polygons(interval_bwp)))
interval_upptime <- filter(intervalfile,Interval==i)$UpperTimeBound
interval_lowtime <- filter(intervalfile,Interval==i)$LowerTimeBound
interval_meantemp <- terra::mask(terra::rast(paste0("E:/Pittas/nichemodelling/krapp_paleoclim_meantemp/meantemp_interval",i,".tif")),total_shapefile)
interval_meantemp <- terra::project(interval_meantemp,"epsg:3857")
interval_meanprecip <- terra::mask(terra::rast(paste0("E:/Pittas/nichemodelling/krapp_paleoclim_meanprecip/meanprecip_interval",i,".tif")),total_shapefile)
interval_meanprecip <- terra::project(interval_meanprecip,"epsg:3857")
temp_ppt <- c(interval_meantemp,interval_meanprecip)
names(temp_ppt) <- c("temp","ppt")
temp_ppt_df <- temp_ppt |>
  as.data.frame(xy=TRUE)
temp_ppt_df$temp_bin <- cut(temp_ppt_df$temp,breaks=c(min(temp_ppt_df$temp),20,24,28,max(temp_ppt_df$temp)),include.lowest=TRUE)
temp_ppt_df$ppt_bin <- cut(temp_ppt_df$ppt,breaks=c(min(temp_ppt_df$ppt),100,300,600,max(temp_ppt_df$ppt)),include.lowest=TRUE)
data <- bi_class(temp_ppt_df,
                 x = temp_bin, 
                 y = ppt_bin, 
                 style = "quantile", dim = 4)
data |> 
  count(bi_class) |> 
  ggplot(aes(x = bi_class, y = n)) +
  geom_col() +  # Create a bar plot to show the count of each bivariate class
  labs(title = "Distribution of Bivariate Classes", x = "Bivariate Class", y = "Frequency")
# Set the color palette for the bivariate map
pallet <- "BlueOr"
landshape <- sf::st_transform(sf::st_read(paste0("E:/Pittas/nichemodelling/shapefile_clean/interval",i,".shp")),crs="epsg:3857")

# Create the bivariate map using ggplot2
temp_precip_map <- ggplot() +
  theme_void(base_size = 14) +  # Set a minimal theme for the map
  # Plot the bivariate raster data with appropriate fill color based on bivariate classes
  geom_raster(data = data, mapping = aes(x = x, y = y, fill = bi_class), color = NA, linewidth = 0.1, show.legend = FALSE) +
  # Apply the bivariate color scale using the selected palette and dimensions
  bi_scale_fill(pal = pallet, dim = 4, flip_axes = FALSE, rotate_pal = FALSE) +
  # Overlay the first administrative level boundaries of the United Kingdom
  # new_scale_fill() +
  # geom_sf(data=interval_bwp,colour=NA,alpha=0.8) +
  geom_sf(data = landshape, fill = NA, color = "black", linewidth = 0.20)


# Create the legend for the bivariate map
legend <- bi_legend(pal = pallet,   
                    flip_axes = FALSE,
                    rotate_pal = FALSE,
                    dim = 4,
                    xlab = "Temperature (\u00B0C)",
                    ylab = "Precipitation (mm)",
                    size = 8)

finalPlot <- ggdraw() +
  draw_plot(temp_precip_map, 0, 0, 1, 1) +  # Draw the main map plot
  draw_plot(legend, 0.05, 0.05, 0.28, 0.28)  # Draw the legend in the specified position

metrics <- data.frame()
#calculate area and environmental metrics
for (i in 0:798) {
  print(paste0("Evaluating interval ",i,"..."))
  interval_temp <- terra::mask(terra::rast(paste0("E:/Pittas/nichemodelling/krapp_paleoclim_meantemp/meantemp_interval",i,".tif")),total_shapefile)
  interval_precip <- terra::mask(terra::rast(paste0("E:/Pittas/nichemodelling/krapp_paleoclim_meanprecip/meanprecip_interval",i,".tif")),total_shapefile)
  bwp_range <- terra::mask(terra::rast(paste0("E:/Pittas/nichemodelling/bwp_predicted_breedingrange/predicted_interval",i,".tif")),total_shapefile)
  bwp_area <- terra::expanse(bwp_range,unit="m",transform=FALSE)
  bwp_earlier <- terra::mask(terra::rast(paste0("E:/Pittas/nichemodelling/bwp_predicted_breedingrange/predicted_interval",i+1,".tif")),total_shapefile)
  bwp_overlap <- terra::clamp(bwp_earlier+bwp_range,lower=1.5,value=FALSE)
  bwp_overlap_area <- terra::expanse(bwp_overlap,unit="m",transform=FALSE)
  bwp_overlap_percent <- (bwp_overlap_area$area/bwp_area$area * 100)
  mean_temp <- as.numeric(terra::global(interval_temp,"mean",na.rm=T))
  mean_precip <- as.numeric(terra::global(interval_precip,"mean",na.rm=T))
  metrics <- rbind(metrics,c(i,bwp_area$area,mean_temp,mean_precip,bwp_overlap_area$area,bwp_overlap_percent))
  gc()
}
names(metrics) <- c("Interval","BWP_Area","Temp","Precip","Overlap_Area","Overlap_Percent")
write.csv(metrics,"E:/Pittas/nichemodelling/climate_metrics.csv")

for (i in 0:799) {
  interval_bwp <- terra::mask(terra::rast(paste0("E:/Pittas/nichemodelling/bwp_predicted_breedingrange/predicted_interval",i,".tif")),total_shapefile)
  interval_bwp <- terra::project(interval_bwp,"epsg:3857")
  interval_upptime <- filter(intervalfile,Interval==i)$UpperTimeBound
  interval_lowtime <- filter(intervalfile,Interval==i)$LowerTimeBound
  interval_meantemp <- terra::mask(terra::rast(paste0("E:/Pittas/nichemodelling/krapp_paleoclim_meantemp/meantemp_interval",i,".tif")),total_shapefile)
  interval_meantemp <- terra::project(interval_meantemp,"epsg:3857")
  interval_meanprecip <- terra::mask(terra::rast(paste0("E:/Pittas/nichemodelling/krapp_paleoclim_meanprecip/meanprecip_interval",i,".tif")),total_shapefile)
  interval_meanprecip <- terra::project(interval_meanprecip,"epsg:3857")
  temp_ppt <- c(interval_meantemp,interval_meanprecip)
  names(temp_ppt) <- c("temp","ppt")
  temp_ppt_df <- temp_ppt |>
    as.data.frame(xy=TRUE)
  temp_ppt_df$temp_bin <- cut(temp_ppt_df$temp,breaks=c(min(temp_ppt_df$temp),20,24,28,max(temp_ppt_df$temp)),include.lowest=TRUE)
  temp_ppt_df$ppt_bin <- cut(temp_ppt_df$ppt,breaks=c(min(temp_ppt_df$ppt),100,300,600,max(temp_ppt_df$ppt)),include.lowest=TRUE)
  data <- bi_class(temp_ppt_df,
                   x = temp_bin, 
                   y = ppt_bin, 
                   style = "quantile", dim = 4)
  landshape <- sf::st_transform(sf::st_read(paste0("E:/Pittas/nichemodelling/shapefile_clean/interval",i,".shp")),crs="epsg:3857")
  
  # Create the bivariate map using ggplot2
  temp_precip_map <- ggplot() +
    theme_void(base_size = 14) +  # Set a minimal theme for the map
    # Plot the bivariate raster data with appropriate fill color based on bivariate classes
    geom_raster(data = data, mapping = aes(x = x, y = y, fill = bi_class), color = NA, linewidth = 0.1, show.legend = FALSE) +
    # Apply the bivariate color scale using the selected palette and dimensions
    bi_scale_fill(pal = pallet, dim = 4, flip_axes = FALSE, rotate_pal = FALSE) +
    # Overlay the first administrative level boundaries of the United Kingdom
    # new_scale_fill() +
    # geom_sf(data=interval_bwp,colour=NA,alpha=0.8) +
    geom_sf(data = landshape, fill = NA, color = "black", linewidth = 0.20)
    

  # Create the legend for the bivariate map
  legend <- bi_legend(pal = pallet,   
                      flip_axes = FALSE,
                      rotate_pal = FALSE,
                      dim = 4,
                      xlab = "Temperature (\u00B0C)",
                      ylab = "Precipitation (mm)",
                      size = 5)
  
  finalPlot <- ggdraw() +
    draw_plot(temp_precip_map, 0, 0, 1, 1) +  # Draw the main map plot
    draw_plot(legend, 0.05, 0.05, 0.28, 0.28)  # Draw the legend in the specified position
  
  map <- cartolight+
    geom_spatvector(data=landshape,fill="green4",alpha=0.4,colour=NA)+
    geom_spatraster(data=interval_bwp,alpha=0.8)+
    theme_void()+
    annotate("text",label=paste0("Time: ",interval_lowtime,"-",interval_upptime," kya\nMean Sea Level: ",filter(intervalfile,Interval==i)$MeanDepth,"m"),x=9600000,y=-0.8E06,hjust=0)
  sealvlgraph <- ggplot(sealvl,aes(x=Time, y=Sealevel_Corrected))+
    geom_line(colour="blue4",linewidth=1.2)+xlim(0,799)+
    xlab("Time (kya)")+
    ylab("Sea level (m)")+
    ggtitle("Sea level from 0 to 799 kya")+
    annotate("rect",xmin=interval_lowtime,xmax=interval_upptime,ymin=-Inf,ymax=Inf,fill="black",alpha=0.95)
  tempgraph <- ggplot(metrics,aes(x=Interval,y=Temp))+
    geom_line(colour="orange2",linewidth=1.2)+xlim(0,799)+
    xlab("Time (kya)")+
    ylab("Mean Monthly Temp (\u00B0C)")+
    ggtitle("Mean breeding season temperature")+
    annotate("rect",xmin=interval_lowtime,xmax=interval_upptime,ymin=-Inf,ymax=Inf,fill="black",alpha=0.95)
  precipgraph <- ggplot(metrics,aes(x=Interval,y=Precip))+
    geom_line(colour="skyblue3",linewidth=1.2)+xlim(0,799)+
    xlab("Time (kya)")+
    ylab("Mean Monthly Precip (mm)")+
    ggtitle("Mean breeding season precipitation")+
    annotate("rect",xmin=interval_lowtime,xmax=interval_upptime,ymin=-Inf,ymax=Inf,fill="black",alpha=0.95)
  bwpareagraph <- ggplot(metrics,aes(x=Interval,y=BWP_Area))+
    geom_line(colour="black",linewidth=1.2)+xlim(0,799)+
    xlab("Time (kya)")+
    ylab("Total Suitable Area (sqm)")+
    ggtitle("Total suitable breeding area")+
    annotate("rect",xmin=interval_lowtime,xmax=interval_upptime,ymin=-Inf,ymax=Inf,fill="black",alpha=0.95)
  metric_graphs <- ggarrange(sealvlgraph,bwpareagraph,tempgraph,precipgraph,nrow=2,ncol=2,align="hv",common.legend = TRUE)
  mapgraphs <- ggarrange(map,finalPlot,ncol=2)
  ggarrange(mapgraphs,metric_graphs,nrow=2)
  ggsave(paste0("E:/Pittas/nichemodelling/bwp_maps/interval",i,".png"),width = 8, height = 8, units = "in",bg="white")
  gc()
}

#calculate average range of BWPitta over Pleistocene time
total_shapefile <- terra::project(terra::vect("E:/Pittas/nichemodelling/Indoburma_Sundaland_hotspots.shp"),"epsg:6933")
bwp_cumulativerange <- terra::rast("E:/Pittas/nichemodelling/bwp_predicted_breedingrange/predicted_interval0.tif")
for (i in 1:799) {
  newmap <- terra::rast(paste0("E:/Pittas/nichemodelling/bwp_predicted_breedingrange/predicted_interval",i,".tif"))
  combinedmaps <- terra::sds(bwp_cumulativerange,newmap)
  bwp_cumulativerange <- terra::app(combinedmaps,fun="sum",na.rm=TRUE)
}
bwp_cumulativerange <- bwp_cumulativerange/800
bwp_cumulativerange <- terra::mask(bwp_cumulativerange,total_shapefile)
terra::writeRaster(bwp_cumulativerange,"E:/Pittas/nichemodelling/bwp_nichemodel/bwp_cumulativerange.tif",overwrite=TRUE)

#calculate average range of BWPitta over Pleistocene time
total_shapefile <- terra::project(terra::vect("E:/Pittas/nichemodelling/Indoburma_Sundaland_hotspots.shp"),"epsg:6933")
bwp_cumulativerange <- terra::rast("E:/Pittas/nichemodelling/bwp_predicted_breedingrange/predicted_interval0.tif")
for (i in 1:600) {
  newmap <- terra::rast(paste0("E:/Pittas/nichemodelling/bwp_predicted_breedingrange/predicted_interval",i,".tif"))
  combinedmaps <- terra::sds(bwp_cumulativerange,newmap)
  bwp_cumulativerange <- terra::app(combinedmaps,fun="sum",na.rm=TRUE)
  gc()
}
bwp_cumulativerange <- bwp_cumulativerange/601
bwp_cumulativerange <- terra::mask(bwp_cumulativerange,total_shapefile)
terra::writeRaster(bwp_cumulativerange,"E:/Pittas/nichemodelling/bwp_nichemodel/bwp_cumulativerange_600kya.tif",overwrite=TRUE)

############################

#Rerun the models with no elevation -- actually don't. The model metrics are dogshit

pred_stack_noelev <- raster::stack(pred_meantemp,pred_meanprecip,pred_temprange,pred_preciprange,pred_slope)
bwp_occ <- read.csv("E:/Pittas/nichemodelling/gbif_bwp_jun_jul/gbif_bwp_jun_jul.csv",header=T,sep=",")
bwp_occ <- bwp_occ %>% dplyr::select(x,y)

#run Maxent
tune.args = list(fc= c("L","LQ","Q","P","LQP","LP"),rm=1:5)
bwp_model_noelev <- ENMeval::ENMevaluate(occs = bwp_occ,
                                  envs = pred_stack_noelev,
                                  bg = bg,
                                  algorithm = "maxent.jar",
                                  partitions = 'block',
                                  tune.args = tune.args)

eval.results(bwp_model_noelev)
write.csv(eval.results(bwp_model_noelev),"E:/Pittas/nichemodelling/bwp_nichemodel/eval_results_noelev.csv")
eval.results.partitions(bwp_model_noelev)
eval.predictions(bwp_model_noelev)
eval.models(bwp_model_noelev)[["fc.LQP_rm.1"]]$betas

bwp.mod.null <- ENMeval::ENMnulls(bwp_model,mod.settings = list(fc="LQP",rm=1),no.iter=100)
null.results(bwp.mod.null) %>% head()
null.emp.results(bwp.mod.null)
write.csv(null.emp.results(bwp.mod.null),"E:/Pittas/nichemodelling/bwp_nichemodel/null_results.csv")
evalplot.nulls(bwp.mod.null, stats = c("or.10p", "auc.val"), plot.type = "histogram")

#empirical results don't seem significantly different from the null model

bwp_results <- eval.results(bwp_model)
#opt.aicc <- mangrove_results %>% filter(delta.AICc == 0)
opt.aicc_LQP1 <- bwp_results %>% dplyr::filter(tune.args == "fc.LQP_rm.1")

opt.aicc_LQP1

opt.seq <- mangrove_results %>%
  dplyr::filter(or.10p.avg == min(or.10p.avg)) %>%
  dplyr::filter(auc.val.avg == max(auc.val.avg))

opt.seq

best_bwp_aicc <- eval.models(bwp_model)[[opt.aicc_LQP1[1,]$tune.args]]
#best_mangrove_aicc_alternative <- eval.models(mangrove_model)[[opt.aicc_LP1[1,]$tune.args]]
best_bwp_aicc$betas
#best_mangrove_aicc_alternative$betas
plot(best_bwp_aicc, type="cloglog")
#plot(best_mangrove_aicc_alternative, type="cloglog")

bestmod = which(bwp_model@results$AICc==min(bwp_model@results$AICc))

#plot model outputs
dev.off()
predicted_bwp <- eval.predictions(bwp_model)[[opt.aicc_LQP1[1,]$tune.args]]
plot(predicted_bwp)
terra::writeRaster(predicted_bwp,filename="E:/Pittas/nichemodelling/predicted_bwp/predicted_interval0.tif",overwrite=TRUE)

#calculate model threshold
pr <- predict(pred_stack,bwp_model@models[[bestmod[1]]],type="cloglog")
est.loc <- raster::extract(pr,bwp_model@occs[1:2])
est.bg <- raster::extract(pr,bwp_model@bg[1:2])
ev <- dismo::evaluate(est.loc,est.bg)
thr <- dismo::threshold(ev)
thr
pr_thr <- pr > 0.3775883
plot(pr_thr)
pr_thr_filt <- raster::clamp(pr_thr,lower=0.5,useValues=FALSE)
#terra::writeRaster(pr_thr,filename="E:/Pittas/nichemodelling/bwp_breedingrange_withland/predicted_interval0.tif",overwrite=TRUE)
pr_thr_filt <- terra::rast(pr_thr_filt)
pr_thr_filt_clean <- terra::sieve(pr_thr_filt,threshold=4,directions=8)
pr_thr_filt_clean <- terra::classify(pr_thr_filt_clean,rcl = matrix(c(0,NA),ncol=2))
terra::writeRaster(pr_thr_filt_clean,filename="E:/Pittas/nichemodelling/bwp_predicted_breedingrange/predicted_interval0.tif",overwrite=TRUE)


##################################################
# Estimating BWP Refugia using Krapp BIOME4 data #
##################################################

#Instead of running niche models, perhaps just use the BIOME4 paleoreconstructions from Krapp et al to see if the results are similar
#Note, however, that resolution will be very coarse, with avg pixel size of approx 50-60km2

biome_raster <- terra::rast(paste0("E:/Pittas/nichemodelling/krapp_paleoclim_biomes/biomes_interval0.tif"))
#reclassify all non-tropical deciduous forest as NA
biome_raster <- terra::classify(biome_raster,rcl=matrix(c(4,28,NA),ncol=3),include.lowest=T)
biome_raster <- terra::classify(biome_raster,rcl=matrix(c(0,1,NA),ncol=3),include.lowest=T)
deciduous <- terra::classify(biome_raster,rcl=matrix(c(2,3,1),ncol=3),include.lowest=T)

for (i in 1:799) {
  biome_raster <- terra::rast(paste0("E:/Pittas/nichemodelling/krapp_paleoclim_biomes/biomes_interval",i,".tif"))
  #reclassify all non-tropical deciduous forest as NA
  biome_raster <- terra::classify(biome_raster,rcl=matrix(c(4,28,NA),ncol=3),include.lowest=T)
  biome_raster <- terra::classify(biome_raster,rcl=matrix(c(0,1,NA),ncol=3),include.lowest=T)
  biome_raster <- terra::classify(biome_raster,rcl=matrix(c(2,3,1),ncol=3),include.lowest=T)
  deciduous <- terra::app(c(deciduous,biome_raster),fun="sum",na.rm=T)
  gc()
}

terra::writeRaster(deciduous,"E:/Pittas/nichemodelling/BWP_biome_cumulative.tif")
terra::writeRaster(deciduous/799,"E:/Pittas/nichemodelling/BWP_biome_percentage.tif")

#maybe try getting data for a specific location along the Isthmus of Kra
#coordinates to sample: 
#IOK: 7.8933,99.3070
#Kaeng Krachan: 12.3909,99.4234
#Chiang Rai: 20.03753, 99.89913
#Laos: 18.02675, 104.79609
#Vietnam: 11.4412500035, 107.3859247969

# locations <- data.frame(
#   name = c("IsthmusKra", "KaengKrachan", "ChiangRai", "Laos", "Vietnam"),
#   longitude = c(99.3070, 99.4234, 99.89913, 104.79609, 107.3859247969), latitude = c(7.8933, 12.3909, 20.03753, 18.02675, 11.4412500035),
#   time_bp = list(min=-600000,max=0)
# )
# 
# locations_ts <- location_series(
#   x = locations,
#   bio_variables = c("bio01","bio05","bio06","bio12","bio13","bio14","temperature_05","temperature_06","temperature_07","temperature_08","precipitation_05","precipitation_06","precipitation_07","precipitation_08"),
#   dataset = "Krapp2021", 
#   nn_interpol = TRUE,
#   time_bp = list(min=-600000,max=0)
# )
# 
# ggplot(data = locations_ts_new, aes(x = time_bp, y = meanbreedtemp, group = name)) +
#   geom_line(aes(col = name))
