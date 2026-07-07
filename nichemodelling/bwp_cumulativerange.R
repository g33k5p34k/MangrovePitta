library(terra)
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