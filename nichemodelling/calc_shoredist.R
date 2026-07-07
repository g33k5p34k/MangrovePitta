calc_shoredist <- function(interval) {
  library(sf)
  library(terra)
  library(smoothr)
  landvect <- sf::st_read(paste0("E:/Pittas/nichemodelling/shapefile/interval",interval,".shp"))
  #close all internal holes so distance is only calculated from shoreline
  landvect_noholes <- smoothr::fill_holes(landvect,threshold=units::set_units(100000, km^2))
  landvect_noholes <- sf::st_union(landvect_noholes)
  sf::st_write(landvect_noholes,paste0("E:/Pittas/nichemodelling/shapefile_clean/interval",interval,".shp"))
  #rasterise the nohole vector
  template <- terra::rast(terra::vect(landvect_noholes),res=1000)
  landraster <- terra::rasterize(vect(landvect_noholes),template)
  m <- c(NA,1,1,NA)
  convmat <- matrix(m,ncol=2,byrow=TRUE)
  landraster_inverse <- terra::classify(landraster,convmat)
  shoredist <- terra::distance(landraster_inverse)
  terra::writeRaster(shoredist,paste0("E:/Pittas/nichemodelling/raster_shoredist/shoredist_interval",interval,".tif"))
}

for (i in 129:350) {
  calc_shoredist(i)
}