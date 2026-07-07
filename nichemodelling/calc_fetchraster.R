library(terra)
library(fetchr)

calc_fetchraster <- function(interval) {
	land_raster <- terra::rast(paste0("/carc/scratch/projects/andersen2016005/pitta/nichemodelling/raster_flat/interval",interval,".tif"))
	fetchraster <- fetchr::get_fetch(r=land_raster, max_dist=300000, in_parallel=TRUE, verbose=TRUE, func = "mean")
	terra::writeRaster(fetchraster,paste0("/carc/scratch/projects/andersen2016005/pitta/nichemodelling/shorefetch/fetchraster_interval",interval,".tif"))
)

calc_fetchraster(0)