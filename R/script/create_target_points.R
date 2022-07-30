# ------------------------------------------------------------------------------
# Create target points for crops, grassland and trees for the area of interest
# ------------------------------------------------------------------------------
generatePoints <- function(lu_path, output_path){

  rm(list = ls())

  message(noquote("Loading necessary packages ..."))
  require(raster)
  require(rgdal)
  require(sp)
  require(dplyr)

  setwd(lu_path)

  message(noquote("Reading landuse raster ..."))
  lu <- raster("Land_Cover_10class_AOI.tif")
  message(noquote("Generating target points ..."))
  point <-
    rasterToPoints(
      x = lu,
      fun = function(x) {
        x == 2 | x == 3 | x == 4
      },
      spatial = TRUE
    )
  setwd(output_path)
  message(noquote("Writing target points ..."))
  writeOGR(point, dsn = ".", layer = "eth_points", driver = "ESRI Shapefile", overwrite_layer = T)
}
