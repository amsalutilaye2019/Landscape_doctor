
# ------------------------------------------------------------------------------
# Script to convert Landuse class to FAO class
# ------------------------------------------------------------------------------

lu_ESA_FAO <- function(landuse_path, output_path){

  message(noquote("Loading necessary packages ..."))
  require(raster)
  require(rgdal)
  require(sp)
  require(dplyr)

  WD_LU <- landuse_path
  WD_OUT <- output_path
  
  setwd(WD_LU)
  message(noquote("Reading landuse raster ..."))
  ESA_LU_AOI <- raster("2016_LULC_Level_I.tif")
  
  #plot(ESA_LU_AOI)
  #Reclassify ESA LAND USE to FAO LAND USE classes
  #ours  FAO
  #0 = 0	  No Data
  #10 = 1 Artificial
  #4 = 2 Croplands
  #5 = 3 Grassland
  #1,2,9 = 4 Tree Covered
  #3  = 5 Shrubs Covered
  #6 = 9 Baresoil
  #8 = 11 Waterbodies
  #7 = 13 Paddy fields(rice/ flooded crops)
  # Create a reclassification matrix. "Is" to "become"
  
  is <- c(0, 10, 4, 5, 1, 2, 9, 3, 6, 8, 7)
  become <- c(0, 1, 2, 3, 4, 4, 4, 5, 9, 11, 13)
  recMat <- matrix(c(is, become), ncol = 2, nrow = 11)
  
  #Reclassify
  message(noquote("reclassifying landuse to FAO class ..."))
  ESA_FAO <- reclassify(ESA_LU_AOI, recMat)
  setwd(WD_SOC)
  
  #Save Landuse raster
  setwd(WD_OUT)
  message(noquote("Writing the reclassed landuse raster ..."))
  writeRaster(ESA_FAO,
              filename = "Land_Cover_10class_AOI.tif",
              format = 'GTiff',
              overwrite = TRUE)
}

