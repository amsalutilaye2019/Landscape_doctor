# ------------------------------------------------------------------------------
# Create points for crops, grassland and trees
# ------------------------------------------------------------------------------

generatePoints <- function(lu_path, lu_raster, output_path, output_lyr){
  setwd(lu_path)
  lu <- raster(lu_raster)
  point <-
    rasterToPoints(
      x = lu,
      fun = function(x) {
        x == 2 | x == 3 | x == 4  # 2 =crops, 3 = grassland, 4 = tree crops
      },
      spatial = TRUE
    )
  setwd(output_path)
  writeOGR(point, dsn = ".", layer = output_lyr, driver = "ESRI Shapefile", overwrite_layer = T)
}

#Example
generatePoints("C:\\Users\\ATilaye\\Documents\\01My_Docs\\01CIAT\\08RothC\\Ethiopia_for_ld3\\data\\final\\woreda_test","woreda_lu.tif", "C:\\Users\\ATilaye\\Documents\\01My_Docs\\01CIAT\\08RothC\\Ethiopia_for_ld3\\data\\final\\woreda_test","lu_woreda2")
