# Create a points for crops, grassland and trees

generatePoints <- function(lu_path, output_path){
  setwd(lu_path)
  lu <- raster("Land_Cover_10class_AOI.tif")
  point <-
    rasterToPoints(
      x = lu,
      fun = function(x) {
        x == 2 | x == 3 | x == 4
      },
      spatial = TRUE
    )
  setwd(output_path)
  writeOGR(point, dsn = ".", layer = "eth_points", driver = "ESRI Shapefile", overwrite_layer = T)
}

#Example
generatePoints("C:\\Users\\ATilaye\\Documents\\01My_Docs\\01CIAT\\08RothC\\Ethiopia_for_ld3\\data\\final\\script9", "C:\\Users\\ATilaye\\Documents\\01My_Docs\\01CIAT\\08RothC\\Ethiopia_for_ld3\\data\\final")
