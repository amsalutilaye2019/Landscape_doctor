# ------------------------------------------------------------------------------
# Script for generating carbon stock maps
# ------------------------------------------------------------------------------

rm(list = ls())
generateMap <-
  function(aoi_path,
           aoi,
           forward_path,
           forward_aoi,
           soc_aoi,
           soc,
           init_year,
           pred_year,
           output_path,
           country_iso) {
    message(noquote("Loading necessary packages ..."))
    require(raster)
    require(rgdal)
    require(sp)
    require(dplyr)
    
    #Open FORWARD vector
    setwd(forward_path)
    message(noquote("Reading forward target points ..."))
    FOWARD <- readOGR(forward_aoi)
    
    #Open SOC MAP (master layer)
    setwd(soc_aoi)
    SOC_MAP <- raster(soc)
    
    #Creates emtpy raster
    message(noquote("create an empty raster ..."))
    empty_raster <- SOC_MAP * 0
    
    # Open the country vector boundaries
    setwd(aoi_path)
    message(noquote("Reading aoi ..."))
    Country <- readOGR(aoi)
    
    # Cut the raster with the country vector
    message(noquote("cropping empty raster by aoi ..."))
    Country_raster <- crop(empty_raster, Country)
    
    # Replace Na values for zero values
    FOWARD@data[is.na(FOWARD@data)] <- 0
    
    # Points to Raster BAU
    setwd(output_path)
    message(noquote("Writing rasters ..."))
    Country_BAU_2023_Map <-
      rasterize(FOWARD, Country_raster , FOWARD$SOC_BAU_20, updateValue = 'all')
    writeRaster(
      Country_BAU_2023_Map,
      filename = paste0(country_iso, "_finalSOC_BAU_",pred_year,".tif"),
      format = "GTiff",
      overwrite = TRUE
    )
    
    # Points to Raster Low Scenario
    Country_Lwr_2030_Map <-
      rasterize(FOWARD, Country_raster , FOWARD$Lw_Sc, updateValue = 'all')
    writeRaster(
      Country_Lwr_2030_Map,
      filename = paste0(country_iso, "_finalSOC_SSM1_",pred_year,".tif"),
      format = "GTiff",
      overwrite = TRUE
    )
    
    # Points to Raster Med Scenario
    Country_Med_2030_Map <-
      rasterize(FOWARD, Country_raster , FOWARD$Md_Sc, updateValue = 'all')
    writeRaster(
      Country_Med_2030_Map,
      filename = paste0(country_iso, "_finalSOC_SSM2_",pred_year,".tif"),
      format = "GTiff",
      overwrite = TRUE
    )
    
    # Points to Raster High Scenario
    Country_Hgh_2030_Map <-
      rasterize(FOWARD, Country_raster , FOWARD$Hgh_S, updateValue = 'all')
    writeRaster(
      Country_Hgh_2030_Map,
      filename = paste0(country_iso, "_finalSOC_SSM3_", pred_year, ".tif"),
      format = "GTiff",
      overwrite = TRUE
    )
    
    # Points to Raster initial SOC (t0) 2018/2020
    Country_SOC_2018_Map <-
      rasterize(FOWARD, Country_raster , FOWARD$SOC_t0, updateValue = 'all')
    writeRaster(
      Country_SOC_2018_Map,
      filename = paste0(country_iso, "_SOC_T0_" , init_year, ".tif"),
      format = "GTiff",
      overwrite = TRUE
    )
  }

#Example
generateMap("C:\\Users\\ATilaye\\Documents\\01My_Docs\\01CIAT\\08RothC\\Ethiopia_for_ld3\\data\\input",
           "siya_debir.shp",
           "C:\\Users\\ATilaye\\Documents\\01My_Docs\\01CIAT\\08RothC\\Ethiopia_for_ld3\\data\\final\\woreda_test\\output",
           "FOWARD_woreda.shp",
           "C:\\Users\\ATilaye\\Documents\\01My_Docs\\01CIAT\\08RothC\\Ethiopia_for_ld3\\data\\final\\woreda_test\\stacks",
           "soc_woreda.tif",
           2018, #initial year @t0
           2023,
           "C:\\Users\\ATilaye\\Documents\\01My_Docs\\01CIAT\\08RothC\\Ethiopia_for_ld3\\data\\final\\woreda_test\\output",
           "siya")
