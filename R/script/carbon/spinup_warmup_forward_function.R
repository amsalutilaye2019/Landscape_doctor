# ------------------------------------------------------------------------------
# Prepare stacks for the spin up, warm up & forward process of the Roth C Model
# ------------------------------------------------------------------------------

spinWarmForward_Stack <- function(aoi_path, raster_path, output_path){

message(noquote("Loading necessary packages ..."))
require(raster)
require(rgdal)
require(dplyr)
require(sp)

WD_AOI <- aoi_path
WD_RASTERS <- raster_path
WD_OUT <- output_path

#Open the shapefile
setwd(WD_AOI)
message(noquote("Reading aoi shapefile ..."))
AOI<-readOGR("eth.shp")

#Open SOC MAP
setwd(WD_RASTERS)
message(noquote("Reading SOC raster ..."))
SOC_MAP_AOI <- raster("SOC_MAP_AOI.tif")

#Clay layer
message(noquote("Reading clay raster ..."))
Clay <- raster("Clay_WA_AOI.tif")
message(noquote("synchronizing clay raster with soc..."))
Clay_AOI <- Clay %>% crop(AOI) %>%
            resample(SOC_MAP_AOI, method = 'bilinear') %>% 
            crop(AOI)

#Precipitation layer 
message(noquote("Reading precipitation 1981-2000 stack ..."))
PREC <- stack("Prec_Stack_81-00_CRU.tif")
message(noquote("synchronizing prec81-00 with soc..."))
PREC_AOI <- PREC %>% crop(AOI) %>%
            resample(SOC_MAP_AOI, method = "bilinear") %>%
            mask(AOI)

#Temperatures layer 
message(noquote("Reading temerature 1981-2000 stack ..."))
TEMP <- stack("Temp_Stack_81-00_CRU.tif")
message(noquote("synchronizing temp81-00 with soc..."))
TEMP_AOI <- TEMP %>% crop(AOI) %>%
            resample(SOC_MAP_AOI, method = "bilinear") %>%
            mask(AOI)

#Potential Evapotranspiration layer 
message(noquote("Reading evapotranspiration 1981-2000 stack ..."))
PET <- stack("PET_Stack_81-00_CRU.tif")
message(noquote("synchronizing pet81-00 with soc..."))
PET_AOI <- PET %>% crop( AOI) %>%
           resample(SOC_MAP_AOI, method = "bilinear") %>%
           mask(AOI)

#Forward stack
message(noquote("Reading precipitation 2000-2018 stack ..."))
PREC_FORWARD <- stack("Prec_Stack_01-18_CRU.tif")
message(noquote("synchronizing prec01-0018 with soc..."))
PREC_FORWARD_AOI <- PREC_FORWARD %>% crop(AOI) %>%
                    resample(SOC_MAP_AOI, method = "bilinear") %>%
                    mask(AOI)

#Open Temperatures layer 
message(noquote("Reading temprature 1981-2000 raster ..."))
TEMP_FORWARD <- stack("Temp_Stack_01-18_CRU.tif")
message(noquote("synchronizing temp01-0018 with soc..."))
TEMP_FORWARD_AOI <- TEMP_FORWARD %>% crop(AOI) %>%
                    resample(SOC_MAP_AOI, method = "bilinear") %>%
                    mask(AOI)

#Potential Evapotranspiration layer 
message(noquote("Reading evapotranspiration 1981-2000 stack ..."))
PET_FORWARD <- stack("PET_Stack_01-18_CRU.tif")
message(noquote("synchronizing pet01-0018 with soc..."))
PET_FORWARD_AOI <- PET_FORWARD %>% crop(AOI) %>%
                  resample(SOC_MAP_AOI, method = "bilinear") %>%
                  mask(AOI)

#Land Use layer
message(noquote("Reading landcover stack ..."))
LU_AOI <- raster("Land_Cover_10class_AOI.tif")
message(noquote("synchronizing landcover with soc..."))
LU_AOI <- LU_AOI %>% crop(AOI) %>%
          resample(SOC_MAP_AOI, method = "ngb") %>%
          mask(AOI)

#Open Vegetation Cover layer 
message(noquote("Reading NDVI stack ..."))
Cov_AOI <- stack('NDVI_Cov_stack_AOI.tif')
message(noquote("synchronizing ndvi with soc..."))
Cov_AOI <- Cov_AOI %>% crop(AOI) %>% 
           resample(SOC_MAP_AOI, method = "bilinear") %>% 
           mask(AOI)

# DR stack for spinup 
message(noquote("Genarating spinup DR layer ..."))
DR_SPINUP <-
  (LU_AOI == 2 |
     LU_AOI == 12 |
     LU_AOI == 13) * 1.44 + (LU_AOI == 4) * 0.25 + (LU_AOI == 3 |
     LU_AOI == 5 | LU_AOI == 6 | LU_AOI == 8) * 0.67

message(noquote("Generating spinup stack ..."))
Stack_SPINUP_AOI <- stack(SOC_MAP_AOI,Clay_AOI,TEMP_AOI,PREC_AOI,PET_AOI,DR_SPINUP,LU_AOI,Cov_AOI)

# warmup
# Set the number of years of the warm up
nWUP<-18
LU_Stack <- stack(replicate(nWUP, LU_AOI))

message(noquote("Genarating warmaup DR layer ..."))
DR_WARMUP <-
  (LU_AOI == 2 |
     LU_AOI == 12 |
     LU_AOI == 13) * 1.44 + (LU_AOI == 4) * 0.25 + (LU_AOI == 3 |
     LU_AOI == 5 | LU_AOI == 6 | LU_AOI == 8) * 0.67
DR_WARMUP_Stack <- LU_Stack

for (i in 1:nlayers(LU_Stack)) {
  DR_WARMUP_Stack[[i]] <-
    (LU_Stack[[i]] == 2 |
       LU_Stack[[i]] == 12 |
       LU_Stack[[i]] == 13) * 1.44 + (LU_Stack[[i]] == 4) * 0.25 + 
       (LU_Stack[[i]] == 3 | LU_Stack[[i]] == 5 | LU_Stack[[i]] == 6 | 
       LU_Stack[[i]] == 8) * 0.67
}

message(noquote("Generating warmup stack ..."))
Stack_WARMUP_AOI <-
  stack(SOC_MAP_AOI, Clay_AOI, Cov_AOI, LU_Stack, DR_WARMUP_Stack)

#Forward
message(noquote("Genarating forward DR layer ..."))
DR_FORWARD <-
  (LU_AOI == 2) * 0.25 + (LU_AOI == 4 | LU_AOI == 7 | LU_AOI == 3) * 0.67

message(noquote("Generating forward stack ..."))
Stack_FORWARD_AOI <-
  stack(
    SOC_MAP_AOI,
    Clay_AOI,
    TEMP_FORWARD_AOI,
    PREC_FORWARD_AOI,
    PET_FORWARD_AOI,
    DR_FORWARD,
    LU_AOI,
    Cov_AOI
  )

setwd(WD_OUT)
message(noquote("Writing spinup, warmup and forward stacks ..."))
writeRaster(
  Stack_SPINUP_AOI,
  filename = ("Stack_Set_SPIN_UP_AOI.tif"),
  format = "GTiff",
  overwrite = TRUE
)
writeRaster(
  Stack_WARMUP_AOI,
  filename = ("Stack_Set_WARM_UP_AOI.tif"),
  format = "GTiff",
  overwrite = TRUE
)
writeRaster(
  Stack_FORWARD_AOI,
  filename = ("Stack_Set_FORWARD_AOI.tif"),
  format = "GTiff",
  overwrite = TRUE
)
}