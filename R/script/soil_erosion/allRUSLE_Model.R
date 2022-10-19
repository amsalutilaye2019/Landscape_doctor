#-------------------------------------------------------------------------------
# Main Function to calculate rusle
#-------------------------------------------------------------------------------

setwd("location of scripts")
source("load_packages.R")
source("RUSLE_Model_R.R")
source("RUSLE_Model_K.R")
source("RUSLE_Model_C.R")

#example
input_wd = "C:\\Users\\ATilaye\\Documents\\01My_Docs\\01CIAT\\01Landscape_Doctor\\RScrips\\test_data"
rf_st_date = "2020-01-01"
rf_en_date = "2020-03-01"
aoi = "woreda.shp"
clay  = "clay.tif"
sand = "sand.tif"
silt = "silt.tif"
orgc = "soc.tif"
nd_en_date = "2020-05-01" 
nd_st_date = "2020-03-15"

#call the functions in the scripts
r <- calculateR(input_wd, rf_st_date, rf_en_date, aoi)
k <- calculateK(input_wd, aoi, sand, silt, clay, orgc)
c <- calculateC(input_wd, aoi, nd_st_date, nd_en_date)

#read ls from file
ls <- rast("ls.tif") 

r <- terra::resample(r, c, method = 'bilinear') %>% crop(c) %>% mask(c)
k <- terra::resample(k, c, method = 'bilinear') %>% crop(c) %>% mask(c)
ls <- terra::resample(ls, c, method = 'bilinear') %>% crop(c) %>% mask(c)

#read p factor from database



