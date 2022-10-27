#-------------------------------------------------------------------------------
# Function to calculate soil cover factor c from hsitorical ndvi
#-------------------------------------------------------------------------------

calculateC <- function(input_wd, aoi, nd_st_date, nd_en_date){
    #Calculate C factor from MODIS data
    setwd(input_wd)
    aoi <- st_read(dsn = getwd(), layer = aoi) %>% st_as_sf()
    
    #find centroid of polygon in long-lat decimal degrees
    aoi.cent <- sf::st_centroid(aoi)
    aoi.proj <- sf::st_transform(aoi, "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=km +no_defs")
    
    #find width and height of country in km
    aoi.km_lr <- diff(st_bbox(aoi.proj)[c(1,3)])
    aoi.km_ab <- diff(st_bbox(aoi.proj)[c(2,4)])
    
    #download the MODIS tiles for the area we defined
    aoi_ndvi <- MODISTools::mt_subset(product = "MOD13Q1",
                          lat = aoi.cent$y, lon =  aoi.cent$x,
                          band = "250m_16_days_NDVI",
                          start = nd_st_date, end = nd_en_date,
                          km_lr = aoi.km_lr, km_ab = aoi.km_ab,
                          site_name = "woreda",
                          internal = TRUE, progress = TRUE)
    
    # convert to a spatial raster
    aoi.rast <-  aoi_ndvi %>% MODISTools::mt_to_raster(reproject = TRUE) 
    aoi.rast  <- aoi.rast %>% terra::rast() %>% terra::crop(terra::vect(aoi)) %>% terra::mask(terra::vect(aoi))
    class(aoi.rast)
    med_ndvi <- terra::median(aoi.rast, na.rm = T) #median ndvi
    alpha <- 2 # as suggested by Knijff 2000
    beta <- 1 # as suggested by Knijff 2000
    c <- min(exp(-alpha * (med_ndvi/(beta - med_ndvi))), 1)
    return(c)
    terra::writeRaster(c, filename = "C.tif", filetype = "GTiff")
}

input_wd = "C:\\Users\\ATilaye\\Documents\\01My_Docs\\01CIAT\\01Landscape_Doctor\\RScrips\\test_data"
aoi = "woreda"
nd_en_date = as.Date("2020-05-01")
nd_st_date = as.Date("2020-01-15")

#call the functions in the scripts
c <- calculateC(input_wd, aoi, nd_st_date, nd_en_date)


