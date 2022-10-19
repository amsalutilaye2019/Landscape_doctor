#-------------------------------------------------------------------------------
# Function to calculate soil erodability K using K_williams
#-------------------------------------------------------------------------------

calculateK <- function(input_wd, aoi, sand, silt, clay, orgc){
    # silt and clay are percent but orgc should be divided by 10 to get percent
    setwd(input_wd)
    woreda <- vect(aoi)
    new_crs <- "epsg:4326"
    woreda <- terra::project(woreda, new_crs)
    
    sndprc <- terra::rast(sand) 
    sndprc <- sndprc %>% terra::project(new_crs) 
    sndprc <- sndprc %>% terra::crop(woreda) %>% terra::mask(woreda)
    sltprc <- terra::rast(silt) 
    sltprc <- sltprc %>% terra::project(new_crs)
    sltprc <- sltprc %>% terra::crop(woreda) %>% terra::mask(woreda)
    clyprc <- terra::rast(clay) 
    clyprc <- clyprc %>% terra::project(new_crs)
    clyprc <- clyprc %>% terra::crop(woreda) %>% terra::mask(woreda)
    orcprc <- (terra::rast(orgc) * 0.1) 
    orcprc <- orcprc %>% terra::project(new_crs) 
    orcprc <- orcprc %>% terra::crop(woreda) %>% terra::mask(woreda)
    
    a <- (0.2 + 0.3 * exp(-0.0256 * sndprc * (1 - sltprc / 100)))
    b <- (sltprc / (clyprc + sltprc)) ^ 0.3
    c <- 1 - (0.25 * orcprc) / (orcprc + exp(3.72 - 2.95 * orcprc))
    sn1 <- 1 - sndprc / 100
    d <- 1 - (0.7 * sn1) / (sn1 + exp(-5.51 + 22.9 * sn1))
    k <- 0.1317 * a * b * c * d 
    return(k)
}