#-------------------------------------------------------------------------------
# Function to calculate soil erosivity r using r_moore
#-------------------------------------------------------------------------------

calculateR <- function(input_wd, st_date, en_date, aoi){
    setwd(input_wd)
    woreda <- terra::vect(aoi)
    dates <- c(st_date, en_date)
    datevec <-
      seq(as.Date(dates[1]),
          by = "day",
          length.out = difftime(en_date, st_date, 'days') + 1)
    rf <-
      get_chirps(
        woreda,
        dates =  c(datevec[1], datevec[length(datevec)]),
        server = "CHC",
        as.raster = TRUE
      )
    rf <- rts(rf, datevec)
    
    # Calculate R factor
    # first calculates the mean monthly from the daily data then
    # calculates the mean for each year and the sums the whole years
    #rf2 <- rf %>% apply.monthly(sum) %>% apply.yearly(mean)
    rf2 <- rf %>% rts::apply.monthly(sum) %>% rts::apply.yearly(sum)
    precp <- sum(rf2@raster) %>% crop(woreda) %>% mask(woreda)
    
    ke <- 11.46 * precp - 2226
    r <- 0.029 * ke - 26
    r_si <- 17.02 * r 
    return(r_si)
}
