suppressPackageStartupMessages({
  library(terra)
  library(tidyterra)
  library(pbapply)
  library(tidyverse)
  })

# 0. Configure functions ----
GetSubDailyPerimeters <- function(FIRE_ID,
                                  PATH_FEATURE_LAYERS = "data/feature_layers"
                                  ){
  # Load fire perimeters
  fire_perimeters <- vect(
    sprintf("%s/feature_layers/fire_atlas/viirs_perimeters_in_cavm_e113.gpkg",DATA_DIR)
  )
  
  # Subset to single perimeter
  selected_fire_perimeter <- fire_perimeters %>% 
    filter(fireid  == FIRE_ID)
  
  # Extract time window of burn
  start_date <- make_date(year = selected_fire_perimeter$tst_year,
                          month = selected_fire_perimeter$tst_month,
                          day = selected_fire_perimeter$tst_day)
  
  end_date <- make_date(year = selected_fire_perimeter$ted_year,
                        month = selected_fire_perimeter$ted_month,
                        day = selected_fire_perimeter$ted_day)
  
  # List sub-daily perimeters for this fire scar
  FLIST_SD <- list.files(path = file.path(PATH_FEATURE_LAYERS,
                                          "fire_atlas/sub_daily",
                                          selected_fire_perimeter$tst_year,
                                          "Snapshot"
                                          ),
                         pattern = "*M.gpkg", full.names = TRUE)
  
  date_list <- lapply(FLIST_SD, function(f) {
    stem <- sub("\\.gpkg$", "", basename(f))
    
    # extract date (YYYYMMDD) and AM/PM
    datestr <- substr(stem, 1, 8)
    ampm    <- substr(stem, 9, 10)
    
    # parse the date
    date <- strptime(datestr, format = "%Y%m%d", tz = "UTC")
    
    # build final datetime using strftime + paste
    datetime_str <- paste(
      strftime(date, "%Y-%m-%d"),
      ifelse(ampm == "AM", "06:00:00", "18:00:00")
    )
    
    as.POSIXct(datetime_str, tz = "UTC")
  })
  
  date_vals <- as.Date(do.call(c, date_list))
  
  # Logical vector of which are in range
  in_range <- date_vals >= start_date & date_vals <= end_date
  
  # Extract matching dates
  FLIST_SD_IN_RANGE <- FLIST_SD[in_range]
  
  vect_perimeters <- do.call(rbind, lapply(FLIST_SD_IN_RANGE, vect))
  
  return(vect_perimeters)
}

ZonalExtraction <- function(optimality_lut,
                          FIRE_ID,
                          PATH_RASTER_LAYERS = "data/raster",
                          burn_severity_index = "dNBR_corr"){
  cat(FIRE_ID)
  
  # Load optimal burn severity raster
  FNAME_SEVERITY_RASTER <- optimality_lut %>% 
    filter(fireid == FIRE_ID,
           severity_index == "dNBR") %>% 
    pull(fname_severity_raster) %>% 
    basename() %>% 
    str_replace(.,"dNBR",burn_severity_index)
  
  rast_severity <- rast(file.path(PATH_RASTER_LAYERS,
                                  "hls/severity_rasters",
                                  FNAME_SEVERITY_RASTER))
  
  
  # Load sub-daily perimeters
  vect_subdaily <- GetSubDailyPerimeters(FIRE_ID = FIRE_ID)
  
  if (is.null(vect_subdaily)) {
    cat("No subdaily perimeters available here.")
    return()
  } 
  
  v <- vect_subdaily %>% 
    project(crs(rast_severity)) %>% 
    crop(ext(rast_severity)) %>%
    filter(fireid == FIRE_ID)
  
  n <- length(v)
  
  # store incremental polygons
  increments <- vector("list", n)
  
  # last layer is the final full perimeter
  increments[[n]] <- v[n]
  
  # loop backward from n-1 → 1
  for (i in (n-1):1) {
    
    current     <- v[i]     # smaller, earlier perimeter
    next_perim  <- v[i + 1] # larger, later perimeter
    
    # growth between step i and i+1
    inc <- erase(next_perim, current)
    
    increments[[i]] <- inc
  }
  
  # combine to SpatVector (skip missing layers)
  increments <- do.call(rbind, increments)
  
  vect_ee <- terra::extract(rast_severity,increments,
                            fun = "mean",bind = TRUE,
                            na.rm = TRUE)
  df_extract <- as.data.frame(vect_ee) 
  return(df_extract)
}

# 1. Apply function ----
OS <- Sys.info()[['sysname']]
DATA_DIR <- ifelse(OS == "Linux","~/data/","data/")

TABLE_DIR <- paste0(DATA_DIR,"tables/")
optimality_lut <- read_csv2(paste0(TABLE_DIR,"optimality_LUT.csv"),
                            show_col_types = FALSE)


results_dnbr_corr <- pblapply(optimality_lut$fireid,
                              ZonalExtraction, 
                              optimality_lut = optimality_lut,
                              burn_severity_index = "dNBR_corr")

df_dnbr_corr <- do.call(rbind,results_dnbr_corr)

results_dgemi <- pblapply(optimality_lut$fireid,
                              ZonalExtraction, 
                              optimality_lut = optimality_lut,
                              burn_severity_index = "dGEMI")


# 3. Analyse ----
df_ee <- df_dnbr_corr %>% 
  filter(n_newpixels != 0) # select perimeters with new burn pixels ()

m1 <-lm(log(meanFRP) ~ dNBR_corr, data = df_ee)
m1 <-lme4::lmer(meanFRP ~ dNBR_corr + (1|fireid), data = df_ee)
summary(m1)
m2 <- lm(log(FRP95) ~ dNBR_corr, data = df_ee)
summary(m2)

ggplot(df_ee,aes(x = dNBR_corr,y = log(meanFRP))) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ fireid)+
  theme_cowplot()

ggplot(df_ee,aes(x = dNBR_corr,y = log(FRP95))) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_cowplot()
