library(sp)
library(ggplot2)
library(broom)
library(dplyr)
library(tidyr)
library(purrr)
library(KernSmooth)
library(sdmTMB)
library(ggplot2)
library(visreg)
library(ggeffects)
library(mgcv)
if (!require("here")) install.packages("here")
library(here)
if (!require("gsw")) install.packages("gsw")
library(gsw)
  
  
# calc o2 solubility, relies on o2 in umol/kg
gsw_O2sol_SP_pt <- function(sal,pt) {
  x = sal
  pt68 = pt*1.00024
  y = log((298.15 - pt68)/(273.15 + pt68))
  
  a0 =  5.80871
  a1 =  3.20291
  a2 =  4.17887
  a3 =  5.10006
  a4 = -9.86643e-2
  a5 =  3.80369
  b0 = -7.01577e-3
  b1 = -7.70028e-3
  b2 = -1.13864e-2
  b3 = -9.51519e-3
  c0 = -2.75915e-7
  
  O2sol = exp(a0 + y*(a1 + y*(a2 + y*(a3 + y*(a4 + a5*y)))) + x*(b0 + y*(b1 + y*(b2 + b3*y)) + c0*x))
  return(O2sol)
}

load_all_hauls <- function() {
  install.packages("remotes")
  remotes::install_github("nwfsc-assess/nwfscSurvey")
  haul = nwfscSurvey::PullHaul.fn(SurveyName = "NWFSC.Combo")
  haul <- plyr::rename(haul, replace=c("salinity_at_gear_psu_der" = "sal", 
                                       "temperature_at_gear_c_der" = "temp", 
                                       "o2_at_gear_ml_per_l_der" = "o2",
                                       "depth_hi_prec_m" = "depth"))
  
  # read in the grid cell data from the survey design
  grid_cells = readxl::read_excel("data/Selection Set 2018 with Cell Corners.xlsx")
  grid_cells = dplyr::mutate(grid_cells,
                             depth_min = as.numeric(unlist(strsplit(grid_cells$Depth.Range,"-"))[1]),
                             depth_max = as.numeric(unlist(strsplit(grid_cells$Depth.Range,"-"))[2]))
  
  # convert grid_cells to sp object
  grid = SpatialPoints(cbind(grid_cells$Cent.Long,grid_cells$Cent.Lat),
                       proj4string = CRS("+proj=longlat +datum=WGS84"))
  r = raster::rasterize(x=grid, y = raster(nrow=length(unique(grid_cells$Cent.Lat)),
                                           ncol=length(unique(grid_cells$Cent.Long))))
  rasterToPoints(r)
  
  raster = aggregate(r, fact = 2)
  raster = projectRaster(raster, crs = "+proj=tmerc +lat_0=31.96 +lon_0=-121.6 +k=1 +x_0=390000 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
  
  # create matrix of point data with coordinates and depth from raster
  grid = as.data.frame(rasterToPoints(raster))
  
  # Figure out the grid cell corresponding to each tow location
  haul$Cent.Lat = NA
  haul$Cent.Lon = NA
  haul$Cent.ID = NA
  for(i in 1:nrow(haul)) {
    indx = which(grid_cells$NW.LAT > haul$latitude_dd[i] &
                   grid_cells$SW.LAT < haul$latitude_dd[i] &
                   grid_cells$NW.LON < haul$longitude_dd[i] &
                   grid_cells$NE.LON > haul$longitude_dd[i])
    if(length(indx) > 0) {
      haul$Cent.ID[i] = grid_cells$Cent.ID[indx]
      haul$Cent.Lat[i] = grid_cells$Cent.Lat[indx]
      haul$Cent.Lon[i] = grid_cells$Cent.Long[indx]
    }
  }
  
  # project lat/lon to UTM, after removing missing values and unsatisfactory hauls
  haul = haul %>% filter(!is.na(Cent.Lon), performance == "Satisfactory")
  
  haul_trans = haul
  coordinates(haul_trans) <- c("Cent.Lon", "Cent.Lat")
  proj4string(haul_trans) <- CRS("+proj=longlat +datum=WGS84")
  newproj = paste("+proj=utm +zone=10 ellps=WGS84")
  haul_trans <- spTransform(haul_trans, CRS(newproj))
  haul_trans = as.data.frame(haul_trans)
  haul_trans$Cent.Lon = haul_trans$Cent.Lon/10000
  haul_trans$Cent.Lat = haul_trans$Cent.Lat/10000
  haul_trans$year = as.numeric(substr(haul_trans$date_yyyymmdd,1,4))
  
  haul$X = haul_trans$Cent.Lon
  haul$Y = haul_trans$Cent.Lat
  haul$year = haul_trans$year
  #haul$year_centered = haul$year - mean(unique(haul$year))
  
  return(haul)
  
  
}

calc_po2_mi <- function(dat) {
  #O2 from trawl data is in ml/l - may need to be converted to umol/kg
  gas_const = 8.31
  partial_molar_vol = 0.000032
  kelvin = 273.15
  boltz = 0.000086173324
  
  #calculate percent saturation for O2 - assumes  units of mL O2/L
  # Input:       S = Salinity (pss-78)
  #              T = Temp (deg C) ! use potential temp
  #depth is in meters
  #[umole/kg] = [ml/L]*44660/(sigmatheta(P=0,theta,S) + 1000)
  dat$SA = gsw_SA_from_SP(dat$sal,dat$depth,dat$longitude_dd,dat$latitude_dd) #absolute salinity for pot T calc
  dat$pt = gsw_pt_from_t(dat$SA,dat$temp,dat$depth) #potential temp at a particular depth
  dat$CT = gsw_CT_from_t(dat$SA,dat$temp,dat$depth) #conservative temp
  dat$sigma0 = gsw_sigma0(dat$SA,dat$CT)
  dat$o2_umolkg = dat$o2*44660/(dat$sigma0+1000)
  
  
  dat$O2_Sat0 = gsw_O2sol_SP_pt(dat$sal,dat$pt)
  
  #= o2satv2a(sal,pt) #uses practical salinity and potential temp - solubity at p =1 atm
  dat$press = exp(dat$depth*10000*partial_molar_vol/gas_const/(dat$temp+kelvin))
  dat$O2_satdepth = dat$O2_Sat0*dat$press
  
  #solubility at p=0
  dat$sol0 = dat$O2_Sat0/0.209
  dat$sol_Dep = dat$sol0*dat$press
  dat$po2 = dat$o2_umolkg/dat$sol_Dep
  dat$po2 <- dat$po2 * 101.325 # convert to kPa
  
  # species-specific parameters, estimated by fitting linear model to all data in Deutsch (fish only)
  #Eo <- 0.4480175
  Eo <- 1.3
  Ao <- 4.4733481
  n <- -0.2201635
  avgbn <- 0.112  # this works for the adult class (0.5 - 6 kg).  for the large adult, adjust
  inv.temp <- (1 / boltz) * ( 1 / (dat$temp + kelvin) - 1 / (15 + kelvin))
  dat$mi = avgbn*Ao*dat$po2 *exp(Eo * inv.temp)
  return(dat)
}

load_data <- function(spc,dat.by.size) {
  dat <- readRDS("data/joined_nwfsc_data.rds")
  dat = dplyr::filter(dat, species == spc, year%in%seq(2010,2015))
  dat <- left_join(dat, dat.by.size, by = "trawl_id")
  # remove tows where there was positive catch but no length measurements
  dat <- dplyr::filter(dat, !is.na(p1))
  # analyze or years and hauls with adequate oxygen and temperature data, within range of occurrence
  
  # get julian day
  dat$julian_day <- rep(NA, nrow(dat))
  for (i in 1:nrow(dat)) dat$julian_day[i] <- as.POSIXlt(dat$date[i], format = "%Y-%b-%d")$yday
  
  
  #O2 from trawl data is in ml/l 
  # just in case, remove any missing or nonsense values from sensors
  dat <- dplyr::filter(dat, !is.na(o2), !is.na(sal), !is.na(temp), is.finite(sal))
  dat <- calc_po2_mi(dat)
  dat <- dplyr::filter(dat, !is.na(temp), !is.na(mi))
  
  # prepare data and models -------------------------------------------------
  
  dat <- dplyr::select(dat, trawl_id, species, year, longitude_dd, latitude_dd, cpue_kg_km2,
                       o2, temp, depth, mi, po2, julian_day, pass, p1, p2, p3, p4)
  
  
  # UTM transformation
  dat_ll = dat
  coordinates(dat_ll) <- c("longitude_dd", "latitude_dd")
  proj4string(dat_ll) <- CRS("+proj=longlat +datum=WGS84")
  # convert to utm with spTransform
  dat_utm = spTransform(dat_ll, 
                        CRS("+proj=utm +zone=10 +datum=WGS84 +units=km"))
  # convert back from sp object to data frame
  dat = as.data.frame(dat_utm)
  dat = dplyr::rename(dat, longitude = coords.x1,
                      latitude = coords.x2)
  return(dat)
}


# Species of interest and max. juvenile lengths (define ontogenetic classes)
length_expand <- function(sci_name) {
# load, clean, and join data
bio = readRDS("data/wcbts_bio_2019-08-01.rds")
haul = readRDS("data/wcbts_haul_2019-08-01.rds")
catch = readRDS("data/wcbts_catch_2019-08-01.rds")
names(catch) = tolower(names(catch))
names(bio) = tolower(names(bio))
names(haul) = tolower(names(haul))

bio$trawl_id = as.character(bio$trawl_id)
haul$trawl_id = as.character(haul$trawl_id)
haul$date_yyyymmdd = as.numeric(haul$date_yyyymmdd)
haul$sampling_end_hhmmss = as.numeric(haul$sampling_end_hhmmss)
haul$sampling_start_hhmmss = as.numeric(haul$sampling_start_hhmmss)

dat = dplyr::left_join(catch[,c("trawl_id","scientific_name","year","subsample_count",
                                "subsample_wt_kg","total_catch_numbers","total_catch_wt_kg","cpue_kg_km2")], haul, relationship = "many-to-many") %>%
  dplyr::left_join(filter(bio, !is.na(length_cm)), relationship = "many-to-many") %>%
  filter(performance == "Satisfactory")  %>%
  mutate(depth_m = depth_hi_prec_m)



# filter out species of interest from joined (catch/haul/bio) dataset
dat_sub = dplyr::filter(dat, scientific_name == sci_name)

# fit length-weight regression by year to predict fish weights that have lengths only.
# note a rank-deficiency warning may indicate there is insufficient data for some year/sex combinations (likely for unsexed group)

fitted = dat_sub %>%
  filter(!is.na(length_cm), !is.na(weight_kg)) %>%
  dplyr::select(trawl_id,year,
                subsample_wt_kg, total_catch_wt_kg, area_swept_ha_der, cpue_kg_km2,
                individual_tracking_id, sex, length_cm, weight_kg) %>%
  group_nest(year)  %>%
  mutate(
    model = map(data, ~ lm(log(weight_kg) ~ log(length_cm), data = .x)),
    tidied = map(model, tidy),
    augmented = map(model, augment),
    predictions = map2(data, model, modelr::add_predictions)
  )

# replace missing weights with predicted weights
dat_pos = fitted %>%
  unnest(predictions) %>%
  dplyr::select(-data, -model, -tidied, -augmented) %>%
  mutate(weight = ifelse(is.na(weight_kg), exp(pred), weight_kg))

trawlids <- unique(dat_pos$trawl_id)
p <- data.frame(trawl_id = trawlids,
                p1 = 0,
                p2 = 0,
                p3 = 0,
                p4 = 0)

sizethresholds <- quantile(dat_pos$weight, c(0.15, 0.5, 0.85, 1), na.rm = T)
for (i in 1:length(trawlids)) {
  haul_sample<- dplyr::filter(dat_pos, trawl_id == trawlids[i])
  if(nrow(haul_sample) > 0 | var(haul_sample$weight >0)) {
    # fit kernel density to weight frequency
    smoothed_w <- bkde(haul_sample$weight, range.x = c(min(dat_pos$weight), max(dat_pos$weight)), bandwidth = 2)
    # make sure smoother predicts positive or zero density
    smoothed_w$y[smoothed_w$y<0] <- 0
    # calculate proportion by biomass and by number
    p_w_byweight <- smoothed_w$y * smoothed_w$x / sum(smoothed_w$x*smoothed_w$y)

    
    p_w_byweight[p_w_byweight<0] <- 0
    #p_w_bynum[p_w_bynum<0] <- 0
    
    p1 <- sum(p_w_byweight[smoothed_w$x<=sizethresholds[1]])
    p2 <- sum(p_w_byweight[smoothed_w$x>sizethresholds[1] & smoothed_w$x <=sizethresholds[2]])
    p3 <- sum(p_w_byweight[smoothed_w$x>sizethresholds[2] & smoothed_w$x <=sizethresholds[3]])
    p4 <- sum(p_w_byweight[smoothed_w$x>sizethresholds[3]])
  
    
    
    p[i,2:5] <- c(p1, p2, p3, p4)
    
  }
  else {
    indx <- which(sizethresholds>haul_sample$weight)
    p[i, min(indx)+1] <- 1
  }
}


# add hauls with zero catch back in
absent = filter(dat_sub, cpue_kg_km2 == 0)
trawlids <- unique(absent$trawl_id)
absent.df <- data.frame(trawl_id = trawlids,
                        p1 = 0,
                        p2 = 0,
                        p3 = 0,
                        p4 = 0)

all_hauls <- rbind(p, absent.df)
all_hauls$trawl_id <- as.numeric(all_hauls$trawl_id)
return(all_hauls)
}

brkptfun <- function(x, b_slope, b_thresh) min(0, b_slope *  (x - b_thresh))
logfun <- function(x, model, mi = F) {
  if (mi) {
    parfit <- model$sd_report
    npars <- length(parfit$value)
    parnames <- names(parfit$value)
    
    s50 <- parfit$value[grep("s50", parnames)]
    delta <- parfit$value[grep("s95", parnames)]
    smax <- parfit$value[grep("s_max", parnames)]
  }
  if (!mi) {
    parfit <- model$sd_report
    npars <- length(parfit$value)
    parnames <- names(parfit$value)
    
    s50 <- parfit$value[grep("s50", parnames)]
    s95 <- (parfit$value[grep("s95", parnames)] )
    delta <- s95 - s50
    smax <- parfit$value[grep("s_max", parnames)]
  }
  a <- log(smax / (log(0.5) + smax) - 1)
  b <- log(smax / (log(0.95) + smax) - 1)
  beta0 <- -a + s50 * (b - a) / delta
  beta1 <- (a - b) / delta
  logmu <- smax * (1 / ( 1 + exp( - beta0 - beta1 * x)) -1)
  
  return(logmu)
} 
getEo <- function(model) {
  parfit <- model$sd_report
  npars <- length(parfit$value)
  parnames <- names(parfit$value)
  Eo <- parfit$value[grep("Eo", parnames)]
  return(Eo)
}

get_inits <- function() {
  # This function stores initial values and lower / upper bounds for breakpoint / logistic models.
  init_vals <- list(species <- c("sablefish", "petralesole"))
  modelnames <- c("m1", "m2", "m2a", "m3")
  init_vals$sablefish <- list(model = modelnames)
  init_vals$petralesole <- list(model = modelnames)
  
  ### Specify sablefish starting values and lower / upper limits #####
  start <- matrix(0, nrow = 2, ncol = 1)
  start[1,1] <- 20
  start[2,1] <- -1.1
  init_vals$sablefish$m1$start <- start
  init_vals$sablefish$m1$lower <- matrix(data=-Inf, ncol=1, nrow=2)
  init_vals$sablefish$m1$upper <- matrix(data=Inf, ncol=1, nrow=2)
  
  start <- matrix(0, ncol = 1, nrow = 4)
  start[1, 1] <- 2 #s50
  start[2, 1] <- 1 #delta
  start[3, 1] <- 20 #smax 
  start[4, 1] <- 1.50 #Eo  #Tim used 1.00 on 6/26
  init_vals$sablefish$m2$start <- start
  init_vals$sablefish$m2$lower <- c(-Inf, .001, 0, 0)
  init_vals$sablefish$m2$upper <- c(Inf, Inf,Inf, Inf)
  
  start <- matrix(0, ncol = 1, nrow = 4)
  start[1, 1] <- 2 #s50
  start[2, 1] <- (1) #delta
  start[3, 1] <- 10 #smax 
  start[4, 1] <- 0.68 #Eo
  init_vals$sablefish$m2a$start <- start
  init_vals$sablefish$m2a$lower <- c(-Inf, 0.01, 0, 0)
  init_vals$sablefish$m2a$upper <- c(Inf, Inf,50, Inf)
  init_vals$sablefish$m2a$prior <- normal(c(NA, NA, NA, 0.448), c(NA, NA, NA, 0.15))
  
  start <- matrix(0, ncol = 1, nrow = 3)
  start[1, 1] <- -1.65 #s50 #Tim used -1.5 on 6/26
  start[2, 1] <- log(0.5) # log delta
  start[3, 1] <- 40 #smax
  init_vals$sablefish$m3$start <- start
  init_vals$sablefish$m3$lower <- matrix(data=-Inf, ncol=1, nrow=3)
  init_vals$sablefish$m3$upper <- matrix(data=Inf, ncol=1, nrow=3)

  
  ### Specify petrale sole starting values and lower / upper limits #####
  start <- matrix(0, nrow = 2, ncol = 1)

  start[1,1] <- 0
  start[2,1] <- -0.2
  start[1,1] <- 200
  start[2,1] <- -1

  
  init_vals$petralesole$m1$start <- start
  init_vals$petralesole$m1$lower <- c(-Inf -Inf, -Inf, 0.01)
  init_vals$petralesole$m1$upper <- NULL
  
  start <- matrix(0, nrow = 4, ncol = 1)
  start[1, 1] <- -1 #s50
  start[2, 1] <- 2 #delta
  start[3, 1] <- 20 #smax 
  start[4, 1] <- 1.00 #Eo
  lower <- c(-5, .001, 0.01, 0.01)
  upper <- c(10, 10, 100, 3)
  
  
  init_vals$petralesole$m2$start <- start
  init_vals$petralesole$m2$lower <- lower
  init_vals$petralesole$m2$upper <- upper

  start <- matrix(0, ncol = 1, nrow = 4)
  start[1, 1] <- 0.67 #s50
  start[2, 1] <- 0.44 #delta
  start[3, 1] <- 150 #smax 
  start[4, 1] <- 0.01 #Eo
  lower <- c(-2, .01, 0.01, 0.01)
  upper <- c(10, 10, 200, 2)
  
  init_vals$petralesole$m2a$start <- start
  init_vals$petralesole$m2a$lower <-  lower
  init_vals$petralesole$m2a$upper <- upper
  init_vals$petralesole$m2a$prior <- normal(c(NA, NA, NA, 0.448), c(NA, NA, NA, 0.15))

  start <- matrix(0, ncol = 1, nrow = 3)
  start[1, 1] <- -1.5 #s50
  start[2, 1] <- log(1.1) # log delta
  start[3, 1] <- 10 #smax
  lower <- c(-4, -Inf, 0.01)
  upper <- c(4, 3,300)
  
  init_vals$petralesole$m3$start <- start
  init_vals$petralesole$m3$lower <-   lower
  init_vals$petralesole$m3$upper <- upper

  return(init_vals)
}
  


