#### Simulating data from sablefish trawl data ####

### Install Packages ###
install.packages("remotes")
library(remotes)
install.packages("devtools")
library(devtools)
install.packages("pkgbuild")
library(pkgbuild)
remotes::install_github("pbs-assess/sdmTMB", dependencies = TRUE, force=TRUE, ref="newlogistic")
library(sdmTMB)
library(here)
library(mvtnorm)
library(mgcv)
library(dplyr)
library(devtools)
library(zoo)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(MASS)
install.packages("ggpubr")
library(ggpubr)
library(scales)
library(visreg)
library(ggeffects)
library(stringr)
library(TMB)

### Set ggplot themes ###
theme_set(theme_bw(base_size = 15))
theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### load helper functions ####
source("Code/util_funs.R")

### Load Data ####

sci_name <- "Anoplopoma fimbria" 
spc <- "sablefish" 
dat.by.size <- length_expand(sci_name)
dat <- load_data(spc = spc, dat.by.size = dat.by.size)

### Set up environmental data ###
kelvin = 273.15
boltz = 0.000086173324
tref <- 12
#Calculate inverse temp
dat$invtemp <- (1 / boltz)  * ( 1 / (dat$temp + 273.15) - 1 / (tref + 273.15))
#Calculate MI for usual case
dat$mi_usual = dat$po2*exp(0.291* dat$invtemp)
#Calculate MI for unusual case (99th percentile)
dat$mi_weird = dat$po2*exp(0.733* dat$invtemp)
#Scale and such
dat$temp_s <- (scale(dat$temp))
dat$po2_s <- (scale(dat$po2))
dat$log_depth_scaled <- scale(log(dat$depth))
dat$log_depth_scaled2 <- with(dat, log_depth_scaled ^ 2)
dat$jday_scaled <- scale(dat$julian_day)
dat$jday_scaled2 <- with(dat, jday_scaled ^ 2)
dat$X <- dat$longitude
dat$Y <- dat$latitude
dat$cpue_kg_km2 <- dat$cpue_kg_km2 * (dat$p2+dat$p3)
dat$year <- as.factor(dat$year)

### Fish density distribution simulation ###
## Simulation parameters ##

mesh <- make_mesh(dat, xy_cols = c("X", "Y"), n_knots=250)

simulate_fish<- function(dat,Eo) {
  x50 <- -1 # same as sablefish; had been 2 in previous data simulation
  delta <- 1 #same as sablefish; had been 2 in previous data simulation
  b_years <- c(4.47,4.53,4.44,4.43,4.68,4.68) #basically same as real sablefish
  beta1 <- 1.5 #depth #same as sablefish
  beta2 <- -1 #depth^2 #same as sablefish
  beta3 <- 50 # maximum effect of MI; had been 4 in previous data simulation
  phi <- 7 #estimated at 16 in real sablefish data
  p <- 1.5 #same as sablefish
  range <- 80 #80 in sablefish; had been 0.3 in previous data simulation
  sigma_O <- 1.5 #1.81 in real sablefish; had been 0.5 in previous data simulation
  seed <- sample(1:1000, 1)
  mesh <- make_mesh(dat, xy_cols = c("X", "Y"), n_knots=250)
  Eo <- Eo
  dat <- dat
  sim <- sdmTMB_simulate(formula=~-1+as.factor(year)+logistic(mi)+log_depth_scaled+log_depth_scaled2,
                      data=dat,
                      family=tweedie(link="log"),
                      tweedie_p=p,
                      phi=phi,
                      range=range,
                      sigma_O=sigma_O,
                      sigma_E=NULL,
                      mesh=mesh,
                      threshold_coefs=c(x50, delta, beta3, Eo),
                      B=c(b_years, beta1, beta2),
                      seed=seed)

  #Merge back with original dataframe
  data <- merge(sim, dat)
  data$seed <- seed
  data_complete <- list(data, mesh)
  return(data_complete)
}

##Run data simulations
n <-10 #Number of simulations
data_sims_usual <- map(seq_len(n), ~simulate_fish(dat, 0.291)) 
data_sims_weird <- map(seq_len(n), ~simulate_fish(dat, 0.733)) 

##Run again to create cross-validation data
data_sims_usual_CV <- map(seq_len(n), ~simulate_fish(dat, 0.291)) 
data_sims_weird_CV <- map(seq_len(n), ~simulate_fish(dat, 0.733)) 

###Sanity checks on simulated data

## Compare one to the real sablefish data ##
sim_test <- data_sims_usual[[1]][[1]]
sim_test2 <- data_sims_weird[[1]][[1]]

ggplot(dat, aes(x=po2,y=cpue_kg_km2))+geom_point()
ggplot(sim_test, aes(x=po2,y=observed))+geom_point()
ggplot(sim_test2, aes(x=po2,y=observed))+geom_point()

#Save simulated data
saveRDS(data_sims_usual, "data_sims_usual.rds")
saveRDS(data_sims_weird, "data_sims_weird.rds")
saveRDS(data_sims_usual_CV, "data_sims_usual_CV.rds")
saveRDS(data_sims_weird_CV, "data_sims_weird_CV.rds")