#### Simulating data from sablefish trawl data ####

### Install Packages ###
#install.packages("remotes")
library(remotes)
#install.packages("devtools")
library(devtools)
#install.packages("pkgbuild")
library(pkgbuild)
remotes::install_github("pbs-assess/sdmTMB", dependencies = TRUE,  ref="newlogistic")
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
#install.packages("ggpubr")
library(ggpubr)
library(scales)
library(visreg)
library(ggeffects)
library(stringr)
library(TMB)
#install.packages("tweedie")
library(tweedie)
#install.packages("ggridges")
library(ggridges)
library(viridis)

### Set ggplot themes ###
theme_set(theme_bw(base_size = 30))
theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### load helper functions ####
source("Code/util_funs.R")
source("Code/sim_funs.R")

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

### Fish density distribution simulation ####
##Set parameter values for data generation and model fitting ###
s50 <-2 # same as sablefish= 0.88; had been 2 in previous data simulation
delta <- 2 #same as sablefish = 0.57; had been 2 in previous data simulation
smax <- 30 # maximum effect of MI; had been 4 in previous data simulation
Eo <- 0.3
Eo2 <- 0.7

b_years <- rnorm(n = 6, mean = 4, sd = 1)
beta1 <- 1.5
beta2 <- -1
phi <- 10 # 16 corresponds to sablefish 
p <- 1.51
range <- 85
sigma_O <- 1.77

## Set how many data sets to produce ##
n <- 250

## Make list of parameter values ##
model.pars <- list(b_years = b_years,
                   beta1 = beta1,
                   beta2 = beta2,
                   phi = phi,
                   p = p,
                   range = range,
                   sigma_O=sigma_O)
## Make mesh ##
mesh <- make_mesh(dat, xy_cols = c("X", "Y"), n_knots=250)

### Simulate data under typical Eo ####
use_previous <- F
if(use_previous){
  simdat <- readRDS("~/Dropbox/GitHub/estimating_mi_from_distribution2/Model Outputs/data_sims_usual.rds")
 simdat2 <- readRDS("~/Dropbox/GitHub/estimating_mi_from_distribution2/Model Outputs/data_sims_weird.rds")
}

if(!use_previous){
simdat <- map(seq_len(n), ~simulate_fish(dat = dat,
                        mesh = mesh,
                        s50 = s50,
                        delta = delta,
                        smax = smax,
                        Eo = Eo,
                        modelpars = model.pars))
### Simulate data undate atypical Eo
simdat2 <- map(seq_len(n), ~simulate_fish(dat = dat,
                            mesh = mesh,
                            s50 = s50,
                            delta = delta,
                            smax = smax,
                            Eo = Eo2,
                            modelpars = model.pars))
}
###Sanity checks on simulated data

## Compare one to the real sablefish data ##
sanity <- F
if(sanity){
sim_test <- simdat[[1]]
sim_test2 <- simdat2[[1]]

ggplot(dat, aes(x=po2_s,y=cpue_kg_km2))+geom_point(size=0.6)+geom_point(sim_test, mapping=aes(x=po2_s, y=sim), color="blue", size=0.6)

#Percent below threshold
sim_test3 <- subset(sim_test, sim_test$mi_usual <s50)
}

#Save simulated data
save <- T
if(save){
saveRDS(simdat, "data_sims_usual.rds")
saveRDS(simdat2, "data_sims_weird.rds")
}


## Fit Model 1: No priors ##
##Set starting parameters
#Correct values
start <- matrix(0, ncol = 1, nrow = 4)
start[1,1] <- s50
start[2,1] <- delta
start[3,1] <- smax
start[4,1] <- Eo

start_unusual <- start
start_unusual[4,1] <- Eo2

#Just above zero
start2 <- matrix(0.08, ncol = 1, nrow = 4)

#Slightly farther
start3 <- matrix(0, ncol = 1, nrow = 4)
start3[1,1] <- s50*0.3
start3[2,1] <- delta*0.3
start3[3,1] <- smax*0.3
start3[4,1] <- Eo*0.3

#With weird Eo
start4 <- matrix(0, ncol = 1, nrow = 4)
start4[1,1] <- s50
start4[2,1] <- delta
start4[3,1] <- smax
start4[4,1] <- Eo2


## Fit model to all simulated datasets ##
if(use_previous){
  load("model_fits.Rdata")
}
if(!use_previous){
fits <- lapply(simdat, run_sdmTMB_noprior, 
                start=start, mesh=mesh)
fits2 <- lapply(simdat, run_sdmTMB_prior, 
                start=start, mesh=mesh)

fits3 <- lapply(simdat2, run_sdmTMB_noprior, 
                start=start_unusual, mesh=mesh)
fits4 <- lapply(simdat2, run_sdmTMB_prior, 
                start=start_unusual, mesh=mesh)
}
#Save models
save <- T
if(save){
save(fits, fits2, fits3, fits4, file="model_fits.Rdata")
}

