#### Simulating data from sablefish trawl data ####

# Specify Parameters ####
s50 <- 2 # same as sablefish= 0.88; had been 2 in previous data simulation
delta <- 2 #same as sablefish = 0.57; had been 2 in previous data simulation
smax <- 30 # maximum effect of MI; had been 4 in previous data simulation
Eo <- 0.3
Eo2 <- 0.4

phi <- 5 #estimated at 16 in real sablefish data
p <- 1.5 #same as sablefish
range <- 90 #80 in sablefish; had been 0.3 in previous data simulation
sigma_O <- 0.8 #1.81 in real sablefish; had been 0.5 in previous data simulation

### Install Packages ###
library(remotes)
library(devtools)
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

simulate_fish<- function(dat,mesh, x50, delta, smax, Eo) {
 # x50 <- 2 # same as sablefish= 0.88; had been 2 in previous data simulation
#  delta <- 2 #same as sablefish = 0.57; had been 2 in previous data simulation
  b_years <- c(4.47,4.53,4.44,4.43,4.68,4.68) #basically same as real sablefish
  beta1 <- 1.5 #depth #same as sablefish
  beta2 <- -1 #depth^2 #same as sablefish
#  smax <- 50 # maximum effect of MI; had been 4 in previous data simulation
  phi <- 16 #estimated at 16 in real sablefish data
  p <- 1.51 #same as sablefish
  range <- 85 #80 in sablefish; had been 0.3 in previous data simulation
  sigma_O <- 1.77 #1.81 in real sablefish; had been 0.5 in previous data simulation
  seed <- sample(1:1000, 1)
  sim <- sdmTMB_simulate(formula=~-1+as.factor(year)+logistic(mi)+log_depth_scaled+log_depth_scaled2,
                         data=dat,
                         family=tweedie(link="log"),
                         tweedie_p=p,
                         phi=phi,
                         range=range,
                         sigma_O=sigma_O,
                         sigma_E=NULL,
                         mesh=mesh,
                         threshold_coefs=c(x50, delta, smax, Eo),
                         B=c(b_years, beta1, beta2),
                         seed=seed)
  
 
  dat$sim <- sim$observed
  return(dat)
}



simdat <- simulate_fish(dat = dat,
                        mesh = mesh,
                        x50 = x50,
                        delta = delta,
                        smax = smax,
                        Eo = Eo)

##Run data simulations
n <-25 #Number of simulations
data_sims_usual <- map(seq_len(n), ~simulate_fish(dat = dat,
                                                            mesh = mesh,
                                                            x50 = s50,
                                                            delta = delta,
                                                            smax = smax,
                                                            Eo = Eo))

data_sims_unusual <- map(seq_len(n), ~simulate_fish(dat = dat,
                                                  mesh = mesh,
                                                  x50 = s50,
                                                  delta = delta,
                                                  smax = smax,
                                                  Eo = Eo2))

##Sanity checks on simulated data
## Compare one to the real sablefish data ##
#sim_test <- data_sims_usual[[1]][[1]]
#sim_test2 <- data_sims_weird[[1]][[1]]

# ggplot(dat, aes(x=po2,y=cpue_kg_km2))+geom_point()
# ggplot(sim_test, aes(x=po2,y=observed))+geom_point()
# ggplot(sim_test2, aes(x=po2,y=observed))+geom_point()
# 
## Specify all parameters from data generation for reference ##
## Simulation parameters ##
#Eo <- 0.291

## Fit Model 1: No priors ##
##Set starting parameters
#Correct values



start1 <- matrix(0, ncol = 1, nrow = 4)
start1[1,1] <- x50
start1[2,1] <- delta
start1[3,1] <- smax
start1[4,1] <- Eo
start2 <- start1
start2[4,1] <- Eo2

##Function to run model and return list of model outputs
run_sdmTMB_1 <- function(simdat, start, mesh) {
  m2 <- try(sdmTMB(sim ~ -1+as.factor(year)+logistic(mi)+log_depth_scaled+log_depth_scaled2, 
                   data = simdat, 
                   spatial = "on",
                   spatiotemporal="off",
                   mesh=mesh,
                   family =tweedie(link="log"),
                   control = sdmTMBcontrol(
                     start = list(b_threshold = start),
                     #lower = list(b_threshold = c(-Inf, -Inf, -Inf, -Inf)), 
                     # upper = list(b_threshold = c(Inf, Inf, 100, Inf)),
                     newton_loops = 2)))
  try(tidy(m2))
  try(return(m2))
}

## Function to run model 2 (with prior) ##
run_sdmTMB_2 <- function(simdat, start, mesh) {
  m2 <- try(sdmTMB(sim ~ -1+as.factor(year)+logistic(mi)+log_depth_scaled+log_depth_scaled2, 
                   data = simdat, 
                   spatial = "on",
                   spatiotemporal="off",
                   mesh=mesh,
                   family =tweedie(link="log"),
                   control = sdmTMBcontrol(
                     start = list(b_threshold = start),
                     newton_loops = 2),
                   priors=sdmTMBpriors(threshold = normal(c(NA, NA, NA, 0.3477), c(NA, NA, NA, 0.1455)))))
  
  try(tidy(m2))
  try(return(m2))
}

## Fit model to all simulated datasets ##
fits <- lapply(data_sims_usual, run_sdmTMB_1, 
               start=start1, mesh=mesh)
fits2 <- lapply(data_sims_usual, run_sdmTMB_2, 
                start=start1, mesh=mesh)

fits3 <- lapply(data_sims_unusual, run_sdmTMB_1, 
                start=start2, mesh=mesh)
fits4 <- lapply(data_sims_unusual, run_sdmTMB_2, 
                start=start2, mesh=mesh)


##Aggregate errors
extract_convergence <- function(x){
  if(!is.character(x)){
    convergence <- as.data.frame(x[["model"]][["convergence"]])
    colnames(convergence) <- "convergence"
    convergence$iterations <- x$model$iterations
    convergence$message <- x$model$message
    convergence$evaluations_function <- x$model$evaluations[1]
    convergence$evaluations_gradient <- x$model$evaluations[2]
  }
  if(is.character(x)){
    convergence <- as.data.frame(matrix(ncol=5, nrow=1, NA))
    colnames(convergence) <- c("convergence", "iterations", "message", "evaluations_function", "evaluations_gradient")
    convergence[,3] <- paste(x[[1]])
  }
  cons <- bind_rows(convergence)
  return(cons)
}

convergence <- lapply(fits, extract_convergence)
convergence <- bind_rows(convergence1)

##Positive definite Hessian
extract_pdHess <- function(x){
  if(!is.character(x)){
    pdh <- as.data.frame(x$pos_def_hessian)
    return(pdh)
  }
}

pdHess1 <- lapply(fits, extract_pdHess)
pdHess <-bind_rows(pdHess1)

extract_pdHess2 <- function(pdHess){
  pdHess <- bind_rows(pdHess)
  pdHess$sim <- as.numeric(row.names(pdHess))
  pdHess$sim <- as.character(pdHess$sim)
  pdHess$neg_def_hessian <- ifelse(pdHess$"x$pos_def_hessian"=="TRUE", 1,0)
  number_negdHess <- print(sum(pdHess$neg_def_hessian))
  return(number_negdHess)
}
hess<- extract_pdHess2(pdHess1)
