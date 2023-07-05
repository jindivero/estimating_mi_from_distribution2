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
<<<<<<< Updated upstream
  phi <- 7 #estimated at 16 in real sablefish data
  p <- 1.5 #same as sablefish
  range <- 80 #80 in sablefish; had been 0.3 in previous data simulation
  sigma_O <- 1.5 #1.81 in real sablefish; had been 0.5 in previous data simulation
  seed <- sample(1:1000, 1)
  mesh <- make_mesh(dat, xy_cols = c("X", "Y"), n_knots=50)
=======
  phi <- 4 #estimated at 16 in real sablefish data
  p <- 1.5 #same as sablefish
  range <- 20 #80 in sablefish; had been 0.3 in previous data simulation
  sigma_O <- 0.8 #1.81 in real sablefish; had been 0.5 in previous data simulation
  seed <- sample(1:1000, 1)
  mesh <- make_mesh(dat, xy_cols = c("X", "Y"), n_knots=250)
>>>>>>> Stashed changes
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

##Sanity checks on simulated data
## Compare one to the real sablefish data ##
sim_test <- data_sims_usual[[1]][[1]]
sim_test2 <- data_sims_weird[[1]][[1]]

ggplot(dat, aes(x=po2,y=cpue_kg_km2))+geom_point()
ggplot(sim_test, aes(x=po2,y=observed))+geom_point()
ggplot(sim_test2, aes(x=po2,y=observed))+geom_point()

## Specify all parameters from data generation for reference ##
## Simulation parameters ##
x50 <- -1 # same as sablefish; had been 2 in previous data simulation
delta <- 1 #same as sablefish; had been 2 in previous data simulation
b_years <- c(4.47,4.53,4.44,4.43,4.68,4.68) #basically same as real sablefish
beta1 <- 1.5 #depth #same as sablefish
beta2 <- -1 #depth^2 #same as sablefish
beta3 <- 50 # maximum effect of MI; had been 4 in previous data simulation
<<<<<<< Updated upstream
phi <- 7 #estimated at 16 in real sablefish data
p <- 1.5 #same as sablefish
range <- 80 #80 in sablefish; had been 0.3 in previous data simulation
sigma_O <- 1.5 #1.81 in real sablefish; had been 0.5 in previous data simulation
=======
phi <- 5 #estimated at 16 in real sablefish data
p <- 1.5 #same as sablefish
range <- 20 #80 in sablefish; had been 0.3 in previous data simulation
sigma_O <- 0.8 #1.81 in real sablefish; had been 0.5 in previous data simulation
>>>>>>> Stashed changes
Eo <- 0.291

## Fit Model 1: No priors ##
##Set starting parameters
#Correct values
start <- matrix(0, ncol = 1, nrow = 4)
start[1,1] <- x50
start[2,1] <- delta
start[3,1] <- beta3
start[4,1] <- Eo

#Just above zero
start2 <- matrix(0.08, ncol = 1, nrow = 4)
#Slightly farther
start3 <- matrix(0, ncol = 1, nrow = 4)
<<<<<<< Updated upstream
start3[1,1] <- x50*0.1
start3[2,1] <- delta*0.1
start3[3,1] <- beta3*0.1
start3[4,1] <- Eo*0.1
=======
start3[1,1] <- x50*0.5
start3[2,1] <- delta*0.5
start3[3,1] <- beta3*0.5
start3[4,1] <- Eo*0.5
>>>>>>> Stashed changes

##Function to run model and return list of model outputs
run_sdmTMB <- function(x, start) {
  m2 <- try(sdmTMB(observed ~ -1+as.factor(year)+logistic(mi)+log_depth_scaled+log_depth_scaled2, 
                   data = x[[1]], 
                   spatial = "on",
                   spatiotemporal="off",
                   mesh=x[[2]],
                   family =tweedie(link="log"),
                   control = sdmTMBcontrol(
                     start = list(b_threshold = start),
<<<<<<< Updated upstream
                     #lower = list(b_threshold = c(0, 0, 0, 0)), 
                     #upper = list(b_threshold = c(Inf, Inf, Inf, Inf)),
                     #priors=sdmTMBpriors(threshold = normal(c(NA, NA, NA, 0.3306), c(NA, NA, NA, 0.173))),
=======
                     lower = list(b_threshold = c(-Inf, -Inf, 0.01, 0.01)), 
                     upper = list(b_threshold = c(Inf, Inf, 100, 1)),
                     priors=sdmTMBpriors(threshold = normal(c(NA, NA, NA, 0.3306), c(NA, NA, NA, 0.173))),
>>>>>>> Stashed changes
                     newton_loops = 2)))
  try(tidy(m2))
  try(return(m2))
}

## Fit model to all simulated datasets ##
fits <- lapply(data_sims_usual, run_sdmTMB, 
<<<<<<< Updated upstream
               start=start)
=======
               start=start3)
>>>>>>> Stashed changes

extract_pars <- function(x){
  if(!is.character(x)){
    par_estimates <- as.data.frame(tidy(x, conf.int = TRUE, effects="fixed"))
    par_estimates_rand <- as.data.frame(tidy(x, conf.int = TRUE, effects="ran_pars"))
    par_estimates <- bind_rows(par_estimates, par_estimates_rand)
    return(par_estimates)
  }
  if(is.character(x)){
    return(NA)
  }
}
##Apply to each sim in each model
#Model starting threshold values at correct parameter values
pars1 <- lapply(fits, extract_pars)
##Function to clean up pars for plotting
clean_pars <- function(pars, fits){
  names(pars) <- c(1:length(fits))
  #Remove models with errors
  pars <- keep(pars, function(x) !is.logical(x))
  #Combine into single dataframe, with column of simulation number
  pars <- bind_rows(pars,.id="id")
  return(pars)
}

##Set true pars vectors of values and name
true_pars2 <- data.frame(term=c("log_depth_scaled", "log_depth_scaled2", "mi-delta", "mi-s50", "mi-smax", "range", "sigma_O", "phi", "tweedie_p", "mi-Eo", "as.factor(year)2010","as.factor(year)2011","as.factor(year)2012", "as.factor(year)2013", "as.factor(year)2014", "as.factor(year)2015"), estimate=c(beta1, beta2, delta, x50, beta3, range, sigma_O, phi, p, Eo, b_years))
pars_names <- c("(Intercept)","as.factor(year)2011","as.factor(year)2012", "as.factor(year)2013", "as.factor(year)2014", "as.factor(year)2015", "log_depth_sc", "log_depth_sc2", "mi-s50", "mi-delta", "mi-smax", "mi-Eo", "range", "phi", "sigma_O", "tweedie_p") 

##Extract pars for the other two models
##Apply to each sim in each model
pars1 <- lapply(fits, extract_pars)
pars1 <- clean_pars(pars1, fits=fits)
#Plot comparison
ggplot(pars1, aes(y=estimate, x=term))+geom_boxplot()+facet_wrap("term", scales="free")+geom_hline(data = true_pars2, aes(yintercept = estimate),linetype="dashed", size=1.2)+theme(legend.position="left", strip.text = element_blank())

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
