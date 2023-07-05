### Fit model in sdmTMB ###
install.packages("remotes")
library(remotes)
install.packages("devtools")
library(devtools)
install.packages("pkgbuild")
library(pkgbuild)
remotes::install_github("pbs-assess/sdmTMB", dependencies = TRUE, ref="newlogistic")
library(sdmTMB)

## Load simulated data if needed ##
readRDS("~/Dropbox/Mac/Documents/GitHub/thresholds_mi_distribution/thresholds_mi_distribution2/data_sims_usual.rds")
readRDS("~/Dropbox/Mac/Documents/GitHub/thresholds_mi_distribution/thresholds_mi_distribution2/data_sims_weird.rds")

## Specify all parameters from data generation for reference ##
## Simulation parameters ##
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
start3[1,1] <- x50*0.1
start3[2,1] <- delta*0.1
start3[3,1] <- beta3*0.1
start3[4,1] <- Eo*0.1

##Function to run model and return list of model outputs
run_sdmTMB_1 <- function(x, start) {
    m2 <- try(sdmTMB(observed ~ -1+as.factor(year)+logistic(mi)+log_depth_scaled+log_depth_scaled2, 
                     data = x[[1]], 
                     spatial = "on",
                     spatiotemporal="off",
                     mesh=x[[2]],
                     family =tweedie(link="log"),
                     control = sdmTMBcontrol(
                      start = list(b_threshold = start),
                     #lower = list(b_threshold = c(0, 0, 0, 0)), 
                     #upper = list(b_threshold = c(Inf, Inf, Inf, Inf)),
                     newton_loops = 2)))
    try(tidy(m2))
    try(return(m2))
}

## Model 2: Prior (correct) ##
run_sdmTMB_2 <- function(x, start) {
  m2 <- try(sdmTMB(observed ~ 1+as.factor(year)+logistic(mi)+log_depth_scaled+log_depth_scaled2, 
                   data = x[[1]], 
                   spatial = "on",
                   spatiotemporal="off",
                   mesh=x[[2]],
                   family =tweedie(link="log"),
                   control = sdmTMBcontrol(
                     start = list(b_threshold = start),
                     newton_loops = 1),
                   priors=sdmTMBpriors(threshold = normal(c(NA, NA, NA, 0.3306), c(NA, NA, NA, 0.173)))))
  try(tidy(m2))
  try(return(m2))
}


## Fit model to all simulated datasets ##
fits <- lapply(data_sims_usual, run_sdmTMB_1, 
               start=start)
fits2 <- lapply(data_sims_usual, run_sdmTMB_2, 
               start=start)
fits3 <- lapply(data_sims_weird, run_sdmTMB_1, 
                start=start)
fits4 <- lapply(data_sims_weird, run_sdmTMB_2, 
                start=start)

#Save models
#save(fits, fits1, fits2, fits3, fits_b, fits_c, file="model_fits.Rdata")
save(fits, fits2, fits3, file="model_fits.Rdata")
