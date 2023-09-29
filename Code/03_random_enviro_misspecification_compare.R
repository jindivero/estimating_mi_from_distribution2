#### Comparing real data to randomly drawn environmental data ####
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

#
### load helper functions ####
source("Code/util_funs.R")
source("Code/sim_funs.R")

### Load Data ####
sci_name <- "Anoplopoma fimbria" 
spc <- "sablefish" 
dat.by.size <- length_expand(sci_name)
dat <- load_data(spc = spc, dat.by.size = dat.by.size)

###Environmental data simulation
simulate_environmental <- function(ndata, dat, po2_lim, temp_lim, depth_lim){
  ##Set seed
  seed <- sample(1:1000, 1)
  set.seed(seed)
  ##Simulate oxygen and temperature data
  #Get variance covariances of depth, temperature and po2
  env.dat <- dplyr::select(dat, "po2", "depth", "temp")
  var.covar <- var(log(env.dat))
  env.bar <- as.matrix(colMeans(log(env.dat)))
  #Simulate from MVN
  log.env <- matrix(nrow=ndata, ncol=3, NA)
  colnames(log.env) <- c("po2", "depth", "temp")
  #Random generation of environmental data, with constraints
  for (i in 1:ndata){
    log.env[i,] <- rmvnorm(n = 1, mean=(env.bar - diag(var.covar) / 2), sigma=(var.covar), method="eigen")
    while(log.env[i,1]>log(po2_lim) | log.env[i,2]>log(depth_lim)| log.env[i,3]>log(temp_lim)){
      log.env[i,] <- rmvnorm(n = 1, mean=(env.bar - diag(var.covar) / 2), sigma=(var.covar), method="eigen")
    }
  }
  #Exponentiate
  env <- as.data.frame(exp(log.env))
  ##Metabolic Index
  #Constants
  #Ao <- 1.81
  #n <- -0.12
  #B = 1000 # size in grams, roughly average
  #avgbn <-B^n  # this works for the adult class (0.5 - 6 kg).  for the large adult, adjust
  kelvin = 273.15
  boltz = 0.000086173324
  tref <- 12
  #Calculate inverse temp
  env$invtemp <- (1 / boltz)  * ( 1 / (env$temp + 273.15) - 1 / (tref + 273.15))
  #Calculate MI
  env$mi_0.3 = env$po2*exp(0.3* env$invtemp)
  env$mi_0.7 = env$po2*exp(0.7* env$invtemp)
  
  ##Add variables to dataset
  #La/lon
  env$lat <- dat$latitude
  env$lon <- dat$longitude
  #Year
  env$year <- dat$year
  ##Log and scale depth
  env$log_depth <- log(env$depth)
  env$log_depth_scaled <- scale(env$log_depth)
  env$log_depth_scaled2 <- env$log_depth_sc ^ 2
  ##Scale mi and oxygen, just in case
  env$mi_0.3_sc <- scale(env$mi_0.3)
  env$mi_0.7_sc <- scale(env$mi_0.7)
  env$po2_sc <- scale(env$po2)
  #Seed
  env$seed <- seed
  return(env)
}

##Run environmental data simulation
n <- 100 #Number of simulations
ndata <- nrow(dat) #Number of data points
#Environmental data limits
po2_lim <- 35
temp_lim <-15
depth_lim <- 1500
#Simulate data
env_sims <- map(seq_len(n), ~simulate_environmental(ndata, dat, po2_lim, temp_lim,depth_lim)) 

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
n <- 100

## Make list of parameter values ##
model.pars <- list(b_years = b_years,
                   beta1 = beta1,
                   beta2 = beta2,
                   phi = phi,
                   p = p,
                   range = range,
                   sigma_O=sigma_O)
## Make mesh ##
mesh <- make_mesh(dat, xy_cols = c("longitude", "latitude"), n_knots=250)
simdat_ran <- lapply(env_sims, simulate_fish,
                 mesh = mesh,
                 s50 = s50,
                 delta = delta,
                 smax = smax,
                 Eo = Eo,
                 modelpars = model.pars)

### Simulate data undate atypical Eo
simdat2_ran <- lapply(env_sims, simulate_fish,
                  mesh = mesh,
                  s50 = s50,
                  delta = delta,
                  smax = smax,
                  Eo = Eo2,
                  modelpars = model.pars)
#Save? #
save <- T
if(save){
  saveRDS(simdat_ran, "data_sims_usual_ran.rds")
  saveRDS(simdat2_ran, "data_sims_weird_ran.rds")
}


# Make list of parameter names"
pars_names <- c("log_depth_scaled", "log_depth_scaled2", "mi-delta", "mi-s50", "mi-smax", "range", "sigma_O", "phi", "tweedie_p", "mi-Eo", "as.factor(year)2010","as.factor(year)2011","as.factor(year)2012", "as.factor(year)2013", "as.factor(year)2014", "as.factor(year)2015")
true_pars <- data.frame(term=c("log_depth_scaled", "log_depth_scaled2", "mi-delta", "mi-s50", "mi-smax", "range", "sigma_O", "phi", "tweedie_p", "mi-Eo", "as.factor(year)2010","as.factor(year)2011","as.factor(year)2012", "as.factor(year)2013", "as.factor(year)2014", "as.factor(year)2015"), 
                        estimate=c(beta1, beta2, delta, s50, smax, range, sigma_O, phi, p, Eo, b_years))
true_pars2 <- data.frame(term=c("log_depth_scaled", "log_depth_scaled2", "mi-delta", "mi-s50", "mi-smax", "range", "sigma_O", "phi", "tweedie_p", "mi-Eo", "as.factor(year)2010","as.factor(year)2011","as.factor(year)2012", "as.factor(year)2013", "as.factor(year)2014", "as.factor(year)2015"), 
                         estimate=c(beta1, beta2, delta, s50, smax, range, sigma_O, phi, p, Eo2, b_years))
# Set model names #
model_names <- c("Typical Case, Unconstrained", "Typical Case, Prior Constrained", "Unusual Case, Unconstrained", "Unusual Case, Prior Constrained")

## Apply to model fits ##
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
fits_ran <- lapply(simdat_ran, run_sdmTMB_noprior, 
                 start=start, mesh=mesh)
fits2_ran <- lapply(simdat_ran, run_sdmTMB_prior, 
                  start=start, mesh=mesh)
  
fits3_ran <- lapply(simdat2_ran, run_sdmTMB_noprior, 
                  start=start_unusual, mesh=mesh)
fits4_ran <- lapply(simdat2_ran, run_sdmTMB_prior, 
                  start=start_unusual, mesh=mesh)
# Save? #
save <- T
if(save){
  save(fits_ran, fits2_ran, fits3_ran, fits4_ran, file="model_fits_ran.Rdata")
}

## Parameter estimates ##
pars1 <- lapply(fits_ran, extract_pars)
pars2 <- lapply(fits2_ran, extract_pars)
pars3 <-lapply(fits3_ran, extract_pars)
pars4 <- lapply(fits4_ran, extract_pars)

pars1 <- clean_pars(pars1, fits=fits_ran)
pars2 <- clean_pars(pars2, fits=fits2_ran)
pars3 <- clean_pars(pars3, fits=fits3_ran)
pars4 <- clean_pars(pars4, fits=fits4_ran)

#Add column to label model
pars1$model <- model_names[1]
pars1$data <- "Typical Case"
pars1$analysis <- "Unconstrained"
pars2$model <- model_names[2]
pars2$data <- "Typical Case"
pars2$analysis <- "Prior Constrained"
pars3$model <- model_names[3]
pars3$data <- "Unusual Case"
pars3$analysis <- "Unconstrained"
pars4$model <- model_names[4]
pars4$data <- "Unusual Case"
pars4$analysis <- "Prior Constrained"

#Merge into one and combine
pars <- rbind(pars1,pars2, pars3, pars4)

### Parameter performance measures ###
## Average ##
avg <- aggregate(estimate ~ term+model, pars, FUN=mean)

# Average standard
sd_avg <- aggregate(estimate ~ term+model, pars, FUN=stats::sd)

## Root mean square error (accuracy) ##
# Calculate error #
pars$error <- case_when(pars$term=="mi-s50"~pars$estimate-s50,
                        pars$term=="mi-delta"~pars$estimate-delta,
                        pars$term=="mi-smax"~pars$estimate-smax,
                        pars$term=="mi-Eo"& pars$model=="Typical Case, Unconstrained"~pars$estimate-Eo,
                        pars$term=="mi-Eo"& pars$model=="Typical Case, Prior Constrained"~pars$estimate-Eo,
                        pars$term=="mi-Eo"& pars$model=="Unusual Case, Unconstrained"~pars$estimate-Eo2,
                        pars$term=="mi-Eo"& pars$model=="Unusual Case, Prior Constrained"~pars$estimate-Eo2,
                        pars$term=="log_depth_scaled"~pars$estimate-beta1,
                        pars$term=="log_depth_scaled2"~pars$estimate-beta2,
                        pars$term=="range"~pars$estimate-range,
                        pars$term=="sigma_O"~pars$estimate-sigma_O,
                        pars$term=="phi"~pars$estimate-phi,
                        pars$term=="tweedie_p"~pars$estimate-p,
                        pars$term=="as.factor(year)2010"~pars$estimate-b_years[1],
                        pars$term=="as.factor(year)2011"~pars$estimate-b_years[2],
                        pars$term=="as.factor(year)2012"~pars$estimate-b_years[3],
                        pars$term=="as.factor(year)2013"~pars$estimate-b_years[4],
                        pars$term=="as.factor(year)2014"~pars$estimate-b_years[5],
                        pars$term=="as.factor(year)2015"~pars$estimate-b_years[6])
pars$error2 <- pars$error^2
rmse <- aggregate(error2 ~ term+model, pars, FUN=sum)
rmse2 <- aggregate(error2 ~ term+model, pars, FUN=length)
rmse$n <- rmse2$error2
rmse$rmse <- sqrt(rmse$error2/rmse$n)
rmse$n <- NULL
rmse$error2 <- NULL

#Create dataframe for plotting

Eo_values_ran <- as.data.frame(matrix(nrow=4))
Eo_values_ran$V1 <- NULL
Eo_values_ran$data <- c("Typical Case", "Typical Case", "Unusual Case", "Unusual Case")
Eo_values_ran$analysis <- c( "Prior Constrained", "Unconstrained", "Prior Constrained", "Unconstrained")
Eo_values_ran$model <- c("Typical Case, Prior Constrained", "Typical Case, Unconstrained","Unusual Case, Prior Constrained", "Unusual Case, Unconstrained")
MLE_avg <- aggregate(estimate~model, subset(pars, term=="mi-Eo"), FUN=mean)
Eo_values_ran$MLE_avg <- MLE_avg$estimate
Eo_values_ran$true <- c(Eo, Eo, Eo2, Eo2)

s50_values_ran <- as.data.frame(matrix(nrow=4))
s50_values_ran$V1 <- NULL
s50_values_ran$data <- c("Typical Case", "Typical Case", "Unusual Case", "Unusual Case")
s50_values_ran$analysis <- c( "Prior Constrained", "Unconstrained", "Prior Constrained", "Unconstrained")
s50_values_ran$model <- c("Typical Case, Prior Constrained", "Typical Case, Unconstrained","Unusual Case, Prior Constrained", "Unusual Case, Unconstrained")
MLE_avg <- aggregate(estimate~model, subset(pars, term=="mi-s50"), FUN=mean)
s50_values_ran$MLE_avg <- MLE_avg$estimate
s50_values_ran$true <- c(s50, s50,s50,s50)

## Make a table ##
rmse$"average" <- avg$estimate
par_performance <- rmse
par_performance$sd <- sd_avg$estimate
colnames(par_performance) <- c("Parameter", "Model", "RMSE", "Average", "Precision")

par_performance$Bias <- case_when(par_performance$Parameter=="mi-s50"~par_performance$Average-s50,
                                  par_performance$Parameter=="mi-delta"~par_performance$Average-delta,
                                  par_performance$Parameter=="mi-smax"~par_performance$Average-smax,
                                  par_performance$Parameter=="mi-Eo"& par_performance$Model=="Typical Case, Unconstrained"~par_performance$Average-Eo,
                                  par_performance$Parameter=="mi-Eo"& par_performance$Model=="Typical Case, Prior Constrained"~par_performance$Average-Eo,
                                  par_performance$Parameter=="mi-Eo"& par_performance$Model=="Unusual Case, Unconstrained"~par_performance$Average-Eo2,
                                  par_performance$Parameter=="mi-Eo"& par_performance$Model=="Unusual Case, Prior Constrained"~par_performance$Average-Eo2,
                                  par_performance$Parameter=="log_depth_scaled"~par_performance$Average-beta1,
                                  par_performance$Parameter=="log_depth_scaled2"~par_performance$Average-beta2,
                                  par_performance$Parameter=="range"~par_performance$Average-range,
                                  par_performance$Parameter=="sigma_O"~par_performance$Average-sigma_O,
                                  par_performance$Parameter=="phi"~par_performance$Average-phi,
                                  par_performance$Parameter=="tweedie_p"~par_performance$Average-p,
                                  par_performance$Parameter=="as.factor(year)2010"~par_performance$Average-b_years[1],
                                  par_performance$Parameter=="as.factor(year)2011"~par_performance$Average-b_years[2],
                                  par_performance$Parameter=="as.factor(year)2012"~par_performance$Average-b_years[3],
                                  par_performance$Parameter=="as.factor(year)2013"~par_performance$Average-b_years[4],
                                  par_performance$Parameter=="as.factor(year)2014"~par_performance$Average-b_years[5],
                                  par_performance$Parameter=="as.factor(year)2015"~par_performance$Average-b_years[6])

Eo_performance_ran <- subset(par_performance, Parameter=="mi-Eo")
s50_performance_ran <- subset(par_performance, Parameter=="mi-s50")

### Plot parameter estimates ###
# Reorder
pars$analysis <- factor(pars$analysis, levels = c("Unconstrained", "Prior Constrained"))
Eo_values_ran$analysis <- factor(Eo_values_ran$analysis, levels = c("Unconstrained", "Prior Constrained"))
#Density plot just Eo #
ggplot(subset(pars, pars$term=="mi-Eo"), aes(x=estimate)) +
  geom_density(fill="lightblue", adjust = 1.5) +
  geom_vline(data = Eo_values, aes(xintercept = MLE_avg),linetype="dashed", size=1.2, color="darkorange", show.legend=T)+
  geom_vline(data = Eo_values, aes(xintercept = true),linetype="dashed", size=1.2)+
  facet_grid(analysis~data)+
  xlab("Eo estimate") + 
  theme(strip.text = element_text(size = 14))

# Density plot s50 #
#Reorder 
s50_values_ran$analysis <- factor(s50_values_ran$analysis, levels = c("Unconstrained", "Prior Constrained"))
ggplot(subset(pars, pars$term=="mi-s50"), aes(x=estimate)) +
  geom_density(fill="lightblue", adjust = 1.5) +
  geom_vline(data = s50_values, aes(xintercept = MLE_avg),linetype="dashed", size=1.2, color="darkorange", show.legend=T)+
  geom_vline(data = s50_values, aes(xintercept = true),linetype="dashed", size=1.2)+
  facet_grid(analysis~data)+
  xlab("s50 estimate") + 
  theme(strip.text = element_text(size = 14))

# Plot accuracy and precision comparison #
# Merge #
Eo_performance$data <- "real"
Eo_performance_ran$data <- "random"
Eo_comparison <- full_join(Eo_performance, Eo_performance_ran)
Eo_comparison$bias2 <- abs(Eo_comparison$Bias)
ggplot(data=Eo_comparison, aes(x=Precision, y=bias2))+geom_point(aes(group=c(data),color=Model, shape=data), size=5)


### Fit mis-specified depth model ###
use_previous <- T
if(use_previous){
  load("model_fits_mis.Rdata")
}

fits_mis <- lapply(simdat, run_sdmTMB_misspecified_noprior, 
                   start=start, mesh=mesh)
fits2_mis <- lapply(simdat, run_sdmTMB_prior_misspecified, 
                    start=start, mesh=mesh)

fits3_mis <- lapply(simdat2, run_sdmTMB_misspecified_noprior, 
                    start=start_unusual, mesh=mesh)
fits4_mis <- lapply(simdat2, run_sdmTMB_prior_misspecified, 
                    start=start_unusual, mesh=mesh)

save <- T
if(save){
  save(fits_mis, fits2_mis, fits3_mis, fits4_mis, file="model_fits_mis.Rdata")
}

## Parameter estimates ##
pars1 <- lapply(fits_mis, extract_pars)
pars2 <- lapply(fits2_mis, extract_pars)
pars3 <-lapply(fits3_mis, extract_pars)
pars4 <- lapply(fits4_mis, extract_pars)

pars1 <- clean_pars(pars1, fits=fits_mis)
pars2 <- clean_pars(pars2, fits=fits2_mis)
pars3 <- clean_pars(pars3, fits=fits3_mis)
pars4 <- clean_pars(pars4, fits=fits4_mis)

#Add column to label model
pars1$model <- model_names[1]
pars1$data <- "Typical Case"
pars1$analysis <- "Unconstrained"
pars2$model <- model_names[2]
pars2$data <- "Typical Case"
pars2$analysis <- "Prior Constrained"
pars3$model <- model_names[3]
pars3$data <- "Unusual Case"
pars3$analysis <- "Unconstrained"
pars4$model <- model_names[4]
pars4$data <- "Unusual Case"
pars4$analysis <- "Prior Constrained"

#Merge into one and combine
pars <- rbind(pars1,pars2, pars3, pars4)

### Parameter performance measures ###
## Average ##
avg <- aggregate(estimate ~ term+model, pars, FUN=mean)

# Average standard
sd_avg <- aggregate(estimate ~ term+model, pars, FUN=stats::sd)

## Root mean square error (accuracy) ##
# Calculate error #
pars$error <- case_when(pars$term=="mi-s50"~pars$estimate-s50,
                        pars$term=="mi-delta"~pars$estimate-delta,
                        pars$term=="mi-smax"~pars$estimate-smax,
                        pars$term=="mi-Eo"& pars$model=="Typical Case, Unconstrained"~pars$estimate-Eo,
                        pars$term=="mi-Eo"& pars$model=="Typical Case, Prior Constrained"~pars$estimate-Eo,
                        pars$term=="mi-Eo"& pars$model=="Unusual Case, Unconstrained"~pars$estimate-Eo2,
                        pars$term=="mi-Eo"& pars$model=="Unusual Case, Prior Constrained"~pars$estimate-Eo2,
                        pars$term=="log_depth_scaled"~pars$estimate-beta1,
                        pars$term=="log_depth_scaled2"~pars$estimate-beta2,
                        pars$term=="range"~pars$estimate-range,
                        pars$term=="sigma_O"~pars$estimate-sigma_O,
                        pars$term=="phi"~pars$estimate-phi,
                        pars$term=="tweedie_p"~pars$estimate-p,
                        pars$term=="as.factor(year)2010"~pars$estimate-b_years[1],
                        pars$term=="as.factor(year)2011"~pars$estimate-b_years[2],
                        pars$term=="as.factor(year)2012"~pars$estimate-b_years[3],
                        pars$term=="as.factor(year)2013"~pars$estimate-b_years[4],
                        pars$term=="as.factor(year)2014"~pars$estimate-b_years[5],
                        pars$term=="as.factor(year)2015"~pars$estimate-b_years[6])
pars$error2 <- pars$error^2
rmse <- aggregate(error2 ~ term+model, pars, FUN=sum)
rmse2 <- aggregate(error2 ~ term+model, pars, FUN=length)
rmse$n <- rmse2$error2
rmse$rmse <- sqrt(rmse$error2/rmse$n)
rmse$n <- NULL
rmse$error2 <- NULL

#Create dataframe for plotting

Eo_values_mis <- as.data.frame(matrix(nrow=4))
Eo_values_mis$V1 <- NULL
Eo_values_mis$data <- c("Typical Case", "Typical Case", "Unusual Case", "Unusual Case")
Eo_values_mis$analysis <- c( "Prior Constrained", "Unconstrained", "Prior Constrained", "Unconstrained")
Eo_values_mis$model <- c("Typical Case, Prior Constrained", "Typical Case, Unconstrained","Unusual Case, Prior Constrained", "Unusual Case, Unconstrained")
MLE_avg <- aggregate(estimate~model, subset(pars, term=="mi-Eo"), FUN=mean)
Eo_values_mis$MLE_avg <- MLE_avg$estimate
Eo_values_mis$true <- c(Eo, Eo, Eo2, Eo2)

s50_values_mis <- as.data.frame(matrix(nrow=4))
s50_values_mis$V1 <- NULL
s50_values_mis$data <- c("Typical Case", "Typical Case", "Unusual Case", "Unusual Case")
s50_values_mis$analysis <- c( "Prior Constrained", "Unconstrained", "Prior Constrained", "Unconstrained")
s50_values_mis$model <- c("Typical Case, Prior Constrained", "Typical Case, Unconstrained","Unusual Case, Prior Constrained", "Unusual Case, Unconstrained")
MLE_avg <- aggregate(estimate~model, subset(pars, term=="mi-s50"), FUN=mean)
s50_values_mis$MLE_avg <- MLE_avg$estimate
s50_values_mis$true <- c(s50, s50,s50,s50)

## Make a table ##
rmse$"average" <- avg$estimate
par_performance <- rmse
par_performance$sd <- sd_avg$estimate
colnames(par_performance) <- c("Parameter", "Model", "RMSE", "Average", "Precision")

par_performance$Bias <- case_when(par_performance$Parameter=="mi-s50"~par_performance$Average-s50,
                                  par_performance$Parameter=="mi-delta"~par_performance$Average-delta,
                                  par_performance$Parameter=="mi-smax"~par_performance$Average-smax,
                                  par_performance$Parameter=="mi-Eo"& par_performance$Model=="Typical Case, Unconstrained"~par_performance$Average-Eo,
                                  par_performance$Parameter=="mi-Eo"& par_performance$Model=="Typical Case, Prior Constrained"~par_performance$Average-Eo,
                                  par_performance$Parameter=="mi-Eo"& par_performance$Model=="Unusual Case, Unconstrained"~par_performance$Average-Eo2,
                                  par_performance$Parameter=="mi-Eo"& par_performance$Model=="Unusual Case, Prior Constrained"~par_performance$Average-Eo2,
                                  par_performance$Parameter=="log_depth_scaled"~par_performance$Average-beta1,
                                  par_performance$Parameter=="log_depth_scaled2"~par_performance$Average-beta2,
                                  par_performance$Parameter=="range"~par_performance$Average-range,
                                  par_performance$Parameter=="sigma_O"~par_performance$Average-sigma_O,
                                  par_performance$Parameter=="phi"~par_performance$Average-phi,
                                  par_performance$Parameter=="tweedie_p"~par_performance$Average-p,
                                  par_performance$Parameter=="as.factor(year)2010"~par_performance$Average-b_years[1],
                                  par_performance$Parameter=="as.factor(year)2011"~par_performance$Average-b_years[2],
                                  par_performance$Parameter=="as.factor(year)2012"~par_performance$Average-b_years[3],
                                  par_performance$Parameter=="as.factor(year)2013"~par_performance$Average-b_years[4],
                                  par_performance$Parameter=="as.factor(year)2014"~par_performance$Average-b_years[5],
                                  par_performance$Parameter=="as.factor(year)2015"~par_performance$Average-b_years[6])

Eo_performance_mis <- subset(par_performance, Parameter=="mi-Eo")
s50_performance_mis <- subset(par_performance, Parameter=="mi-s50")

### Plot parameter estimates ###
# Reorder
pars$analysis <- factor(pars$analysis, levels = c("Unconstrained", "Prior Constrained"))
Eo_values_mis$analysis <- factor(Eo_values_mis$analysis, levels = c("Unconstrained", "Prior Constrained"))
#Density plot just Eo #
ggplot(subset(pars, pars$term=="mi-Eo"), aes(x=estimate)) +
  geom_density(fill="lightblue", adjust = 1.5) +
  geom_vline(data = Eo_values, aes(xintercept = MLE_avg),linetype="dashed", size=1.2, color="darkorange", show.legend=T)+
  geom_vline(data = Eo_values, aes(xintercept = true),linetype="dashed", size=1.2)+
  facet_grid(analysis~data)+
  xlab("Eo estimate") + 
  theme(strip.text = element_text(size = 14))

# Density plot s50 #
#Reorder 
s50_values_mis$analysis <- factor(s50_values_ran$analysis, levels = c("Unconstrained", "Prior Constrained"))
ggplot(subset(pars, pars$term=="mi-s50"), aes(x=estimate)) +
  geom_density(fill="lightblue", adjust = 1.5) +
  geom_vline(data = s50_values, aes(xintercept = MLE_avg),linetype="dashed", size=1.2, color="darkorange", show.legend=T)+
  geom_vline(data = s50_values, aes(xintercept = true),linetype="dashed", size=1.2)+
  facet_grid(analysis~data)+
  xlab("s50 estimate") + 
  theme(strip.text = element_text(size = 14))

#### Covariance of parameters
## Make dataframe of all parameter estimates wide ##
pars_wide <- pivot_wider(pars, id_cols=c(id, model), names_from=term, values_from=estimate)
## Plot Eo vs s50##
pars_wide2 <- subset(pars_wide, model=="Typical Case, Prior Constrained")
ggplot(pars_wide2, aes(x=pars_wide2$"mi-Eo", y=pars_wide2$"mi-s50"))+ geom_point(size=5)+xlab("Eo estimate")+ylab("s50 estimate")
## Plot vs covariance of real data
pars_wide <- pivot_wider(pars, id_cols=c(id, model), names_from=term, values_from=estimate)
pars_wide3 <- subset(pars_wide, model=="Typical Case, Prior Constrained")
ggplot(pars_wide2, aes(x=pars_wide2$"mi-Eo", y=pars_wide2$"mi-s50"))+ geom_point(size=5)+geom_point(data=pars_wide3, aes(x=pars_wide3$"mi-Eo", y=pars_wide3$"mi-s50"), size=5, color='blue')+xlab("Eo estimate")+ylab("s50 estimate")


# Plot accuracy and precision comparison #
# Merge #
Eo_performance$data <- "real"
Eo_performance_mis$data <- "mis-specified"
Eo_comparison <- full_join(Eo_comparison, Eo_performance_mis)
Eo_comparison$bias2 <- abs(Eo_comparison$Bias)
ggplot(data=Eo_comparison, aes(x=Precision, y=bias2))+geom_point(aes(group=c(data),color=Model, shape=data), size=5)

#### Comparing po2' effect from parameter estimates ####
#calculate the po2' estimated from the Eo, use that to calculate f(po2')
simdats1 <- mapply(FUN=calculate_po2_prime, simdat,fits2_ran, SIMPLIFY=F)
simdats2 <- mapply(FUN=calculate_po2_prime, simdat2,fits4_ran, SIMPLIFY=F)

#plot true po2 vs exp(f(po2'estimated) (either as points or lines)
simdats1 <- mapply(FUN=logfun2, simdats1,"po2_prime", fits2_ran, SIMPLIFY=F)
simdats2 <- mapply(FUN=logfun2, simdats2,"po2_prime", fits4_ran, SIMPLIFY=F)

#Calculate true po2' effect
true_effect <- as.data.frame(logfun_basic(dat$mi_usual, smax, s50, delta))
colnames(true_effect) <- "mi_effect"
true_effect$mi <- dat$mi_usual

##Remove datasets with NA and then bind together, adding a column with the data frame number ##
simdats1 <- keep(.x=simdats1, .p=is.data.frame)
simdats1 <- bind_rows(simdats1, .id="df")
simdats2 <- keep(.x=simdats2, .p=is.data.frame)
simdats2 <- bind_rows(simdats2, .id="df")

ggplot(simdats1, aes(mi_usual, logmu, colour=Eo)) +
  geom_point(size=0.1)+
  scale_colour_viridis()+
  xlab("True pO2'")+
  ylab("Estimated pO2' effect")+
  theme(legend.position=c(0.9,0.2))+
  geom_line(data=true_effect, aes(x=mi, y=mi_effect), color="black", linetype="dashed", size=2)
