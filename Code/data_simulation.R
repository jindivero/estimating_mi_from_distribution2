#### Simulating data from sablefish trawl data ####

##Set parameter values for data generation and model fitting ####
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

## How many data sets to produce
n <- 100

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
## Simulation parameters ##
model.pars <- list(b_years = b_years,
                   beta1 = beta1,
                   beta2 = beta2,
                   phi = phi,
                   p = p,
                   range = range,
                   sigma_O=sigma_O)
mesh <- make_mesh(dat, xy_cols = c("X", "Y"), n_knots=250)




### Simulate data under typical Eo ####
use_previous <- T
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
save <- F
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
save <- F
if(save){
save(fits, fits2, fits3, fits4, file="model_fits.Rdata")
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
## Parameter estimates ##
pars1 <- lapply(fits, extract_pars)
pars2 <- lapply(fits2, extract_pars)
pars3 <-lapply(fits3, extract_pars)
pars4 <- lapply(fits4, extract_pars)

pars1 <- clean_pars(pars1, fits=fits)
pars2 <- clean_pars(pars2, fits=fits2)
pars3 <- clean_pars(pars3, fits=fits3)
pars4 <- clean_pars(pars4, fits=fits4)

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

Eo_values <- as.data.frame(matrix(nrow=4))
Eo_values$V1 <- NULL
Eo_values$data <- c("Typical Case", "Typical Case", "Unusual Case", "Unusual Case")
Eo_values$analysis <- c( "Prior Constrained", "Unconstrained", "Prior Constrained", "Unconstrained")
Eo_values$model <- c("Typical Case, Prior Constrained", "Typical Case, Unconstrained","Unusual Case, Prior Constrained", "Unusual Case, Unconstrained")
MLE_avg <- aggregate(estimate~model, subset(pars, term=="mi-Eo"), FUN=mean)
Eo_values$MLE_avg <- MLE_avg$estimate
Eo_values$true <- c(Eo, Eo, Eo2, Eo2)

s50_values <- as.data.frame(matrix(nrow=4))
s50_values$V1 <- NULL
s50_values$data <- c("Typical Case", "Typical Case", "Unusual Case", "Unusual Case")
s50_values$analysis <- c( "Prior Constrained", "Unconstrained", "Prior Constrained", "Unconstrained")
s50_values$model <- c("Typical Case, Prior Constrained", "Typical Case, Unconstrained","Unusual Case, Prior Constrained", "Unusual Case, Unconstrained")
MLE_avg <- aggregate(estimate~model, subset(pars, term=="mi-s50"), FUN=mean)
s50_values$MLE_avg <- MLE_avg$estimate
s50_values$true <- c(s50, s50,s50,s50)

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

Eo_performance <- subset(par_performance, Parameter=="mi-Eo")
s50_performance <- subset(par_performance, Parameter=="mi-s50")

### Plot parameter estimates ###
# Reorder
pars$analysis <- factor(pars$analysis, levels = c("Unconstrained", "Prior Constrained"))
Eo_values$analysis <- factor(Eo_values$analysis, levels = c("Unconstrained", "Prior Constrained"))
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
s50_values$analysis <- factor(s50_values$analysis, levels = c("Unconstrained", "Prior Constrained"))
ggplot(subset(pars, pars$term=="mi-s50"), aes(x=estimate)) +
  geom_density(fill="lightblue", adjust = 1.5) +
  geom_vline(data = s50_values, aes(xintercept = MLE_avg),linetype="dashed", size=1.2, color="darkorange", show.legend=T)+
  geom_vline(data = s50_values, aes(xintercept = true),linetype="dashed", size=1.2)+
  facet_grid(analysis~data)+
  xlab("s50 estimate") + 
  theme(strip.text = element_text(size = 14))


#### Cross-validation with model predictions ####
### Simulate new cross-validation data ###
simdat_cv <- map(seq_len(n), ~simulate_fish(dat = dat,
                                         mesh = mesh,
                                         s50 = s50,
                                         delta = delta,
                                         smax = smax,
                                         Eo = Eo,
                                         modelpars = model.pars))

simdat2_cv <- map(seq_len(n), ~simulate_fish(dat = dat,
                                          mesh = mesh,
                                          s50 = s50,
                                          delta = delta,
                                          smax = smax,
                                          Eo = Eo2,
                                          modelpars = model.pars))

### Predictions ###
## Function to predict from each model and each dataset ##

## Make predictions from existing model fits ##
preds1<- mapply(FUN=predict_sims, fits, simdat_cv, p, phi, SIMPLIFY=F)
preds2<- mapply(FUN=predict_sims, fits2, simdat_cv, p, phi, SIMPLIFY=F)
preds3<- mapply(FUN=predict_sims, fits3, simdat2_cv, p, phi, SIMPLIFY=F)
preds4<- mapply(FUN=predict_sims, fits4, simdat2_cv, p, phi, SIMPLIFY=F)



## Apply to each ##
nll1 <- mapply(FUN=calculate_nll, simdat_cv, preds1, "sim", "est", p, phi, SIMPLIFY=F)
nll2 <- mapply(FUN=calculate_nll, simdat_cv, preds2, "sim", "est", p, phi, SIMPLIFY=F)
nll3 <- mapply(FUN=calculate_nll, simdat2_cv, preds3, "sim", "est", p, phi, SIMPLIFY=F)
nll4 <- mapply(FUN=calculate_nll, simdat2_cv, preds4, "sim", "est", p, phi, SIMPLIFY=F)


nll_sum1 <- mapply(FUN=sum_nll, nll1, "nll", SIMPLIFY=F)
nll_sum2 <- mapply(FUN=sum_nll, nll2, "nll", SIMPLIFY=F)
nll_sum3 <- mapply(FUN=sum_nll, nll3, "nll", SIMPLIFY=F)
nll_sum4 <- mapply(FUN=sum_nll, nll4, "nll", SIMPLIFY=F)


nll1_thresh <- mapply(FUN=sum_nll_thresh, nll1, "nll", "mi_usual", SIMPLIFY=F)
nll2_thresh <- mapply(FUN=sum_nll_thresh, nll2, "nll", "mi_usual", SIMPLIFY=F)
nll3_thresh <- mapply(FUN=sum_nll_thresh, nll3, "nll", "mi_weird", SIMPLIFY=F)
nll4_thresh <- mapply(FUN=sum_nll_thresh, nll4, "nll", "mi_weird", SIMPLIFY=F)


nll1_above <- mapply(FUN=sum_nll_above, nll1, "nll", "mi_usual", SIMPLIFY=F)
nll2_above <- mapply(FUN=sum_nll_above, nll2, "nll", "mi_usual", SIMPLIFY=F)
nll3_above <- mapply(FUN=sum_nll_above, nll3, "nll", "mi_weird", SIMPLIFY=F)
nll4_above <- mapply(FUN=sum_nll_above, nll4, "nll", "mi_weird", SIMPLIFY=F)


#Subset to where observations have no fish
#nll1_thresh <- lapply(nll1, subset, pred2==0)
#nll2_thresh <- lapply(nll2, subset, pred2==0)
#nll3_thresh <- lapply(nll3, subset, pred2==0)

## Combine into one dataset for plotting##
# Overall #
nll_all <- as.data.frame(unlist(nll_sum1))
nll_all <- cbind(nll_all, unlist(nll_sum2))
nll_all <- cbind(nll_all, unlist(nll_sum3))
nll_all <- cbind(nll_all, unlist(nll_sum4))
colnames(nll_all) <- c( "Typical Case, Unconstrained", "Typical Case, Prior Constrained", "Unusual Case, Unconstrained", "Unusual Case, Prior Constrained")

# Above s95 #
nll_95 <- as.data.frame(unlist(nll1_above))
nll_95$model2 <- as.vector(unlist(nll2_above))
nll_95$model3 <- as.vector(unlist(nll3_above))
nll_95$model4 <- as.vector(unlist(nll4_above))
colnames(nll_95) <- c( "Typical Case, Unconstrained", "Typical Case, Prior Constrained", "Unusual Case, Unconstrained", "Unusual Case, Prior Constrained")

# Below s95 #
nll_below <- as.data.frame(unlist(nll1_thresh))
nll_below$model2 <- as.vector(unlist(nll2_thresh))
nll_below$model3 <- as.vector(unlist(nll3_thresh))
nll_below$model4 <- as.vector(unlist(nll4_thresh))
colnames(nll_below) <- c( "Typical Case, Unconstrained", "Typical Case, Prior Constrained", "Unusual Case, Unconstrained", "Unusual Case, Prior Constrained")

# Pivot long #
nll_all <- pivot_longer(nll_all, c(1:4), names_to="Model")
nll_below <- pivot_longer(nll_below, c(1:4), names_to="Model")
nll_95 <- pivot_longer(nll_95, c(1:4), names_to="Model")

# Relabel #
colnames(nll_all)[2] <- "Overall"
colnames(nll_below)[2] <- "Below_s50"
colnames(nll_95)[2] <- "Above_s95"

# Combine #
nll_combined <- cbind(nll_all,nll_95, nll_below)
nll_combined <- nll_combined[, c(1,2,4,6)]

### Run logistic pO2 equations and run predictions ###
## Function to run model 2 (with prior) ##
start <- matrix(0, ncol = 1, nrow = 3)
start[1,1] <-  -1.1
start[2,1] <- -1.1
start[3,1] <- 15


## Fit model to all simulated datasets ##
fits5 <- lapply(simdat, run_sdmTMB_3, 
                start=start, mesh=mesh)

start <- matrix(0, ncol = 1, nrow = 3)
start[1,1] <-  -1.1
start[2,1] <- -1
start[3,1] <- 1000

fits6 <- lapply(simdat2, run_sdmTMB_3, 
                start=start, mesh=mesh)

### Make predictions from logistic pO2 model fits ###
preds5<- mapply(FUN=predict_sims, fits5, simdat_cv, p, phi, SIMPLIFY=F)
preds6<- mapply(FUN=predict_sims, fits6, simdat2_cv, p, phi, SIMPLIFY=F)

nll5 <- mapply(FUN=calculate_nll, simdat_cv, preds5, "sim", "est", p, phi, SIMPLIFY=F)
nll6 <- mapply(FUN=calculate_nll, simdat2_cv, preds6, "sim", "est", p, phi, SIMPLIFY=F)
nll_sum5 <- mapply(FUN=sum_nll, nll5, "nll", SIMPLIFY=F)
nll_sum6 <- mapply(FUN=sum_nll, nll6, "nll", SIMPLIFY=F)

# Subset logistic(po2) into thresholds #
# Above s95 #
nll5_thresh <- mapply(FUN=sum_nll_thresh, nll5, "nll", "mi_usual", SIMPLIFY=F)
nll6_thresh <- mapply(FUN=sum_nll_thresh, nll6, "nll", "mi_unusual", SIMPLIFY=F)

nll5_above <- mapply(FUN=sum_nll_above, nll5, "nll", "mi_usual", SIMPLIFY=F)
nll6_above <- mapply(FUN=sum_nll_above, nll6, "nll", "mi_unusual", SIMPLIFY=F)

## Combine into one dataset for plotting##
# Above s95 #
nll_log_95 <- as.data.frame(unlist(nll5_above))
nll_log_95$model2 <- as.vector(unlist(nll6_above))
colnames(nll_log_95) <- c( "Typical Case, Logistic pO2", "Unusual Case, Logistic pO2")

# Below s95 #
nll_log_below <- as.data.frame(unlist(nll5_thresh))
nll_log_below$model2 <- as.vector(unlist(nll6_thresh))
colnames(nll_log_below) <- c( "Typical Case, Logistic pO2", "Unusual Case, Logistic pO2")

# Pivot long #
nll_log_below <- pivot_longer(nll_log_below, c(1:2), names_to="Model")
nll_log_95 <- pivot_longer(nll_log_95, c(1:2), names_to="Model")
nll_log <- pivot_longer(nll_log, c(1:2), names_to="Model")

# Relabel #
colnames(nll_log_below)[2] <- "Below_s50"
colnames(nll_log_95)[2] <- "Above_s95"
colnames(nll_log)[2] <- "Overall"

# Combine #
nll_log_combined <- cbind(nll_log_below, nll_log_95, nll_log)
nll_log_combined <- nll_log_combined[,c(1,2,4,6)]

nll_combined <- bind_rows(nll_combined, nll_log_combined)

#Average NLL for plotting
LL_values <- as.data.frame(matrix(nrow=6))
LL_values$V1 <- NULL
LL_values$data <- c("Typical Case", "Typical Case", "Typical Case", "Unusual Case", "Unusual Case", "Unusual Case")
LL_values$analysis <- c( "Logistic", "Prior", "Unconstrained", "Logistic", "Prior", "Unconstrained")
LL_values$model <- c("Typical Case, Logistic pO2", "Typical Case, Prior Constrained", "Typical Case, Unconstrained","Unusual Case, Logistic pO2", "Unusual Case, Prior Constrained", "Unusual Case, Unconstrained")
LL_avg <- aggregate(Below_s50~Model, nll_combined, FUN=mean)
LL_values$LL_avg <- LL_avg$Below_s50

# Add separate model and data columns (and remove the comma from the data type) #
nll_combined$data <- word(nll_combined$Model, 1,2, sep=" ")
nll_combined$data <- str_sub(nll_combined$data, 1, str_length(nll_combined$data)-1)
nll_combined$analysis <- word(nll_combined$Model, 3)

ggplot(nll_combined, aes(x=Below_s50))+
  geom_density(fill="lightblue", adjust = 1.5)+
  facet_grid(analysis~data)+
  geom_vline(data = LL_values, aes(xintercept = LL_avg),linetype="dashed", size=1.2, color="darkorange", show.legend=T)+
  ylab("Density")+
  xlab("Sum Log-Likelihood for Observations Below True s50 of pO2'")

ggplot(nll_combined, aes(x=Overall))+
  geom_density(fill="lightblue", adjust = 1.5)+
  facet_grid(analysis~data)+
  ylab("Density")+
  xlab("Sum Log-Likelihood for All Observations'")

#### Comparing po2' and predictions for all simulations (to see how predictions compare to parameters) ####
#calculate the po2' estimated from the Eo, use that to calculate f(po2')
simdats1 <- mapply(FUN=calculate_po2_prime, simdat,fits, SIMPLIFY=F)
simdats2 <- mapply(FUN=calculate_po2_prime, simdat2,fits3, SIMPLIFY=F)

#plot true po2 vs exp(f(po2'estimated) (either as points or lines)
simdats1 <- mapply(FUN=logfun2, simdats1,"po2_prime", fits, SIMPLIFY=F)
simdats2 <- mapply(FUN=logfun2, simdats2,"po2_prime", fits3, SIMPLIFY=F)

##Remove datasets with NA and then bind together, adding a column with the data frame number ##
simdats1 <- keep(.x=simdats1, .p=is.data.frame)
simdats1 <- bind_rows(simdats1, .id="df")
simdats2 <- keep(.x=simdats2, .p=is.data.frame)
simdats2 <- bind_rows(simdats2, .id="df")

ggplot(simdats1, aes(mi_usual, logmu, colour=as.factor(df))) +
  geom_point(size=0.1)+
  scale_colour_viridis(discrete=TRUE)+
  xlab("True pO2'")+
  ylab("Estimated pO2' effect")+
  theme(legend.position="none")

ggplot(simdats1, aes(mi_usual, logmu, colour=Eo)) +
  geom_point(size=0.1)+
  scale_colour_viridis()+
  xlab("True pO2'")+
  ylab("Estimated pO2' effect")

ggplot(simdats1, aes(mi_usual, logmu, colour=Eo)) +
  geom_point(size=0.1)+
  scale_colour_viridis()+
  xlab("True pO2'")+
  ylab("Estimated pO2' effect")

##### Extra stuff---ignore #####
###Comparing po2' and predictions for all simulations (to see how predictions compare to parameters) ###
## Predictions of fish density on original data for unconstrained models ##
pred1<- mapply(FUN=predict_sims, fits, simdat, p, phi, SIMPLIFY=F)
pred2<- mapply(FUN=predict_sims, fits3, simdat2, p, phi, SIMPLIFY=F)
## Calculate po2 prime from estimated parameters for each ##
simdats1 <- mapply(FUN=calculate_po2_prime, pred1,fits, SIMPLIFY=F)
simdats2 <- mapply(FUN=calculate_po2_prime, pred2,fits3, SIMPLIFY=F)

##Remove datasets with NA and then bind together, adding a column with the data frame number ##
simdats1 <- keep(.x=simdats1, .p=is.data.frame)
simdats1 <- bind_rows(simdats1, .id="df")
simdats2 <- keep(.x=simdats2, .p=is.data.frame)
simdats2 <- bind_rows(simdats2, .id="df")
## Plot all ##
# Points #
ggplot(simdats1, aes(po2_prime, pred2, colour=as.factor(df))) +
  geom_point(size=0.1)+
  scale_colour_viridis(discrete=TRUE)+
  xlab("pO2'")+
  ylab("Predicted density")+
  theme(legend.position="none")
# Density plot #
ggplot(simdats1, aes(po2_prime, pred2, colour=as.factor(df))) +
  geom_density_2d()+
  scale_colour_viridis(discrete=TRUE)+
  xlab("pO2'")+
  ylab("Predicted density")+
  theme(legend.position="none")

## Just a couple
# Points #
ggplot(subset(simdats1, simdats1$df=="1"|simdats1$df=="19"), aes(po2_prime, pred2, colour=as.factor(df))) +
  geom_point(size=2)+
 scale_colour_viridis(discrete=TRUE, labels=c('E0=0.68, s50=2.4, smax=129', "E0=0 -0.056, s50=1.4, smax=22"))+
  xlab("pO2'")+
  ylab("Predicted density")+
  theme(legend.position=c(0.7,0.8))

## Just a couple
# Density # (this looks weird)
ggplot(subset(simdats1, simdats1$df=="1"|simdats1$df=="19"), aes(po2_prime, pred2, colour=as.factor(df))) +
  stat_density_2d(
    geom = "raster",
    aes(fill = after_stat(density)),
    contour = FALSE
  ) +
  scale_colour_viridis(discrete=TRUE)+
  xlab("pO2'")+
  ylab("Predicted density")

#### Covariance of parameters
# Make dataframe of all parameter estimates wide #
pars_wide <- pivot_wider(pars, id_cols=c(id, model), names_from=term, values_from=estimate)
### Plot Eo vs s50 ###
## Plot all together ##
ggplot(pars_wide, aes(x=pars_wide$"mi-Eo", y=pars_wide$"mi-s50"))+ geom_point(aes(group=model, color=model), size=5)+xlab("Eo estimate")+ylab("s50 estimate")+
  theme(legend.position=c(0.3,0.8))
# Plot all together in facet grid # 
ggplot(pars_wide, aes(x=pars_wide$"mi-Eo", y=pars_wide$"mi-s50"))+ geom_point(size=5)+xlab("Eo estimate")+facet_wrap("model", scales="free")+ylab("s50 estimate")

# Subset just one model # 
pars_wide1 <- subset(pars_wide, pars_wide$model=="Typical Case, Unconstrained")
ggplot(subset(pars_wide1), aes(x=pars_wide1$"mi-Eo", y=pars_wide1$"mi-s50"))+ geom_point(size=5)

ggplot(pars_wide, aes(x=pars_wide$"mi-Eo", y=pars_wide$"mi-smax"))+ geom_point(aes(group=model, color=model), size=5)+xlab("Eo estimate")+ylab("s50 estimate")+
  theme(legend.position=c(0.2,0.8))
# Apply #
simdats1 <- mapply(FUN=calculate_po2_prime, simdat,fits, SIMPLIFY=F)

# Apply #
simdats1 <- mapply(FUN=logfun_all, simdats1, fits)
ggplot(bind_rows(simdats2, .id="df"), aes(po2_prime, logfun, colour=as.factor(df))) +
  geom_point(size=0.1)+
  scale_colour_viridis(discrete=TRUE)+
  xlab("pO2'")+
  ylab("Predicted density")


# Plot #
simdats1 <- keep(.x=simdats1, .p=is.data.frame)
simdats2 <- simdats1[1:5]
ggplot(bind_rows(simdats2, .id="df"), aes(po2_prime, logfun, colour=as.factor(smax))) +
  geom_point(size=0.1)+
  scale_colour_viridis(discrete=TRUE)+
  xlab("pO2'")+
  ylab("pO2' effect")

# simulation version #
smax_test <- as.data.frame(dat$mi_usual)
colnames(smax_test) <- "po2_prime"

smax_test$logmu1 <- logfun_basic(smax_test$po2_prime, smax=5, s50=s50, delta)
smax_test$logmu2 <- logfun_basic(smax_test$po2_prime, smax=20, s50=s50, delta)
smax_test$logmu3 <- logfun_basic(smax_test$po2_prime, smax=50, s50=s50, delta)
smax_test$logmu4 <- logfun_basic(smax_test$po2_prime, smax=100, s50=s50, delta)
smax_test$logmu5 <- logfun_basic(smax_test$po2_prime, smax=500, s50=s50, delta)
smax_test$logmu6 <- logfun_basic(smax_test$po2_prime, smax=1000, s50=s50, delta)
smax_test$logmu7 <- logfun_basic(smax_test$po2_prime, smax=10000, s50=s50, delta)

smax_test <- pivot_longer(smax_test, 2:8)

ggplot(smax_test, aes(x=po2_prime, y=value))+geom_point(aes(group=name, color=name))+
  scale_colour_viridis(discrete=T)+
  xlab("pO2'")+
  ylab("pO2' effect")+
  scale_colour_discrete(type="viridis", "smax value", labels=c("5", "20", "50", "100", "500", "1000", "10000"))+theme(legend.position=c(0.8,0.2))


##Correlation of Eo and number of zero observations in dataset #
# Make wide #
pars_wide <- pivot_wider(pars, id_cols=c(id, model), names_from=term, values_from=estimate)
# Plots #
ggplot(pars_wide, aes(y=pars_wide$"mi-smax", x=pars_wide$"mi-Eo"))+geom_point()+
  facet_wrap("model", scales="free")+
  geom_vline(data = subset(true_pars, true_pars$term=="mi-Eo"), aes(xintercept = estimate),linetype="dashed", size=1.2)+
  geom_hline(data = subset(true_pars,  true_pars$term=='mi-smax'), aes(yintercept = estimate),linetype="dashed", size=1.2)

ggplot(pars_wide, aes(y=pars_wide$"mi-s50", x=pars_wide$"mi-Eo"))+geom_point()+
  facet_wrap("model", scales="free")+
  geom_vline(data = subset(true_pars, true_pars$term=="mi-Eo"), aes(xintercept = estimate),linetype="dashed", size=1.2)+
  geom_hline(data = subset(true_pars,  true_pars$term=='mi-s50'), aes(yintercept = estimate),linetype="dashed", size=1.2)

ggplot(pars_wide, aes(y=pars_wide$"mi-delta", x=pars_wide$"mi-Eo"))+geom_point()+
  facet_wrap("model", scales="free")+
  geom_vline(data = subset(true_pars, true_pars$term=="mi-Eo"), aes(xintercept = estimate),linetype="dashed", size=1.2)+
  geom_hline(data = subset(true_pars,  true_pars$term=='mi-delta'), aes(yintercept = estimate),linetype="dashed", size=1.2)

ggplot(pars_wide, aes(y=pars_wide$"log_depth_scaled2", x=pars_wide$"mi-Eo"))+geom_point()+
  facet_wrap("model", scales="free")


counts_below_zero <- lapply(simdat, count_below_zero,threshold=(s50+delta))
counts_below_zero <- as.data.frame(unlist(counts_below_zero))
counts_below_zero$id <- as.character(1:100)
counts_below_zero2 <- lapply(simdat2, count_below_zero,threshold=(s50+delta))
counts_below_zero2 <- as.data.frame(unlist(counts_below_zero2))
counts_below_zero2$id <- as.character(1:100)

# Add to pars_wide #
pars_wide1 <- subset(pars_wide, model=='Unusual Case, Unconstrained')
pars_wide1 <- left_join(pars_wide1, counts_below_zero2, by="id")

ggplot(pars_wide1, aes(y=pars_wide1$"unlist(counts_below_zero2)", x=pars_wide1$"mi-Eo"))+geom_point()+
  facet_wrap("model", scales="free")


### Other plots of parameter estimates###
#plot boxplots of a single model #
ggplot(pars1, aes(y=estimate, x=term))+geom_boxplot()+facet_wrap("term", scales="free")+geom_hline(data = true_pars, aes(yintercept = estimate),linetype="dashed", size=1.2)+theme(legend.position="left", strip.text = element_blank())
ggplot(pars2, aes(y=estimate, x=term))+geom_boxplot()+facet_wrap("term", scales="free")+geom_hline(data = true_pars, aes(yintercept = estimate),linetype="dashed", size=1.2)+theme(legend.position="left", strip.text = element_blank())
ggplot(pars3, aes(y=estimate, x=term))+geom_boxplot()+facet_wrap("term", scales="free")+geom_hline(data = true_pars2, aes(yintercept = estimate),linetype="dashed", size=1.2)+theme(legend.position="left", strip.text = element_blank())
ggplot(pars4, aes(y=estimate, x=term))+geom_boxplot()+facet_wrap("term", scales="free")+geom_hline(data = true_pars2, aes(yintercept = estimate),linetype="dashed", size=1.2)+theme(legend.position="left", strip.text = element_blank())

#Plot boxplots of all #
ggplot(pars, aes(y=estimate, x=model))+geom_boxplot(aes(group=model, fill=model))+
  facet_wrap("term", scales="free")+  
  geom_hline(data = true_pars2, aes(yintercept = estimate),linetype="dashed", size=1.2)+
  geom_hline(data = true_pars, aes(yintercept = estimate),linetype="dashed", size=1.2)+
  theme(legend.position="left")+
  scale_x_discrete(labels=c("1", "2", "3", "4"))

#Plot boxplots of just Eo and logistic parameters #
ggplot(subset(pars, pars$term=="mi-Eo"|pars$term=="mi-smax"|pars$term=="mi-s50"|pars$term=="mi-delta"), aes(y=estimate, x=model))+geom_boxplot(aes(group=model, fill=model))+
  facet_wrap("term", scales="free")+  
  geom_hline(data = subset(true_pars, true_pars2$term=="mi-Eo"|true_pars$term=="mi-smax"|true_pars$term=="mi-s50"|true_pars$term=="mi-delta"), aes(yintercept = estimate),linetype="dashed", size=1.2)+
  geom_hline(data = subset(true_pars2, true_pars2$term=="mi-Eo"), aes(yintercept = estimate),linetype="dashed", size=1.2)+
  theme(legend.position="left")+
  scale_x_discrete(labels=c("1", "2", "3", "4"))

#Restrict smax #
pars$estimate2 <- ifelse(pars$estimate>1000 & pars$term=="mi-smax", 1000, pars$estimate)

# Plot just Eo and smax #
ggplot(subset(pars, pars$term=="mi-Eo"|pars$term=="mi-smax"|pars$term=="mi-s50"|pars$term=="mi-delta"), aes(y=estimate2, x=model))+geom_boxplot(aes(group=model, fill=model))+
  facet_wrap("term", scales="free")+  
  geom_hline(data = subset(true_pars, true_pars2$term=="mi-Eo"|true_pars$term=="mi-smax"|true_pars$term=="mi-s50"|true_pars$term=="mi-delta"), aes(yintercept = estimate),linetype="dashed", size=1.2)+
  geom_hline(data = subset(true_pars2, true_pars2$term=="mi-Eo"), aes(yintercept = estimate),linetype="dashed", size=1.2)+
  theme(legend.position="left")+
  scale_x_discrete(labels=c("1", "2", "3", "4"))

### Other Performance Evlauation Options ###

# SD of average estimate (precision) #
sd_avg <- aggregate(estimate ~ term+model, pars, FUN=stats::sd)
# Range of average estimates (precision) #
range_pars <- aggregate(estimate ~ term+model, pars, FUN=fivenum)
# Range of standard error (precision) #
sd_range <- aggregate(std.error~ term+model, pars, FUN=fivenum)
# Standard deviation of standard error #
sd_sd <- aggregate(std.error ~ term+model, pars, FUN=stats::sd)
# Precision #
#95% range (upper 95% CI--lower 95% CI) #
pars$con.range <- pars$conf.high-pars$conf.low
#Average 95% range
avg_conf.range <- aggregate(con.range ~ term+model, pars, FUN=mean)

### Simulation diagnostics ###
## Convergence/fitting errors ##
convergence1 <- lapply(fits, extract_convergence)
convergence1 <- bind_rows(convergence1)
convergence2 <- lapply(fits2, extract_convergence)
convergence2 <- bind_rows(convergence2)
convergence3 <- lapply(fits3, extract_convergence)
convergence3 <- bind_rows(convergence2)
convergence4 <- lapply(fits4, extract_convergence)
convergence4 <- bind_rows(convergence2)

# Add column for data simulation and model
convergence1$sim <- 1:100
convergence1$model <- model_names[1]
convergence2$sim <- 1:100
convergence2$model <- model_names[2]
convergence3$sim <- 1:100
convergence3$model <- model_names[3]
convergence4$sim <- 1:100
convergence4$model <- model_names[4]

#Bind all together
convergence <-bind_rows(convergence1, convergence2,convergence3, convergence4)
#Truncate message
convergence$message <- substr(convergence$message, 1, 100)
##Summarize number of each type of convergence for each model
cons <- convergence %>% group_by(message, model) %>% summarise(count=n()) 
##Flip to wide
cons <-pivot_wider(cons, names_from=model, values_from=count)
colnames(cons)[1] <- "metric"

## Hessian matrix ##
pdHess1 <- lapply(fits, extract_pdHess)
pdHess <-bind_rows(pdHess1)
hess<- extract_pdHess2(pdHess)

pdHess2 <- lapply(fits2, extract_pdHess)
pdHess2 <-bind_rows(pdHess2)
hess2<- extract_pdHess2(pdHess2)

pdHess3 <- lapply(fits3, extract_pdHess)
pdHess3 <-bind_rows(pdHess3)
hess3<- extract_pdHess2(pdHess3)

pdHess4 <- lapply(fits4, extract_pdHess)
pdHess4 <-bind_rows(pdHess4)
hess4<- extract_pdHess2(pdHess4)

hess <- as.data.frame(matrix(nrow=1, ncol=3))
colnames(hess) <- model_names
hess[1,1] <- "Positive definite hessian matrix"
hess[1,2] <- extract_pdHess2(pdHess1)
hess[1,3] <- extract_pdHess2(pdHess2)
hess[1,4] <- extract_pdHess2(pdHess3)
hess[1,5] <- extract_pdHess2(pdHess4)
# Add to convergence table #
diagnostics <- rbind(cons,setNames(hess, names(cons)))

## Gradients ##

# Apply to all models #
grad1<- lapply(fits, extract_grad)
grad2<- lapply(fits2, extract_grad)
grad3<- lapply(fits3, extract_grad)
grad4<- lapply(fits4, extract_grad)

grads1<- clean_grad(grad1, fits)
grads2<- clean_grad(grad2, fits2)
grads3<- clean_grad(grad3, fits3)
grads4<- clean_grad(grad4, fits4)

# Bind together #
grads1$model2 <-grads2$big_grad
grads1$model3 <-grads3$big_grad
grads1$model4 <-grads4$big_grad
colnames(grads1) <-c("metric", model_names[1], model_names[2], model_names[3], model_names[4])
grads1$metric <- paste0("Low gradient ", grads1$metric)

## NA in SD ##
pars$NAs <- ifelse(pars$std.error=="NaN"|is.na(pars$std.error), 0,1)
NAs <- aggregate(pars$NAs~pars$term+pars$model, FUN=sum)
colnames(NAs) <- c("metric", "model", "value")
NAs$metric <- paste0("SE estimated", NAs$metric)
# Flip to wide #
NAs <-pivot_wider(NAs, names_from=model, values_from=value)

## Combine into table summarizing all diagnostics ##
diagnostics <- bind_rows(diagnostics, grads1, NAs)

#Plot RMSE
#ggplot(rmse, aes(x=model, y=rmse))+geom_point(aes(group=model, color=model))

# Average standard error (precision)--also to compare to SD of prior #
#std_error_within_model <- aggregate(std.error ~ term+model, pars, FUN=mean) 
