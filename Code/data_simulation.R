#### Simulating data from sablefish trawl data ####

### Install Packages ###
install.packages("remotes")
library(remotes)
install.packages("devtools")
library(devtools)
install.packages("pkgbuild")
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
install.packages("ggpubr")
library(ggpubr)
library(scales)
library(visreg)
library(ggeffects)
library(stringr)
library(TMB)
install.packages("tweedie")
library(tweedie)
install.packages("ggridges")
library(ggridges)
library(viridis)

### Set ggplot themes ###
theme_set(theme_bw(base_size = 25))
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

##Set parameter values for data generation and model fitting
x50 <-0.88 # same as sablefish= 0.88; had been 2 in previous data simulation
delta <- 0.56 #same as sablefish = 0.57; had been 2 in previous data simulation
smax <- 50 # maximum effect of MI; had been 4 in previous data simulation
Eo <- 0.291
Eo2 <- 0.733

#Set number of iterations
n <-100 #Number of simulations
simdat <- map(seq_len(n), ~simulate_fish(dat = dat,
                        mesh = mesh,
                        x50 = x50,
                        delta = delta,
                        smax = smax,
                        Eo = Eo))

simdat2 <- map(seq_len(n), ~simulate_fish(dat = dat,
                            mesh = mesh,
                            x50 = x50,
                            delta = delta,
                            smax = smax,
                            Eo = Eo2))

###Sanity checks on simulated data

## Compare one to the real sablefish data ##
sim_test <- simdat[[1]]
sim_test2 <- simdat2[[1]]

ggplot(dat, aes(x=po2_s,y=cpue_kg_km2))+geom_point(size=0.6)+geom_point(sim_test, mapping=aes(x=po2_s, y=sim), color="blue", size=0.6)
ggplot(sim_test, aes(x=po2,y=observed))+geom_point()
ggplot(sim_test2, aes(x=po2,y=observed))+geom_point()

#Save simulated data
saveRDS(data_sims_usual, "data_sims_usual.rds")
saveRDS(data_sims_weird, "data_sims_weird.rds")
saveRDS(data_sims_usual_CV, "data_sims_usual_CV.rds")
saveRDS(data_sims_weird_CV, "data_sims_weird_CV.rds")

## Fit Model 1: No priors ##
##Set starting parameters
#Correct values
start <- matrix(0, ncol = 1, nrow = 4)
start[1,1] <- x50
start[2,1] <- delta
start[3,1] <- smax
start[4,1] <- Eo

#Just above zero
start2 <- matrix(0.08, ncol = 1, nrow = 4)

#Slightly farther
start3 <- matrix(0, ncol = 1, nrow = 4)
start3[1,1] <- x50*0.3
start3[2,1] <- delta*0.3
start3[3,1] <- smax*0.3
start3[4,1] <- Eo*0.3

#With weird Eo
start4 <- matrix(0, ncol = 1, nrow = 4)
start4[1,1] <- x50
start4[2,1] <- delta
start4[3,1] <- smax
start4[4,1] <- Eo2

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
fits <- lapply(simdat, run_sdmTMB_1, 
                start=start, mesh=mesh)
fits2 <- lapply(simdat, run_sdmTMB_2, 
                start=start, mesh=mesh)

fits3 <- lapply(simdat2, run_sdmTMB_1, 
                start=start4, mesh=mesh)
fits4 <- lapply(simdat2, run_sdmTMB_2, 
                start=start4, mesh=mesh)

#Save models
#save(fits, fits1, fits2, fits3, fits_b, fits_c, file="model_fits.Rdata")
save(fits, fits2, fits3, fits4, file="model_fits.Rdata")

### Evaluate model performance ###
## Set true pars vectors of values and name ##
x50 <- 0.88
delta <-0.56
b_years <- c(4.47,4.53,4.44,4.43,4.68,4.68)
beta1 <- 1.5
beta2 <- -1
smax <- 50
phi <- 16 
p <- 1.51
range <- 85
sigma_O <- 1.77
Eo <- 0.291
Eo2 <- 0.733
# Make list of parameter names"
pars_names <- c("log_depth_scaled", "log_depth_scaled2", "mi-delta", "mi-s50", "mi-smax", "range", "sigma_O", "phi", "tweedie_p", "mi-Eo", "as.factor(year)2010","as.factor(year)2011","as.factor(year)2012", "as.factor(year)2013", "as.factor(year)2014", "as.factor(year)2015")
true_pars <- data.frame(term=c("log_depth_scaled", "log_depth_scaled2", "mi-delta", "mi-s50", "mi-smax", "range", "sigma_O", "phi", "tweedie_p", "mi-Eo", "as.factor(year)2010","as.factor(year)2011","as.factor(year)2012", "as.factor(year)2013", "as.factor(year)2014", "as.factor(year)2015"), 
                         estimate=c(beta1, beta2, delta, x50, smax, range, sigma_O, phi, p, Eo, b_years))
true_pars2 <- data.frame(term=c("log_depth_scaled", "log_depth_scaled2", "mi-delta", "mi-s50", "mi-smax", "range", "sigma_O", "phi", "tweedie_p", "mi-Eo", "as.factor(year)2010","as.factor(year)2011","as.factor(year)2012", "as.factor(year)2013", "as.factor(year)2014", "as.factor(year)2015"), 
                         estimate=c(beta1, beta2, delta, x50, smax, range, sigma_O, phi, p, Eo2, b_years))
# Set model names #
model_names <- c("Typical Case, Unconstrained", "Typical Case, Prior Constrained", "Unusual Case, Unconstrained", "Unusual Case, Prior Constrained")


## Functions for extracting parameter estimates and diagnostics ##
# Extract parameter estimates #
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

# Function to clean up pars for plotting #
clean_pars <- function(pars, fits){
  names(pars) <- c(1:length(fits))
  #Remove models with errors
  pars <- keep(pars, function(x) !is.logical(x))
  #Combine into single dataframe, with column of simulation number
  pars <- bind_rows(pars,.id="id")
  return(pars)
}

# Aggregate errors #
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

# Extract if positive definite hessian #

# Extract hessian matrix T/F value #
extract_pdHess <- function(x){
  if(!is.character(x)){
    pdh <- as.data.frame(x$pos_def_hessian)
    return(pdh)
  }
}

# Sum up across all iterations #
extract_pdHess2 <- function(pdHess){
  pdHess <- bind_rows(pdHess)
  pdHess$sim <- as.numeric(row.names(pdHess))
  pdHess$sim <- as.character(pdHess$sim)
  pdHess$neg_def_hessian <- ifelse(pdHess$"x$pos_def_hessian"=="TRUE", 1,0)
  number_negdHess <- print(sum(pdHess$neg_def_hessian))
  return(number_negdHess)
}

# Extract gradients for each parameter #
extract_grad <- function(x){
  if(!is.character(x)){
    grad <- as.data.frame(x$gradients)
    return(grad)
  }
}

# Clean up gradients into dataframe #
clean_grad <- function(grad, fits){
  if(!is.null(grad)){
    grad <- bind_cols(grad)
  }
  grad <- as.data.frame(t(grad))
  names <- names(fits[[14]]$sd_report$par.fixed)
  colnames(grad) <- pars_names
  grad$sim <- 1:nrow(grad)
  grad <- pivot_longer(grad,cols=1:11, names_to="term", values_to="gradient")
  grad$big_grad <- ifelse(grad$gradient<0.001, 1,0)
  large_grad <- aggregate(big_grad ~ term, grad, FUN=sum)
  #large_grad$proportion <- large_grad$big_grad/length(fits)
  return(large_grad)
}

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
pars2$model <- model_names[2]
pars3$model <- model_names[3]
pars4$model <- model_names[4]

#plot a single model
ggplot(pars1, aes(y=estimate, x=term))+geom_boxplot()+facet_wrap("term", scales="free")+geom_hline(data = true_pars, aes(yintercept = estimate),linetype="dashed", size=1.2)+theme(legend.position="left", strip.text = element_blank())
ggplot(pars2, aes(y=estimate, x=term))+geom_boxplot()+facet_wrap("term", scales="free")+geom_hline(data = true_pars, aes(yintercept = estimate),linetype="dashed", size=1.2)+theme(legend.position="left", strip.text = element_blank())
ggplot(pars3, aes(y=estimate, x=term))+geom_boxplot()+facet_wrap("term", scales="free")+geom_hline(data = true_pars2, aes(yintercept = estimate),linetype="dashed", size=1.2)+theme(legend.position="left", strip.text = element_blank())
ggplot(pars4, aes(y=estimate, x=term))+geom_boxplot()+facet_wrap("term", scales="free")+geom_hline(data = true_pars2, aes(yintercept = estimate),linetype="dashed", size=1.2)+theme(legend.position="left", strip.text = element_blank())

#Merge into one and combine
pars <- rbind(pars1,pars2, pars3, pars4)
ggplot(pars, aes(y=estimate, x=model))+geom_boxplot(aes(group=model, fill=model))+
facet_wrap("term", scales="free")+  
  geom_hline(data = true_pars2, aes(yintercept = estimate),linetype="dashed", size=1.2)+
  geom_hline(data = true_pars, aes(yintercept = estimate),linetype="dashed", size=1.2)+
  theme(legend.position="left")+
  scale_x_discrete(labels=c("1", "2", "3", "4"))

#Just Eo 
ggplot(subset(pars, pars$term=="mi-Eo"), aes(y=estimate, x=model))+geom_boxplot(aes(group=model, fill=model))+
  geom_hline(data = subset(true_pars, true_pars$term=="mi-Eo"), aes(yintercept = estimate),linetype="dashed", size=1.2)+
  geom_hline(data = subset(true_pars2, true_pars2$term=="mi-Eo"), aes(yintercept = estimate),linetype="dashed", size=1.2)+
  theme(legend.position="left")+
  scale_x_discrete(labels=c("1", "2", "3", "4"))
#Eo and logistic parameters
ggplot(subset(pars, pars$term=="mi-Eo"|pars$term=="mi-smax"|pars$term=="mi-s50"|pars$term=="mi-delta"), aes(y=estimate, x=model))+geom_boxplot(aes(group=model, fill=model))+
  facet_wrap("term", scales="free")+  
  geom_hline(data = subset(true_pars, true_pars2$term=="mi-Eo"|true_pars$term=="mi-smax"|true_pars$term=="mi-s50"|true_pars$term=="mi-delta"), aes(yintercept = estimate),linetype="dashed", size=1.2)+
  geom_hline(data = subset(true_pars2, true_pars2$term=="mi-Eo"), aes(yintercept = estimate),linetype="dashed", size=1.2)+
  theme(legend.position="left")+
  scale_x_discrete(labels=c("1", "2", "3", "4"))

#Restrict smax
pars$estimate2 <- ifelse(pars$estimate>1000 & pars$term=="mi-smax", 1000, pars$estimate)

#Eo and smax
ggplot(subset(pars, pars$term=="mi-Eo"|pars$term=="mi-smax"|pars$term=="mi-s50"|pars$term=="mi-delta"), aes(y=estimate2, x=model))+geom_boxplot(aes(group=model, fill=model))+
  facet_wrap("term", scales="free")+  
  geom_hline(data = subset(true_pars, true_pars2$term=="mi-Eo"|true_pars$term=="mi-smax"|true_pars$term=="mi-s50"|true_pars$term=="mi-delta"), aes(yintercept = estimate),linetype="dashed", size=1.2)+
  geom_hline(data = subset(true_pars2, true_pars2$term=="mi-Eo"), aes(yintercept = estimate),linetype="dashed", size=1.2)+
  theme(legend.position="left")+
  scale_x_discrete(labels=c("1", "2", "3", "4"))

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

### Parameter performance measures ###
## Average ##
avg <- aggregate(estimate ~ term+model, pars, FUN=mean)
# Average standard
sd_avg <- aggregate(estimate ~ term+model, pars, FUN=stats::sd)

## Root mean square error (accuracy) ##
# Calculate error #
pars$error <- case_when(pars$term=="mi-s50"~pars$estimate-x50,
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

#Plot RMSE
ggplot(rmse, aes(x=model, y=rmse))+geom_point(aes(group=model, color=model))

# Average standard error (precision)--also to compare to SD of prior #
std_error_within_model <- aggregate(std.error ~ term+model, pars, FUN=mean) 

## Make a table ##
rmse$"average" <- avg$estimate
rmse$sd <- sd_avg$estimate
rmse <- left_join(rmse, std_error_within_model)
colnames(rmse) <- c("Parameter", "Model", "RMSE of inter-iteration average", "Inter-iteration average", "Inter-iteration SD of average estimate", "Average intra-iteration SE")

## Other Options ##

# SD of average estimate (precision) #
#sd_avg <- aggregate(estimate ~ term+model, pars, FUN=stats::sd)
# Range of average estimates (precision) #
#range_pars <- aggregate(estimate ~ term+model, pars, FUN=fivenum)
# Range of standard error (precision) #
#sd_range <- aggregate(std.error~ term+model, pars, FUN=fivenum)
# Standard deviation of standard error #
#sd_sd <- aggregate(std.error ~ term+model, pars, FUN=stats::sd)
# Precision #
# 95% range (upper 95% CI--lower 95% CI) #
#pars$con.range <- pars$conf.high-pars$conf.low
#Average 95% range
#avg_conf.range <- aggregate(con.range ~ term+model, pars, FUN=mean)

#### Cross-validation with model predictions ####
### Simulate new cross-validation data ###
n <-100 #Number of simulations
simdat_cv <- map(seq_len(n), ~simulate_fish(dat = dat,
                                         mesh = mesh,
                                         x50 = x50,
                                         delta = delta,
                                         smax = smax,
                                         Eo = Eo))

simdat2_cv <- map(seq_len(n), ~simulate_fish(dat = dat,
                                          mesh = mesh,
                                          x50 = x50,
                                          delta = delta,
                                          smax = smax,
                                          Eo = Eo2))

### Predictions ###
## Function to predict from each model and each dataset ##
predict_sims <- function(x, new_data, p, phi){
  if(!is.character(x)){
    preds <- predict(x, newdata=new_data, type="response", return_tmb_object=F)
    #Add observation error from predictions with Tweedie parameters from model fi
    preds$pred2 <- rTweedie(preds$est, p = p, phi = phi)
  }
  if(is.character(x)){
    preds <- NA
  }
    return(preds)
  }
 
## Make predictions from existing model fits ##
preds1<- mapply(FUN=predict_sims, fits, simdat_cv, p, phi, SIMPLIFY=F)
preds2<- mapply(FUN=predict_sims, fits2, simdat_cv, p, phi, SIMPLIFY=F)
preds3<- mapply(FUN=predict_sims, fits3, simdat2_cv, p, phi, SIMPLIFY=F)
preds4<- mapply(FUN=predict_sims, fits4, simdat2_cv, p, phi, SIMPLIFY=F)

### Negative log probability of observations given model fit ###
calculate_nll <- function(data_observed, data_predicted, column_obs, column_preds, ps, phis){
  if(is.data.frame(data_predicted)){
  observed <- as.numeric(data_observed[ , column_obs])
  predicted <- as.numeric(data_predicted[, column_preds])
  nll <- dtweedie(y=observed, mu=predicted,power=p, phi=phi)
  predicted2 <- cbind(data_predicted, nll)}
  if(!is.data.frame(data_predicted)){
    predicted2 <- NA
  }
  return(predicted2)
  }

## Apply to each ##
nll1 <- mapply(FUN=calculate_nll, simdat_cv, preds1, "sim", "est", p, phi, SIMPLIFY=F)
nll2 <- mapply(FUN=calculate_nll, simdat_cv, preds2, "sim", "est", p, phi, SIMPLIFY=F)
nll3 <- mapply(FUN=calculate_nll, simdat2_cv, preds3, "sim", "est", p, phi, SIMPLIFY=F)
nll4 <- mapply(FUN=calculate_nll, simdat2_cv, preds4, "sim", "est", p, phi, SIMPLIFY=F)

# Sum overall #
sum_nll <- function(nll, column_nll){
  if(is.data.frame(nll)){
  nlls <- as.numeric(nll[ , column_nll])
  sum_nll <- sum(nlls)
  }
  if(!is.data.frame(nll)){
    sum_nll <- NA
  }
  return(sum_nll)
}

nll_sum1 <- mapply(FUN=sum_nll, nll1, "nll", SIMPLIFY=F)
nll_sum2 <- mapply(FUN=sum_nll, nll2, "nll", SIMPLIFY=F)
nll_sum3 <- mapply(FUN=sum_nll, nll3, "nll", SIMPLIFY=F)
nll_sum4 <- mapply(FUN=sum_nll, nll4, "nll", SIMPLIFY=F)

## Below simulated true threshold (=s95) ##
s95 <- x50+delta
sum_nll_thresh <- function(nll, column_nll, mi){
  if(is.data.frame(nll)){
    nll <- filter(nll, mi <x50) 
    nlls <- as.numeric(nll[ , column_nll])
    sum_nll <- sum(nlls)
  }
  if(!is.data.frame(nll)){
    sum_nll <- NA
  }
  return(sum_nll)
}

nll1_thresh <- mapply(FUN=sum_nll_thresh, nll1, "nll", "mi_usual", SIMPLIFY=F)
nll2_thresh <- mapply(FUN=sum_nll_thresh, nll2, "nll", "mi_usual", SIMPLIFY=F)
nll3_thresh <- mapply(FUN=sum_nll_thresh, nll3, "nll", "mi_weird", SIMPLIFY=F)
nll4_thresh <- mapply(FUN=sum_nll_thresh, nll4, "nll", "mi_weird", SIMPLIFY=F)

## Above simulated true threshold (=s95) ##
sum_nll_above <- function(nll, column_nll, mi){
  if(is.data.frame(nll)){
    nll <- filter(nll, mi >s95) 
    nlls <- as.numeric(nll[ , column_nll])
    sum_nll <- sum(nlls)
  }
  if(!is.data.frame(nll)){
    sum_nll <- NA
  }
  return(sum_nll)
}
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

# Plot #
ggplot(nll_combined, aes(x=Above_s95, y=Below_s50))+
  stat_density_2d(geom = "raster",aes(fill = after_stat(density)),contour = FALSE)+
  scale_fill_viridis_c()+geom_point(color="white", size=3.5, shape="*")+
  facet_grid("Model")+
  xlab("Above s95")+
  ylab("Below s50")

ggplot(nll_combined, aes(x=Overall, y=Below_s50))+
  stat_density_2d(geom = "raster",aes(fill = after_stat(density)),contour = FALSE)+
  scale_fill_viridis_c()+geom_point(color="white", size=3.5, shape="*")+
  facet_grid("Model")+
  xlab("All Data")+
  ylab("Below s50")

### Run logistic pO2 equations and run predictions ###
## Function to run model 2 (with prior) ##
start <- matrix(0, ncol = 1, nrow = 3)
start[1,1] <-  -1.1
start[2,1] <- -1.1
start[3,1] <- 15

run_sdmTMB_3 <- function(simdat, start, mesh) {
  m2 <- try(sdmTMB(sim ~ -1+as.factor(year)+logistic(po2_s)+log_depth_scaled+log_depth_scaled2, 
                   data = simdat, 
                   spatial = "on",
                   spatiotemporal="off",
                   mesh=mesh,
                   family =tweedie(link="log"),
                   control = sdmTMBcontrol(
                     start = list(b_threshold = start),
                     newton_loops = 2)))
  
  try(tidy(m2))
  try(return(m2))
}

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

nll_log <- as.data.frame(unlist(nll_sum5))
nll_log <- cbind(nll_log, unlist(nll_sum6))
colnames(nll_log) <- c( "Typical Case, Logistic pO2", "Unusual Case, Logistic pO2")

# Combine with other #
nll_combined$logistic_typical <- nll_log[,1]
nll_combined$logistic_unusual <- nll_log[,2]

# Plot #
ggplot(subset(nll_combined, Model=="Typical Case, Unconstrained"|Model=="Typical Case, Prior Constrained"), aes(x=logistic_typical, y=Overall))+
  stat_density_2d(geom = "raster",aes(fill = after_stat(density)),contour = FALSE)+
  scale_fill_viridis_c()+geom_point(color="white", size=3.5, shape="*")+
  facet_grid("Model")+
  xlab("Logistic Typical")+
  ylab("pO2'")

ggplot(subset(nll_combined, Model=="Unusual Case, Unconstrained"|Model=="Unusual Case, Prior Constrained"), aes(x=logistic_unusual, y=Overall))+
  stat_density_2d(geom = "raster",aes(fill = after_stat(density)),contour = FALSE)+
  scale_fill_viridis_c()+geom_point(color="white", size=3.5, shape="*")+
  facet_grid("Model")+
  xlab("Logistic Unusual")+
  ylab("pO2'")

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

# Reorder factor # 
nll_combined$Model <-factor(nll_combined$Model, levels=c("Typical Case, Unconstrained", "Unusual Case, Unconstrained","Typical Case, Prior Constrained", "Unusual Case, Prior Constrained","Typical Case, Logistic pO2","Unusual Case, Logistic pO2"))
# Plot #
ggplot(nll_combined, aes(x=Above_s95, y=Below_s50))+
  stat_density_2d(geom = "raster",aes(fill = after_stat(density)),contour = FALSE)+
  scale_fill_viridis_c()+geom_point(color="white", size=3.5, shape="*")+
  facet_wrap("Model", ncol=2)+
  xlab("Above s95")+
  ylab("Below s50")

#1D kernel density
ggplot(nll_combined, aes(y=Below_s50, x=Model, group=Model, fill = Model))+
  geom_violin(aes(group=Model))+
  geom_jitter(shape=16, position=position_jitter(0.2))+
 # facet_wrap("Model", ncol=2)+
  xlab("Model")+
  ylab("Log-Likelihood for Observations Below true s50 of pO2'")

ggplot(nll_combined, aes(x=Below_s50, y=Model, group=Model, fill = Model))+
  geom_density_ridges(alpha=0.1)+
  #geom_point(color="white", size=3.5, shape="*")+
  #scale_fill_viridis_c()+
  #facet_wrap("Model", ncol=2)+
  ylab("Density")+
  xlab("Sum Log-Likelihood for Observations Below True s50 of pO2'")

ggplot(nll_combined, aes(y=Below_s50,x=Model, group=Model, fill = Model))+
  geom_boxplot(aes(group=Model))+
 # geom_jitter(shape=16, position=position_jitter(0.2))+
  # facet_wrap("Model", ncol=2)+
  xlab("Model")+
  ylab("Log-Likelihood for Observations Below true s50 of pO2'")

### Logistic shape comparison ###
## Calculate po2_prime effect for each model ##
# Function #
calculate_po2_prime <- function(dat, model) {
  if(!is.character(model)){
  parfit <- model$sd_report
  npars <- length(parfit$value)
  parnames <- names(parfit$value)

  Eo <- parfit$value[grep("Eo", parnames)]
  dat$po2_prime <-  dat$po2 * exp(Eo * dat$invtemp)
  dat$smax <- parfit$value[grep("s_max", parnames)]
  }
  if(is.character(model)){
    dat <- NA
}
  return(dat)
}
# Apply #
simdats1 <- mapply(FUN=calculate_po2_prime, simdat,fits, SIMPLIFY=F)


## Calculate po2 prime effect ##
# Function #
logfun_all <- function(dat, model) {
  if(!is.character(model)){
    dat$logfun <- exp(logfun(dat$po2_prime, model = model, mi = T))
  }
  if(is.character(model)){
    dat <- NA
  }
  return(dat)
}

# Apply #
simdats1 <- mapply(FUN=logfun_all, simdats1, fits)


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
logfun_basic <- function(mi, smax, s50, delta){
  a <- log(smax / (log(0.5) + smax) - 1)
  b <- log(smax / (log(0.95) + smax) - 1)
  beta0 <- -a + s50 * (b - a) / delta
  beta1 <- (a - b) / delta
  logmu <- exp(smax * (1 / ( 1 + exp( - beta0 - beta1 * mi)) -1))
}
smax_test$logmu1 <- logfun_basic(smax_test$po2_prime, smax=5, s50=x50, delta)
smax_test$logmu2 <- logfun_basic(smax_test$po2_prime, smax=20, s50=x50, delta)
smax_test$logmu3 <- logfun_basic(smax_test$po2_prime, smax=50, s50=x50, delta)
smax_test$logmu4 <- logfun_basic(smax_test$po2_prime, smax=100, s50=x50, delta)
smax_test$logmu5 <- logfun_basic(smax_test$po2_prime, smax=500, s50=x50, delta)
smax_test$logmu6 <- logfun_basic(smax_test$po2_prime, smax=1000, s50=x50, delta)
smax_test$logmu7 <- logfun_basic(smax_test$po2_prime, smax=10000, s50=x50, delta)

smax_test <- pivot_longer(smax_test, 2:8)

ggplot(smax_test, aes(x=po2_prime, y=value))+geom_point(aes(group=name, color=name))+
  scale_colour_viridis(discrete=T)+
  xlab("pO2'")+
  ylab("pO2' effect")+
  scale_colour_discrete(type="viridis", "smax value", labels=c("5", "20", "50", "100", "500", "1000", "10000"))+theme(legend.position=c(0.8,0.2))

