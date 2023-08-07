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

## Calculate logistic threshold effect ##
logfun2 <- function(smax, s50, delta, mi) {
  a <- log(smax / (log(0.5) + smax) - 1)
  b <- log(smax / (log(0.95) + smax) - 1)
  beta0 <- -a + s50 * (b - a) / delta
  beta1 <- (a - b) / delta
  logmu <- smax * (1 / ( 1 + exp( - beta0 - beta1 * mi)) -1)
}

## Simulate fish ##
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
  dat$po2_prime <- dat$po2 * exp(Eo * dat$invtemp)
  dat$sim <- sim$observed
  a <- log(smax / (log(0.5) + smax) - 1)
  b <- log(smax / (log(0.95) + smax) - 1)
  beta0 <- -a + x50 * (b - a) / delta
  beta1 <- (a - b) / delta
  dat$logmu <- exp(smax * (1 / ( 1 + exp( - beta0 - beta1 * dat$po2_prime)) -1))
  return(dat)
}

## Set parameter values for data generation and model fitting ##
Eo <- 0.291
smax <- 50

# For plotting #
n <-1 #Number of simulations
simdat <- map(seq_len(n), ~simulate_fish(dat = dat,
                                         mesh = mesh,
                                         x50 = 0.5,
                                         delta = 0.5,
                                         smax = smax,
                                         Eo = Eo))
simdat2 <- map(seq_len(n), ~simulate_fish(dat = dat,
                                         mesh = mesh,
                                         x50 = 0.8,
                                         delta =0.8,
                                         smax = smax,
                                         Eo = Eo))
simdat3 <- map(seq_len(n), ~simulate_fish(dat = dat,
                                         mesh = mesh,
                                         x50 = 1,
                                         delta = 1,
                                         smax = smax,
                                         Eo = Eo))
simdat4 <- map(seq_len(n), ~simulate_fish(dat = dat,
                                         mesh = mesh,
                                         x50 = 2,
                                         delta = 2,
                                         smax = smax,
                                         Eo = Eo))
simdat5 <- map(seq_len(n), ~simulate_fish(dat = dat,
                                          mesh = mesh,
                                          x50 = 3,
                                          delta = 3,
                                          smax = smax,
                                          Eo = Eo))
simdat6 <- map(seq_len(n), ~simulate_fish(dat = dat,
                                          mesh = mesh,
                                          x50 = 4,
                                          delta = 4,
                                          smax = smax,
                                          Eo = Eo))

# Plot #
ggplot(simdat[[1]], aes(po2_prime, logmu))+geom_line(color="red")+geom_line(simdat2[[1]], mapping=aes(po2_prime, logmu), color="orange")+geom_line(simdat3[[1]], mapping=aes(po2_prime, logmu), color="yellow")+geom_line(simdat4[[1]], mapping=aes(po2_prime, logmu), color="green")+geom_line(simdat5[[1]], mapping=aes(po2_prime, logmu), color="blue")+geom_line(simdat6[[1]], mapping=aes(po2_prime, logmu), color="purple")


## Simulate a bunch for model fitting ##
n <-100 #Number of simulations
simdat <- map(seq_len(n), ~simulate_fish(dat = dat,
                                         mesh = mesh,
                                         x50 = 0.56,
                                         delta = 0.88,
                                         smax = smax,
                                         Eo = Eo))
simdat2 <- map(seq_len(n), ~simulate_fish(dat = dat,
                                          mesh = mesh,
                                          x50 = 1,
                                          delta =1,
                                          smax = smax,
                                          Eo = Eo))
simdat3 <- map(seq_len(n), ~simulate_fish(dat = dat,
                                          mesh = mesh,
                                          x50 = 2,
                                          delta = 2,
                                          smax = smax,
                                          Eo = Eo))
simdat4 <- map(seq_len(n), ~simulate_fish(dat = dat,
                                          mesh = mesh,
                                          x50 = 5,
                                          delta = 5,
                                          smax = smax,
                                          Eo = Eo))

#Plot of different logistic shapes

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

start <- matrix(0, ncol = 1, nrow = 4)
start[1,1] <- 0.56
start[2,1] <- 0.88
start[3,1] <- smax
start[4,1] <- Eo

start2 <- matrix(0, ncol = 1, nrow = 4)
start2[1,1] <- 1
start2[2,1] <-1
start2[3,1] <- smax
start2[4,1] <- Eo

start3 <- matrix(0, ncol = 1, nrow = 4)
start3[1,1] <- 2
start3[2,1] <-2
start3[3,1] <- smax
start3[4,1] <- Eo

start4 <- matrix(0, ncol = 1, nrow = 4)
start4[1,1] <- 4
start4[2,1] <-4
start4[3,1] <- smax
start4[4,1] <- Eo

fits <- lapply(simdat, run_sdmTMB_1, 
               start=start, mesh=mesh)
fits2 <- lapply(simdat2, run_sdmTMB_1, 
               start=start2, mesh=mesh)
fits3 <- lapply(simdat3, run_sdmTMB_1, 
               start=start3, mesh=mesh)
fits4 <- lapply(simdat4, run_sdmTMB_1, 
               start=start4, mesh=mesh)


#Save models
#save(fits, fits1, fits2, fits3, fits_b, fits_c, file="model_fits.Rdata")
save(fits, fits2, fits3, fits4, file="model_fits.Rdata")

### Evaluate model performance ###
## Set true pars vectors of values and name ##
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
model_names <- c("0.56,0.88", "1,1", "2,2", "4,4")


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
  names <- names(fits[[3]]$sd_report$par.fixed)
  colnames(grad) <- pars_names
  grad$sim <- 1:nrow(grad)
  grad <- pivot_longer(grad,cols=1:11, names_to="term", values_to="gradient")
  grad$big_grad <- ifelse(grad$gradient<0.001, 1,0)
  large_grad <- aggregate(big_grad ~ term, grad, FUN=sum)
  #large_grad$proportion <- large_grad$big_grad/length(fits)
  return(large_grad)
}

### Apply to model fits ###
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
pars$estimate <- ifelse(pars$estimate>100 & pars$term=="mi-smax", 100, pars$estimate)

#Eo and smax
ggplot(subset(pars, pars$term=="mi-Eo"|pars$term=="mi-smax"|pars$term=="mi-s50"|pars$term=="mi-delta"), aes(y=estimate, x=model))+geom_boxplot(aes(group=model, fill=model))+
  facet_wrap("term", scales="free")+  
 # geom_hline(data = subset(true_pars, true_pars2$term=="mi-Eo"|true_pars$term=="mi-smax"|true_pars$term=="mi-s50"|true_pars$term=="mi-delta"), aes(yintercept = estimate),linetype="dashed", size=1.2)+
  #geom_hline(data = subset(true_pars2, true_pars2$term=="mi-Eo"), aes(yintercept = estimate),linetype="dashed", size=1.2)+
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

# Apply to all models
grad1<- lapply(fits, extract_grad)
grad2<- lapply(fits2, extract_grad)
grad3<- lapply(fits3, extract_grad)
grad4<- lapply(fits4, extract_grad)

grads1<- clean_grad(grad1, fits)
grads2<- clean_grad(grad2, fits2)
grads3<- clean_grad(grad3, fits3)
grads4<- clean_grad(grad4, fits4)

# Bind together
grads1$model2 <-grads2$big_grad
grads1$model3 <-grads3$big_grad
grads1$model4 <-grads4$big_grad
colnames(grads1) <-c("metric", model_names[1], model_names[2], model_names[3], model_names[4])
grads1$metric <- paste0("Low gradient ", grads1$metric)

##NA in SD ##
pars$NAs <- ifelse(pars$std.error=="NaN"|is.na(pars$std.error), 0,1)
NAs <- aggregate(pars$NAs~pars$term+pars$model, FUN=sum)
colnames(NAs) <- c("metric", "model", "value")
NAs$metric <- paste0("SE estimated", NAs$metric)
##Flip to wide
NAs <-pivot_wider(NAs, names_from=model, values_from=value)

### Combine into table summarizing all diagnostics ##
diagnostics <- bind_rows(diagnostics, grads1, NAs)

### Parameter performance measures ###
## Root mean square error (accuracy) ##
# Calculate error #
pars$error <- case_when(pars$term=="mi-s50"~pars$estimate-x50,
                        pars$term=="mi-delta"~pars$estimate-delta,
                        pars$term=="mi-smax"~pars$estimate-smax,
                        pars$term=="mi-Eo"~pars$estimate-Eo,
                        pars$term=="log_depth_sc"~pars$estimate-beta1,
                        pars$term=="log_depth_sc2"~pars$estimate-beta2,
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

# Average standard error (precision)--also to compare to SD of prior #
sd_avg <- aggregate(std.error ~ term+model, pars, FUN=mean) 

# Coefficient of variation of estimate #
pars$cv <- pars$std.error/pars$estimate

