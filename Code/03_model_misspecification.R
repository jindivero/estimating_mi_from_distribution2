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
library(ggpubr)

### Set ggplot themes ###
theme_set(theme_bw(base_size = 35))
theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### load helper functions ####
source("Code/util_funs.R")
source("Code/sim_funs.R")

#Load data and models if needed
use_previous <- T
if(use_previous){
  simdat <- readRDS("~/Dropbox/GitHub/estimating_mi_from_distribution2/Model Outputs/data_sims_usual.rds")
  simdat2 <- readRDS("~/Dropbox/GitHub/estimating_mi_from_distribution2/Model Outputs/data_sims_weird.rds")
}

# Set parameter values for data generation and model fitting #
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

## Load data to make mesh ##
sci_name <- "Anoplopoma fimbria" 
spc <- "sablefish" 
dat.by.size <- length_expand(sci_name)
dat <- load_data(spc = spc, dat.by.size = dat.by.size)
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
## Make mesh ##
mesh <- make_mesh(dat, xy_cols = c("X", "Y"), n_knots=250)
#Set starting values
start <- matrix(0, ncol = 1, nrow = 4)
start[1,1] <- s50
start[2,1] <- delta
start[3,1] <- smax
start[4,1] <- Eo

start_unusual <- start
start_unusual[4,1] <- Eo2

### Fit mis-specified depth model ###
fits2_mis <- lapply(simdat, run_sdmTMB_prior_misspecified, 
                    start=start, mesh=mesh)
fits4_mis <- lapply(simdat2, run_sdmTMB_prior_misspecified, 
                    start=start_unusual, mesh=mesh)

save <- T
if(save){
  save(fits2_mis, fits4_mis, file="model_fits_mis.Rdata")
}
## Set how many data sets to produce ##
n <- 250

pars_names <- c("log_depth_scaled", "log_depth_scaled2", "mi-delta", "mi-s50", "mi-smax", "range", "sigma_O", "phi", "tweedie_p", "mi-Eo", "as.factor(year)2010","as.factor(year)2011","as.factor(year)2012", "as.factor(year)2013", "as.factor(year)2014", "as.factor(year)2015")
true_pars <- data.frame(term=c("log_depth_scaled", "log_depth_scaled2", "mi-delta", "mi-s50", "mi-smax", "range", "sigma_O", "phi", "tweedie_p", "mi-Eo", "as.factor(year)2010","as.factor(year)2011","as.factor(year)2012", "as.factor(year)2013", "as.factor(year)2014", "as.factor(year)2015"), 
                        estimate=c(beta1, beta2, delta, s50, smax, range, sigma_O, phi, p, Eo, b_years))
true_pars2 <- data.frame(term=c("log_depth_scaled", "log_depth_scaled2", "mi-delta", "mi-s50", "mi-smax", "range", "sigma_O", "phi", "tweedie_p", "mi-Eo", "as.factor(year)2010","as.factor(year)2011","as.factor(year)2012", "as.factor(year)2013", "as.factor(year)2014", "as.factor(year)2015"), 
                         estimate=c(beta1, beta2, delta, s50, smax, range, sigma_O, phi, p, Eo2, b_years))

## Apply to model fits ##
## Parameter estimates ##
pars1 <- lapply(fits2_mis, extract_pars)
pars2 <- lapply(fits4_mis, extract_pars)

pars1 <- clean_pars(pars1, fits=fits2_mis)
pars2 <- clean_pars(pars2, fits=fits4_mis)

#Add column to label model
pars1$model <- "Typical Case, Misspecified"
pars1$data <- "Typical Case"
pars1$analysis <- "Misspecified"
pars2$model <- "Unusual Case, Misspecified"
pars2$data <- "Unusual Case"
pars2$analysis <- "Misspecified"

#Merge into one and combine
pars <- rbind(pars1,pars2)

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
                        pars$term=="mi-Eo"& pars$model=="Typical Case, Misspecified"~pars$estimate-Eo,
                       # pars$term=="mi-Eo"& pars$model=="Typical Case, Prior Constrained"~pars$estimate-Eo,
                        pars$term=="mi-Eo"& pars$model=="Unusual Case, Misspecified"~pars$estimate-Eo2,
                      #  pars$term=="mi-Eo"& pars$model=="Unusual Case, Prior Constrained"~pars$estimate-Eo2,
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

## Root mean square error (accuracy) ##
# Calculate error #
pars$error <- case_when(pars$term=="mi-s50"~pars$estimate-s50,
                        pars$term=="mi-delta"~pars$estimate-delta,
                        pars$term=="mi-smax"~pars$estimate-smax,
                        pars$term=="mi-Eo"& pars$model=="Typical Case, Misspecified"~pars$estimate-Eo,
                       # pars$term=="mi-Eo"& pars$model=="Typical Case, Prior Constrained"~pars$estimate-Eo,
                        pars$term=="mi-Eo"& pars$model=="Unusual Case, Misspecified"~pars$estimate-Eo2,
                       # pars$term=="mi-Eo"& pars$model=="Unusual Case, Prior Constrained"~pars$estimate-Eo2,
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

Eo_values <- as.data.frame(matrix(nrow=2))
Eo_values$V1 <- NULL
Eo_values$data <- c("Typical Case", "Unusual Case")
Eo_values$analysis <- c( "Misspecified", "Misspecified")
Eo_values$model <- c("Typical Case, Misspecified", "Unusual Case, Misspecified")
MLE_avg <- aggregate(estimate~model, subset(pars, term=="mi-Eo"), FUN=mean)
Eo_values$MLE_avg <- MLE_avg$estimate
Eo_values$true <- c(Eo, Eo2)

s50_values <- as.data.frame(matrix(nrow=2))
s50_values$V1 <- NULL
s50_values$data <- c("Typical Case", "Unusual Case")
s50_values$analysis <- c( "Misspecified", "Misspecified")
s50_values$model <- c("Typical Case, Misspecified", "Unusual Case, Misspecified")
MLE_avg <- aggregate(estimate~model, subset(pars, term=="mi-s50"), FUN=mean)
s50_values$MLE_avg <- MLE_avg$estimate
s50_values$true <- c(s50, s50)

## Make a table ##
rmse$"average" <- avg$estimate
par_performance <- rmse
par_performance$sd <- sd_avg$estimate
colnames(par_performance) <- c("Parameter", "Model", "RMSE", "Average", "Precision")

par_performance$Bias <- case_when(par_performance$Parameter=="mi-s50"~par_performance$Average-s50,
                                  par_performance$Parameter=="mi-delta"~par_performance$Average-delta,
                                  par_performance$Parameter=="mi-smax"~par_performance$Average-smax,
                                  par_performance$Parameter=="mi-Eo"& par_performance$Model=="Typical Case, Misspecified"~par_performance$Average-Eo,
                                 # par_performance$Parameter=="mi-Eo"& par_performance$Model=="Typical Case, Prior Constrained"~par_performance$Average-Eo,
                                  par_performance$Parameter=="mi-Eo"& par_performance$Model=="Unusual Case, Misspecified"~par_performance$Average-Eo2,
                                  #par_performance$Parameter=="mi-Eo"& par_performance$Model=="Unusual Case, Prior Constrained"~par_performance$Average-Eo2,
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
#Density plot just Eo #
ggplot(subset(pars, pars$term=="mi-Eo"), aes(x=estimate)) +
  geom_density(fill="lightblue", adjust = 1.5) +
  geom_vline(data = Eo_values, aes(xintercept = MLE_avg),linetype="dashed", size=1.2, color="darkorange", show.legend=T)+
  geom_vline(data = Eo_values, aes(xintercept = true),linetype="dashed", size=1.2)+
  facet_grid(analysis~data)+
  xlab("Eo estimate") + 
  theme(strip.text = element_text(size = 14))

#Add prior distribution
ggplot(subset(pars, pars$term=="mi-Eo"), aes(x=estimate)) +
  geom_density(fill="lightblue", adjust = 1.5) +
  geom_vline(data = Eo_values, aes(xintercept = MLE_avg),linetype="dashed", size=1.2, color="darkorange", show.legend=T)+
  geom_vline(data = Eo_values, aes(xintercept = true),linetype="dashed", size=1.2)+
  facet_grid(analysis~data)+
  xlab("Eo estimate") + 
  theme(strip.text = element_text(size = 14))+
  stat_function(fun = dnorm, n = n, args = list(mean = 0.3477, sd = 0.1455), linetype="dashed", geom="area", alpha=0.2)+
  stat_function(fun = dnorm, n = n, args = list(mean = 0.3477, sd = 0.1455), linetype="dashed")+
  ylab("Density of data simulations")

# Density plot s50 #
ggplot(subset(pars, pars$term=="mi-s50"), aes(x=estimate)) +
  geom_density(fill="lightblue", adjust = 1.5) +
  geom_vline(data = s50_values, aes(xintercept = MLE_avg),linetype="dashed", size=1.2, color="darkorange", show.legend=T)+
  geom_vline(data = s50_values, aes(xintercept = true),linetype="dashed", size=1.2)+
  facet_grid(analysis~data)+
  xlab("s50 Maximum Likelihood Estimate")+
  ylab("Density of iterations")

#### Comparing po2' effect from parameter estimates ####
#calculate the po2' estimated from the Eo, use that to calculate f(po2')
simdats1 <- mapply(FUN=calculate_po2_prime, simdat,fits2_mis, SIMPLIFY=F)
simdats2 <- mapply(FUN=calculate_po2_prime, simdat2,fits4_mis, SIMPLIFY=F)

#plot true po2 vs exp(f(po2'estimated) (either as points or lines)
simdats1 <- mapply(FUN=logfun2, simdats1,"po2_prime", fits2_mis, SIMPLIFY=F)
simdats2 <- mapply(FUN=logfun2, simdats2,"po2_prime", fits4_mis, SIMPLIFY=F)

#Calculate true po2' effect
dat <- simdats1[[3]]
true_effect <- as.data.frame(logfun_basic(dat$mi_usual, smax, s50, delta))
true_effect2<- as.data.frame(logfun_basic(dat$mi_weird, smax, s50, delta))
colnames(true_effect) <- "mi_effect"
colnames(true_effect2) <- "mi_effect"
true_effect$mi <- dat$mi_usual
true_effect2$mi <- dat$mi_weird

##Remove datasets with NA and then bind together, adding a column with the data frame number ##
simdats1 <- keep(.x=simdats1, .p=is.data.frame)
simdats1 <- bind_rows(simdats1, .id="df")
simdats2 <- keep(.x=simdats2, .p=is.data.frame)
simdats2 <- bind_rows(simdats2, .id="df")
simdats1$title <- "Typical Case"
simdats2$title <- "Unusual Case"
simdats2$side <- "Misspecified"

q1 <- ggplot(simdats1, aes(mi_usual, logmu, colour=Eo)) +
  geom_point(size=0.1)+
  scale_colour_viridis()+
  xlab("True pO2'")+
  ylab("Estimated pO2' effect")+
  theme(legend.position=c(0.9,0.2))+
  geom_line(data=true_effect, aes(x=mi, y=mi_effect), color="black", linetype="dashed", size=2)+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())+
  facet_wrap("title")
q2 <- ggplot(simdats2, aes(mi_weird, logmu, colour=Eo)) +
  geom_point(size=0.1)+
  scale_colour_viridis()+
  xlab("True pO2'")+
  ylab("Estimated pO2' effect")+
  theme(legend.position=c(0.9,0.2))+
  geom_line(data=true_effect, aes(x=mi, y=mi_effect), color="black", linetype="dashed", size=2)+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())+
  facet_grid(side~title)


figure2 <- ggarrange(q1,q2,
                     common.legend = TRUE, legend = "right")

annotate_figure(figure2, left=text_grob(expression(paste("pO"[2], "'", " Effect")), size=30, rot=90),bottom=text_grob(expression(paste("pO"[2], "'", " True Value")), size=30))

###Calculate RMSE ####
simdats1$true <- true_effect$mi_effect
simdats2$true <- true_effect2$mi_effect

simdats1$error <- simdats1$logmu-simdats1$true
simdats1$error2 <- simdats1$error^2
rmse <- aggregate(error2 ~ df, simdats1, FUN=sum)
rmse2 <- aggregate(error2 ~ df, simdats1, FUN=length)
rmse$n <- rmse2$error2
rmse$rmse <- sqrt(rmse$error2/rmse$n)
rmse$n <- NULL
rmse$error2 <- NULL
rmse1 <- mean(rmse$rmse)

simdats2$error <- simdats2$logmu-simdats2$true
simdats2$error2 <- simdats2$error^2
rmse <- aggregate(error2 ~ df, simdats2, FUN=sum)
rmse2 <- aggregate(error2 ~ df, simdats2, FUN=length)
rmse$n <- rmse2$error2
rmse$rmse <- sqrt(rmse$error2/rmse$n)
rmse$n <- NULL
rmse$error2 <- NULL
rmse2a <- mean(rmse$rmse)


