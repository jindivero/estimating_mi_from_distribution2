### Install packages ####
#remotes::install_github("pbs-assess/sdmTMB", dependencies = TRUE, ref="mi")


#install.packages("devtools")
#library(devtools)
#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)
#library(INLA)
devtools::install_local("/Users/jindiv/Dropbox/GitHub/sdmTMB-mi.zip")
library(sdmTMB)

library(ggplot2)
library(visreg)
library(ggeffects)
library(dplyr)
library(tidyr)
library(mgcv)
library(here)

brkptfun <- function(x, b_slope, b_thresh) min(0, b_slope *  (x - b_thresh))
logfun <- function(x, s50, delta, smax) smax * ((1 + exp(-log(19) * (x - s50)/(delta)))^(-1) - 1)

### Load Data ####
#here("thresholds_mi_distribution")
setwd("/Users/jindiv/Dropbox/GitHub/estimating_mi_from_distribution2")
dat <- readRDS("data_sablefish2.rds")

### Calculate inverse temp ####
kelvin = 273.15 #To convert to Kelvin
boltz = 0.000086173324 #Boltzman's constant
tref <- 7 #Reference temperature in celsius
dat$invtemp <- (1 / boltz)  * ( 1 / (dat$temp + 273.15) - 1 / (tref + 273.15)) #invtemp #I used z because "psi" was too confusing

### Make mesh ####
mesh <- make_mesh(dat, xy_cols = c("X", "Y"), n_knots = 250)

### Fit Breakpoint model to po2 ####
start <- matrix(0, nrow = 2, ncol = 1)
start[1,1] <- 20
start[2,1] <- -1.1

m1 <- sdmTMB(cpue_kg_km2 ~ -1+year+breakpt(po2_sc)+log_depth_scaled+log_depth_scaled2, 
             data = dat,
             time = NULL,
             reml = F,
             anisotropy = TRUE,
             spatiotemporal = FALSE,
             mesh=mesh,
             family =tweedie(link="log"),
             control = sdmTMBcontrol(
               start = list(b_threshold = start),
               newton_loops = 1))

summary(m1)
AIC(m1)

#### Plot fitted relationship ####

##### extract estimates ####
m1_pars <- m1$sd_report$par.fixed
m1_parnames <- names(m1_pars)
b_threshold <- m1_pars[grep("b_threshold", m1_parnames)]
##### plot ####
plot(dat$po2, exp(sapply(X = dat$po2_sc, FUN = brkptfun, b_slope = b_threshold[1], b_thresh = b_threshold[2])),
     ylab = "po2 marginal effect", 
     main="Breakpoint-pO2",
     xlab="pO2",
     xlim = c(0,5),)

## Identify number below threshold ##
dat$po2_prime <- po2_prime
dat_below2 <- subset(dat, po2 < exp(-1.09))

### Fit Eo estimation - po2 prime model ####
#Set starting parameters: 

start <- matrix(0, ncol = 1, nrow = 4)
start[1, 1] <- 2 #s50
start[2, 1] <- (1) #delta
start[3, 1] <- 10 #smax (ie beta_3)
start[4, 1] <- 0.68 #Eo
m2 <- sdmTMB(cpue_kg_km2 ~ -1+year+logistic(mi)+log_depth_scaled+log_depth_scaled2, 
             data = dat, 
             time = NULL,
             reml = F,
             anisotropy = TRUE,
             spatiotemporal = FALSE,
             mesh=mesh,
             family =tweedie(link="log"),
             control = sdmTMBcontrol(
               start = list(b_threshold = start),
               lower = list(b_threshold = c(-Inf, -Inf, 0, 0)), upper = list(b_threshold = c(Inf, Inf,50, Inf)),
               newton_loops = 1))

summary(m2)
AIC(m2)

#### Plot fitted relationshop ####
##### extract estimates ####

re2 <- m2$tmb_random
parfit <- m2$sd_report
npars <- length(parfit$value)
parnames <- names(parfit$value)

s50 <- parfit$value[grep("s50", parnames)]
delta <- parfit$value[grep("s95", parnames)]
smax <- parfit$value[grep("s_max", parnames)]
Eo <- parfit$value[grep("Eo", parnames)]
betas <-  m2$sd_report$par.fixed[grep("b_j", names( m2$sd_report$par.fixed))]

beta_depth <- betas[(length(betas)-1):length(betas)]

##### plot ####
po2_prime <- dat$po2 * exp(Eo * dat$invtemp)

plot(po2_prime, exp(logfun(po2_prime, s50, delta, smax)),
     xlim = c(0,5),
     main="Eo estimation and logistic pO2'",
     ylab = "po2-prime marginal effect")

## Identify number below threshold ##
dat$po2_prime <- po2_prime
dat_below <- subset(dat, po2_prime<1)

### Fit Eo estimation - po2 prime model (with prior) ####
#Set starting parameters: 

start <- matrix(0, ncol = 1, nrow = 4)
start[1, 1] <- 2 #s50
start[2, 1] <- (1) #delta
start[3, 1] <- 10 #smax (ie beta_3)
start[4, 1] <- 0.68 #Eo
m2a <- sdmTMB(cpue_kg_km2 ~ -1+year+logistic(mi)+log_depth_scaled+log_depth_scaled2, 
             data = dat, 
             time = NULL,
             reml = F,
             anisotropy = TRUE,
             spatiotemporal = FALSE,
             mesh=mesh,
             family =tweedie(link="log"),
             priors=sdmTMBpriors(threshold = normal(c(NA, NA, NA, 0.448), c(NA, NA, NA, 0.15))),
             control = sdmTMBcontrol(
               start = list(b_threshold = start),
               lower = list(b_threshold = c(-Inf, -Inf, 0, 0)), upper = list(b_threshold = c(Inf, Inf,50, Inf)),
               newton_loops = 1))

summary(m2a)
AIC(m2a)

#### Plot fitted relationshop ####
##### extract estimates ####

re2a <- m2a$tmb_random
parfit <- m2a$sd_report
npars <- length(parfit$value)
parnames <- names(parfit$value)

s50 <- parfit$value[grep("s50", parnames)]
delta <- parfit$value[grep("s95", parnames)]
smax <- parfit$value[grep("s_max", parnames)]
Eo <- parfit$value[grep("Eo", parnames)]
betas <-  m2a$sd_report$par.fixed[grep("b_j", names( m2a$sd_report$par.fixed))]

beta_depth <- betas[(length(betas)-1):length(betas)]

##### plot ####
## Calculate po2 prime from parameter ##
po2_prime <- dat$po2 * exp(Eo * dat$invtemp)

plot(po2_prime, exp(logfun(po2_prime, s50, delta, smax)),
     xlim = c(0,5),
     main="Eo estimation and logistic pO2 (with prior)'",
     ylab = "po2-prime marginal effect")

### Fit logistic po2 model ####
start <- matrix(0, ncol = 1, nrow = 3)
start[1, 1] <- -1.5 #s50
start[2, 1] <- log(0.5) # log delta
start[3, 1] <- 40 #smax
m3 <- sdmTMB(cpue_kg_km2 ~ -1+year+logistic(po2_sc)+log_depth_scaled+log_depth_scaled2, 
             data = dat, 
             spatial = "on",
             mesh=mesh,
             anisotropy=T,
             reml=F,
             time=NULL,
             family =tweedie(link="log"),
             control = sdmTMBcontrol(
               start = list(b_threshold=start),
               newton_loops = 2
             ))
summary(m3)


#### Plot fitted relationship ####
##### extract estimates ####
parfit <- m3$sd_report
npars <- length(parfit$value)
parnames <- names(parfit$value)

s50 <- parfit$value[grep("s50", parnames)]
s95 <- (parfit$value[grep("s95", parnames)] )
delta <- s95 - s50
smax <- parfit$value[grep("s_max", parnames)]

##### plot ####
plot(dat$po2, exp(logfun(dat$po2_sc, s50, delta, smax)),
     main="Logistic-pO2",
     xlim = c(0,5),
     ylab = "po2 marginal effect",
     xlab="pO2"
     )

### Fit null model ####
m4 <- sdmTMB(cpue_kg_km2 ~ -1+year+log_depth_scaled+log_depth_scaled2, 
             data = dat, 
             spatial = "on",
             mesh=mesh,
             anisotropy=T,
             reml=F,
             time=NULL,
             family =tweedie(link="log"),
             control = sdmTMBcontrol(newton_loops = 2,
               eval.max = 1000000L,
               iter.max = 1000000L,
               nlminb_loops = 100L
               )
             )
summary(m4)
AIC(m4)

### Create dAIC table ###
## Make list of model names##
models <- c("breakpt-pO2", "Eo estimation and logistic po2' (no prior)", "Eo estimation and logistic po2' (prior)","logistic-pO2", "Null")
## Create table and add AIC for each ##
AIC <- as.data.frame(matrix(NA, ncol = 1, nrow =5, dimnames = list(models)))
AIC[1,] <- AIC(m1)
AIC[2,] <- AIC(m2)
AIC[3,] <- AIC(m2a)
AIC[4,] <- AIC(m3)
AIC[5,] <- AIC(m4)
## Calculate delta-AIC ##
AIC$dAIC <- abs(min(AIC$V1)-(AIC$V1))

### Plot depth effects ###
## Set ggplot theme ##
theme_set(theme_bw(base_size = 30))
theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
## Internal Eo mode ##
# Extract parameters another way for internal Eo model #
par_estimates <- as.data.frame(tidy(m2, conf.int = TRUE, effects="fixed"))
par_estimates_rand <- as.data.frame(tidy(m2, conf.int = TRUE, effects="ran_pars"))
par_estimates <- bind_rows(par_estimates, par_estimates_rand)
pars <- par_estimates[1:2]
pars <- pivot_wider(pars, names_from=term, values_from=estimate)
# Calculate temp-corrected po2 from Eo value #
dat$mi_pred <- dat$po2*exp(pars$"mi-Eo"* dat$invtemp)
# Calculate effect of MI #
dat$mi_effect <- pars$"mi-smax" * (1 / (1 + exp(-log(19) * (dat$mi_pred - pars$`mi-s50`) / pars$"mi-delta")) - 1)
# Calculate combined depth effect #
dat$depth_effect_combined <- (dat$log_depth_scaled*pars$log_depth_scaled)+(dat$log_depth_scaled2*pars$log_depth_scaled2)

## For breakpt(o2) model ##
par_estimates2 <- as.data.frame(tidy(m1, conf.int = TRUE, effects="fixed"))
par_estimates_rand2 <- as.data.frame(tidy(m1, conf.int = TRUE, effects="ran_pars"))
par_estimates2 <- bind_rows(par_estimates2, par_estimates_rand2)
pars2 <- par_estimates2[1:2]
pars2 <- pivot_wider(pars2, names_from=term, values_from=estimate)
# alculate depth effects #
dat$depth_effect_combined2 <- (dat$log_depth_scaled*pars2$log_depth_scaled)+(dat$log_depth_scaled2*pars2$log_depth_scaled2)
# Calculate po2 effect #
dat$po2_effect <- ifelse(dat$po2_sc>pars2$`po2_sc-breakpt`,pars2$`po2_sc-breakpt`*pars2$`po2_sc-slope`,dat$`po2_sc`*pars2$`po2_sc-slope`)

# Plot depth effects compared between model #
ggplot(dat, aes(y=depth_effect_combined, x=depth))+geom_line(size=1.3)+geom_line(dat, mapping=aes(y=depth_effect_combined2, x=depth), color="red", size=1.3)+ylab("Depth Effect Combined")
ggplot(subset(dat, depth>450), aes(y=depth_effect_combined, x=depth))+geom_point()+geom_point(subset(dat, depth>450), mapping=aes(y=depth_effect_combined2, x=depth), color="red")+ylab("Depth Effect Combined")

# Plot metabolic index vs depth #
ggplot(dat, aes(y=depth, x=mi_pred))+geom_point(aes(color=log(cpue_kg_km2)))+xlab("pO2'")+theme(legend.position=c(0.8, 0.8))

