### Install packages ####
remotes::install_github("pbs-assess/sdmTMB", dependencies = TRUE, ref="mi")
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
#setwd("~/Dropbox/Mac/Documents/GitHub/thresholds_mi_distribution/thresholds_mi_distribution")
dat <- readRDS("data/data_sablefish2.rds")

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
     xlim = c(0,5),)


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
     xlim = c(0,5),
     ylab = "po2 marginal effect")

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
