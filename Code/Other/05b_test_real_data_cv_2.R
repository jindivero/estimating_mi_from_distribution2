### Install packages ####
install_local <- F
library(devtools)

if (install_local) devtools::install_local("/Users/juliaindivero/Library/CloudStorage/Dropbox/sdmTMB-newlogistic.zip")
if (!install_local) remotes::install_github("pbs-assess/sdmTMB", dependencies = TRUE, ref="newlogistic")
library(sdmTMB)

### load helper functions ####
source("Code/util_funs.R")

### Set ggplot themes ###
theme_set(theme_bw(base_size = 25))
theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### Load Data ####
sci_name <-"Anoplopoma fimbria" # "Sebastolobus altivelis"   #"Eopsetta jordani" 
spc <- "sablefish" #"longspine thornyhead" # "petrale sole" #
dat.by.size <- length_expand(sci_name)
dat <- load_data(spc = spc, dat.by.size = dat.by.size)

#Constrain depth?
constrain_depth <- F
if(constrain_depth) dat <- subset(dat, depth<600)

#Constrain depth for petrale
constrain_depth <- F
if(constrain_depth) dat <- subset(dat, depth<500)

## Scale, rename, and calculate variables ##
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

# Remove outliers = catch > 10 sd above the mean
dat$cpue_s <- scale(dat$cpue_kg_km2)
dat <- dplyr::filter(dat, cpue_s <=20)
### Calculate inverse temp ####
kelvin = 273.15 #To convert to Kelvin
boltz = 0.000086173324 #Boltzman's constant
tref <- 7 #Reference temperature in celsius
dat$invtemp <- (1 / boltz)  * ( 1 / (dat$temp + 273.15) - 1 / (tref + 273.15)) #invtemp 

### Make mesh ####
mesh <- make_mesh(dat, xy_cols = c("X", "Y"), n_knots = 250)

### Fit Breakpoint model to po2 ####
seed <- sample(1929)
set.seed(seed)
k_folds <- 5
dat$fold_ids <- sample(seq_len(k_folds), nrow(dat), replace = TRUE)
future::plan(future::multisession)
start <- matrix(c(20, 1.1), ncol = 1L)

#Breakpoint model works:
m1 <- sdmTMB_cv(cpue_kg_km2 ~ -1+year+breakpt(po2_s)+log_depth_scaled+log_depth_scaled2, 
             data = dat,
             time = NULL,
             reml = F,
             anisotropy = TRUE,
             spatiotemporal = FALSE,
             mesh=mesh,
             family =tweedie(link="log"),
            future_globals = c("start"),
             control = sdmTMBcontrol(
               start = list(b_threshold = start)))

m1$models
m1$sum_loglik
m1$fold_loglik

# Logistic model (on po2) works:
#Starting values (work for both all depths and depth-constraine)
start <- matrix(c(-1.65, 0.5,200), ncol = 1L)
m3 <- sdmTMB_cv(cpue_kg_km2 ~ -1+year+logistic(po2_s)+log_depth_scaled+log_depth_scaled2,
                data = dat, 
                spatial = "on",
                mesh=mesh,
                anisotropy=T,
                reml=F,
                time=NULL,
                family =tweedie(link="log"),
                future_globals = c("start"),
                control = sdmTMBcontrol(
                  start = list(b_threshold=start),
                )
)


m3$models
m3$sum_loglik
m3$fold_loglik

#Logistic metabolic index model gets an error:
#Set starting parameters: 
start <- matrix(c(1,1,1500,1), ncol = 1L)
#This causes an error:
m2 <- sdmTMB_cv(cpue_kg_km2 ~ -1+year+logistic(mi)+log_depth_scaled+log_depth_scaled2, 
             data = dat, 
             time = NULL,
             reml = F,
             anisotropy = TRUE,
             spatiotemporal = FALSE,
             mesh=mesh,
             family =tweedie(link="log"),
           future_globals = c("start"),
             control = sdmTMBcontrol(
              start = list(b_threshold = start)))
m2$models
m2$sum_loglik
m2$fold_loglik

#And same with a prior:
prior <- normal(c(NA, NA, NA, 0.331), c(NA, NA, NA, 0.176))
m2a <- sdmTMB_cv(cpue_kg_km2 ~ -1+year+logistic(mi)+log_depth_scaled+log_depth_scaled2,
             data = dat, 
             time = NULL,
             reml = F,
             anisotropy = TRUE,
             spatiotemporal = FALSE,
             mesh=mesh,
             family =tweedie(link="log"),
             priors=sdmTMBpriors(threshold = prior),
             future_globals = c("prior", "start"),
             control = sdmTMBcontrol(
              start = list(b_threshold = start),
               newton_loops = 2))

m2a$models
m2a$sum_loglik
m2a$fold_loglik