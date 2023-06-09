### Install packages ####
install_local <- F
library(devtools)

if (install_local) devtools::install_local("/Users/juliaindivero/Library/CloudStorage/Dropbox/sdmTMB-mi.zip")
if (!install_local) remotes::install_github("pbs-assess/sdmTMB", dependencies = TRUE, ref="newlogistic")
library(sdmTMB)

### load helper functions ####
source("Code/util_funs.R")

### Set ggplot themes ###
theme_set(theme_bw(base_size = 15))
theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### Load Data ####
sci_name <- "Sebastolobus altivelis" #"Eopsetta jordani" #"Anoplopoma fimbria"
spc <- "longspine thornyhead" # "petrale sole"#"longspine thornyhead"  #"sablefish"
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

## Get initial values ##                                                                                                
init_vals <- get_inits()

### Plot data vs. depth and po2 ####
ggplot(data = dat, aes(x = po2, y = -depth, col = log(cpue_kg_km2+1))) +
  geom_point() +
  scale_colour_viridis_c(limits = c(0, 6), oob = scales::squish) +
  xlim(c(0, 16))

### Fit Breakpoint model to po2 ####
start <- init_vals$"Sebastolobus altivelis"$m1$start
#lower<- init_vals$"Sebastolobus altivelis"$m1$lower
#upper <- init_vals$"Sebastolobus altivelis"$m1$upper
m1 <- sdmTMB(cpue_kg_km2 ~ 1+year+breakpt(po2_s)+log_depth_scaled+log_depth_scaled2, 
             data = dat,
             time = NULL,
             reml = F,
             anisotropy = TRUE,
             spatiotemporal = FALSE,
             mesh=mesh,
             family =tweedie(link="log"),
             control = sdmTMBcontrol(
               start = list(b_threshold = start),
               #lower = list(b_threshold = lower), 
               #upper = list(b_threshold = upper),
               newton_loops = 2))

summary(m1)
AIC(m1)

#### Plot fitted relationship ####

##### extract estimates ####
m1_pars <- m1$sd_report$par.fixed
m1_parnames <- names(m1_pars)
b_threshold <- m1_pars[grep("b_threshold", m1_parnames)]
##### plot ####
dat$thresh_pred <- exp(brkptfun(x = dat$po2_s, b_slope = b_threshold[1], b_thresh = b_threshold[2]))
ggplot(dat, aes(x = po2, y = thresh_pred, col = log(cpue_kg_km2 + 1))) +
  geom_jitter(width = 0.1, height = 0.01) +
  scale_colour_viridis_c(limits = c(0, 5), oob = scales::squish)
  xlim(c(0,5))
  

### Fit Eo estimation - po2 prime model ####
#Set starting parameters: 
start <- init_vals$"Sebastolobus altivelis"$m2$start

start <- matrix(c(-1,1,100,0.4))

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
           # lower = list(b_threshold = lower),
              # upper = list(b_threshold = upper),
               newton_loops = 2,
               nlminb_loops=2))

summary(m2)
AIC(m2)

#### Plot fitted relationshop ####
##### extract estimates ####

Eo <- getEo(model = m2)

##### plot ####
dat$po2_prime <- dat$po2 * exp(Eo * dat$invtemp)
dat$Eo_unconstrained <- exp(logfun(dat$po2_prime, model = m2, mi = T))
ggplot(dat, aes(x = po2_prime, y = Eo_unconstrained, col = log(cpue_kg_km2 + 1))) +
  geom_jitter(width = 0.1, height = 0.01) +
  scale_colour_viridis_c(limits = c(0, 6), oob = scales::squish) + 
  xlim(c(0,5))
  


### Fit Eo estimation - po2 prime model (with prior) ####
## Set starting parameters:

start <- init_vals$petralesole$m2a$start
lower <- init_vals$petralesole$m2a$lower
upper <- init_vals$petralesole$m2a$upper

prior <- normal(c(NA, NA, NA, 0.331), c(NA, NA, NA, 0.176))

m2a <- sdmTMB(cpue_kg_km2 ~ -1+year+logistic(mi)+log_depth_scaled+log_depth_scaled2,
             data = dat, 
             time = NULL,
             reml = F,
             anisotropy = TRUE,
             spatiotemporal = FALSE,
             mesh=mesh,
             family =tweedie(link="log"),
             priors=sdmTMBpriors(threshold = prior),
             control = sdmTMBcontrol(
               start = list(b_threshold = start),
            #   upper = list(b_threshold = upper),
             #  lower = list(b_threshold = lower),
               newton_loops = 2))

summary(m2a)
AIC(m2a)

#### Plot fitted relationship ####
##### extract estimates ####
Eo <- getEo(model = m2a)

##### plot ####
## Calculate po2 prime from parameter ##
dat$po2_prime <- dat$po2 * exp(Eo * dat$invtemp)

##### plot ####

dat$Eo_constrained <- exp(logfun(dat$po2_prime, model = m2a, mi = T))
ggplot(dat, aes(x = po2_prime, y = Eo_constrained, col = log(cpue_kg_km2 + 1))) +
  geom_jitter(width = 0.1, height = 0.01) +
  scale_colour_viridis_c(limits = c(0, 6), oob = scales::squish) +
  xlim(c(0,5))
  


### Fit logistic po2 model ####
#Starting values (work for both all depths and depth-constraine)
start <- init_vals$petralesole$m3$start
upper <-  matrix(init_vals$petralesole$m3$upper)
lower <- matrix(init_vals$petralesole$m3$lower)

m3 <- sdmTMB(cpue_kg_km2 ~ -1+year+logistic(po2_s)+log_depth_scaled+log_depth_scaled2,
             data = dat, 
             spatial = "on",
             mesh=mesh,
             anisotropy=T,
             reml=F,
             time=NULL,
             family =tweedie(link="log"),
             control = sdmTMBcontrol(
               start = list(b_threshold=start),
               #lower = list(b_threshold = lower),
               #upper = list(b_threshold = upper)
               )
        )


summary(m3)


#### Plot fitted relationship ####
##### extract estimates ####
dat$logistic_pos <- exp(logfun(dat$po2_s, model = m3, mi = F))
##### plot ####
ggplot(dat, aes(x = po2, y = logistic_pos, col = log(cpue_kg_km2 + 1))) +
  geom_jitter(width = 0.1, height = 0.01) +
  scale_colour_viridis_c(limits = c(0, 6), oob = scales::squish)+
  xlim(c(0,5))

### Fit null model ####
m4 <- sdmTMB(cpue_kg_km2 ~ -1+year+log_depth_scaled+log_depth_scaled2, 
             data = dat, 
             spatial = "on",
             mesh=mesh,
             anisotropy=T,
             reml=F,
             time=NULL,
             family =tweedie(link="log"),
             control = sdmTMBcontrol(
               newton_loops = 2,
               )
             )
summary(m4)
AIC(m4)

### Fit other alternative models ###
## Temperature ##

m5 <- sdmTMB(cpue_kg_km2 ~ -1+year+log_depth_scaled+log_depth_scaled2+temp_s,
             data = dat, 
             spatial = "on",
             mesh=mesh,
             anisotropy=T,
             reml=F,
             time=NULL,
             family =tweedie(link="log"),
             control = sdmTMBcontrol(
               newton_loops = 1,
             )
)

summary(m5)
AIC(m5)

##Oxygen only##
m6 <- sdmTMB(cpue_kg_km2 ~ -1+year+log_depth_scaled+log_depth_scaled2+po2_s,
             data = dat, 
             spatial = "on",
             mesh=mesh,
             anisotropy=T,
             reml=F,
             time=NULL,
             family =tweedie(link="log"),
             control = sdmTMBcontrol(newton_loops = 1,
             )
)

summary(m6)
AIC(m6)

## Temp and o2 ##
m7 <- sdmTMB(cpue_kg_km2 ~ -1+year+log_depth_scaled+log_depth_scaled2+temp_s + po2_s,
             data = dat, 
             spatial = "on",
             mesh=mesh,
             anisotropy=T,
             reml=F,
             time=NULL,
             family =tweedie(link="log"),
             control = sdmTMBcontrol(newton_loops = 1,
             )
)

summary(m7)
AIC(m7)

##Temp and o2 interaction ##
m8 <- sdmTMB(cpue_kg_km2 ~ -1+year+log_depth_scaled+log_depth_scaled2+temp_s * po2_s,
             data = dat, 
             spatial = "on",
             mesh=mesh,
             anisotropy=T,
             reml=F,
             time=NULL,
             family =tweedie(link="log"),
             control = sdmTMBcontrol(newton_loops = 1,
             )
)

summary(m8)
AIC(m8)


### Create dAIC table ###
## Make list of model names##
models <- c("breakpt-pO2", "Eo estimation and logistic po2' (no prior)", "Eo estimation and logistic po2' (prior)","logistic-pO2", "Null", "temp", "po2", "temp+po2", "temp * po2")
## Create table and add AIC for each ##
AIC <- as.data.frame(matrix(NA, ncol = 1, nrow =length(models), dimnames = list(models)))
AIC[1,] <- AIC(m1)
AIC[2,] <- AIC(m2)
AIC[3,] <- AIC(m2a)
AIC[4,] <- AIC(m3)
AIC[5,] <- AIC(m4)
AIC[6,] <- AIC(m5)
AIC[7,] <- AIC(m6)
AIC[8,] <- AIC(m7)
AIC[9,] <- AIC(m8)

## Calculate delta-AIC ##
AIC$dAIC <- abs(min(AIC$V1)-(AIC$V1))

### Plot depth effects  ####
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
# Calculate effect of MI # Tim note: probably better to use the logfun function above to calculate the log effect size from po2prime or po2.  
dat$mi_effect <- pars$"mi-smax" * (1 / (1 + exp(-log(19) * (dat$mi_pred - pars$`mi-s50`) / pars$"mi-delta")) - 1) # this is exponentiated. 
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
dat$po2_effect <- ifelse(dat$po2_s>pars2$`po2_s-breakpt`,pars2$`po2_s-breakpt`*pars2$`po2_s-slope`,dat$`po2_s`*pars2$`po2_s-slope`)

# Plot depth effects compared between model #
ggplot(dat, aes(y=depth_effect_combined, x=depth))+geom_line(size=1.3)+geom_line(dat, mapping=aes(y=depth_effect_combined2, x=depth), color="red", size=1.3)+ylab("Depth Effect Combined")
ggplot(subset(dat, depth>450), aes(y=depth_effect_combined, x=depth))+geom_point()+geom_point(subset(dat, depth>450), mapping=aes(y=depth_effect_combined2, x=depth), color="red")+ylab("Depth Effect Combined")

# Plot metabolic index vs depth #
ggplot(dat, aes(y=depth, x=mi_pred))+geom_point(aes(color=log(cpue_kg_km2)))+xlab("pO2'")+theme(legend.position=c(0.8, 0.8))

#Plot oxygen vs depth
ggplot(dat, aes(y=-depth, x=po2))+geom_point(aes(color=log(cpue_kg_km2)))+xlab("pO2")+theme(legend.position=c(0.7, 0.3))+scale_colour_viridis_c()

#Save models
save(m1,m2,m2a,m3,m4,m5,m6,m7,m8, file="models_petralesole.RDS")

#Calculate percent below
below <- subset(dat, dat$po2_prime < 1.45 & dat$cpue_kg_km2>0)
below2 <- subset(dat,  dat$cpue_kg_km2>0)
below <- subset(dat, po2_s < -1.09 & dat$cpue_kg_km2>0)
