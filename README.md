
## Purpose
This vignette is to provide guidance on using an optional estimation feature +logistic(mi) added to model equation using  sdmTMB() developed and used in Indivero et al. (in prep). This feature is specifically designed for using spatial data on oxygen, temperature, and marine species density to estimate both:

1. The metabolic index, with an estimated parameter Eo
2. A threshold function of the effect of the metabolic index on fish density (with estimated parameters $s50$, $s95$, and scaling parameter $\psi$ (called "smax" in sdmTMB)) 
all within the estimation model of a generalized linear mixed model, i.e. simultaneously estimating over covariates, spatial random variation, etc. 

For further details on the equations of the metabolic index and the form of the threshold function, see Indivero et al. (in prep).

Here, we will use as an example the data for sablefish (*Anoplopoma fimbria*) from Indivero et al. This data is from 2010-2015 of the National Oceanic and Atmospheric Administration West Coast Bottom Trawl Survey ([Keller et al. 2017])(https://repository.library.noaa.gov/view/noaa/14179), spanning Washington, Oregon, and California of the California Current. 

## Downloading the branch
Currently, this feature is included in a separate branch of sdmTMB, called "sdmTMB-newlogistic". To install the package in R:

```{r install package, echo=TRUE}
remotes::install_github("pbs-assess/sdmTMB", dependencies = TRUE,  ref="newlogistic")
library(sdmTMB)
library(ggplot2)
library(dplyr)
```

## Data Inputs
To download the example used here, load the rds
```{r data, echo=TRUE}
dat <- readRDS("example_data.rds")
```

This model feature requires three columns in the dataframe provided to sdmTMB (in addition to latitude and longitude, if a spatial model): 

1. **species density** (e.g. catch-per-unit-effort, or other response variable) (column can be named anything)
2. **oxygen** (must be in kiloPascals (kPa), and column must be named *po2*
3. **temperature** (must be in inverse temperature and Kelvin, see below) and column must be named *invtemp*. 

Any additional time variable or environmental covariates can be included in the data at the same coordinates as the catch, temperature, and oxygen data.

###Oxygen: 
Oxygen needs to be in saturation of kPa. For example and functions of how to convert from mL per L to kPa, see [util_funs.R](https://github.com/jindivero/estimating_mi_from_distribution2/blob/f696e18bb8bf5b20d2127924eff033822d6fc2e4/Code/util_funs.R).

### Inverse temperature
Temperature needs to be inverse temperature, in Kelvin, and in reference to a set temperature and corrected with Boltzmann's constant in an Arrhenius equation. Below is an example of converting temperature in Celsius to the correct inverse temperature format.
```{r inverse temperature, echo=TRUE}
# Constants
  kelvin = 273.15
  boltz = 0.000086173324
  tref <- 12
# Calculate inverse temp
  dat$invtemp <- (1 / boltz)  * ( 1 / (dat$temp + 273.15) - 1 / (tref + 273.15))
```

### Species density
In this example, species density is catch-per-unit-effort, in biomass kg km^-2^. We have not tested with other response variables, though in theory this model could work with other measurements of density (e.g. counts) or other response variables (e.g. size, condition, etc.).

### Other covariates
The metabolic index threshold feature only requires oxygen, temperature, and species density. However, spatial data, time data (e.g. year), and other environmental covariates can also be included. In this example, we also include year and depth as covariates. 

```{r view data, echo=TRUE}
head(dat)
```

## Model Estimation
The Metabolic Index estimation is implemented using the standard sdmTMB() or sdmTMB_cv() function in sdmTMB, and adding  +logistic(mi) to the model equation. can be added in combination with most other model structure options in sdmTMB. It can be used with or without temporal, spatial, and spatio-temporal variation, and with other fixed effects in the model equation. However, it cannot currently be used with a threshold function or a smoother on an additional covariate. 

###Make mesh
```{r mesh, echo=TRUE}
mesh <- make_mesh(dat, xy_cols = c("X", "Y"), n_knots=250)
```

###Fit model
In this example, we will fit a model of fish density to the metabolic index in addition to year, and a quadratic effect of depth. We will include spatial random variation, but not spatio-temporal, and we will use a Tweedie distribution.

Additionally, to help with estimation, we are going to use starting values and use 2 newton loops. We will not use lower and upper bounds here, but the code below shows where these can be added. These can be specified in a matrix that is fed to the sdmTMB() function (e.g. as in the start matrix in the example below), or added in a list in the the function itself (as i the lower and upper bounds in the example below).

Note that there are four estimated parameters of the metabolic index threshold feature (in order: s50, s95, smax, Eo), so if you want to use any values for starting or bounds, it needs to be a list of four, with NAs for any parameter you do not want to include a value for.

```{r fit model, echo=TRUE}
# Make starting values
  start <- matrix(0, ncol = 1, nrow = 4)
  start[1,1] <- 2 #s50
  start[2,1] <- 2 #s95
  start[3,1] <- 30 #smax
  start[4,1] <- 0.3 #Eo
  
# Fit model
  fit <- sdmTMB(cpue_kg_km2 ~ -1+as.factor(year)+logistic(mi)+log_depth_scaled+log_depth_scaled2, 
                   data = dat, 
                   spatial = "on",
                   spatiotemporal="off",
                   mesh=mesh,
                   family =tweedie(link="log"),
                   control = sdmTMBcontrol(
                     start = list(b_threshold = start),
                     lower = list(b_threshold = c(-Inf, -Inf, -Inf, -Inf)), 
                     upper = list(b_threshold = c(Inf, Inf, Inf, Inf)),
                     newton_loops = 2))
```

### Add a prior
A prior can be added using priors=sdmTMBpriors(threshold=). Here, we will add a prior on Eo only.

```{r prior, echo=TRUE}
fit <- sdmTMB(cpue_kg_km2 ~ -1+as.factor(year)+logistic(mi)+log_depth_scaled+log_depth_scaled2, 
                   data = dat, 
                   spatial = "on",
                   spatiotemporal="off",
                   mesh=mesh,
                   family =tweedie(link="log"),
                   control = sdmTMBcontrol(
                     start = list(b_threshold = start),
                     lower = list(b_threshold = c(-Inf, -Inf, -Inf, -Inf)), 
                     upper = list(b_threshold = c(Inf, Inf, Inf, Inf)),
                     newton_loops = 2),
                   priors=sdmTMBpriors(threshold = normal(c(NA, NA, NA, 0.3477), c(NA, NA, NA, 0.1455))))
```

## Interpreting Model Output
Standard wrapper functions can be used to evaluate the model fit. The parameters estimated from the metabolic index and threshold function start with the prefix "mi-" in the parameter list.

NOTE: The "s95" parameter value reported is actually a delta value (difference between s95 and s50). To calculate the value of s95 itself (i.e. the value of the metabolic index at which the logged response on fish is a 95% reduction on density), add the reported values of s50+s95.

```{r interpret, echo=TRUE}
#View parameter estimates
  summary(fit)

#Get AIC of model fit
  AIC(fit)

```

### Marginal effects
We can pull out the parameter estimates for the metabolic index threshold feature to calculate the metabolic index from observed temperature and oxygen and the marginal effect of the metabolic index on fish density.

To evaluate how the response is impacted across the range of estimated Eo values (e.g. the 95% confidence interval of the maximum likelihood estimate), the lower and upper interval values can be used in the equations.
```{r marginal effect, echo=TRUE}
## Extract parameters
  #Using tidy
    par_estimates <- as.data.frame(tidy(fit, conf.int = TRUE, effects="fixed"))
    par_estimates_rand <- as.data.frame(tidy(fit, conf.int = TRUE, effects="ran_pars"))
    par_estimates <- bind_rows(par_estimates, par_estimates_rand)
    
  #Directly from model object
    parfit <- fit$sd_report
    npars <- length(parfit$value)
    parnames <- names(parfit$value)
    s50 <- parfit$value[grep("s50", parnames)]
    delta <- parfit$value[grep("s95", parnames)]
    smax <- parfit$value[grep("s_max", parnames)]
    Eo <- parfit$value[grep("Eo", parnames)]

## Calculate metabolic index
  dat$metabolic_index = dat$po2*exp(Eo* dat$invtemp)
  
## Calculate marginal effect of threshold function
# Re-arrange parameters
  a <- log(smax / (log(0.5) + smax) - 1)
  b <- log(smax / (log(0.95) + smax) - 1)
  beta0 <- -a + s50 * (b - a) / delta
  beta1 <- (a - b) / delta
# Calculate
  dat$effect <- exp(smax * (1 / ( 1 + exp( - beta0 - beta1 * dat$metabolic_index)) -1))

## Plot
ggplot(dat, aes(x=metabolic_index, y=effect))+geom_point()

```

### Threshold oxygen level across temperatures
The parameter estimates can also be used to calculate the values of oxygen across a range of temperatures that correspond to the limiting metabolic index threshold (either s50 or s95).

```{r threshold, echo=TRUE}
#For the s50 threshold
dat$po2_s50 <- s50/exp(Eo*dat$invtemp)
#For the s95 threshold
s95 <- s50+delta
dat$po2_s95 <- s95/exp(Eo*dat$invtemp)

ggplot(dat, aes(x=invtemp, y=po2_s50))+geom_point()

#Or with temperature in Celsius
ggplot(dat, aes(x=temp, y=po2_s50))+geom_point()
```
