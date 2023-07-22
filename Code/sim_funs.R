### Functions for simulating and fitting data #####

simulate_fish<- function(dat,mesh, s50, delta, smax, Eo, modelpars) {
  # extract model pars
  parnames <- names(modelpars)
  for (i in 1:length(parnames)) eval(parse(text = paste0(parnames[i],"<-", modelpars[i])))
  
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
                         threshold_coefs=c(s50, delta, smax, Eo),
                         B=c(b_years, beta1, beta2),
                         seed=seed)
  dat$sim <- sim$observed
  return(dat)
}

##Function to run model and return list of model outputs
run_sdmTMB_noprior <- function(simdat, start, mesh) {
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
run_sdmTMB_prior <- function(simdat, start, mesh) {
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

## Below simulated true threshold (=s95) ##
s95 <- s50+delta
sum_nll_thresh <- function(nll, column_nll, mi){
  if(is.data.frame(nll)){
    nll <- filter(nll, mi <s50) 
    nlls <- as.numeric(nll[ , column_nll])
    sum_nll <- sum(nlls)
  }
  if(!is.data.frame(nll)){
    sum_nll <- NA
  }
  return(sum_nll)
}

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

logfun_basic <- function(mi, smax, s50, delta){
  a <- log(smax / (log(0.5) + smax) - 1)
  b <- log(smax / (log(0.95) + smax) - 1)
  beta0 <- -a + s50 * (b - a) / delta
  beta1 <- (a - b) / delta
  logmu <- exp(smax * (1 / ( 1 + exp( - beta0 - beta1 * mi)) -1))
}

## Number of zero observations below threshold for each data simulation ##
count_below_zero <- function(dat, threshold) {
  count <- subset(dat, dat$sim==0 & dat$mi_weird < threshold)
  count <- nrow(count)
  return(count)
}
