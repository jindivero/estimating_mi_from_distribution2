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

run_sdmTMB_prior_misspecified <- function(simdat, start, mesh) {
  m2 <- try(sdmTMB(sim ~ -1+as.factor(year)+logistic(mi)+log_depth_scaled, 
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

run_sdmTMB_misspecified_noprior <- function(simdat, start, mesh) {
  m2 <- try(sdmTMB(sim ~ -1+as.factor(year)+logistic(mi)+log_depth_scaled, 
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

### Run alternative models for each data simulation ###
run_alt_models <- function(dat, start, mesh) {
  m1 <- try(sdmTMB(sim ~ 1+year+breakpt(po2_s)+log_depth_scaled+log_depth_scaled2, 
               data = dat,
               time = NULL,
               reml = F,
               anisotropy = TRUE,
               spatiotemporal = FALSE,
               mesh=mesh,
               family =tweedie(link="log"),
              # control = sdmTMBcontrol(
                 #start = list(b_threshold = start),
                 #lower = list(b_threshold = lower), 
                 #upper = list(b_threshold = upper),
                 #newton_loops = 2
  ))
  
  m2 <- try(sdmTMB(sim ~ -1+year+logistic(mi)+log_depth_scaled+log_depth_scaled2, 
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
                 # upper = list(b_threshold = upper),
                 newton_loops = 2,
                 nlminb_loops=2)))
  
  m2a <- try(sdmTMB(sim ~ -1+year+logistic(mi)+log_depth_scaled+log_depth_scaled2,
                data = dat, 
                time = NULL,
                reml = F,
                anisotropy = TRUE,
                spatiotemporal = FALSE,
                mesh=mesh,
                family =tweedie(link="log"),
                control = sdmTMBcontrol(
                  start = list(b_threshold = start),
                  #   upper = list(b_threshold = upper),
                  #  lower = list(b_threshold = lower),
                  newton_loops = 2),
               priors=sdmTMBpriors(threshold = normal(c(NA, NA, NA, 0.3477), c(NA, NA, NA, 0.1455)))))
  
  m3 <- try(sdmTMB(sim ~ -1+year+logistic(po2_s)+log_depth_scaled+log_depth_scaled2,
               data = dat, 
               spatial = "on",
               mesh=mesh,
               anisotropy=T,
               reml=F,
               time=NULL,
               family =tweedie(link="log"),
               control = sdmTMBcontrol(
                 start = list(b_threshold=start),
                 #  lower = list(b_threshold = lower),
                 #  upper = list(b_threshold = upper)
               )))
  m4 <- try(sdmTMB(sim ~ -1+year+log_depth_scaled+log_depth_scaled2, 
               data = dat, 
               spatial = "on",
               mesh=mesh,
               anisotropy=T,
               reml=F,
               time=NULL,
               family =tweedie(link="log"),
               control = sdmTMBcontrol(
                 newton_loops = 2,
               )))

  m5 <- try(sdmTMB(sim ~ -1+year+log_depth_scaled+log_depth_scaled2+temp_s,
               data = dat, 
               spatial = "on",
               mesh=mesh,
               anisotropy=T,
               reml=F,
               time=NULL,
               family =tweedie(link="log"),
               control = sdmTMBcontrol(
                 newton_loops = 1,
               )))
  m6 <- try(sdmTMB(sim ~ -1+year+log_depth_scaled+log_depth_scaled2+po2_s,
               data = dat, 
               spatial = "on",
               mesh=mesh,
               anisotropy=T,
               reml=F,
               time=NULL,
               family =tweedie(link="log"),
               control = sdmTMBcontrol(newton_loops = 1,
               )))
  m7 <- try(sdmTMB(sim ~ -1+year+log_depth_scaled+log_depth_scaled2+temp_s + po2_s,
               data = dat, 
               spatial = "on",
               mesh=mesh,
               anisotropy=T,
               reml=F,
               time=NULL,
               family =tweedie(link="log"),
               control = sdmTMBcontrol(newton_loops = 1,
               )))
  m8 <- try(sdmTMB(sim ~ -1+year+log_depth_scaled+log_depth_scaled2+temp_s * po2_s,
               data = dat, 
               spatial = "on",
               mesh=mesh,
               anisotropy=T,
               reml=F,
               time=NULL,
               family =tweedie(link="log"),
               control = sdmTMBcontrol(newton_loops = 1,
               )))
  ### Create dAIC table ###
  ## Make list of model names##
  models <- c("breakpt-pO2", "Eo estimation and logistic po2' (no prior)","Eo estimation and logistic po2' (prior)","logistic-pO2", "Null", "temp", "po2", "temp+po2", "temp * po2")
  ## Create table and add AIC for each ##
  AIC <- as.data.frame(matrix(NA, ncol = 1, nrow =length(models), dimnames = list(models)))
  AIC[1,] <- try(AIC(m1))
  AIC[2,] <- try(AIC(m2))
  AIC[3,] <- try(AIC(m2a))
  AIC[4,] <- try(AIC(m3))
  AIC[5,] <- try(AIC(m4))
  AIC[6,] <- try(AIC(m5))
  AIC[7,] <- try(AIC(m6))
  AIC[8,] <- try(AIC(m7))
  AIC[9,] <- try(AIC(m8))
  
  ## Calculate delta-AIC ##
  AIC$V1 <- as.numeric(AIC$V1)
  AIC$dAIC <- try(abs(min(AIC$V1, na.rm=T)-(AIC$V1)))
  AIC$model <- models
  try(return(AIC))
}

run_alt_models_mis <- function(dat, start, mesh) {
  m1 <- try(sdmTMB(sim ~ 1+year+breakpt(po2_s)+log_depth_scaled, 
                   data = dat,
                   time = NULL,
                   reml = F,
                   anisotropy = TRUE,
                   spatiotemporal = FALSE,
                   mesh=mesh,
                   family =tweedie(link="log"),
                   # control = sdmTMBcontrol(
                   #start = list(b_threshold = start),
                   #lower = list(b_threshold = lower), 
                   #upper = list(b_threshold = upper),
                   #newton_loops = 2
  ))
  
  m2a <- try(sdmTMB(sim ~ -1+year+logistic(mi)+log_depth_scaled,
                    data = dat, 
                    time = NULL,
                    reml = F,
                    anisotropy = TRUE,
                    spatiotemporal = FALSE,
                    mesh=mesh,
                    family =tweedie(link="log"),
                    control = sdmTMBcontrol(
                      start = list(b_threshold = start),
                      #   upper = list(b_threshold = upper),
                      #  lower = list(b_threshold = lower),
                      newton_loops = 2),
                    priors=sdmTMBpriors(threshold = normal(c(NA, NA, NA, 0.3477), c(NA, NA, NA, 0.1455)))))
  
  m3 <- try(sdmTMB(sim ~ -1+year+logistic(po2_s)+log_depth_scaled,
                   data = dat, 
                   spatial = "on",
                   mesh=mesh,
                   anisotropy=T,
                   reml=F,
                   time=NULL,
                   family =tweedie(link="log"),
                   control = sdmTMBcontrol(
                     start = list(b_threshold=start),
                     #  lower = list(b_threshold = lower),
                     #  upper = list(b_threshold = upper)
                   )))
  m4 <- try(sdmTMB(sim ~ -1+year+log_depth_scaled, 
                   data = dat, 
                   spatial = "on",
                   mesh=mesh,
                   anisotropy=T,
                   reml=F,
                   time=NULL,
                   family =tweedie(link="log"),
                   control = sdmTMBcontrol(
                     newton_loops = 2,
                   )))
  
  m5 <- try(sdmTMB(sim ~ -1+year+log_depth_scaled+temp_s,
                   data = dat, 
                   spatial = "on",
                   mesh=mesh,
                   anisotropy=T,
                   reml=F,
                   time=NULL,
                   family =tweedie(link="log"),
                   control = sdmTMBcontrol(
                     newton_loops = 1,
                   )))
  m6 <- try(sdmTMB(cpue_kg_km2 ~ -1+year+log_depth_scaled+po2_s,
                   data = dat, 
                   spatial = "on",
                   mesh=mesh,
                   anisotropy=T,
                   reml=F,
                   time=NULL,
                   family =tweedie(link="log"),
                   control = sdmTMBcontrol(newton_loops = 1,
                   )))
  m7 <- try(sdmTMB(sim ~ -1+year+log_depth_scaled+temp_s + po2_s,
                   data = dat, 
                   spatial = "on",
                   mesh=mesh,
                   anisotropy=T,
                   reml=F,
                   time=NULL,
                   family =tweedie(link="log"),
                   control = sdmTMBcontrol(newton_loops = 1,
                   )))
  m8 <- try(sdmTMB(sim ~ -1+year+log_depth_scaled+temp_s * po2_s,
                   data = dat, 
                   spatial = "on",
                   mesh=mesh,
                   anisotropy=T,
                   reml=F,
                   time=NULL,
                   family =tweedie(link="log"),
                   control = sdmTMBcontrol(newton_loops = 1,
                   )))
  ### Create dAIC table ###
  ## Make list of model names##
  models <- c("breakpt-pO2", "Eo estimation and logistic po2' (prior)","logistic-pO2", "Null", "temp", "po2", "temp+po2", "temp * po2")
  ## Create table and add AIC for each ##
  AIC <- as.data.frame(matrix(NA, ncol = 1, nrow =length(models), dimnames = list(models)))
  AIC[1,] <- try(AIC(m1))
  AIC[2,] <- try(AIC(m2a))
  AIC[3,] <- try(AIC(m3))
  AIC[4,] <- try(AIC(m4))
  AIC[5,] <- try(AIC(m5))
  AIC[6,] <- try(AIC(m6))
  AIC[7,] <- try(AIC(m7))
  AIC[8,] <- try(AIC(m8))
  
  ## Calculate delta-AIC ##
  AIC$V1 <- as.numeric(AIC$V1)
  AIC$dAIC <- try(abs(min(AIC$V1, na.rm=T)-(AIC$V1)))
  AIC$model <- models
  try(return(AIC))
}
### Run alternative models for each data simulation ###
run_alt_models_cv <- function(dat, start, mesh) {
  seed <- sample(1:2000, 1)
  set.seed(seed)
  k_folds <- 5
  dat$fold_ids <- sample(seq_len(k_folds), nrow(dat), replace = TRUE)
  future::plan(future::multisession)
  m1 <- try(sdmTMB_cv(sim ~ 1+year+breakpt(po2_s)+log_depth_scaled+log_depth_scaled2, 
                   data = dat,
                   time = NULL,
                   reml = F,
                   anisotropy = TRUE,
                   spatiotemporal = FALSE,
                   mesh=mesh,
                   k_folds=k_folds,
                   fold_ids = dat$fold_ids,
                   family =tweedie(link="log"),
                   # control = sdmTMBcontrol(
                   #start = list(b_threshold = start),
                   #lower = list(b_threshold = lower), 
                   #upper = list(b_threshold = upper),
                   #newton_loops = 2
  ))
  
  m2 <- try(sdmTMB_cv(sim ~ -1+year+logistic(mi)+log_depth_scaled+log_depth_scaled2, 
                   data = dat, 
                   time = NULL,
                   reml = F,
                   anisotropy = TRUE,
                   spatiotemporal = FALSE,
                   mesh=mesh,
                   k_folds=k_folds,
                   fold_ids = dat$fold_ids,
                   family =tweedie(link="log"),
                   control = sdmTMBcontrol(
                     start = list(b_threshold = start),
                     #lower = list(b_threshold = lower),
                     # upper = list(b_threshold = upper),
                     newton_loops = 2,
                     nlminb_loops=2)))
  
  m2a <- try(sdmTMB_cv(sim ~ -1+year+logistic(mi)+log_depth_scaled+log_depth_scaled2,
                    data = dat, 
                    time = NULL,
                    reml = F,
                    anisotropy = TRUE,
                    spatiotemporal = FALSE,
                    k_folds=k_folds,
                    fold_ids = dat$fold_ids,
                    mesh=mesh,
                    family =tweedie(link="log"),
                    control = sdmTMBcontrol(
                      start = list(b_threshold = start),
                      #   upper = list(b_threshold = upper),
                      #  lower = list(b_threshold = lower),
                      newton_loops = 2),
                    priors=sdmTMBpriors(threshold = normal(c(NA, NA, NA, 0.3477), c(NA, NA, NA, 0.1455)))))
  
  m3 <- try(sdmTMB_cv(sim ~ -1+year+logistic(po2_s)+log_depth_scaled+log_depth_scaled2,
                   data = dat, 
                   spatial = "on",
                   mesh=mesh,
                   anisotropy=T,
                   reml=F,
                   time=NULL,
                   family =tweedie(link="log"),
                   k_folds=k_folds,
                   fold_ids = dat$fold_ids,
                   control = sdmTMBcontrol(
                     start = list(b_threshold=start),
                     #  lower = list(b_threshold = lower),
                     #  upper = list(b_threshold = upper)
                   )))
  m4 <- try(sdmTMB_cv(sim ~ -1+year+log_depth_scaled+log_depth_scaled2, 
                   data = dat, 
                   spatial = "on",
                   mesh=mesh,
                   anisotropy=T,
                   reml=F,
                   time=NULL,
                   k_folds=k_folds,
                   fold_ids = dat$fold_ids,
                   family =tweedie(link="log"),
                   control = sdmTMBcontrol(
                     newton_loops = 2,
                   )))
  
  m5 <- try(sdmTMB_cv(sim ~ -1+year+log_depth_scaled+log_depth_scaled2+temp_s,
                   data = dat, 
                   spatial = "on",
                   mesh=mesh,
                   anisotropy=T,
                   reml=F,
                   k_folds,
                   fold_ids = dat$fold_ids,
                   time=NULL,
                   family =tweedie(link="log"),
                   control = sdmTMBcontrol(
                     newton_loops = 1,
                   )))
  m6 <- try(sdmTMB_cv(sim ~ -1+year+log_depth_scaled+log_depth_scaled2+po2_s,
                   data = dat, 
                   spatial = "on",
                   mesh=mesh,
                   anisotropy=T,
                   reml=F,
                   time=NULL,
                   k_folds,
                   fold_ids = dat$fold_ids,
                   family =tweedie(link="log"),
                   control = sdmTMBcontrol(newton_loops = 1,
                   )))
  m7 <- try(sdmTMB_cv(sim ~ -1+year+log_depth_scaled+log_depth_scaled2+temp_s + po2_s,
                   data = dat, 
                   spatial = "on",
                   mesh=mesh,
                   anisotropy=T,
                   reml=F,
                   time=NULL,
                   k_folds,
                   fold_ids = dat$fold_ids,
                   family =tweedie(link="log"),
                   control = sdmTMBcontrol(newton_loops = 1,
                   )))
  m8 <- try(sdmTMB_cv(sim ~ -1+year+log_depth_scaled+log_depth_scaled2+temp_s * po2_s,
                   data = dat, 
                   spatial = "on",
                   mesh=mesh,
                   anisotropy=T,
                   reml=F,
                   time=NULL,
                   k_folds,
                   fold_ids = dat$fold_ids,
                   family =tweedie(link="log"),
                   control = sdmTMBcontrol(newton_loops = 1,
                   )))
  ### Create dAIC table ###
  ## Make list of model names##
  models <- c("breakpt-pO2", "Eo estimation and logistic po2' (no prior)","Eo estimation and logistic po2' (prior)","logistic-pO2", "Null", "temp", "po2", "temp+po2", "temp * po2")
  ## Create table and add LL for each ##
  LL <- as.data.frame(matrix(NA, ncol = 3, nrow =length(models), dimnames = list(models)))
  LL$model <- models
  LL$model_id <- 1:9
  LL[1,3] <-try(m1$sum_loglik)
  LL[2,3] <-try(m2$sum_loglik)
  LL[3,3] <-try(m2a$sum_loglik)
  LL[4,3] <-try(m3$sum_loglik)
  LL[5,3] <-try(m4$sum_loglik)
  LL[6,3] <-try(m5$sum_loglik)
  LL[7,3] <-try(m6$sum_loglik)
  LL[8,3] <-try(m7$sum_loglik)
  LL[9,3] <-try(m8$sum_loglik)
  LL$V1 <- as.numeric(LL$V1)
  try(return(LL))
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

logfun2 <- function(data, po2_prime, model) {
  if(!is.character(model)){
  parfit <- model$sd_report
  npars <- length(parfit$value)
  parnames <- names(parfit$value)
  
  s50 <- parfit$value[grep("s50", parnames)]
  delta <- parfit$value[grep("s95", parnames)]
  smax <- parfit$value[grep("s_max", parnames)]
  Eo <- parfit$value[grep("Eo", parnames)]
  
  a <- log(smax / (log(0.5) + smax) - 1)
  b <- log(smax / (log(0.95) + smax) - 1)
  beta0 <- -a + s50 * (b - a) / delta
  beta1 <- (a - b) / delta
  x <- as.numeric(data[ , po2_prime])
  data$logmu <- exp(smax * (1 / ( 1 + exp( - beta0 - beta1 * x)) -1))
  data$Eo <- Eo
  data$s50 <- s50
  }
  if(is.character(model)){
   data <- NA
  }
  return(data)
} 

## Number of zero observations below threshold for each data simulation ##
count_below_zero <- function(dat, threshold) {
  count <- subset(dat, dat$sim==0 & dat$mi_weird < threshold)
  count <- nrow(count)
  return(count)
}
