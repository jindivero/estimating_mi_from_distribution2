#Load packages
library(tmbstan)

#Load model fits
load("model_fits.Rdata")

#Isolate one to test
fit_mle <- fits[[1]]
#Pass to tmbstan
fit_stan <- tmbstan::tmbstan(
  fit_mle$tmb_obj,
  iter = 1000, chains = 3,
  seed = 8217 # ensures repeatability
  #control = list(adapt_delta = 0.9, max_treedepth = 12)
)

fit_stan

#Look at rhat (<1.05) and n_ff (>100)

pairs(fit_stan)

plot(fit_stan)

pars_plot <- c("b_j[1]", "b_j[2]", "ln_tau_O", "omega_s[1]")

bayesplot::mcmc_trace(fit_stan, pars = pars_plot)

bayesplot::mcmc_pairs(fit_stan, pars = pars_plot)