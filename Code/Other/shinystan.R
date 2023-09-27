library(tmstan)
library(shinystan)
mcmc <-tmbstan(fits_test)

# Run to get parallel chains
init.fn <- function()
  list(u=rnorm(114), beta=rnorm(2), logsdu=runif(1,0,10), logsd0=runif(1,0,1))
fit <- tmbstan(obj, chains=5, cores=5, open_progress=FALSE, init=init.fn)

init='last.par.best'
chains=5, cores=5 to get parallel chains.
mcmc = tmbstan( tmbobj )
launch_shinystan(as.shinystan( mcmc ))
