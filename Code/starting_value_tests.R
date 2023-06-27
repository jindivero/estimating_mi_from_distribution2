##For model 2, with full data no depth constraints
#converges and no NAs, and a more gradual logistic, but a very negative Eo (-7.7)
start <- matrix(c(1.5,0.3,1,0))
lower <- matrix(c(-Inf, -Inf, -Inf,-Inf))
upper <- matrix(c(Inf, Inf, Inf, Inf))

#non-positive-definite-Hessian, steep logistic, NA on Eo, and Eo at 0.01
start <- matrix(c(1.5,0.3,1,0))
lower <- matrix(c(0.01, 0.01, 1,0.01))
upper <- matrix(c(Inf, Inf, Inf, Inf))

#non-positive-definite-hessian, breakpoint-looking plot (vertical slope, with s50=0 and smax=0) and a hugely positive Eo (25)
start <- matrix(c(1.5,0.3,125,2)) 
lower <- matrix(c(-Inf, -Inf, -Inf,-Inf))
upper <- matrix(c(Inf, Inf, Inf, Inf))

#non-positive definite hessian, but get Eo=4.79, smax=68.04, an NA in the delta and a pretty sharp logistic (s50 and delta above 0), but still more gradual than above
start <- matrix(c(1.5,0.3,125,2))
lower <- matrix(c(0.1, 0.1,1, 0))
upper <- matrix(c(Inf, Inf, Inf, 5))

# non-positive-definite hessian, 0.10 for Eo, but NAs on SD for s50, delta, and Eo, and a two-stage logistic curve (really flat at low levels, than sharp jump up to the threshold)
start <- matrix(c(1.5,0.2,0.01,0.2))
lower <- matrix(c(0.01, 0.1,0.01, 0.1))
upper <- matrix(c(Inf, Inf, Inf, 5))


### For model m2a ###
#constraining depth
prior <- normal(c(NA, NA, NA, 0.448), c(NA, NA, NA, 0.15))
start <- matrix(c(3,1,0.5,0.1))
lower <- matrix(c(0.01, 0.01,0.01, 0.01))
upper <- matrix(c(3, 1, 3, 5))

###For model 3 ###
#w/out constraining depth
start <- matrix(c(1.5,0.3,120,0.3))
lower <- matrix(c(0.1, 0.1,1, 0.01))
upper <- matrix(c(Inf, Inf, Inf, 5))