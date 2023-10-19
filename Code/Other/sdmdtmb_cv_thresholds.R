remotes::install_github("pbs-assess/sdmTMB", dependencies = TRUE, ref="newlogistic")
library(sdmTMB)

dat <- pcod
mesh <- make_mesh(dat, c("X", "Y"), cutoff = 10)

#Set CV folds
seed <- sample(1:2000, 1)
set.seed(seed)
k_folds <- 5
dat$fold_ids <- sample(seq_len(k_folds), nrow(dat), replace = TRUE)
future::plan(future::multisession)
#Starting values
start <- c(1,1)

m1 <- sdmTMB_cv(density ~ 1+as.factor(year)+breakpt(depth_scaled), 
                data = dat,
                time = NULL,
                reml = F,
                anisotropy = TRUE,
                spatiotemporal = FALSE,
                mesh=mesh,
                family =tweedie(link="log"),
                control = sdmTMBcontrol(
                start = list(b_threshold = start)))

prior <-normal(c(NA, 0.331), c(NA, 0.176))
m1 <- sdmTMB_cv(density ~ 1+as.factor(year)+breakpt(depth_scaled), 
                data = dat,
                time = NULL,
                reml = F,
                anisotropy = TRUE,
                spatiotemporal = FALSE,
                mesh=mesh,
                family =tweedie(link="log"),
                priors=sdmTMBpriors(threshold = prior))
