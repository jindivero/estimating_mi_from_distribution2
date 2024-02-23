library(sdmTMB)
library(stringr)
install.packages("ggbeeswarm")
library(ggbeeswarm)

### Set ggplot themes ###
theme_set(theme_bw(base_size = 35))
theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Set parameter values for data generation and model fitting #
s50 <-2 # same as sablefish= 0.88; had been 2 in previous data simulation
delta <- 2 #same as sablefish = 0.57; had been 2 in previous data simulation
smax <- 30 # maximum effect of MI; had been 4 in previous data simulation
Eo <- 0.3
Eo2 <- 0.7

b_years <- rnorm(n = 6, mean = 4, sd = 1)
beta1 <- 1.5
beta2 <- -1
phi <- 10 # 16 corresponds to sablefish 
p <- 1.51
range <- 85
sigma_O <- 1.77

### load helper functions ####
source("Code/util_funs.R")
source("Code/sim_funs.R")

use_previous <- T
if(use_previous==T){
  simdat <- readRDS("Model Outputs/data_sims_usual.rds")
  simdat2 <- readRDS("Model Outputs/data_sims_weird.rds")
}

#Load data, start, mesh
mesh <- make_mesh(simdat[[1]], xy_cols = c("X", "Y"), n_knots=250)

start <- matrix(0, ncol = 1, nrow = 4)
start[1,1] <- s50
start[2,1] <- delta
start[3,1] <- smax
start[4,1] <- Eo

start_unusual <- start
start_unusual[4,1] <- Eo2

compare1 <- lapply(simdat, run_alt_models, start, mesh)
compare2 <- lapply(simdat2, run_alt_models, start_unusual, mesh)

if(save){
  save(compare1, compare2, file="AIC.Rdata")
}
if(use_previous){
  load("Model Outputs/simulation_AIC.R")
}
# Combine #
AIC1 <- keep(.x=compare1, .p=is.data.frame)
AIC1 <- bind_rows(compare1, .id="df")

AIC2 <- keep(.x=compare2, .p=is.data.frame)
AIC2 <- bind_rows(compare2, .id="df")

#Average and SD AIC
AIC_summary <- aggregate(dAIC~model, data=AIC1, FUN=mean)
AIC_summary1a <- aggregate(dAIC~model, data=AIC1, FUN=sd)
AIC_summary1b <- aggregate(dAIC~model, data=AIC1, FUN=fivenum)

AIC_summary2 <- aggregate(dAIC~model, data=AIC2, FUN=mean)
AIC_summary2a <- aggregate(dAIC~model, data=AIC2, FUN=sd)

#Count best-fitting in each model type
AIC_count1 <- subset(AIC1, AIC1$dAIC==0)
AIC_count2 <- subset(AIC2, AIC2$dAIC==0)
AIC_count1 <- aggregate(df~model, data=AIC_count1, FUN=length)
AIC_count2 <- aggregate(df~model, data=AIC_count2, FUN=length)

#Get list of data iterations that had estimation issues
AIC_na <- subset(AIC1, AIC1$model=="Eo estimation and logistic po2' (no prior)"&is.na(AIC1$dAIC))
nas <- as.vector(as.numeric(AIC_na$df))
AIC_na2 <- subset(AIC2, AIC2$model=="Eo estimation and logistic po2' (no prior)"&is.na(AIC2$dAIC))
nas2 <- as.vector(as.numeric(AIC_na2$df))

#Isolate data iterations 
simdat_na <- simdat[nas]
simdat2_na <- simdat[nas2]

#Run alternative models including prior
#Try different start
starta <- start*0.9
start_unusuala <- start_unusual*0.9
compare1a <- lapply(simdat_na, run_alt_models_prior, start, mesh)
compare2a <- lapply(simdat2_na, run_alt_models_prior, start_unusual, mesh)

#Get average, sd, counts like above
AIC1a <- keep(.x=compare1a, .p=is.data.frame)
AIC1a <- bind_rows(compare1a, .id="df")

AIC2a <- keep(.x=compare2a, .p=is.data.frame)
AIC2a <- bind_rows(compare2a, .id="df")

#Average and SD AIC
AIC_summarya <- aggregate(dAIC~model, data=AIC1a, FUN=mean)
AIC_summary1aa <- aggregate(dAIC~model, data=AIC1a, FUN=sd)

AIC_summary2a <- aggregate(dAIC~model, data=AIC2a, FUN=mean)
AIC_summary2aa <- aggregate(dAIC~model, data=AIC2a, FUN=sd)

#Count best-fitting in each model type
AIC_count1a <- subset(AIC1a, AIC1a$dAIC==0)
AIC_count2a <- subset(AIC2a, AIC2a$dAIC==0)
AIC_count1a <- aggregate(df~model, data=AIC_count1a, FUN=length)
AIC_count2a <- aggregate(df~model, data=AIC_count2a, FUN=length)

### Mis-specified
#Count best-fitting model of each type
start <- matrix(c(1,1,1000,-0.3), ncol = 1, nrow = 4)
compare3 <- lapply(simdat, run_alt_models_mis, start, mesh)
compare4 <- lapply(simdat2, run_alt_models_mis, start, mesh)

save(compare3, compare4, file="AIC_mis.Rdata")

AIC3 <- keep(.x=compare3, .p=is.data.frame)
AIC3 <- bind_rows(compare3, .id="df")

AIC4 <- keep(.x=compare4, .p=is.data.frame)
AIC4 <- bind_rows(compare4, .id="df")

AIC_summary3 <- aggregate(dAIC~model, data=AIC3, FUN=mean)
AIC_summary4 <- aggregate(dAIC~model, data=AIC4, FUN=mean)
AIC_summary3a <- aggregate(dAIC~model, data=AIC3, FUN=sd)
AIC_summary4a <- aggregate(dAIC~model, data=AIC4, FUN=sd)

AIC_count3 <- subset(AIC3, AIC3$dAIC==0)
AIC_count4 <- subset(AIC4, AIC4$dAIC==0)
AIC_count3 <- aggregate(df~model, data=AIC_count3, FUN=length)
AIC_count4 <- aggregate(df~model, data=AIC_count4, FUN=length)







###Cross-validation ###

compare1 <- lapply(test, run_alt_models_cv, start, mesh)




## Make predictions from existing model fits ##
preds1<- mapply(FUN=predict_sims, fits, simdat_cv, p, phi, SIMPLIFY=F)
preds2<- mapply(FUN=predict_sims, fits2, simdat_cv, p, phi, SIMPLIFY=F)
preds3<- mapply(FUN=predict_sims, fits3, simdat2_cv, p, phi, SIMPLIFY=F)
preds4<- mapply(FUN=predict_sims, fits4, simdat2_cv, p, phi, SIMPLIFY=F)

## Apply to each ##
nll1 <- mapply(FUN=calculate_nll, simdat_cv, preds1, "sim", "est", p, phi, SIMPLIFY=F)
nll2 <- mapply(FUN=calculate_nll, simdat_cv, preds2, "sim", "est", p, phi, SIMPLIFY=F)
nll3 <- mapply(FUN=calculate_nll, simdat2_cv, preds3, "sim", "est", p, phi, SIMPLIFY=F)
nll4 <- mapply(FUN=calculate_nll, simdat2_cv, preds4, "sim", "est", p, phi, SIMPLIFY=F)


nll_sum1 <- mapply(FUN=sum_nll, nll1, "nll", SIMPLIFY=F)
nll_sum2 <- mapply(FUN=sum_nll, nll2, "nll", SIMPLIFY=F)
nll_sum3 <- mapply(FUN=sum_nll, nll3, "nll", SIMPLIFY=F)
nll_sum4 <- mapply(FUN=sum_nll, nll4, "nll", SIMPLIFY=F)


nll1_thresh <- mapply(FUN=sum_nll_thresh, nll1, "nll", "mi_usual", SIMPLIFY=F)
nll2_thresh <- mapply(FUN=sum_nll_thresh, nll2, "nll", "mi_usual", SIMPLIFY=F)
nll3_thresh <- mapply(FUN=sum_nll_thresh, nll3, "nll", "mi_weird", SIMPLIFY=F)
nll4_thresh <- mapply(FUN=sum_nll_thresh, nll4, "nll", "mi_weird", SIMPLIFY=F)


nll1_above <- mapply(FUN=sum_nll_above, nll1, "nll", "mi_usual", SIMPLIFY=F)
nll2_above <- mapply(FUN=sum_nll_above, nll2, "nll", "mi_usual", SIMPLIFY=F)
nll3_above <- mapply(FUN=sum_nll_above, nll3, "nll", "mi_weird", SIMPLIFY=F)
nll4_above <- mapply(FUN=sum_nll_above, nll4, "nll", "mi_weird", SIMPLIFY=F)


#Subset to where observations have no fish
#nll1_thresh <- lapply(nll1, subset, pred2==0)
#nll2_thresh <- lapply(nll2, subset, pred2==0)
#nll3_thresh <- lapply(nll3, subset, pred2==0)

## Combine into one dataset for plotting##
# Overall #
nll_all <- as.data.frame(unlist(nll_sum1))
nll_all <- cbind(nll_all, unlist(nll_sum2))
nll_all <- cbind(nll_all, unlist(nll_sum3))
nll_all <- cbind(nll_all, unlist(nll_sum4))
colnames(nll_all) <- c( "Typical Case, Unconstrained", "Typical Case, Prior Constrained", "Unusual Case, Unconstrained", "Unusual Case, Prior Constrained")

# Above s95 #
nll_95 <- as.data.frame(unlist(nll1_above))
nll_95$model2 <- as.vector(unlist(nll2_above))
nll_95$model3 <- as.vector(unlist(nll3_above))
nll_95$model4 <- as.vector(unlist(nll4_above))
colnames(nll_95) <- c( "Typical Case, Unconstrained", "Typical Case, Prior Constrained", "Unusual Case, Unconstrained", "Unusual Case, Prior Constrained")

# Below s95 #
nll_below <- as.data.frame(unlist(nll1_thresh))
nll_below$model2 <- as.vector(unlist(nll2_thresh))
nll_below$model3 <- as.vector(unlist(nll3_thresh))
nll_below$model4 <- as.vector(unlist(nll4_thresh))
colnames(nll_below) <- c( "Typical Case, Unconstrained", "Typical Case, Prior Constrained", "Unusual Case, Unconstrained", "Unusual Case, Prior Constrained")

# Pivot long #
nll_all <- pivot_longer(nll_all, c(1:4), names_to="Model")
nll_below <- pivot_longer(nll_below, c(1:4), names_to="Model")
nll_95 <- pivot_longer(nll_95, c(1:4), names_to="Model")

# Relabel #
colnames(nll_all)[2] <- "Overall"
colnames(nll_below)[2] <- "Below_s50"
colnames(nll_95)[2] <- "Above_s95"

# Combine #
nll_combined <- cbind(nll_all,nll_95, nll_below)
nll_combined <- nll_combined[, c(1,2,4,6)]

### Run logistic pO2 equations and run predictions ###
## Function to run model 2 (with prior) ##
start <- matrix(0, ncol = 1, nrow = 3)
start[1,1] <-  -1.1
start[2,1] <- -1.1
start[3,1] <- 15


## Fit model to all simulated datasets ##
fits5 <- lapply(simdat, run_sdmTMB_3, 
                start=start, mesh=mesh)

start <- matrix(0, ncol = 1, nrow = 3)
start[1,1] <-  -1.1
start[2,1] <- -1
start[3,1] <- 1000

fits6 <- lapply(simdat2, run_sdmTMB_3, 
                start=start, mesh=mesh)

### Make predictions from logistic pO2 model fits ###
preds5<- mapply(FUN=predict_sims, fits5, simdat_cv, p, phi, SIMPLIFY=F)
preds6<- mapply(FUN=predict_sims, fits6, simdat2_cv, p, phi, SIMPLIFY=F)

nll5 <- mapply(FUN=calculate_nll, simdat_cv, preds5, "sim", "est", p, phi, SIMPLIFY=F)
nll6 <- mapply(FUN=calculate_nll, simdat2_cv, preds6, "sim", "est", p, phi, SIMPLIFY=F)
nll_sum5 <- mapply(FUN=sum_nll, nll5, "nll", SIMPLIFY=F)
nll_sum6 <- mapply(FUN=sum_nll, nll6, "nll", SIMPLIFY=F)

# Subset logistic(po2) into thresholds #
# Above s95 #
nll5_thresh <- mapply(FUN=sum_nll_thresh, nll5, "nll", "mi_usual", SIMPLIFY=F)
nll6_thresh <- mapply(FUN=sum_nll_thresh, nll6, "nll", "mi_unusual", SIMPLIFY=F)

nll5_above <- mapply(FUN=sum_nll_above, nll5, "nll", "mi_usual", SIMPLIFY=F)
nll6_above <- mapply(FUN=sum_nll_above, nll6, "nll", "mi_unusual", SIMPLIFY=F)

## Combine into one dataset for plotting##
# Above s95 #
nll_log_95 <- as.data.frame(unlist(nll5_above))
nll_log_95$model2 <- as.vector(unlist(nll6_above))
colnames(nll_log_95) <- c( "Typical Case, Logistic pO2", "Unusual Case, Logistic pO2")

# Below s95 #
nll_log_below <- as.data.frame(unlist(nll5_thresh))
nll_log_below$model2 <- as.vector(unlist(nll6_thresh))
colnames(nll_log_below) <- c( "Typical Case, Logistic pO2", "Unusual Case, Logistic pO2")

# Pivot long #
nll_log_below <- pivot_longer(nll_log_below, c(1:2), names_to="Model")
nll_log_95 <- pivot_longer(nll_log_95, c(1:2), names_to="Model")
nll_log <- pivot_longer(nll_log, c(1:2), names_to="Model")

# Relabel #
colnames(nll_log_below)[2] <- "Below_s50"
colnames(nll_log_95)[2] <- "Above_s95"
colnames(nll_log)[2] <- "Overall"

# Combine #
nll_log_combined <- cbind(nll_log_below, nll_log_95, nll_log)
nll_log_combined <- nll_log_combined[,c(1,2,4,6)]

nll_combined <- bind_rows(nll_combined, nll_log_combined)

#Average NLL for plotting
LL_values <- as.data.frame(matrix(nrow=6))
LL_values$V1 <- NULL
LL_values$data <- c("Typical Case", "Typical Case", "Typical Case", "Unusual Case", "Unusual Case", "Unusual Case")
LL_values$analysis <- c( "Logistic", "Prior", "Unconstrained", "Logistic", "Prior", "Unconstrained")
LL_values$model <- c("Typical Case, Logistic pO2", "Typical Case, Prior Constrained", "Typical Case, Unconstrained","Unusual Case, Logistic pO2", "Unusual Case, Prior Constrained", "Unusual Case, Unconstrained")
LL_avg <- aggregate(Below_s50~Model, nll_combined, FUN=mean)
LL_values$LL_avg <- LL_avg$Below_s50

# Add separate model and data columns (and remove the comma from the data type) #
nll_combined$data <- word(nll_combined$Model, 1,2, sep=" ")
nll_combined$data <- str_sub(nll_combined$data, 1, str_length(nll_combined$data)-1)
nll_combined$analysis <- word(nll_combined$Model, 3)

ggplot(nll_combined, aes(x=Below_s50))+
  geom_density(fill="lightblue", adjust = 1.5)+
  facet_grid(analysis~data)+
  geom_vline(data = LL_values, aes(xintercept = LL_avg),linetype="dashed", size=1.2, color="darkorange", show.legend=T)+
  ylab("Density")+
  xlab("Sum Log-Likelihood for Observations Below True s50 of pO2'")

ggplot(nll_combined, aes(x=Overall))+
  geom_density(fill="lightblue", adjust = 1.5)+
  facet_grid(analysis~data)+
  ylab("Density")+
  xlab("Sum Log-Likelihood for All Observations'")

##### Extra stuff---ignore #####
###Comparing po2' and predictions for all simulations (to see how predictions compare to parameters) ###
## Predictions of fish density on original data for unconstrained models ##
pred1<- mapply(FUN=predict_sims, fits, simdat, p, phi, SIMPLIFY=F)
pred2<- mapply(FUN=predict_sims, fits3, simdat2, p, phi, SIMPLIFY=F)
## Calculate po2 prime from estimated parameters for each ##
simdats1 <- mapply(FUN=calculate_po2_prime, pred1,fits, SIMPLIFY=F)
simdats2 <- mapply(FUN=calculate_po2_prime, pred2,fits3, SIMPLIFY=F)

##Remove datasets with NA and then bind together, adding a column with the data frame number ##
simdats1 <- keep(.x=simdats1, .p=is.data.frame)
simdats1 <- bind_rows(simdats1, .id="df")
simdats2 <- keep(.x=simdats2, .p=is.data.frame)
simdats2 <- bind_rows(simdats2, .id="df")
## Plot all ##
# Points #
ggplot(simdats1, aes(po2_prime, pred2, colour=as.factor(df))) +
  geom_point(size=0.1)+
  scale_colour_viridis(discrete=TRUE)+
  xlab("pO2'")+
  ylab("Predicted density")+
  theme(legend.position="none")
# Density plot #
ggplot(simdats1, aes(po2_prime, pred2, colour=as.factor(df))) +
  geom_density_2d()+
  scale_colour_viridis(discrete=TRUE)+
  xlab("pO2'")+
  ylab("Predicted density")+
  theme(legend.position="none")

## Just a couple
# Points #
ggplot(subset(simdats1, simdats1$df=="1"|simdats1$df=="19"), aes(po2_prime, pred2, colour=as.factor(df))) +
  geom_point(size=2)+
  scale_colour_viridis(discrete=TRUE, labels=c('E0=0.68, s50=2.4, smax=129', "E0=0 -0.056, s50=1.4, smax=22"))+
  xlab("pO2'")+
  ylab("Predicted density")+
  theme(legend.position=c(0.7,0.8))

## Just a couple
# Density # (this looks weird)
ggplot(subset(simdats1, simdats1$df=="1"|simdats1$df=="19"), aes(po2_prime, pred2, colour=as.factor(df))) +
  stat_density_2d(
    geom = "raster",
    aes(fill = after_stat(density)),
    contour = FALSE
  ) +
  scale_colour_viridis(discrete=TRUE)+
  xlab("pO2'")+
  ylab("Predicted density")



# simulation version #
smax_test <- as.data.frame(dat$mi_usual)
colnames(smax_test) <- "po2_prime"

smax_test$logmu1 <- logfun_basic(smax_test$po2_prime, smax=5, s50=s50, delta)
smax_test$logmu2 <- logfun_basic(smax_test$po2_prime, smax=20, s50=s50, delta)
smax_test$logmu3 <- logfun_basic(smax_test$po2_prime, smax=50, s50=s50, delta)
smax_test$logmu4 <- logfun_basic(smax_test$po2_prime, smax=100, s50=s50, delta)
smax_test$logmu5 <- logfun_basic(smax_test$po2_prime, smax=500, s50=s50, delta)
smax_test$logmu6 <- logfun_basic(smax_test$po2_prime, smax=1000, s50=s50, delta)
smax_test$logmu7 <- logfun_basic(smax_test$po2_prime, smax=10000, s50=s50, delta)

smax_test <- pivot_longer(smax_test, 2:8)

ggplot(smax_test, aes(x=po2_prime, y=value))+geom_point(aes(group=name, color=name))+
  scale_colour_viridis(discrete=T)+
  xlab("pO2'")+
  ylab("pO2' effect")+
  scale_colour_discrete(type="viridis", "smax value", labels=c("5", "20", "50", "100", "500", "1000", "10000"))+theme(legend.position=c(0.8,0.2))


##Correlation of Eo and number of zero observations in dataset #
# Make wide #
pars_wide <- pivot_wider(pars, id_cols=c(id, model), names_from=term, values_from=estimate)
# Plots #
ggplot(pars_wide, aes(y=pars_wide$"mi-smax", x=pars_wide$"mi-Eo"))+geom_point()+
  facet_wrap("model", scales="free")+
  geom_vline(data = subset(true_pars, true_pars$term=="mi-Eo"), aes(xintercept = estimate),linetype="dashed", size=1.2)+
  geom_hline(data = subset(true_pars,  true_pars$term=='mi-smax'), aes(yintercept = estimate),linetype="dashed", size=1.2)

ggplot(pars_wide, aes(y=pars_wide$"mi-s50", x=pars_wide$"mi-Eo"))+geom_point()+
  facet_wrap("model", scales="free")+
  geom_vline(data = subset(true_pars, true_pars$term=="mi-Eo"), aes(xintercept = estimate),linetype="dashed", size=1.2)+
  geom_hline(data = subset(true_pars,  true_pars$term=='mi-s50'), aes(yintercept = estimate),linetype="dashed", size=1.2)

ggplot(pars_wide, aes(y=pars_wide$"mi-delta", x=pars_wide$"mi-Eo"))+geom_point()+
  facet_wrap("model", scales="free")+
  geom_vline(data = subset(true_pars, true_pars$term=="mi-Eo"), aes(xintercept = estimate),linetype="dashed", size=1.2)+
  geom_hline(data = subset(true_pars,  true_pars$term=='mi-delta'), aes(yintercept = estimate),linetype="dashed", size=1.2)

ggplot(pars_wide, aes(y=pars_wide$"log_depth_scaled2", x=pars_wide$"mi-Eo"))+geom_point()+
  facet_wrap("model", scales="free")


counts_below_zero <- lapply(simdat, count_below_zero,threshold=(s50+delta))
counts_below_zero <- as.data.frame(unlist(counts_below_zero))
counts_below_zero$id <- as.character(1:100)
counts_below_zero2 <- lapply(simdat2, count_below_zero,threshold=(s50+delta))
counts_below_zero2 <- as.data.frame(unlist(counts_below_zero2))
counts_below_zero2$id <- as.character(1:100)

# Add to pars_wide #
pars_wide1 <- subset(pars_wide, model=='Unusual Case, Unconstrained')
pars_wide1 <- left_join(pars_wide1, counts_below_zero2, by="id")

ggplot(pars_wide1, aes(y=pars_wide1$"unlist(counts_below_zero2)", x=pars_wide1$"mi-Eo"))+geom_point()+
  facet_wrap("model", scales="free")


### Other plots of parameter estimates###
#plot boxplots of a single model #
ggplot(pars1, aes(y=estimate, x=term))+geom_boxplot()+facet_wrap("term", scales="free")+geom_hline(data = true_pars, aes(yintercept = estimate),linetype="dashed", size=1.2)+theme(legend.position="left", strip.text = element_blank())
ggplot(pars2, aes(y=estimate, x=term))+geom_boxplot()+facet_wrap("term", scales="free")+geom_hline(data = true_pars, aes(yintercept = estimate),linetype="dashed", size=1.2)+theme(legend.position="left", strip.text = element_blank())
ggplot(pars3, aes(y=estimate, x=term))+geom_boxplot()+facet_wrap("term", scales="free")+geom_hline(data = true_pars2, aes(yintercept = estimate),linetype="dashed", size=1.2)+theme(legend.position="left", strip.text = element_blank())
ggplot(pars4, aes(y=estimate, x=term))+geom_boxplot()+facet_wrap("term", scales="free")+geom_hline(data = true_pars2, aes(yintercept = estimate),linetype="dashed", size=1.2)+theme(legend.position="left", strip.text = element_blank())

#Plot boxplots of all #
ggplot(pars, aes(y=estimate, x=model))+geom_boxplot(aes(group=model, fill=model))+
  facet_wrap("term", scales="free")+  
  geom_hline(data = true_pars2, aes(yintercept = estimate),linetype="dashed", size=1.2)+
  geom_hline(data = true_pars, aes(yintercept = estimate),linetype="dashed", size=1.2)+
  theme(legend.position="left")+
  scale_x_discrete(labels=c("1", "2", "3", "4"))

#Plot boxplots of just Eo and logistic parameters #
ggplot(subset(pars, pars$term=="mi-Eo"|pars$term=="mi-smax"|pars$term=="mi-s50"|pars$term=="mi-delta"), aes(y=estimate, x=model))+geom_boxplot(aes(group=model, fill=model))+
  facet_wrap("term", scales="free")+  
  geom_hline(data = subset(true_pars, true_pars2$term=="mi-Eo"|true_pars$term=="mi-smax"|true_pars$term=="mi-s50"|true_pars$term=="mi-delta"), aes(yintercept = estimate),linetype="dashed", size=1.2)+
  geom_hline(data = subset(true_pars2, true_pars2$term=="mi-Eo"), aes(yintercept = estimate),linetype="dashed", size=1.2)+
  theme(legend.position="left")+
  scale_x_discrete(labels=c("1", "2", "3", "4"))

#Restrict smax #
pars$estimate2 <- ifelse(pars$estimate>1000 & pars$term=="mi-smax", 1000, pars$estimate)

# Plot just Eo and smax #
ggplot(subset(pars, pars$term=="mi-Eo"|pars$term=="mi-smax"|pars$term=="mi-s50"|pars$term=="mi-delta"), aes(y=estimate2, x=model))+geom_boxplot(aes(group=model, fill=model))+
  facet_wrap("term", scales="free")+  
  geom_hline(data = subset(true_pars, true_pars2$term=="mi-Eo"|true_pars$term=="mi-smax"|true_pars$term=="mi-s50"|true_pars$term=="mi-delta"), aes(yintercept = estimate),linetype="dashed", size=1.2)+
  geom_hline(data = subset(true_pars2, true_pars2$term=="mi-Eo"), aes(yintercept = estimate),linetype="dashed", size=1.2)+
  theme(legend.position="left")+
  scale_x_discrete(labels=c("1", "2", "3", "4"))

### Other Performance Evlauation Options ###

# SD of average estimate (precision) #
sd_avg <- aggregate(estimate ~ term+model, pars, FUN=stats::sd)
# Range of average estimates (precision) #
range_pars <- aggregate(estimate ~ term+model, pars, FUN=fivenum)
# Range of standard error (precision) #
sd_range <- aggregate(std.error~ term+model, pars, FUN=fivenum)
# Standard deviation of standard error #
sd_sd <- aggregate(std.error ~ term+model, pars, FUN=stats::sd)
# Precision #
#95% range (upper 95% CI--lower 95% CI) #
pars$con.range <- pars$conf.high-pars$conf.low
#Average 95% range
avg_conf.range <- aggregate(con.range ~ term+model, pars, FUN=mean)

### Simulation diagnostics ###
## Convergence/fitting errors ##
convergence1 <- lapply(fits, extract_convergence)
convergence1 <- bind_rows(convergence1)
convergence2 <- lapply(fits2, extract_convergence)
convergence2 <- bind_rows(convergence2)
convergence3 <- lapply(fits3, extract_convergence)
convergence3 <- bind_rows(convergence2)
convergence4 <- lapply(fits4, extract_convergence)
convergence4 <- bind_rows(convergence2)

# Add column for data simulation and model
convergence1$sim <- 1:100
convergence1$model <- model_names[1]
convergence2$sim <- 1:100
convergence2$model <- model_names[2]
convergence3$sim <- 1:100
convergence3$model <- model_names[3]
convergence4$sim <- 1:100
convergence4$model <- model_names[4]

#Bind all together
convergence <-bind_rows(convergence1, convergence2,convergence3, convergence4)
#Truncate message
convergence$message <- substr(convergence$message, 1, 100)
##Summarize number of each type of convergence for each model
cons <- convergence %>% group_by(message, model) %>% summarise(count=n()) 
##Flip to wide
cons <-pivot_wider(cons, names_from=model, values_from=count)
colnames(cons)[1] <- "metric"

## Hessian matrix ##
pdHess1 <- lapply(fits, extract_pdHess)
pdHess <-bind_rows(pdHess1)
hess<- extract_pdHess2(pdHess)

pdHess2 <- lapply(fits2, extract_pdHess)
pdHess2 <-bind_rows(pdHess2)
hess2<- extract_pdHess2(pdHess2)

pdHess3 <- lapply(fits3, extract_pdHess)
pdHess3 <-bind_rows(pdHess3)
hess3<- extract_pdHess2(pdHess3)

pdHess4 <- lapply(fits4, extract_pdHess)
pdHess4 <-bind_rows(pdHess4)
hess4<- extract_pdHess2(pdHess4)

hess <- as.data.frame(matrix(nrow=1, ncol=3))
colnames(hess) <- model_names
hess[1,1] <- "Positive definite hessian matrix"
hess[1,2] <- extract_pdHess2(pdHess1)
hess[1,3] <- extract_pdHess2(pdHess2)
hess[1,4] <- extract_pdHess2(pdHess3)
hess[1,5] <- extract_pdHess2(pdHess4)
# Add to convergence table #
diagnostics <- rbind(cons,setNames(hess, names(cons)))

## Gradients ##

# Apply to all models #
grad1<- lapply(fits, extract_grad)
grad2<- lapply(fits2, extract_grad)
grad3<- lapply(fits3, extract_grad)
grad4<- lapply(fits4, extract_grad)

grads1<- clean_grad(grad1, fits)
grads2<- clean_grad(grad2, fits2)
grads3<- clean_grad(grad3, fits3)
grads4<- clean_grad(grad4, fits4)

# Bind together #
grads1$model2 <-grads2$big_grad
grads1$model3 <-grads3$big_grad
grads1$model4 <-grads4$big_grad
colnames(grads1) <-c("metric", model_names[1], model_names[2], model_names[3], model_names[4])
grads1$metric <- paste0("Low gradient ", grads1$metric)

## NA in SD ##
pars$NAs <- ifelse(pars$std.error=="NaN"|is.na(pars$std.error), 0,1)
NAs <- aggregate(pars$NAs~pars$term+pars$model, FUN=sum)
colnames(NAs) <- c("metric", "model", "value")
NAs$metric <- paste0("SE estimated", NAs$metric)
# Flip to wide #
NAs <-pivot_wider(NAs, names_from=model, values_from=value)

## Combine into table summarizing all diagnostics ##
diagnostics <- bind_rows(diagnostics, grads1, NAs)

#Plot RMSE
#ggplot(rmse, aes(x=model, y=rmse))+geom_point(aes(group=model, color=model))

# Average standard error (precision)--also to compare to SD of prior #
#std_error_within_model <- aggregate(std.error ~ term+model, pars, FUN=mean) 
