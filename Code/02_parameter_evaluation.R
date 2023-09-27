### Install Packages ###
#install.packages("remotes")
library(remotes)
#install.packages("devtools")
library(devtools)
#install.packages("pkgbuild")
library(pkgbuild)
remotes::install_github("pbs-assess/sdmTMB", dependencies = TRUE,  ref="newlogistic")
library(sdmTMB)
library(here)
library(mvtnorm)
library(mgcv)
library(dplyr)
library(devtools)
library(zoo)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(MASS)
#install.packages("ggpubr")
library(ggpubr)
library(scales)
library(visreg)
library(ggeffects)
library(stringr)
library(TMB)
#install.packages("tweedie")
library(tweedie)
#install.packages("ggridges")
library(ggridges)
library(viridis)

### Set ggplot themes ###
theme_set(theme_bw(base_size = 30))
theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### load helper functions ####
source("Code/util_funs.R")
source("Code/sim_funs.R")

# Make list of parameter names"
pars_names <- c("log_depth_scaled", "log_depth_scaled2", "mi-delta", "mi-s50", "mi-smax", "range", "sigma_O", "phi", "tweedie_p", "mi-Eo", "as.factor(year)2010","as.factor(year)2011","as.factor(year)2012", "as.factor(year)2013", "as.factor(year)2014", "as.factor(year)2015")
true_pars <- data.frame(term=c("log_depth_scaled", "log_depth_scaled2", "mi-delta", "mi-s50", "mi-smax", "range", "sigma_O", "phi", "tweedie_p", "mi-Eo", "as.factor(year)2010","as.factor(year)2011","as.factor(year)2012", "as.factor(year)2013", "as.factor(year)2014", "as.factor(year)2015"), 
                        estimate=c(beta1, beta2, delta, s50, smax, range, sigma_O, phi, p, Eo, b_years))
true_pars2 <- data.frame(term=c("log_depth_scaled", "log_depth_scaled2", "mi-delta", "mi-s50", "mi-smax", "range", "sigma_O", "phi", "tweedie_p", "mi-Eo", "as.factor(year)2010","as.factor(year)2011","as.factor(year)2012", "as.factor(year)2013", "as.factor(year)2014", "as.factor(year)2015"), 
                         estimate=c(beta1, beta2, delta, s50, smax, range, sigma_O, phi, p, Eo2, b_years))
# Set model names #
model_names <- c("Typical Case, Unconstrained", "Typical Case, Prior Constrained", "Unusual Case, Unconstrained", "Unusual Case, Prior Constrained")

## Apply to model fits ##
## Parameter estimates ##
pars1 <- lapply(fits, extract_pars)
pars2 <- lapply(fits2, extract_pars)
pars3 <-lapply(fits3, extract_pars)
pars4 <- lapply(fits4, extract_pars)

pars1 <- clean_pars(pars1, fits=fits)
pars2 <- clean_pars(pars2, fits=fits2)
pars3 <- clean_pars(pars3, fits=fits3)
pars4 <- clean_pars(pars4, fits=fits4)

#Add column to label model
pars1$model <- model_names[1]
pars1$data <- "Typical Case"
pars1$analysis <- "Unconstrained"
pars2$model <- model_names[2]
pars2$data <- "Typical Case"
pars2$analysis <- "Prior Constrained"
pars3$model <- model_names[3]
pars3$data <- "Unusual Case"
pars3$analysis <- "Unconstrained"
pars4$model <- model_names[4]
pars4$data <- "Unusual Case"
pars4$analysis <- "Prior Constrained"

#Merge into one and combine
pars <- rbind(pars1,pars2, pars3, pars4)

### Parameter performance measures ###
## Average ##
avg <- aggregate(estimate ~ term+model, pars, FUN=mean)

# Average standard
sd_avg <- aggregate(estimate ~ term+model, pars, FUN=stats::sd)

## Root mean square error (accuracy) ##
# Calculate error #
pars$error <- case_when(pars$term=="mi-s50"~pars$estimate-s50,
                        pars$term=="mi-delta"~pars$estimate-delta,
                        pars$term=="mi-smax"~pars$estimate-smax,
                        pars$term=="mi-Eo"& pars$model=="Typical Case, Unconstrained"~pars$estimate-Eo,
                        pars$term=="mi-Eo"& pars$model=="Typical Case, Prior Constrained"~pars$estimate-Eo,
                        pars$term=="mi-Eo"& pars$model=="Unusual Case, Unconstrained"~pars$estimate-Eo2,
                        pars$term=="mi-Eo"& pars$model=="Unusual Case, Prior Constrained"~pars$estimate-Eo2,
                        pars$term=="log_depth_scaled"~pars$estimate-beta1,
                        pars$term=="log_depth_scaled2"~pars$estimate-beta2,
                        pars$term=="range"~pars$estimate-range,
                        pars$term=="sigma_O"~pars$estimate-sigma_O,
                        pars$term=="phi"~pars$estimate-phi,
                        pars$term=="tweedie_p"~pars$estimate-p,
                        pars$term=="as.factor(year)2010"~pars$estimate-b_years[1],
                        pars$term=="as.factor(year)2011"~pars$estimate-b_years[2],
                        pars$term=="as.factor(year)2012"~pars$estimate-b_years[3],
                        pars$term=="as.factor(year)2013"~pars$estimate-b_years[4],
                        pars$term=="as.factor(year)2014"~pars$estimate-b_years[5],
                        pars$term=="as.factor(year)2015"~pars$estimate-b_years[6])
pars$error2 <- pars$error^2
rmse <- aggregate(error2 ~ term+model, pars, FUN=sum)
rmse2 <- aggregate(error2 ~ term+model, pars, FUN=length)
rmse$n <- rmse2$error2
rmse$rmse <- sqrt(rmse$error2/rmse$n)
rmse$n <- NULL
rmse$error2 <- NULL

#Create dataframe for plotting

Eo_values <- as.data.frame(matrix(nrow=4))
Eo_values$V1 <- NULL
Eo_values$data <- c("Typical Case", "Typical Case", "Unusual Case", "Unusual Case")
Eo_values$analysis <- c( "Prior Constrained", "Unconstrained", "Prior Constrained", "Unconstrained")
Eo_values$model <- c("Typical Case, Prior Constrained", "Typical Case, Unconstrained","Unusual Case, Prior Constrained", "Unusual Case, Unconstrained")
MLE_avg <- aggregate(estimate~model, subset(pars, term=="mi-Eo"), FUN=mean)
Eo_values$MLE_avg <- MLE_avg$estimate
Eo_values$true <- c(Eo, Eo, Eo2, Eo2)

s50_values <- as.data.frame(matrix(nrow=4))
s50_values$V1 <- NULL
s50_values$data <- c("Typical Case", "Typical Case", "Unusual Case", "Unusual Case")
s50_values$analysis <- c( "Prior Constrained", "Unconstrained", "Prior Constrained", "Unconstrained")
s50_values$model <- c("Typical Case, Prior Constrained", "Typical Case, Unconstrained","Unusual Case, Prior Constrained", "Unusual Case, Unconstrained")
MLE_avg <- aggregate(estimate~model, subset(pars, term=="mi-s50"), FUN=mean)
s50_values$MLE_avg <- MLE_avg$estimate
s50_values$true <- c(s50, s50,s50,s50)

## Make a table ##
rmse$"average" <- avg$estimate
par_performance <- rmse
par_performance$sd <- sd_avg$estimate
colnames(par_performance) <- c("Parameter", "Model", "RMSE", "Average", "Precision")

par_performance$Bias <- case_when(par_performance$Parameter=="mi-s50"~par_performance$Average-s50,
                                  par_performance$Parameter=="mi-delta"~par_performance$Average-delta,
                                  par_performance$Parameter=="mi-smax"~par_performance$Average-smax,
                                  par_performance$Parameter=="mi-Eo"& par_performance$Model=="Typical Case, Unconstrained"~par_performance$Average-Eo,
                                  par_performance$Parameter=="mi-Eo"& par_performance$Model=="Typical Case, Prior Constrained"~par_performance$Average-Eo,
                                  par_performance$Parameter=="mi-Eo"& par_performance$Model=="Unusual Case, Unconstrained"~par_performance$Average-Eo2,
                                  par_performance$Parameter=="mi-Eo"& par_performance$Model=="Unusual Case, Prior Constrained"~par_performance$Average-Eo2,
                                  par_performance$Parameter=="log_depth_scaled"~par_performance$Average-beta1,
                                  par_performance$Parameter=="log_depth_scaled2"~par_performance$Average-beta2,
                                  par_performance$Parameter=="range"~par_performance$Average-range,
                                  par_performance$Parameter=="sigma_O"~par_performance$Average-sigma_O,
                                  par_performance$Parameter=="phi"~par_performance$Average-phi,
                                  par_performance$Parameter=="tweedie_p"~par_performance$Average-p,
                                  par_performance$Parameter=="as.factor(year)2010"~par_performance$Average-b_years[1],
                                  par_performance$Parameter=="as.factor(year)2011"~par_performance$Average-b_years[2],
                                  par_performance$Parameter=="as.factor(year)2012"~par_performance$Average-b_years[3],
                                  par_performance$Parameter=="as.factor(year)2013"~par_performance$Average-b_years[4],
                                  par_performance$Parameter=="as.factor(year)2014"~par_performance$Average-b_years[5],
                                  par_performance$Parameter=="as.factor(year)2015"~par_performance$Average-b_years[6])

Eo_performance <- subset(par_performance, Parameter=="mi-Eo")
s50_performance <- subset(par_performance, Parameter=="mi-s50")

### Plot parameter estimates ###
# Reorder
pars$analysis <- factor(pars$analysis, levels = c("Unconstrained", "Prior Constrained"))
Eo_values$analysis <- factor(Eo_values$analysis, levels = c("Unconstrained", "Prior Constrained"))
#Density plot just Eo #
ggplot(subset(pars, pars$term=="mi-Eo"), aes(x=estimate)) +
  geom_density(fill="lightblue", adjust = 1.5) +
  geom_vline(data = Eo_values, aes(xintercept = MLE_avg),linetype="dashed", size=1.2, color="darkorange", show.legend=T)+
  geom_vline(data = Eo_values, aes(xintercept = true),linetype="dashed", size=1.2)+
  facet_grid(analysis~data)+
  xlab("Eo estimate") + 
  theme(strip.text = element_text(size = 14))

# Density plot s50 #
#Reorder 
s50_values$analysis <- factor(s50_values$analysis, levels = c("Unconstrained", "Prior Constrained"))
ggplot(subset(pars, pars$term=="mi-s50"), aes(x=estimate)) +
  geom_density(fill="lightblue", adjust = 1.5) +
  geom_vline(data = s50_values, aes(xintercept = MLE_avg),linetype="dashed", size=1.2, color="darkorange", show.legend=T)+
  geom_vline(data = s50_values, aes(xintercept = true),linetype="dashed", size=1.2)+
  facet_grid(analysis~data)+
  xlab("s50 estimate")
  #To do: add RMSE, precision, accuracy as text boxes directly to plot

# Check cumulative RMSE # 
ggplot(subset(pars, pars$term=="mi-Eo"), aes(x=estimate))+
  stat_ecdf(geom = "step")+
  facet_grid(analysis~data)

# Alternative plot: dot and whiskers #
ggplot(subset(pars, pars$term=="mi-Eo"), aes(x=id)) +
  geom_errorbar(aes(ymin = conf.low, ymax =conf.high)) +  
  geom_point(aes(y=estimate))+
  facet_grid(analysis~data)+
  ylab("Eo estimate")+
  xlab("Simulation")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  geom_hline(data = Eo_values, aes(yintercept = MLE_avg),linetype="dashed", size=1.2, color="darkorange", show.legend=T)+
  geom_hline(data = Eo_values, aes(yintercept = true),linetype="dashed", size=1.2)
  
ggplot(subset(pars, pars$term=="mi-s50"), aes(x=id)) +
  geom_errorbar(aes(ymin = conf.low, ymax =conf.high)) +  
  geom_point(aes(y=estimate))+
  facet_grid(analysis~data)+
  ylab("s50 estimate") + 
  geom_density(fill="lightblue", adjust = 1.5) +
  geom_vline(data = s50_values, aes(xintercept = MLE_avg),linetype="dashed", size=1.2, color="darkorange", show.legend=T)+
  geom_vline(data = s50_values, aes(xintercept = true),linetype="dashed", size=1.2)+
  facet_grid(analysis~data)+
  xlab("s50 estimate") + 
  theme(strip.text = element_text(size = 14))

# To Do: Supplemental figures: other parameters # 

#Plot accuracy and precision
ggplot(data=Eo_performance, aes(x=Precision, y=Bias))+geom_point(aes(group=Model,color=Model), size=5)

#### Comparing po2' effect from parameter estimates ####
#calculate the po2' estimated from the Eo, use that to calculate f(po2')
simdats1 <- mapply(FUN=calculate_po2_prime, simdat,fits, SIMPLIFY=F)
simdats2 <- mapply(FUN=calculate_po2_prime, simdat2,fits3, SIMPLIFY=F)

#plot true po2 vs exp(f(po2'estimated) (either as points or lines)
simdats1 <- mapply(FUN=logfun2, simdats1,"po2_prime", fits, SIMPLIFY=F)
simdats2 <- mapply(FUN=logfun2, simdats2,"po2_prime", fits3, SIMPLIFY=F)

#Calculate true po2' effect
true_effect <- as.data.frame(logfun_basic(dat$mi_usual, smax, s50, delta))
colnames(true_effect) <- "mi_effect"
true_effect$mi <- dat$mi_usual

##Remove datasets with NA and then bind together, adding a column with the data frame number ##
simdats1 <- keep(.x=simdats1, .p=is.data.frame)
simdats1 <- bind_rows(simdats1, .id="df")
simdats2 <- keep(.x=simdats2, .p=is.data.frame)
simdats2 <- bind_rows(simdats2, .id="df")

ggplot(simdats1, aes(mi_usual, logmu, colour=Eo)) +
  geom_point(size=0.1)+
  scale_colour_viridis()+
  xlab("True pO2'")+
  ylab("Estimated pO2' effect")+
  theme(legend.position=c(0.9,0.2))+
  geom_line(data=true_effect, aes(x=mi, y=mi_effect), color="black", linetype="dashed", size=2)
  
#### Covariance of parameters
## Make dataframe of all parameter estimates wide ##
pars_wide <- pivot_wider(pars, id_cols=c(id, model), names_from=term, values_from=estimate)
## Plot Eo vs s50 ##
ggplot(pars_wide, aes(x=pars_wide$"mi-Eo", y=pars_wide$"mi-s50"))+ geom_point(aes(group=model, color=model), size=5)+xlab("Eo estimate")+ylab("s50 estimate")+
  theme(legend.position=c(0.3,0.8))