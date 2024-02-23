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
library(ggpubr)

### Set ggplot themes ###
theme_set(theme_bw(base_size = 25))
theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Make list of parameter names
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

### load helper functions  ###
source("Code/util_funs.R")
source("Code/sim_funs.R")

### Load data and models if needed ###
use_previous <- T
if(use_previous){
load("Model Outputs/model_fits.Rdata")
simdat <- readRDS("~/Dropbox/GitHub/estimating_mi_from_distribution2/Model Outputs/data_sims_usual.rds")
simdat2 <- readRDS("~/Dropbox/GitHub/estimating_mi_from_distribution2/Model Outputs/data_sims_weird.rds")
}

## Set how many data sets produced ##
n <- 250

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
pars1$data <- "Typical Species"
pars1$analysis <- "Unconstrained"
pars2$model <- model_names[2]
pars2$data <- "Typical Species"
pars2$analysis <- "Prior Constrained"
pars3$model <- model_names[3]
pars3$data <- "Unusual Species"
pars3$analysis <- "Unconstrained"
pars4$model <- model_names[4]
pars4$data <- "Unusual Species"
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

#Cumulative RMSE for determining if enough simulations
Eo_all <- subset(pars, term=="mi-s50")
Eo_all1 <- subset(Eo_all, model=="Typical Case, Unconstrained")

RMSE_sums <- matrix(nrow=250, ncol=1, NA)
for(i in 1:nrow(Eo_all1)){
  x <- Eo_all1[1:i,]
  total <- sum(x$error2)
  cum_rmse <- total/nrow(x)
  RMSE_sums[i] <- cum_rmse
}

plot(RMSE_sums, type="b")

#Create dataframe for plotting

Eo_values <- as.data.frame(matrix(nrow=4))
Eo_values$V1 <- NULL
Eo_values$data <- c("Typical Species", "Typical Species", "Unusual Species", "Unusual Species")
Eo_values$analysis <- c( "Prior Constrained", "Unconstrained", "Prior Constrained", "Unconstrained")
Eo_values$model <- c("Typical Case, Prior Constrained", "Typical Case, Unconstrained","Unusual Case, Prior Constrained", "Unusual Case, Unconstrained")
MLE_avg <- aggregate(estimate~model, subset(pars, term=="mi-Eo"), FUN=mean)
Eo_values$MLE_avg <- MLE_avg$estimate
Eo_values$true <- c(Eo, Eo, Eo2, Eo2)

s50_values <- as.data.frame(matrix(nrow=4))
s50_values$V1 <- NULL
s50_values$data <- c("Typical Species", "Typical Species", "Unusual Species", "Unusual Species")
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
## Figure 3 in paper (Unconstrained; Eo and s50 MLEs for typical & unusual)
# Reorder
pars$analysis <- factor(pars$analysis, levels = c("Unconstrained", "Prior Constrained"))
Eo_values$analysis <- factor(Eo_values$analysis, levels = c("Unconstrained", "Prior Constrained"))

## Plot just unconstrained
#Relabel species 
new_labels <- c("Typical Species" = "Typical species", "Unusual Species" = "Unusual species")
p1 <- ggplot(subset(pars, pars$term=="mi-Eo"&pars$analysis=="Unconstrained"), aes(x=estimate)) +
  geom_density(fill="lightblue", adjust = 1.5) +
  geom_vline(data=subset(Eo_values, Eo_values$analysis=="Unconstrained"), aes(xintercept = MLE_avg),linetype="dashed", size=1.2, color="darkorange", show.legend=T)+
  geom_vline(data=subset(Eo_values, Eo_values$analysis=="Unconstrained"), aes(xintercept = true),linetype="dashed", size=1.2)+
  facet_wrap("data", labeller=labeller(data = new_labels))+
  scale_y_continuous(limits=c(0,6))+
  scale_x_continuous(limits=c(-0.7,1.4))+
  xlab(bquote(E[0]~Maximum~Likelihood~Estimate))+
  theme(axis.title.y=element_blank(),
    axis.title.x=element_text(size=25),
    axis.text=element_text(size=20),
    strip.background = element_rect(fill="white", color="white"), 
    strip.text= element_text(face = "bold", size=25))
  

p2<- ggplot(subset(pars, pars$term=="mi-s50"&pars$analysis=="Unconstrained"), aes(x=estimate)) +
  geom_density(fill="lightblue", adjust = 1.5) +
  geom_vline(data=subset(s50_values, s50_values$analysis=="Unconstrained"), aes(xintercept = MLE_avg),linetype="dashed", size=1.2, color="darkorange", show.legend=T)+
  geom_vline(data=subset(s50_values, s50_values$analysis=="Unconstrained"), aes(xintercept = true),linetype="dashed", size=1.2)+
  facet_wrap("data")+
  xlab("s50 Threshold Maximum Likelihood Estimate")+
  theme(strip.text.x = element_blank(), axis.title.y=element_blank(),
        axis.title.x=element_text(size=25),
        axis.text=element_text(size=20))

figure <- ggarrange(p1, p2, labels=c("A", "B"),
                    ncol = 1, nrow = 2)

annotate_figure(figure, left=text_grob("Estimated probability density", size=25, rot=90))

ggsave(
  "figure_3_par_comparison.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8,
  height = 8,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE
)
## Figure 4 in paper (Prior constrained)
p1<- ggplot(subset(pars, pars$term=="mi-Eo"&pars$analysis=="Prior Constrained"), aes(x=estimate)) +
  geom_density(fill="lightblue", adjust = 1.5) +
  geom_vline(data=subset(Eo_values, Eo_values$analysis=="Prior Constrained"), aes(xintercept = MLE_avg),linetype="dashed", size=1.2, color="darkorange", show.legend=T)+
  geom_vline(data=subset(Eo_values, Eo_values$analysis=="Prior Constrained"), aes(xintercept = true),linetype="dashed", size=1.2)+
  facet_wrap("data", labeller=labeller(data = new_labels))+
  scale_y_continuous(limits=c(0,6))+
  scale_x_continuous(limits=c(-0.7,1.4))+
  stat_function(fun = dnorm, n = n, args = list(mean = 0.3477, sd = 0.1455), linetype="dashed", geom="area", alpha=0.2)+
  stat_function(fun = dnorm, n = n, args = list(mean = 0.3477, sd = 0.1455), linetype="dashed")+
  xlab(bquote(E[0]~Maximum~Likelihood~Estimate))+
  theme(axis.title.y=element_blank(),
        axis.title.x=element_text(size=25),
        axis.text=element_text(size=20),
        strip.background = element_rect(fill="white", color="white"), 
        strip.text= element_text(face = "bold", size=25))

p2<- ggplot(subset(pars, pars$term=="mi-s50"&pars$analysis=="Prior Constrained"), aes(x=estimate)) +
  geom_density(fill="lightblue", adjust = 1.5) +
  geom_vline(data=subset(s50_values, s50_values$analysis=="Prior Constrained"), aes(xintercept = MLE_avg),linetype="dashed", size=1.2, color="darkorange", show.legend=T)+
  geom_vline(data=subset(s50_values, s50_values$analysis=="Prior Constrained"), aes(xintercept = true),linetype="dashed", size=1.2)+
  facet_wrap("data")+
  xlab("s50 Threshold Maximum Likelihood Estimate")+
  scale_x_continuous(limits=c(1,3))+
  theme(strip.text.x = element_blank(), axis.title.y=element_blank(),
        axis.title.x=element_text(size=25),
        axis.text=element_text(size=20))

figure <- ggarrange(p1, p2, labels=c("A", "B"),
                    ncol = 1, nrow = 2)

annotate_figure(figure, left=text_grob("Estimated probability density", size=25, rot=90))

ggsave(
  "figure_4_par_comparison_prior.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8,
  height = 8,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE
)


#### Comparing po2' effect from parameter estimates ####
#calculate the po2' estimated from the Eo, use that to calculate f(po2')
simdats1 <- mapply(FUN=calculate_po2_prime, simdat,fits, SIMPLIFY=F)
simdats2 <- mapply(FUN=calculate_po2_prime, simdat,fits2, SIMPLIFY=F)
simdats3 <- mapply(FUN=calculate_po2_prime, simdat2,fits3, SIMPLIFY=F)
simdats4 <- mapply(FUN=calculate_po2_prime, simdat2,fits4, SIMPLIFY=F)

#plot true po2 vs exp(f(po2'estimated) (either as points or lines)
simdats1 <- mapply(FUN=logfun2, simdats1,"po2_prime", fits, SIMPLIFY=F)
simdats2 <- mapply(FUN=logfun2, simdats2,"po2_prime", fits2, SIMPLIFY=F)
simdats3 <- mapply(FUN=logfun2, simdats3,"po2_prime", fits3, SIMPLIFY=F)
simdats4 <- mapply(FUN=logfun2, simdats4,"po2_prime", fits4, SIMPLIFY=F)

##Remove datasets with NA and then bind together, adding a column with the data frame number ##
simdats1 <- keep(.x=simdats1, .p=is.data.frame)
simdats1 <- bind_rows(simdats1, .id="df")
simdats2 <- keep(.x=simdats2, .p=is.data.frame)
simdats2 <- bind_rows(simdats2, .id="df")
simdats3<- keep(.x=simdats3, .p=is.data.frame)
simdats3 <- bind_rows(simdats3, .id="df")
simdats4 <- keep(.x=simdats4, .p=is.data.frame)
simdats4 <- bind_rows(simdats4, .id="df")
simdats1$title <- "Typical Species"
simdats3$title <- "Unusual Species"
simdats3$side <- "Unconstrained"
simdats4$side <- "Prior Constrained"

#Calculate true po2' effect
dat <- subset(simdats1, df==1)
true_effect <- as.data.frame(logfun_basic(dat$mi_usual, smax, s50, delta))
true_effect2<- as.data.frame(logfun_basic(dat$mi_weird, smax, s50, delta))
colnames(true_effect) <- "mi_effect"
colnames(true_effect2) <- "mi_effect"
true_effect$mi <- dat$mi_usual
true_effect2$mi <- dat$mi_weird

###Calculate RMSE ####
simdats1$true <- true_effect$mi_effect
simdats2$true <- true_effect$mi_effect
simdats3$true <- true_effect2$mi_effect
simdats4$true <- true_effect2$mi_effect

simdats1$error <- simdats1$logmu-simdats1$true
simdats1$error2 <- simdats1$error^2
rmse <- aggregate(error2 ~ df, simdats1, FUN=sum)
rmse2 <- aggregate(error2 ~ df, simdats1, FUN=length)
rmse$n <- rmse2$error2
rmse$rmse <- sqrt(rmse$error2/rmse$n)
rmse$n <- NULL
rmse$error2 <- NULL
rmse1 <- mean(rmse$rmse)
rmse1sd <- sd(rmse$rmse)

simdats2$error <- simdats2$logmu-simdats2$true
simdats2$error2 <- simdats2$error^2
rmse <- aggregate(error2 ~ df, simdats2, FUN=sum)
rmse2 <- aggregate(error2 ~ df, simdats2, FUN=length)
rmse$n <- rmse2$error2
rmse$rmse <- sqrt(rmse$error2/rmse$n)
rmse$n <- NULL
rmse$error2 <- NULL
rmse2a <- mean(rmse$rmse)
rmse2sd <- sd(rmse$rmse)

simdats3$error <- simdats3$logmu-simdats3$true
simdats3$error2 <- simdats3$error^2
rmse <- aggregate(error2 ~ df, simdats3, FUN=sum)
rmse2 <- aggregate(error2 ~ df, simdats3, FUN=length)
rmse$n <- rmse2$error2
rmse$rmse <- sqrt(rmse$error2/rmse$n)
rmse$n <- NULL
rmse$error2 <- NULL
rmse3 <- mean(rmse$rmse)
rmse3sd <- sd(rmse$rmse)

simdats4$error <- simdats4$logmu-simdats4$true
simdats4$error2 <- simdats4$error^2
rmse <- aggregate(error2 ~ df, simdats4, FUN=sum)
rmse2 <- aggregate(error2 ~ df, simdats4, FUN=length)
rmse$n <- rmse2$error2
rmse$rmse <- sqrt(rmse$error2/rmse$n)
rmse$n <- NULL
rmse$error2 <- NULL
rmse4 <- mean(rmse$rmse)
rmse4sd <- sd(rmse$rmse)

#### Covariance of parameters ####
# Make dataframe of all parameter estimates wide ##
pars_wide <- pivot_wider(pars, id_cols=c(id, model), names_from=term, values_from=estimate)

#Just unconstrained and typical species for plotting
pars_wide2 <- subset(pars_wide, model=="Typical Case, Unconstrained")

#### Calculate for simulated data, combos of po2 and temp that are MI=s50 (i.e. 2) ####
dat$pO2_s50 <- 2/exp(0.3*dat$invtemp)


#### Plot Figure 2 in paper ####
p1 <- ggplot(pars_wide2, aes(x=pars_wide2$"mi-Eo", y=pars_wide2$"mi-s50"))+ 
            geom_point(size=5)+
          xlab(bquote(Estimated~E[0]))+
          ylab("Estimated s50")
 # theme(axis.title=element_text(size=25),
   #     axis.text=element_text(size=20))
  

p2 <- ggplot(simdats1, aes(mi_usual, logmu, colour=Eo)) +
  geom_point(size=0.1)+
  scale_colour_viridis()+
  xlab(expression("Observed"~phi[eco]))+
  ylab(expression("Estimated"~f(phi[eco])))+
  theme(legend.position=c(0.85,0.4))+
  geom_line(data=true_effect, aes(x=mi, y=mi_effect), color="black", linetype="dashed", size=2)+
  xlim(0,7)+
  annotate("text", label =expression("RMSE="~"0.035"%+-%"0.016"), size = 6, x = 5, y = 0.05)

p3 <- ggplot(dat, aes(x=temp, y=log(po2)))+
  theme(legend.position="none")+
  geom_ribbon(ymin=-Inf, aes(ymax=pO2_s50), fill='grey88')+
  geom_point()+
  geom_line(dat, mapping=aes(x=temp, y=pO2_s50), size=2)+
  xlab("Observed Temperature (C)")+
  ylab("Observed log Partial Pressure Oxygen (kPa)")+
  theme(axis.title.y=element_text(size=20))


### Impact on predictions: Evaluate impact of different combos on the limiting value of oxygen at each temperature, and then if temperatures increased
##Choose three representative parameter combos (highest and lowest Eo)
predict1 <- pars_wide2[which.max(pars_wide2$"mi-Eo"),]
predict2 <- subset(pars_wide2, pars_wide2$"mi-Eo">0)
predict2 <- pars_wide2[which.min(predict2$"mi-s50"),]

#Pull out parameter estimates
s50_1 <- as.numeric(predict1[,11])
Eo_1 <- as.numeric(predict1[,14])
delta_1 <- as.numeric(predict1[,12])
smax_1 <- as.numeric(predict1[,11])
s50_2 <- as.numeric(predict2[,11])
Eo_2 <- as.numeric(predict2[,14])
delta_2 <- as.numeric(predict2[,12])
smax_2 <- as.numeric(predict2[,11])

##Impact on fish response
#Calculate MI at two different temperatures
kelvin = 273.15
boltz = 0.000086173324
tref <- 12
invtemp1 <- (1 / boltz)  * ( 1 / (6 + 273.15) - 1 / (tref + 273.15))
invtemp2 <- (1 / boltz)  * ( 1 / (7.5 + 273.15) - 1 / (tref + 273.15))
#Oxygen
po2 <- seq(0,10.5,by=0.1)
mi1_high <- po2*exp(Eo_1* invtemp2)
mi2_high <- po2*exp(Eo_2* invtemp2)
mi_true1<- po2*exp(Eo* invtemp1)

#Response at each point
#Extract a high s50
logmu1_high <- logfun_basic(mi1_high, smax_1, s50_1, delta_1)
logmu2_high <- logfun_basic(mi2_high, smax_2, s50_2, delta_2)
logtrue1 <- logfun_basic(mi_true1, smax, s50, delta)

#Combine
preds <- as.data.frame(cbind(po2,mi1_high, mi2_high, logmu1_high, logmu2_high, logtrue1))
preds2 <- pivot_longer(preds, 4:6)

p4 <- ggplot(preds2, aes(x=po2, y=value, group=name,colour=name, linetype=name))+geom_line(size=1.5)+
  xlab("Partial Pressure Oxygen (kPa)")+
  ylab(expression("Estimated"~"f"(phi[eco])))+
  xlim(0,7)+
  scale_colour_manual(labels=c("+1.5 C, High Eo & s50", "+1.5 C, Low Eo & s50", "6 C Simulation"), values=c("gold", "purple3", "black"))+
  scale_linetype_manual(labels=c("+1.5 C, High Eo & s50", "+1.5 C, Low Eo & s50", "6 C Simulation"), values=c("solid", "solid", "dashed"))+
  theme(legend.position=c(0.6,0.2), legend.title=element_blank())

figure <- ggarrange(p1, p2, p3, p4, labels=c("A", "B", "C", "D"),
                    ncol = 2, nrow = 2)

ggsave(
  "figure_2_pars_effect_obs.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 15,
  height = 15,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE
)

