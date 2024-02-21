#### Plot fitted relationship ####
##### extract estimates ####
Eo <- getEo(model = m2a)

##### plot ####
## Calculate po2 prime from parameter ##
dat$po2_prime <- dat$po2 * exp(Eo * dat$invtemp)

##### plot ####

dat$Eo_constrained <- exp(logfun(dat$po2_prime, model = m2a, mi = T))
p1 <- ggplot(dat, aes(x = po2_prime, y = Eo_constrained, col = log(cpue_kg_km2 + 1))) +
  geom_jitter(width = 0.1, height = 0.01) +
  scale_colour_viridis_c(limits = c(0, 6), oob = scales::squish) +
  xlim(c(0,5))+
  xlab(expression("Observed"~phi[eco]))+
  ylab(expression("Estimated"~f(phi[eco])))


#### Plot marginal effects w/ +- SE using Tim's sdm_o2 code ####


### Method 1: Plot effect size ###
##Breakpoint##
m1_pars <- m1$sd_report$par.fixed
m1_sd <- m1$sd_report$sd
m1_parnames <- names(m1_pars)
b_threshold1 <- m1_pars[grep("b_threshold", m1_parnames)]
b_thresh_sd1<- m1_sd[3:4]
z <- 1.96

dat$thresh_pred <- exp(brkptfun(x = dat$po2_s, b_slope = b_threshold1[1], b_thresh = b_threshold1[2]))
dat$thresh_pred_high <- exp(brkptfun(x = dat$po2_s, b_slope = (b_threshold1[1]+z*b_thresh_sd1[1]), b_thresh = (b_threshold1[2]+z*b_thresh_sd1[2])))
dat$thresh_pred_low <- exp(brkptfun(x = dat$po2_s, b_slope = (b_threshold1[1]-z*b_thresh_sd1[1]), b_thresh = (b_threshold1[2]-z*b_thresh_sd1[2])))

# Plot
ggplot(dat, aes(x = po2, y = thresh_pred)) +
  geom_ribbon(aes(x=po2,y=thresh_pred, ymin = thresh_pred_low, ymax = thresh_pred_high), alpha=0.4)+
  geom_point(aes(colour = log(cpue_kg_km2 + 1)))+
  scale_colour_viridis_c(limits = c(0, 5), oob = scales::squish)+
  geom_line(aes(x=po2, y=thresh_pred))+
  xlim(c(0,5))+
  ylab("Marginal Effect")

##pO2' ##
Eo <- getEo(model = m2)
m2_pars <- m2$sd_report$par.fixed
m2_sd <- m2$sd_report$sd
m2_parnames <- names(m2_pars)
b_threshold2 <- m2_pars[grep("b_threshold", m2_parnames)]
b_thresh_sd2 <- m2_sd[3:6]
z <- 1.96
dat$po2_prime <- dat$po2 * exp(Eo * dat$invtemp)
dat$po2_prime_min <- dat$po2 * exp((Eo-b_thresh_sd2[4]) * dat$invtemp)
dat$po2_prime_max <- dat$po2 * exp((Eo+b_thresh_sd2[4]) * dat$invtemp)
dat$thresh_pred2 <- exp(logfun(dat$po2_prime, model = m2, mi = T))

#Get lower and upper bound#
dat$thresh_pred_high2 <- logfun_basic(dat$po2_prime_max, s50=(b_threshold2[1]+z*b_thresh_sd2[1]), s95=(b_threshold2[2]+z*b_thresh_sd2[2]), smax=(b_threshold2[2]+z*b_thresh_sd2[2]))
dat$thresh_pred_low2 <- logfun_basic(dat$po2_prime_min, s50=(b_threshold2[1]-z*b_thresh_sd2[1]), s95=(b_threshold2[2]-z*b_thresh_sd2[2]), smax=(b_threshold2[2]-z*b_thresh_sd2[2]))

ggplot(dat, aes(x = po2_prime, y = thresh_pred2)) +
  geom_ribbon(aes(x=po2,y=thresh_pred2, ymax = thresh_pred_low2, ymin = thresh_pred_high2), alpha=0.4)+
  geom_point(aes(colour = log(cpue_kg_km2 + 1)))+
  scale_colour_viridis_c(limits = c(0, 5), oob = scales::squish)+
  #geom_line(aes(x=po2, y=thresh_pred2))+
  xlim(c(0,5))+
  ylab("Marginal Effect")

ggplot(dat, aes(x = po2_prime, y = thresh_pred2)) +
  geom_ribbon(aes(y=thresh_pred2, xmin= po2_prime_min, ymin = po2_prime_max), alpha=0.4)+
  geom_point(aes(colour = log(cpue_kg_km2 + 1)))+
  scale_colour_viridis_c(limits = c(0, 5), oob = scales::squish)+
  #geom_line(aes(x=po2, y=thresh_pred2))+
  xlim(c(0,5))+
  ylab("Marginal Effect")



### Method 2: Plot marginal effect (in cpue kg km2) ###

## Breakpoint ##
# Make model output tidy #
tidy(m1,"ran_pars",conf.int = TRUE)
tidy(m1,"fixed", conf.int = TRUE)

# Create new data with everything set to 0 except for pO2 #
nd_po2 <- data.frame(po2_s = seq(min(dat$po2_s), max(dat$po2_s), length.out = 300), 
                     temp = 0,
                     log_depth_scaled = 0,
                     log_depth_scaled2 = 0,
                     year = as.factor(2010L)
)
nd_po2 <- convert_class(nd_po2)

# predict to new data #
p_o2 <- predict(m1, newdata = nd_po2, se_fit = TRUE, re_form = NA)

# plot predictions with uncertainty
z <- 1.96 # for 90% CI
p_o2$po2 <-back.convert(p_o2$po2_s, attr(dat$po2_s, "scaled:center"), attr(dat$po2_s, "scaled:scale"))
ggplot(p_o2, aes(x=po2, y=exp(est), 
                 ymin = exp(est - z *est_se), ymax = exp(est + z * est_se))) +
  geom_line() + geom_ribbon(alpha = 0.4) + 
  scale_x_continuous(limits = c(0, 3), expand = expansion(mult = c(0, 0.0))) +
  labs(x = "Partial Pressure of Oxygen (kPa)", y = bquote('Population Density'~(kg~km^-2)))


## Constrained Eo ##
tidy(m2,"ran_pars",conf.int = TRUE)
tidy(m2,"fixed", conf.int = TRUE)

# predict on full dataframe # 
nd_po2 <- data.frame(po2 = seq(min(dat$po2), max(dat$po2), length.out = 300), 
                     invtemp = mean(dat$invtemp),
                     log_depth_scaled = 0,
                     log_depth_scaled2 = 0,
                     year = as.factor(2010L)
)
nd_po2 <- convert_class(nd_po2)
p_o2 <- predict(m2, newdata = nd_po2, se_fit = TRUE, re_form = NA)

#Calculate po2 prime #
Eo <- getEo(model = m2)
p_o2$prime <- p_o2$po2 * exp(Eo * p_o2$invtemp)

# Plot--ends up super squiggly #
ggplot(p_o2, aes(x=prime, y=exp(est),
                 ymin = exp(est - z *est_se), ymax = exp(est + z * est_se))) +
  geom_line() + geom_ribbon(alpha = 0.4) + 
  scale_y_continuous(limits = c(0, 75), expand = expansion(mult = c(0, 0.0))) +
  scale_x_continuous(limits = c(0, 1.5), expand = expansion(mult = c(0, 0.0))) +
  labs(x = "Partial Pressure of Oxygen (kPa)", y = bquote('Population Density'~(kg~km^-2)))

#Try to manually thin it to create a smoother line--this still makes a squiggly plot #
p_o2 <- arrange(p_o2, p_o2$prime)
ggplot(p_o2[seq(1, nrow(p_o2), 10), ], aes(x=prime, y=exp(est),
                                           ymin = exp(est - z *est_se), ymax = exp(est + z * est_se))) +
  geom_line() + geom_ribbon(alpha = 0.4) + 
  scale_y_continuous(limits = c(0, 75), expand = expansion(mult = c(0, 0.0))) +
  scale_x_continuous(limits = c(0, 1.5), expand = expansion(mult = c(0, 0.0))) +
  labs(x = "Partial Pressure of Oxygen (kPa)", y = bquote('Population Density'~(kg~km^-2)))

##Option 2: creating new dataframe with smoothed po2 and invtemp #
# Back-calculate pO2 prime, po2, and invtemp values #
Eo <- getEo(model = m2)
dat$po2_prime <- dat$po2 * exp(Eo * dat$invtemp)
prime_preds <- data.frame(prime =seq(min(dat$po2_prime), max(dat$po2_prime), length.out=1000))
prime_preds$po2 <- seq(min(dat$po2), max(dat$po2), length.out=1000)
prime_preds$invtemp <- log(prime_preds$prime/prime_preds$po2)/Eo
prime_preds$test <- prime_preds$po2 * exp(Eo * prime_preds$invtemp)

nd_po2 <- data.frame(po2 = prime_preds$po2, 
                     invtemp = prime_preds$invtemp,
                     log_depth_scaled = 0,
                     log_depth_scaled2 = 0,
                     year = as.factor(2010L)
)
nd_po2 <- convert_class(nd_po2)

# predict to new data--this Rbombs#
p_o2 <- predict(m2, newdata = nd_po2, se_fit = TRUE, re_form = NA)

ggplot(p_o2, aes(x=po2, y=exp(est),ymin = exp(est - z *est_se), ymax = exp(est + z * est_se))) +
  geom_line()+geom_ribbon(alpha = 0.4) + 
  scale_y_continuous(limits = c(0, 75), expand = expansion(mult = c(0, 0.0))) +
  scale_x_continuous(limits = c(0, 1.5), expand = expansion(mult = c(0, 0.0))) +
  labs(x = "Partial Pressure of Oxygen (kPa)", y = bquote('Population Density'~(kg~km^-2)))

## Plot logistic ##
# Make model output tidy #
tidy(m3,"ran_pars",conf.int = TRUE)
tidy(m3,"fixed", conf.int = TRUE)

# Create new data with everything set to 0 except for pO2 #
nd_po2 <- data.frame(po2_s = seq(min(dat$po2_s), max(dat$po2_s), length.out = 300), 
                     temp = 0,
                     log_depth_scaled = 0,
                     log_depth_scaled2 = 0,
                     year = as.factor(2010L)
)
nd_po2 <- convert_class(nd_po2)

# predict to new data #
p_o2 <- predict(m3, newdata = nd_po2, se_fit = TRUE, re_form = NA)

# plot predictions with uncertainty
z <- 1.645 # for 90% CI
p_o2$po2 <-back.convert(p_o2$po2_s, attr(dat$po2_s, "scaled:center"), attr(dat$po2_s, "scaled:scale"))
plot_logistic <- ggplot(p_o2, aes(x=po2, y=exp(est), 
                                  ymin = exp(est - z *est_se), ymax = exp(est + z * est_se))) +
  geom_line() + geom_ribbon(alpha = 0.4) + 
  scale_y_continuous(limits = c(0, 250), expand = expansion(mult = c(0, 0.0))) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.0))) +
  labs(x = "Partial Pressure of Oxygen (kPa)", y = bquote('Population Density'~(kg~km^-2)))

plot_logistic

multipanel <- gridExtra::grid.arrange(breakpoint, logistic, Eo_constrained, nrow = 3, ncol = 1)
ggsave("plots/marginal_effects.png", multipanel, width = 6, height = 3, units = "in", device = "png")

##pO2 and Temp data and po2/temp curve
#Calculate for simulated data, combos of po2 and temp that are MI=s50 (i.e. 2)
dat$pO2_s50 <-  0.88/exp(Eo*dat$invtemp)

ggplot(dat, aes(x=temp, y=log(po2)))+
  #stat_density_2d_filled(alpha=0.8)+
  #scale_fill_distiller(palette= "Spectral", direction=-1)+
  theme(legend.position="none")+
  geom_ribbon(ymin=-Inf, aes(ymax=pO2_s50), fill='grey88')+
  #geom_ribbon(aes(ymin=pO2_s50), ymax=Inf, fill='grey88', alpha=0.4)+
  #geom_hex()+
  # scale_fill_grey()+
  geom_point(aes(colour=log(cpue_kg_km2+1)))+
  scale_colour_viridis_c(limits = c(0, 6), oob = scales::squish) +
  geom_line(dat, mapping=aes(x=temp, y=pO2_s50), size=2, colour="black")+
  xlab("Observed Temperature (C)")+
  ylab("Observed log Dissolved Oxygen (kPa)")




##### Older versions, ignore ####
### Plot depth effects  ####
## Set ggplot theme ##
theme_set(theme_bw(base_size = 30))
theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
## Internal Eo mode ##
# Extract parameters another way for internal Eo model #
par_estimates <- as.data.frame(tidy(m2, conf.int = TRUE, effects="fixed"))
par_estimates_rand <- as.data.frame(tidy(m2, conf.int = TRUE, effects="ran_pars"))
par_estimates <- bind_rows(par_estimates, par_estimates_rand)
pars <- par_estimates[1:2]
pars <- pivot_wider(pars, names_from=term, values_from=estimate)
# Calculate temp-corrected po2 from Eo value #
dat$mi_pred <- dat$po2*exp(pars$"mi-Eo"* dat$invtemp)
# Calculate effect of MI # Tim note: probably better to use the logfun function above to calculate the log effect size from po2prime or po2.  
dat$mi_effect <- pars$"mi-smax" * (1 / (1 + exp(-log(19) * (dat$mi_pred - pars$`mi-s50`) / pars$"mi-delta")) - 1) # this is exponentiated. 
# Calculate combined depth effect #
dat$depth_effect_combined <- (dat$log_depth_scaled*pars$log_depth_scaled)+(dat$log_depth_scaled2*pars$log_depth_scaled2)

## For breakpt(o2) model ##
par_estimates2 <- as.data.frame(tidy(m1, conf.int = TRUE, effects="fixed"))
par_estimates_rand2 <- as.data.frame(tidy(m1, conf.int = TRUE, effects="ran_pars"))
par_estimates2 <- bind_rows(par_estimates2, par_estimates_rand2)
pars2 <- par_estimates2[1:2]
pars2 <- pivot_wider(pars2, names_from=term, values_from=estimate)
# alculate depth effects #
dat$depth_effect_combined2 <- (dat$log_depth_scaled*pars2$log_depth_scaled)+(dat$log_depth_scaled2*pars2$log_depth_scaled2)
# Calculate po2 effect #
dat$po2_effect <- ifelse(dat$po2_s>pars2$`po2_s-breakpt`,pars2$`po2_s-breakpt`*pars2$`po2_s-slope`,dat$`po2_s`*pars2$`po2_s-slope`)

# Plot depth effects compared between model #
ggplot(dat, aes(y=depth_effect_combined, x=depth))+geom_line(size=1.3)+geom_line(dat, mapping=aes(y=depth_effect_combined2, x=depth), color="red", size=1.3)+ylab("Depth Effect Combined")
ggplot(subset(dat, depth>450), aes(y=depth_effect_combined, x=depth))+geom_point()+geom_point(subset(dat, depth>450), mapping=aes(y=depth_effect_combined2, x=depth), color="red")+ylab("Depth Effect Combined")

# Plot metabolic index vs depth #
ggplot(dat, aes(y=depth, x=mi_pred))+geom_point(aes(color=log(cpue_kg_km2)))+xlab("pO2'")+theme(legend.position=c(0.8, 0.8))

#Plot oxygen vs depth
ggplot(dat, aes(y=-depth, x=po2))+geom_point(aes(color=log(cpue_kg_km2)))+xlab("pO2")+theme(legend.position=c(0.7, 0.3))+scale_colour_viridis_c()

