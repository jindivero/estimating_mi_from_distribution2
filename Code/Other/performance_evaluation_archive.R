#### Cross-validation with model predictions ####
### Simulate new cross-validation data ###
simdat_cv <- map(seq_len(n), ~simulate_fish(dat = dat,
                                            mesh = mesh,
                                            s50 = s50,
                                            delta = delta,
                                            smax = smax,
                                            Eo = Eo,
                                            modelpars = model.pars))

simdat2_cv <- map(seq_len(n), ~simulate_fish(dat = dat,
                                             mesh = mesh,
                                             s50 = s50,
                                             delta = delta,
                                             smax = smax,
                                             Eo = Eo2,
                                             modelpars = model.pars))

### Predictions ###
## Function to predict from each model and each dataset ##

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

#### Comparing po2' and predictions for all simulations (to see how predictions compare to parameters) ####
#calculate the po2' estimated from the Eo, use that to calculate f(po2')
simdats1 <- mapply(FUN=calculate_po2_prime, simdat,fits, SIMPLIFY=F)
simdats2 <- mapply(FUN=calculate_po2_prime, simdat2,fits3, SIMPLIFY=F)

#plot true po2 vs exp(f(po2'estimated) (either as points or lines)
simdats1 <- mapply(FUN=logfun2, simdats1,"po2_prime", fits, SIMPLIFY=F)
simdats2 <- mapply(FUN=logfun2, simdats2,"po2_prime", fits3, SIMPLIFY=F)

##Remove datasets with NA and then bind together, adding a column with the data frame number ##
simdats1 <- keep(.x=simdats1, .p=is.data.frame)
simdats1 <- bind_rows(simdats1, .id="df")
simdats2 <- keep(.x=simdats2, .p=is.data.frame)
simdats2 <- bind_rows(simdats2, .id="df")

ggplot(simdats1, aes(mi_usual, logmu, colour=as.factor(df))) +
  geom_point(size=0.1)+
  scale_colour_viridis(discrete=TRUE)+
  xlab("True pO2'")+
  ylab("Estimated pO2' effect")+
  theme(legend.position="none")

ggplot(simdats1, aes(mi_usual, logmu, colour=Eo)) +
  geom_point(size=0.1)+
  scale_colour_viridis()+
  xlab("True pO2'")+
  ylab("Estimated pO2' effect")

ggplot(simdats1, aes(mi_usual, logmu, colour=Eo)) +
  geom_point(size=0.1)+
  scale_colour_viridis()+
  xlab("True pO2'")+
  ylab("Estimated pO2' effect")

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

#### Covariance of parameters
# Make dataframe of all parameter estimates wide #
pars_wide <- pivot_wider(pars, id_cols=c(id, model), names_from=term, values_from=estimate)
### Plot Eo vs s50 ###
## Plot all together ##
ggplot(pars_wide, aes(x=pars_wide$"mi-Eo", y=pars_wide$"mi-s50"))+ geom_point(aes(group=model, color=model), size=5)+xlab("Eo estimate")+ylab("s50 estimate")+
  theme(legend.position=c(0.3,0.8))
# Plot all together in facet grid # 
ggplot(pars_wide, aes(x=pars_wide$"mi-Eo", y=pars_wide$"mi-s50"))+ geom_point(size=5)+xlab("Eo estimate")+facet_wrap("model", scales="free")+ylab("s50 estimate")

# Subset just one model # 
pars_wide1 <- subset(pars_wide, pars_wide$model=="Typical Case, Unconstrained")
ggplot(subset(pars_wide1), aes(x=pars_wide1$"mi-Eo", y=pars_wide1$"mi-s50"))+ geom_point(size=5)

ggplot(pars_wide, aes(x=pars_wide$"mi-Eo", y=pars_wide$"mi-smax"))+ geom_point(aes(group=model, color=model), size=5)+xlab("Eo estimate")+ylab("s50 estimate")+
  theme(legend.position=c(0.2,0.8))
# Apply #
simdats1 <- mapply(FUN=calculate_po2_prime, simdat,fits, SIMPLIFY=F)

# Apply #
simdats1 <- mapply(FUN=logfun_all, simdats1, fits)
ggplot(bind_rows(simdats2, .id="df"), aes(po2_prime, logfun, colour=as.factor(df))) +
  geom_point(size=0.1)+
  scale_colour_viridis(discrete=TRUE)+
  xlab("pO2'")+
  ylab("Predicted density")


# Plot #
simdats1 <- keep(.x=simdats1, .p=is.data.frame)
simdats2 <- simdats1[1:5]
ggplot(bind_rows(simdats2, .id="df"), aes(po2_prime, logfun, colour=as.factor(smax))) +
  geom_point(size=0.1)+
  scale_colour_viridis(discrete=TRUE)+
  xlab("pO2'")+
  ylab("pO2' effect")

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

### Impact on predictions: Evaluate impact of different combos on the limiting value of oxygen at each temperature, and then if temperatures increased
##Choose three representative parameter combos (highest and lowest Eo)
predict1 <- pars_wide2[which.max(pars_wide2$"mi-Eo"),]
predict2 <- pars_wide2[which.min(pars_wide2$"mi-Eo"),]
predict3 <- subset(pars_wide2, pars_wide2$"mi-Eo">0)
predict3 <- pars_wide2[which.min(predict3$"mi-Eo"),]

##Impact on critical pO2 level
#Pull out the Eo and s50 parameters
s50_1 <- as.numeric(predict1[,c("mi-s50")])
Eo_1 <- as.numeric(predict1[,c("mi-Eo")])
delta_1 <- as.numeric(predict1[,c("mi-delta")])
smax_1 <- as.numeric(predict1[,c("mi-smax")])
s50_2 <- as.numeric(predict2[,c("mi-s50")])
Eo_2 <- as.numeric(predict2[,c("mi-Eo")])
delta_2 <- as.numeric(predict2[,c("mi-delta")])
smax_2 <- as.numeric(predict2[,c("mi-smax")])
s50_3 <- as.numeric(predict3[,c("mi-s50")])
Eo_3 <- as.numeric(predict3[,c("mi-Eo")])
delta_3 <- as.numeric(predict3[,c("mi-delta")])
smax_3 <- as.numeric(predict3[,c("mi-smax")])

#Calculate po2-crit for each (with observed data)
dat$po2_1 <- s50_1/exp(Eo_1*dat$invtemp)
dat$po2_2 <- s50_2/exp(Eo_2*dat$invtemp)
dat$po2_3 <- s50_3/exp(Eo_3*dat$invtemp)


#Increase the temperature 2 degrees
kelvin = 273.15
boltz = 0.000086173324
tref <- 12
dat$invtemp2 <- (1 / boltz)  * ( 1 / ((dat$temp+2) + 273.15) - 1 / (tref + 273.15))
#Re-calculate critical pO2 values
dat$po2_2a <- 2/exp(Eo_2*dat$invtemp2)
dat$po2_3a <- 2/exp(Eo_3*dat$invtemp2)
dat$pO2_s50_predict <- 2/exp(0.3*dat$invtemp2)
dat$po2_2b <- s50_2/exp(0.3*dat$invtemp2)
dat$po2_3b <- s50_3/exp(0.3*dat$invtemp2)

#Plot
p1 <- ggplot(dat, aes(x=temp+2, y=po2_1a))+
  #geom_line(color="grey", size=2)+
  geom_line(aes(x=temp+2, y=po2_2a),colour="gold",  size=2)+
  geom_line(aes(x=temp+2, y=po2_3a), colour="purple4",  size=2)+
  geom_line(aes(x=temp+2, y=pO2_s50_predict), linetype="dashed",  size=2)+
  geom_line(aes(x=temp+2, y=po2_2b),colour="yellow",  size=2)+
  geom_line(aes(x=temp+2, y=po2_3b), colour="purple",  size=2)+
  xlab("Temperature (C)")+
  ylab((expression("Oxygen (kPa) at"~phi["eco,s50"])))+
  annotate("text", label=bquote(Low~E[0]), size = 10, x =17, y = 1)+
  annotate("text", label=bquote(High~E[0]), size = 10, x = 17, y =1.5)

p1 <- ggplot(dat, aes(x=temp+2, y=po2_2a))+
  geom_line(color="purple", size=2)+
  geom_line(aes(x=temp+2, y=po2_3a), colour="purple4",  size=2)+
  geom_line(aes(x=temp+2, y=pO2_s50_predict), linetype="dashed",  size=2)+
  geom_line(aes(x=temp+2, y=po2_2b),colour="gold",  size=2)+
  geom_line(aes(x=temp+2, y=po2_3b), colour="gold3",  size=2)+
  xlab("Temperature (C)")+
  ylab((expression("Oxygen (kPa) at"~phi["eco,s50"])))+
  annotate("text", label=bquote(Low~E[0]), size = 10, x =17, y = 1)+
  annotate("text", label=bquote(High~E[0]), size = 10, x = 17, y =1.5)

ggplot(preds, aes(x=temp, y=logmu1b))+geom_line(colour="purple", size=1.5)+
  geom_line(mapping=aes(x=temp, logmu1c),colour="purple4", size=1.5)+
  geom_line(mapping=aes(x=temp, logmu2b),colour="gold",size=1.5)+
  geom_line(mapping=aes(x=temp, logmu2c),colour="gold3", size=1.5)+
  geom_line(aes(x=temp, y=logtrue), linetype="dashed",  size=2)+
  xlab("Temperature (C)")+
  
  ggsave(
    "figure_temp_po2crit.png",
    plot = last_plot(),
    device = NULL,
    path = NULL,
    scale = 1,
    width = 8,
    height =6,
    units = c("in"),
    dpi = 600,
    limitsize = TRUE
  )

#Only changing one parameter
dat$po2_2a <- 2/exp(Eo_2*dat$invtemp2)
dat$po2_3a <- 2/exp(Eo_3*dat$invtemp2)
dat$pO2_s50_predict <- 2/exp(0.3*dat$invtemp2)
dat$po2_2b <- 1.2/exp(0.3*dat$invtemp2)
dat$po2_3b <- 3/exp(0.3*dat$invtemp2)

##Impact on fish response
#Calculate MI at two different temperatures
invtemp1 <- (1 / boltz)  * ( 1 / (8 + 273.15) - 1 / (tref + 273.15))
invtemp2 <- (1 / boltz)  * ( 1 / (15 + 273.15) - 1 / (tref + 273.15))
#Oxygen
po2 <- seq(0,10.5,by=0.1)
mi1a <- po2*exp(Eo_1* invtemp1) 
mi1b<- po2*exp(Eo_2* invtemp1)
mi1c <- po2*exp(Eo_3* invtemp1)
mi2a <- po2*exp(Eo_1* invtemp2) 
mi2b<- po2*exp(Eo_2* invtemp2)
mi2c <- po2*exp(Eo_3* invtemp2)
#Temperature
temp <- seq(5,17.5,by=0.5)
invtemp <- (1 / boltz)  * ( 1 / (temp + 273.15) - 1 / (tref + 273.15))
mi1b<- 1.7*exp(Eo_2* invtemp)
mi1c <- 1.7*exp(Eo_3* invtemp)
mi_true <- 1.7*exp(0.3*invtemp)

#Response at each point
#Extract a high s50
predict3 <- subset(pars_wide2, pars_wide2$"mi-Eo">0)
high_s50 <- predict3[which.max(predict3$"mi-s50"),]
low_s50 <- pars_wide2[which.min(pars_wide2$"mi-s50"),]
logmu1b <- logfun_basic(mi1b, smax, s50, delta)
logmu1c <- logfun_basic(mi1c, smax, s50, delta)
logmu2b <- logfun_basic(mi_true, smax, 1.2, delta)
logmu2c <- logfun_basic(mi_true, smax, 3, delta)
logtrue <- logfun_basic(mi_true, smax, s50, delta)


#Combine
preds <- as.data.frame(cbind(po2, mi2a, mi2b, mi2c, mi1a, mi1b, mi1c, logmu2a, logmu2b, logmu2c, logmu1a, logmu1b, logmu1c))
preds <- as.data.frame(cbind(temp, mi2b, mi2c, mi1b, mi1c, logmu2b, logmu2c, logmu1b, logmu1c, logtrue))


#Plot
p2 <-ggplot(preds, aes(x=po2, y=logmu1b))+geom_line(colour="purple", size=1.5)+
  geom_line(mapping=aes(x=po2, logmu1c),colour="gold", size=1.5)+
  geom_line(mapping=aes(x=po2, logmu2b),colour="purple4",size=1.5)+
  geom_line(mapping=aes(x=po2, logmu2c),colour="gold3", size=1.5)+
  xlim(0,4.5)+
  xlab("Dissolved Oxygen (kPa)")+
  ylab(expression("Estimated"~f(phi[eco])))

#Only change one parameter
p2 <-ggplot(preds, aes(x=temp, y=logmu1b))+geom_line(colour="purple", size=1.5)+
  geom_line(mapping=aes(x=temp, logmu1c),colour="purple4", size=1.5)+
  geom_line(mapping=aes(x=temp, logmu2b),colour="gold",size=1.5)+
  geom_line(mapping=aes(x=temp, logmu2c),colour="gold3", size=1.5)+
  geom_line(aes(x=temp, y=logtrue), linetype="dashed",  size=2)+
  xlab("Temperature (C)")+
  ylab(expression("Estimated"~f(phi[eco])))
figure <- ggarrange(p1, p2,labels=c("A", "B"),
                    ncol = 2, nrow = 1)

#Keep pairs of parameter
logmu1b <- logfun_basic(mi1b, smax_2, s50_2, delta_2)
logmu1c <- logfun_basic(mi1c, smax_3, s50_2, delta_3)
logtrue <- logfun_basic(mi1c, smax, s50, delta)

preds <- as.data.frame(cbind(temp, mi2b, mi2c, mi1b, mi1c, logmu2b, logmu2c, logmu1b, logmu1c, logtrue))

ggplot(preds, aes(x=temp, y=logmu1b))+geom_line(colour="purple", size=1.5)+
  geom_line(mapping=aes(x=temp, logmu1c),colour="purple4", size=1.5)+
  geom_line(mapping=aes(x=temp, logmu2b),colour="gold",size=1.5)+
  geom_line(mapping=aes(x=temp, logmu2c),colour="gold3", size=1.5)+
  geom_line(aes(x=temp, y=logtrue), linetype="dashed",  size=2)+
  xlab("Temperature (C)")+
  ylab(expression("Estimated"~f(phi[eco])))









ggsave(
  "figure_4.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 30,
  height = 8,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE
)

###Response on fish density

#Make new oxygen and temperature data
dat$temp2 <- runif(3304,10,16)
dat$po22 <- runif(3304,0,10)
invtemp2 <- (1 / boltz)  * ( 1 / ((dat$temp2) + 273.15) - 1 / (tref + 273.15))

#Plot
ggplot(dat, aes(x=temp, y=po2, colour=mi1a))+geom_point()+geom_point(data=dat, mapping=aes(x=temp2, y=po22, colour=mi2a))

#Calculate metabolic index for original values and new values
mi2a <- dat$po22*exp(Eo_1* invtemp2) 
mi2b<- dat$po22*exp(Eo_2* invtemp2)
mi2c <- dat$po22*exp(Eo_3* invtemp2)
mi1a <- dat$po2*exp(Eo_1* invtemp)
mi1b <- dat$po2*exp(Eo_2*invtemp)
mi1c <- dat$po2*exp(Eo_3* invtemp)

#Calculate response at new temperature
logmu2a <- logfun_basic(mi2a, smax_1, s50_1, delta_1)
logmu2b <- logfun_basic(mi2b, smax_2, s50_2, delta_2)
logmu2c <- logfun_basic(mi2c, smax_3, s50_3, delta_3)
logmu1a <- logfun_basic(mi1a, smax_1, s50_1, delta_1)
logmu1b <- logfun_basic(mi1b, smax_2, s50_2, delta_2)
logmu1c <- logfun_basic(mi1c, smax_3, s50_3, delta_3)

#Plot
preds <- as.data.frame(cbind(invtemp, invtemp2, mi2a, mi2b, mi2c, mi1a, mi1b, mi1c, logmu2a, logmu2b, logmu2c, logmu1a, logmu1b, logmu1c))
ggplot(preds, aes(x=mi2a, y=logmu2a))+
  geom_point(colour="orange")+
  geom_point(aes(x=mi2b, y=logmu2b), colour="orange")+
  geom_point(aes(x=mi2c, y=logmu2c), colour="orange")

ggplot(preds, aes(x=mi1a, y=logmu1a))+
  geom_point(colour="darkblue")+
  geom_point(aes(x=mi1b, y=logmu1b),colour="darkblue")+
  geom_point(aes(x=mi1c, y=logmu1c), colour="darkblue")+
  xlim(0,12.5)

ggplot(preds, aes(x=mi2a1, y=logmu3a))+
  geom_line(colour="orange", linetype="dashed")+
  geom_line(aes(x=mi2b1, y=logmu3b), colour="orange", linetype="dotted")+
  geom_line(aes(x=mi2c1, y=logmu3c), colour="orange")+
  geom_line(aes(x=mi1a, y=logmu1a), colour="darkblue", linetype="dashed")+
  geom_line(aes(x=mi1b, y=logmu1b),colour="darkblue", linetype="dotted")+
  geom_line(aes(x=mi1c, y=logmu1c), colour="darkblue")


# How does uncertainty in s50 and Eo compare to uncertainty in other variables?

#What is the proportion of iterations that the true value falls within the 95% confidence interval?
##Coverage
#How often does true parameter value fall within the range of the 95% CI
#1 if true value is within confidence interval, 0 if not
pars$within_95 <- case_when(pars$term=="mi-s50"~ifelse(pars$conf.low <s50 & pars$conf.high >s50, 1, 0),
                            pars$term=="mi-delta"~ifelse(pars$conf.low <delta & pars$conf.high >delta, 1, 0),
                            pars$term=="mi-smax"~ifelse(pars$conf.low <smax & pars$conf.high >smax, 1, 0),
                            pars$term=="mi-Eo"& pars$model=="Typical Case, Unconstrained"~ifelse(pars$conf.low <Eo & pars$conf.high>Eo, 1, 0),
                            pars$term=="mi-Eo"& pars$model=="Unusual Case, Unconstrained"~ifelse(pars$conf.low <Eo2 & pars$conf.high>Eo2, 1, 0),
                            pars$term=="mi-Eo"& pars$model=="Typical Case, Prior Constrained"~ifelse(pars$conf.low <Eo & pars$conf.high>Eo, 1, 0),
                            pars$term=="mi-Eo"& pars$model=="Unusual Case, Prior Constrained"~ifelse(pars$conf.low <Eo2 & pars$conf.high>Eo2, 1, 0),
                            pars$term=="log_depth_sc"~ifelse(pars$conf.low <beta1 & pars$conf.high>beta1, 1, 0),
                            pars$term=="log_depth_sc2"~ifelse(pars$conf.low <beta2 & pars$conf.high>beta2, 1, 0),
                            pars$term=="range"~ifelse(pars$conf.low <range & pars$conf.high>range, 1, 0),
                            pars$term=="sigma_O"~ifelse(pars$conf.low <sigma_O & pars$conf.high>sigma_O, 1, 0),
                            pars$term=="phi"~ifelse(pars$conf.low <phi & pars$conf.high>phi, 1, 0),
                            pars$term=="tweedie_p"~ifelse(pars$conf.low <p & pars$conf.high>p, 1, 0))
#pars$term=="as.factor(year)2011"~ifelse(pars$conf.low <b_years[1] & pars$conf.high>b_years[1], 1, 0),
#pars$term=="as.factor(year)2012"~ifelse(pars$conf.low <b_years[2] & pars$conf.high>b_years[2], 1, 0),
#pars$term=="as.factor(year)2013"~ifelse(pars$conf.low <b_years[3] & pars$conf.high>b_years[3], 1, 0),
#pars$term=="as.factor(year)2014"~ifelse(pars$conf.low <b_years[4] & pars$conf.high>b_years[4], 1, 0),
#pars$term=="as.factor(year)2015"~ifelse(pars$conf.low <b_years[5] & pars$conf.high>b_years[5], 1, 0))

count_within_95 <- aggregate(within_95 ~ term+model, pars, FUN=sum)
count_within_95$proportion <- count_within_95$within_95/length(fits)

#Average 95% range
#95% range (upper 95% CI--lower 95% CI)
pars$con.range <- pars$conf.high-pars$conf.low
avg_conf.range <- aggregate(con.range ~ term+model, pars, FUN=mean)

#Average standard error of each estimate (how does looking at internal model uncertainty serve as proxy?)
sd_avg <- aggregate(std.error ~ term+model, pars, FUN=mean)

###What is the estimated response across a range of pO2 values, given a low, medium, and high temperature







###For the three extreme Eo examples, what is the range of uncertainty in predictions if the 95% confidence interval range is considered?
pars_test <- subset(pars1, pars1$id==20)

#Pull out lower and upper bounds
s50_1 <- as.numeric(pars_test[9,3])
Eo_1 <- as.numeric(pars_test[12,3])
delta_1 <- as.numeric(pars_test[9,3])
smax_1 <- as.numeric(pars_test[9,3])
s50_2 <- as.numeric(pars_test[9,5])
Eo_2 <-as.numeric(pars_test[12,5])
s50_3 <- as.numeric(pars_test[9,6])
Eo_3 <-as.numeric(pars_test[9,6])

#Calculate metabolic index for original values and new values
mi2a <- dat$po22*exp(Eo_1* invtemp2) 
mi2b<- dat$po22*exp(Eo_2* invtemp2)
mi2c <- dat$po22*exp(Eo_3* invtemp2)
mi1a <- dat$po2*exp(Eo_1* invtemp)
mi1b <- dat$po2*exp(Eo_2*invtemp)
mi1c <- dat$po2*exp(Eo_3* invtemp)

#Calculate response at original environmental conditions
#Low s50, regular Eo
logmu_s50a <- logfun_basic(mi1a, smax_1, s50_2, delta_1)
#High s50, regular
logmu_s50b <- logfun_basic(mi1a, smax_1, s50_3, delta_1)
#Regular s50, low Eo
logmu_Eoa <- logfun_basic(mi1b, smax_1, s50_1, delta_1)
#Regular s50, high Eo
logmu_Eob <- logfun_basic(mi1c, smax_1, s50_1, delta_1)

preds <- as.data.frame(cbind(mi1a, mi1b, mi1c, logmu_s50a, logmu_s50b, logmu_Eoa,logmu_Eob))

ggplot(preds, aes(x=mi1a, y=logmu_s50a))+geom_point(colour="blue")+
  geom_point(mapping=aes(x=mi1a, y=logmu_s50b), colour="blue")+
  geom_point(mapping=aes(x=mi1b, y=logmu_Eoa), colour="red")+
  geom_point(mapping=aes(x=mi1c, y=logmu_Eob), colour="red")


#Regular s50, high MI
logmu2b <- logfun_basic(mi2b, smax_2, s50_2, delta_2)
logmu2c <- logfun_basic(mi2c, smax_3, s50_3, delta_3)
logmu1a <- logfun_basic(mi1a, smax_1, s50_1, delta_1)
logmu1b <- logfun_basic(mi1b, smax_2, s50_2, delta_2)
logmu1c <- logfun_basic(mi1c, smax_3, s50_3, delta_3)



#Flip to make plotting easier
preds2 <- pivot_longer(preds,2:5)

#Plot
p2 <- ggplot(preds2, aes(x=po2, y=value))+geom_line(aes(colour=name), size=1.5)+
  xlab("Dissolved Oxygen (kPa)")+
  ylab(expression("Estimated"~f(phi[eco])))+
  scale_colour_manual(values=c("blue2", "skyblue2", "orange2", "gold"), labels=c("High Eo&s50, High Temp", "High Eo&s50, Low Temp", "Low Eo&s50, High Temp", "Low Eo&s50, Low Temp"))+
  geom_line(preds2, mapping=aes(x=po2, logtrue1),linetype="dotted",size=1.5)+
  geom_line(preds2, mapping=aes(x=po2, logtrue2),linetype="dashed",size=1.5)+
  theme(legend.position=c(0.8, 0.2))

ggsave(
  "figure_2_pars_effect_obs.png",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 24,
  height = 6,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE
)

#Calculate percent differences
preds2$diff <- case_when(preds2$name=="logmu1_low"|preds2$name=="logmu2_low"~preds2$value-preds2$logtrue1,
                         preds2$name=="logmu1_high"|preds2$name=="logmu2_high"~preds2$value-preds2$logtrue2)



p2 <- ggplot(preds2, aes(x=po2, y=diff))+geom_line(aes(colour=name), size=1.5)+
  xlab("Dissolved Oxygen (kPa)")+
  ylab(expression("Difference in"~f(phi[eco])))+
  scale_colour_manual(values=c("blue2", "skyblue2", "orange2", "gold"), labels=c("High Eo&s50, High Temp", "High Eo&s50, Low Temp", "Low Eo&s50, High Temp", "Low Eo&s50, Low Temp"))+
  #geom_line(preds2, mapping=aes(x=po2, logtrue),linetype="dashed",size=1.5)+
  theme(legend.position=c(0.8, 0.5))
