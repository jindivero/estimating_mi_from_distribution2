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

#Plot accuracy and precision
ggplot(data=Eo_performance, aes(x=Precision, y=Bias))+geom_point(aes(group=Model,color=Model), size=5)

#Density plot just Eo #
ggplot(subset(pars, pars$term=="mi-Eo"), aes(x=estimate)) +
  geom_density(fill="lightblue", adjust = 1.5) +
  geom_vline(data = Eo_values, aes(xintercept = MLE_avg),linetype="dashed", size=1.2, color="darkorange", show.legend=T)+
  geom_vline(data = Eo_values, aes(xintercept = true),linetype="dashed", size=1.2)+
  facet_grid(analysis~data)+
  xlab("Eo estimate") + 
  theme(strip.text = element_text(size = 14))




#### Comparing po2' effect from parameter estimates ####
#All combined into one figure
q1 <- ggplot(simdats1, aes(mi_usual, logmu, colour=Eo)) +
  geom_point(size=0.1)+
  scale_colour_viridis()+
  xlab("True pO2'")+
  ylab("Estimated pO2' effect")+
  theme(legend.position=c(0.9,0.2))+
  geom_line(data=true_effect, aes(x=mi, y=mi_effect), color="black", linetype="dashed", size=2)+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank())+
  facet_wrap("title")
q2 <- ggplot(simdats2, aes(mi_usual, logmu, colour=Eo)) +
  geom_point(size=0.1)+
  scale_colour_viridis()+
  xlab("True pO2'")+
  ylab("Estimated pO2' effect")+
  theme(legend.position=c(0.9,0.2))+
  geom_line(data=true_effect, aes(x=mi, y=mi_effect), color="black", linetype="dashed", size=2)+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())
q3 <- ggplot(simdats3, aes(mi_weird, logmu, colour=Eo)) +
  geom_point(size=0.1)+
  scale_colour_viridis()+
  xlab("True pO2'")+
  ylab("Estimated pO2' effect")+
  theme(legend.position=c(0.9,0.2))+
  geom_line(data=true_effect2, aes(x=mi, y=mi_effect), color="black", linetype="dashed", size=2)+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text=element_blank())+
  facet_grid(side~title)
q4 <- ggplot(simdats4, aes(mi_weird, logmu, colour=Eo)) +
  geom_point(size=0.1)+
  scale_colour_viridis()+
  xlab("True pO2'")+
  ylab("Estimated pO2' effect")+
  theme(legend.position=c(0.9,0.2))+
  geom_line(data=true_effect2, aes(x=mi, y=mi_effect), color="black", linetype="dashed", size=2)+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank())+
  facet_wrap("side", strip.position="right")

figure2 <- ggarrange(q1, q3,q2,q4,
                     common.legend = TRUE, legend = "right")

annotate_figure(figure2, left=text_grob(expression(paste("pO"[2], "'", " Effect")), size=30, rot=90),bottom=text_grob(expression(paste("pO"[2], "'", " True Value")), size=30))

#Prior constrained
q2 <- ggplot(simdats2, aes(mi_usual, logmu, colour=Eo)) +
  geom_point(size=0.1)+
  scale_colour_viridis()+
  xlab("True pO2'")+
  ylab("Estimated pO2' effect")+
  theme(legend.position=c(0.9,0.2))+
  geom_line(data=true_effect, aes(x=mi, y=mi_effect), color="black", linetype="dashed", size=2)+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())
q4 <- ggplot(simdats4, aes(mi_weird, logmu, colour=Eo)) +
  geom_point(size=0.1)+
  scale_colour_viridis()+
  xlab("True pO2'")+
  ylab("Estimated pO2' effect")+
  theme(legend.position=c(0.9,0.2))+
  geom_line(data=true_effect2, aes(x=mi, y=mi_effect), color="black", linetype="dashed", size=2)+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank())+
  facet_wrap("side", strip.position="right")

figure2 <- ggarrange(q2, q4,
                     common.legend = TRUE, legend = "right")

#unconstrained and typical species only
ggplot(simdats1, aes(mi_usual, logmu, colour=Eo)) +
  geom_point(size=0.1)+
  scale_colour_viridis()+
  xlab("True g(pO2,T)")+
  ylab("Estimated f(pO2, T)")+
  theme(legend.position=c(0.8,0.3))+
  geom_line(data=true_effect, aes(x=mi, y=mi_effect), color="black", linetype="dashed", size=2)

#unconstrained
q1 <- ggplot(simdats1, aes(mi_usual, logmu, colour=Eo)) +
  geom_point(size=0.1)+
  scale_colour_viridis()+
  xlab("True pO2'")+
  ylab("Estimated pO2' effect")+
  theme(legend.position=c(0.9,0.2))+
  geom_line(data=true_effect, aes(x=mi, y=mi_effect), color="black", linetype="dashed", size=2)+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())+
  facet_wrap("title")
q3 <- ggplot(simdats3, aes(mi_weird, logmu, colour=Eo)) +
  geom_point(size=0.1)+
  scale_colour_viridis()+
  xlab("True pO2'")+
  ylab("Estimated pO2' effect")+
  theme(legend.position=c(0.9,0.2))+
  geom_line(data=true_effect2, aes(x=mi, y=mi_effect), color="black", linetype="dashed", size=2)+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())+
  facet_grid(side~title)
figure2 <- ggarrange(q1, q3,
                     common.legend = TRUE, legend = "right")

annotate_figure(figure2, left=text_grob(expression(paste("pO"[2], "'", " Effect")), size=30, rot=90),bottom=text_grob(expression(paste("pO"[2], "'", " True Value")), size=30))

annotate_figure(figure2, left=text_grob(expression(paste("pO"[2], "'", " Effect")), size=30, rot=90),bottom=text_grob(expression(paste("pO"[2], "'", " True Value")), size=30))


##pO2 and Temp data and po2/temp curve


ggplot(dat, aes(x=log(invtemp), y=log(po2)))+geom_point()+geom_line(dat, mapping=aes(x=log(invtemp), y=log(pO2_s50)), linetype="dashed", size=2, colour="purple")
ggplot(dat, aes(x=invtemp, y=log(po2)))+geom_point()+geom_line(dat, mapping=aes(x=invtemp, y=log(pO2_s50)), linetype="dashed", size=2, colour="purple")

ggplot(dat, aes(x=temp, y=po2))+geom_point()+geom_line(dat, mapping=aes(x=temp, y=pO2_s50), linetype="dashed", size=2, colour="purple")

# Plot Eo vs s50 ##
#All combined
ggplot(pars_wide, aes(x=pars_wide$"mi-Eo", y=pars_wide$"mi-s50"))+ geom_point(aes(group=model, color=model), size=5)+xlab("Eo estimate")+ylab("s50 estimate")+
  theme(legend.position=c(0.3,0.8))

ggplot(pars_wide2, aes(x=pars_wide2$"mi-Eo", y=pars_wide2$"mi-s50"))+ 
  geom_point(size=5)+
  theme(legend.position=c(0.3,0.8))+
  xlab("Eo Maximum Likelihood Estimate")+
  ylab("s50 Threshold Maximum Likelihood Estimate")+
  stat_density2d(aes(fill = stat(level)), geom = "polygon")
#With density
ggplot(pars_wide2, aes(x=pars_wide2$"mi-Eo", y=pars_wide2$"mi-s50"))+
  theme(legend.position=c(0.2,0.7))+scale_shape_manual(labels = c("Typical Species", "Unusual Species"), values = c(16,17), name=NULL)+
  xlab("Eo Maximum Likelihood Estimate")+
  ylab("s50 Threshold Maximum Likelihood Estimate")+stat_density2d(aes(fill = stat(level)), geom = "polygon")+
  geom_point(aes(group=model, shape=model), size=5)