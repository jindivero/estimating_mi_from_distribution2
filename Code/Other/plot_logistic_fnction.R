x <- seq(0, 5, length.out = 100)
s50 <- 2
delta <-2
s95 <- s50 + delta
smaxlist <- 2^(0:6)


plotlogistic <- function(x, smax, s50, delta) {
  a <- log(smax / (log(0.5) + smax) - 1)
  b <- log(smax / (log(0.95) + smax) - 1)
  beta0 <- -a + s50 * (b - a) / delta
  beta1 <- (a - b) / delta
  mu <- exp(smax * (1 / ( 1 + exp( - beta0 - beta1 * x)) -1))
  
}


plot.df <- tibble(x = NULL,
                  y = NULL,
                  smax = NULL)


for (i in 1:length(smaxlist)) {
  smax <- smaxlist[i]
  y <- plotlogistic(x, smax, s50, delta)
  plot.df <- rbind(plot.df, cbind(x, y, smax = rep(smax, length(x))))
}

plot.df$smax <- as.factor(plot.df$smax)
ggplot(plot.df) +
  geom_line( aes(x = x, y = y, group = smax, color = smax), size = 1.5) + 
  scale_color_viridis_d(direction = -1, begin = 0, end = 0.75, option = "plasma") +
  xlab("po2'") +
  ylab("marginal effect")
