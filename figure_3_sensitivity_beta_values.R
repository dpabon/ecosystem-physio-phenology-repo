library(circular)
library(glmnet)
library(ggplot2)
library(plotly)
library(ggpubr)
library(scales)
library(utils)
library(R.utils)
library(tictoc)
library(tidyverse)
library(furrr)

abso <- function(x){
  if(x < 0){
    return(-1)
  } 
  if(x > 0){
    return(1)
  }
}


data <- 100
n_simulations <- 100

start_value <- 0.01
end_value <- 3





beta2 <- seq(from =start_value, to = end_value, by = 0.1)

beta1 <- beta2

output <- data.frame(beta1 = rep(beta1, each = length(beta2)), 
                     beta2 = rep(beta2, times = length(beta2)), 
                     linear.beta1 = NA,
                     linear.beta2 = NA,
                     circular.beta1 = NA,
                     circular.beta2 = NA,
                     linear.beta1.se = NA,
                     linear.beta2.se = NA,
                     circular.beta1.se = NA,
                     circular.beta2.se = NA)




mu <- circular(0)

for (i in 1:n_simulations){
  
}

x <- cbind(rnorm(data, mean = 10), rnorm(data, mean = 15, sd = 2))


sd.model <- rep(NA, length(beta2))

for (i in 1:nrow(output)) {
    y <- mu + circular(2*atan(c(x%*%c(output$beta1[i],output$beta2[i])))) + rvonmises(data, mu = mu, kappa = 10)
    sd.model[i] <- sd.circular(y)
    mean.circular(y)
}

hist(sd.model)
min(sd.model)
max(sd.model)

library(foreach)
library(doParallel)

#setup parallel backend to use many processors
cores=30
cl <- makeCluster(cores) #not to overload your computer
registerDoParallel(cl)
tic()
finalMatrix <- foreach(i=1:nrow(output), .combine=rbind, .packages= c("circular", "R.utils")) %dopar% {
  
  y <- mu + circular(2*atan(c(x%*%c(output$beta1[i],output$beta2[i])))) + rvonmises(data, mu = mu, kappa = 40)
  coef <- withTimeout(try(lm.circular(y = y, x = x, init = c(0,0), type = 'c-l', verbose = F, tol = 1e-10)), timeout = 10, onTimeout = "warning") 
  if ((class(coef) == "try-error") || is.null(coef)) {
    tempMatrix = c(rep(NA, 8))
  }else{
    coef.lm <- lm(as.numeric(y)~x[,1] + x[,2])
    tempMatrix = c(coef.lm$coefficients[2:3], 
                   coef$coefficients, 
                   summary(coef.lm)$coefficients[2:3,2],
                   coef$se.coef)
  }
  
  tempMatrix
}
toc()
#stop cluster
stopCluster(cl)

output[,c(3:10)] <- finalMatrix

## differences between the distances

output$linear.beta1.dist <- abs((output$linear.beta1 - output$beta1) / output$beta1) 
output$linear.beta2.dist <- abs((output$linear.beta2 - output$beta2) / output$beta2) 
output$circular.beta1.dist <- abs((output$circular.beta1 - output$beta1) / output$beta1) 
output$circular.beta2.dist <- abs((output$circular.beta2 - output$beta2) / output$beta2) 

## comparison circular vs linear

output$comp.beta1 <-  output$linear.beta1.dist - output$circular.beta1.dist
output$comp.beta2 <-  output$linear.beta2.dist - output$circular.beta2.dist

## summary both coefficients

output$comp.summary <- (output$comp.beta1 + output$comp.beta2) / 2


ggplot(data = output, aes(factor(beta1), factor(beta2), fill = linear.beta1.dist))+
  geom_tile() + 
  scale_fill_gradient2() 




ggplot(data = output, aes(factor(beta1), factor(beta2), fill = comp.beta1))+
  geom_tile() + 
  scale_fill_gradient2() 


ggplot(data = output, aes(factor(beta1), factor(beta2), fill = comp.beta2))+
  geom_tile() + 
  scale_fill_gradient2()

ggplot(data = output, aes(factor(beta1), factor(beta2), fill = comp.summary))+
  geom_tile() +
  scale_fill_gradient2() +
  scale_x_discrete(breaks = 50)

# scaling the compendium

output$comp.summary.scale <- NA
output$comp.summary.class <- NA

for (i in 1:nrow(output)){
  if(is.na(output$comp.summary[i]) == F){
    if(output$comp.summary[i] > 0){
      output$comp.summary.scale[i] <- 1
      output$comp.summary.class[i] <- "circular"
      
    }
    if(output$comp.summary[i] < 0){
      output$comp.summary.scale[i] <- -1
      output$comp.summary.class[i] <- "linear"
    }
    if(output$comp.summary[i] == 0){
      output$comp.summary.scale[i] <- 0
      output$comp.summary.class[i] <- "both"
    }
  }
}

ggplot(data = output, aes(factor(beta1), factor(beta2), fill = comp.summary.scale))+
  geom_tile() +
  scale_fill_viridis_c()

par(mfrow = c(2,1))

hist(abs(output$comp.summary[which(output$comp.summary.class == "linear")]), main = "Distribution of the differences when linear is better")

hist(abs(output$comp.summary[which(output$comp.summary.class == "circular")]), main = "Distribution of the differences when circular is better")



min(abs(output$comp.summary[which(output$comp.summary.class == "linear")]))
max(abs(output$comp.summary[which(output$comp.summary.class == "linear")]))
l1 <- length(abs(output$comp.summary[which(output$comp.summary.class == "linear")]))



min(abs(output$comp.summary[which(output$comp.summary.class == "circular")]))
max(abs(output$comp.summary[which(output$comp.summary.class == "circular")]))
l2 <- length(abs(output$comp.summary[which(output$comp.summary.class == "circular")]))
barplot(c(l1,l2),names.arg = c("linear", "circular"), main = "number of times when the model perfom better")
## comparison by histograms

ggplot(output, aes(x = abs(comp.summary), color = comp.summary.class))+
  geom_histogram(fill = "white", alpha = 0.5) + xlab("Performance distance")

ggplot(output, aes(y = abs(comp.summary), x = comp.summary.class))+
  geom_boxplot()


# plots SE

min.plot <- min(apply(cbind(output$linear.beta1.se, output$linear.beta2, output$circular.beta1.se, output$circular.beta2.se), MARGIN = 2, FUN = min, na.rm = T))

max.plot <- min(apply(cbind(output$linear.beta1.se, output$linear.beta2, output$circular.beta1.se, output$circular.beta2.se), MARGIN = 2, FUN = max, na.rm = T))


ggplot(data = output, aes(factor(beta1), factor(beta2), fill = linear.beta1.se))+
  geom_tile() +
  scale_fill_viridis_c() + ggtitle("Linear beta1")+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())




linear_beta1 <- ggplot(data = output, aes(factor(beta1), factor(beta2), fill = linear.beta1.se))+
  geom_tile() +
  scale_fill_viridis_c(limits = c(min.plot, max.plot)) + ggtitle("Linear beta1")+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())


linear_beta2 <- ggplot(data = output, aes(factor(beta1), factor(beta2), fill = linear.beta2.se))+
  geom_tile() +
  scale_fill_viridis_c()+ ggtitle("Linear beta2")+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())




circular_beta1 <- ggplot(data = output, aes(factor(beta1), factor(beta2), fill = circular.beta1.se))+
  geom_tile() +
  scale_fill_viridis_c()+ ggtitle("circular beta1")+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())


circular_beta2 <- ggplot(data = output, aes(factor(beta1), factor(beta2), fill = circular.beta2.se))+
  geom_tile() +
  scale_fill_viridis_c()+ ggtitle("circular beta2")+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())


ggarrange(linear_beta1, linear_beta2, circular_beta1, circular_beta2)





########################


#  (furrr)
circular.fun <- function(beta1, beta2, x, mu, vonmisesNoise){
  y <- as.numeric(mu) + circular(2*atan(c(x%*%c(beta1,beta2)))) + vonmisesNoise
  coef <- withTimeout(try(lm.circular(y = y, x = x, init = c(0,0), type = 'c-l', verbose = F, tol = 1e-6)), timeout = 50, onTimeout = "warning")
  coef.lm <- lm(as.numeric(y)~x[,1] + x[,2])
  if ((class(coef) == "try-error") || is.null(coef)) {
    return(data.frame(l.beta1 = NA, 
                      l.beta2 = NA,
                      c.beta1 = NA,
                      c.beta2 = NA,
                      sd.circular  = NA
                      ))
  }else{
    return(data.frame(l.beta1 = coef.lm$coefficients[2], 
                      l.beta2 = coef.lm$coefficients[3],
                      c.beta1 = coef$coefficients[1],
                      c.beta2 = coef$coefficients[2],
                      sd.circular  = sd.circular(y)
                      ))
  }
}


# beta1.in <- output$beta1
# 
# names(beta1.in) <- seq(1, length(beta1.in))

plan(multisession(workers = 30))
tic()
test <- future_map2_dfr(.x = output$beta1, .y = output$beta2, .f =  circular.fun, x = x, mu = mu, .progress = F)
toc()
plan(sequential)

circular.fun(beta1 = 0.1, beta2 = 0.1, x = x, mu = mu)

mapply(circular.fun, mu = mu, x = list(x), beta1 = output$beta1, beta2 = output$beta2)

# running bootstrapping -----

output.coef <- data.frame(l.beta1 = NA, 
                          l.beta2 = NA,
                          c.beta1 = NA,
                          c.beta2 = NA,
                          sd.circular  = NA
)




mu <- circular(0)
vonmises <- rvonmises(data, mu = mu, kappa = 30)
plan(multisession(workers = 30))
tic()
for (i in 1:n_simulations) {
  x <- cbind(rnorm(data, mean = 10), rnorm(data, mean = 15, sd = 2))
  output.temp <- future_map2_dfr(.x = output$beta1, .y = output$beta2, .f =  circular.fun, x = x, mu = mu, vonmisesNoise = vonmises,  .progress = F)
  output.coef <- rbind(output.coef, output.temp)
  print(i)
}
toc()
plan(sequential)

saveRDS(output.coef, file = "output_coef_100sim_mu0.RDS")

## mu = pi

mu <- circular(pi)

vonmises <- rvonmises(data, mu = mu, kappa = 30)

plan(multisession(workers = 30))

output.coef <- data.frame(l.beta1 = NA,
                          l.beta2 = NA,
                          c.beta1 = NA,
                          c.beta2 = NA,
                          sd.circular  = NA
)
tic()
for (i in 1:n_simulations) {
  x <- cbind(rnorm(data, mean = 10), rnorm(data, mean = 15, sd = 2))
  output.temp <- future_map2_dfr(.x = output$beta1, .y = output$beta2, .f =  circular.fun, x = x, mu = mu, vonmisesNoise = vonmises,  .progress = F)
  output.coef <- rbind(output.coef, output.temp)
  print(i)
}
toc()
plan(sequential)

saveRDS(output.coef, file = "output_coef_100sim_mupi.RDS")


output.coef.0 <- readRDS("output_coef_100sim_mu0.RDS")
output.coef.pi <- readRDS("output_coef_100sim_mupi.RDS")

output.coef.0 <- output.coef.0[-1,]
output.coef.pi <- output.coef.pi[-1,]

hist(output.coef.0$sd.circular)
hist(output.coef.pi$sd.circular)

# simulation index
output.coef.0$simulation <- as.factor(rep(1, n_simulations, each = length(beta1)))
output.coef.pi$simulation <- as.factor(rep(1, n_simulations, each = length(beta1)))


# mu = 0 ----

# original coefficients
output.coef.0$beta1 <- output$beta1
output.coef.0$beta2 <- output$beta2
output.coef.pi$beta1 <- output$beta1
output.coef.pi$beta2 <- output$beta2
View(output.coef.0)

## differences between the distances

output.coef.0$linear.beta1.dist <- abs((output.coef.0$l.beta1 - output.coef.0$beta1)) # / output.coef.0$beta1)
output.coef.0$linear.beta2.dist <- abs((output.coef.0$l.beta2 - output.coef.0$beta2)) # / output.coef.0$beta2) 
output.coef.0$circular.beta1.dist <- abs((output.coef.0$c.beta1 - output.coef.0$beta1)) # / output.coef.0$beta1) 
output.coef.0$circular.beta2.dist <- abs((output.coef.0$c.beta2 - output.coef.0$beta2)) # / output.coef.0$beta2) 

## comparison circular vs linear

output.coef.0$comp.beta1 <-  output.coef.0$linear.beta1.dist - output.coef.0$circular.beta1.dist
output.coef.0$comp.beta2 <-  output.coef.0$linear.beta2.dist - output.coef.0$circular.beta2.dist

## summary both coefficients

output.coef.0$comp.summary <- (output.coef.0$comp.beta1 + output.coef.0$comp.beta2) / 2


output.coef.0$index <- paste(output.coef.0$beta1, output.coef.0$beta2, sep = "_") 

test <- tapply(output.coef.0$comp.summary, output.coef.0$index, sd)

equals(as.numeric(mean.com$comp.summary), as.numeric(test))
head(test)
head(mean.com$comp.summary)

mean.com.0 <- output.coef.0 %>%
  mutate(index = paste(beta1, beta2, sep = "_")) %>% 
  group_by(simulation, index) %>% 
  summarise(comp.summary.com = mean(comp.summary),
            comp.summary.sd = stats::sd(comp.summary, na.rm = T),
            n = n(),
            mean.linear.beta1.dist = mean(linear.beta1.dist),
            mean.linear.beta2.dist = mean(linear.beta2.dist),
            mean.circular.beta1.dist = mean(circular.beta1.dist),
            mean.circular.beta2.dist = mean(circular.beta2.dist)
            )

mean.com.0 <- mean.com.0 %>% 
  mutate(beta1 = output$beta1, beta2 = output$beta2)


accuracy.0 <- ggplot(data = mean.com.0, aes(factor(beta1), factor(beta2), fill = comp.summary.com))+
  geom_tile() +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  scale_fill_gradient2(low = muted("red"), mid = "gray",
                       high = muted("blue"), midpoint = 0, space = "Lab",
                       na.value = "grey50", guide = "colourbar", aesthetics = "fill") +
  xlab("beta 1") +
  ylab("beta 2") +
  labs(fill = "Differences in \n Accuracy")

accuracy.0

mean.com.0 <- mean.com.0 %>% 
  mutate(abs.com.summary = modify(mean.com.0$comp.summary.com, abso))

ggplot(data = mean.com.0, aes(factor(beta1), factor(beta2), fill = abs.com.summary.com))+
  geom_tile() +
  scale_fill_gradient2(low = muted("red"), mid = "gray",
                       high = muted("blue"), midpoint = 0, space = "Lab",
                       na.value = "grey50", guide = "colourbar", aesthetics = "fill")


# precision metric


output.coef.0$comp.beta1 <-  output.coef.0$linear.beta1.dist - output.coef.0$circular.beta1.dist
output.coef.0$comp.beta2 <-  output.coef.0$linear.beta2.dist - output.coef.0$circular.beta2.dist

precision.metrics.0 <- output.coef.0 %>% 
  mutate(index = paste(beta1, beta2, sep = "_")) %>% 
  group_by(simulation, index) %>% 
  summarise(sd.beta1.com = sd(linear.beta1.dist) - sd(circular.beta1.dist),
            sd.beta2.com = sd(linear.beta2.dist) - sd(circular.beta2.dist)) %>% 
  mutate(sd.com = (sd.beta1.com / sd.beta2.com) / 2)

precision.metrics.0 <- precision.metrics.0 %>% 
  mutate(beta1 = output$beta1, beta2 = output$beta2)


precision.0 <- ggplot(data = precision.metrics.0, aes(factor(beta1), factor(beta2), fill = sd.com))+
  geom_tile() +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  scale_fill_gradient2(low = muted("red"), mid = "gray",
                       high = muted("blue"), midpoint = 0, space = "Lab",
                       na.value = "grey50", guide = "colourbar", aesthetics = "fill") +
  xlab("beta 1") +
  ylab("beta 2") +
  labs(fill = "Differences in \n Precision")

precision.0

# mu = pi ----

# original coefficients
output.coef.pi$beta1 <- output$beta1
output.coef.pi$beta2 <- output$beta2
output.coef.pi$beta1 <- output$beta1
output.coef.pi$beta2 <- output$beta2
View(output.coef.pi)

## differences between the distances

output.coef.pi$linear.beta1.dist <- abs((output.coef.pi$l.beta1 - output.coef.pi$beta1)) # / output.coef.pi$beta1)
output.coef.pi$linear.beta2.dist <- abs((output.coef.pi$l.beta2 - output.coef.pi$beta2)) # / output.coef.pi$beta2) 
output.coef.pi$circular.beta1.dist <- abs((output.coef.pi$c.beta1 - output.coef.pi$beta1)) # / output.coef.pi$beta1) 
output.coef.pi$circular.beta2.dist <- abs((output.coef.pi$c.beta2 - output.coef.pi$beta2)) # / output.coef.pi$beta2) 

## comparison circular vs linear

output.coef.pi$comp.beta1 <-  output.coef.pi$linear.beta1.dist - output.coef.pi$circular.beta1.dist
output.coef.pi$comp.beta2 <-  output.coef.pi$linear.beta2.dist - output.coef.pi$circular.beta2.dist

## summary both coefficients

output.coef.pi$comp.summary <- (output.coef.pi$comp.beta1 + output.coef.pi$comp.beta2) / 2


output.coef.pi



mean.com.pi <- output.coef.pi %>%
  mutate(index = paste(beta1, beta2, sep = "_")) %>% 
  group_by(simulation, index) %>% 
  summarise(comp.summary.mean = mean(comp.summary),
            comp.summary.sd = sd(comp.summary, na.rm = T),
            n = n())

mean.com.pi <- mean.com.pi %>% 
  mutate(beta1 = output$beta1, beta2 = output$beta2)


accuracy.pi <- ggplot(data = mean.com.pi, aes(factor(beta1), factor(beta2), fill = comp.summary.mean))+
  geom_tile() +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  scale_fill_gradient2(low = muted("red"), mid = "gray",
                       high = muted("blue"), midpoint = 0, space = "Lab",
                       na.value = "grey50", guide = "colourbar", aesthetics = "fill") +
  xlab("beta 1") +
  ylab("beta 2") +
  labs(fill = "Differences in \n Accuracy")



mean.com <- mean.com %>% 
  mutate(abs.com.summary = modify(mean.com$comp.summary.mean,abso))

ggplot(data = mean.com, aes(factor(beta1), factor(beta2), fill = abs.com.summary))+
  geom_tile() +
  scale_fill_gradient2(low = muted("red"), mid = "gray",
                       high = muted("blue"), midpoint = 0, space = "Lab",
                       na.value = "grey50", guide = "colourbar", aesthetics = "fill")

# precision metric


output.coef.0$comp.beta1 <-  output.coef.0$linear.beta1.dist - output.coef.0$circular.beta1.dist
output.coef.0$comp.beta2 <-  output.coef.0$linear.beta2.dist - output.coef.0$circular.beta2.dist

precision.metrics.pi <- output.coef.pi %>% 
  mutate(index = paste(beta1, beta2, sep = "_")) %>% 
  group_by(simulation, index) %>% 
  summarise(sd.beta1.com = sd(linear.beta1.dist) - sd(circular.beta1.dist),
            sd.beta2.com = sd(linear.beta2.dist) - sd(circular.beta2.dist)) %>% 
  mutate(sd.com = (sd.beta1.com / sd.beta2.com) / 2)

precision.metrics.pi <- precision.metrics.pi %>% 
  mutate(beta1 = output$beta1, beta2 = output$beta2)


precision.pi <- ggplot(data = precision.metrics.pi, aes(factor(beta1), factor(beta2), fill = sd.com))+
  geom_tile() +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  scale_fill_gradient2(low = muted("red"), mid = "gray",
                       high = muted("blue"), midpoint = 0, space = "Lab",
                       na.value = "grey50", guide = "colourbar", aesthetics = "fill") +
  xlab("beta 1") +
  ylab("beta 2") +
  labs(fill = "Differences in \n Precision")

hist(precision.metrics.0$sd.com)
max(mean.com.pi$comp.summary.mean)
max(accura)
## all plots together

accuracy.0 <- accuracy.0 + ggtitle(expression(paste("a)  ", mu, "=", 0, "(Beggining of the year)", sep = " ")))
accuracy.pi <- accuracy.pi + ggtitle(expression(paste("b)  ", mu, "=", pi, "(Mid-year)", sep = " ")))

precision.0 <- precision.0 + ggtitle(expression(paste("c)  ", mu, "=", 0, "(Beggining of the year)", sep = " ")))
precision.pi <- precision.pi + ggtitle(expression(paste("d)  ", mu, "=", pi, "(Mid-year)", sep = " ")))

all.plots <- ggarrange(accuracy.0, accuracy.pi, precision.0, precision.pi)

ggsave("figure3.pdf" ,all.plots, device = "pdf",  width=12, height=8)
ggsave("figure3.png" ,all.plots, device = "png",  width=12, height=8)
