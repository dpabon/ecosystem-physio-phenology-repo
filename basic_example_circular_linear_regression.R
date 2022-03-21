# basic example of circular-linear regression (circular: timing, lineal: climate conditions)
# created by: Daniel E. Pabon-Moreno
# MPI-BGC
# email: dpabon@bgc-jena.mpg.de


library(circular)
library(lubridate)

# based on the code example presented on the circular package

set.seed(1234)

# tair and precip
x <- cbind(scale(rnorm(100, mean = 25, sd = 2)), scale(rnorm(100, mean = 200, sd = 10)))

# timing + noise
y <- circular(2*atan(c(x%*%c(0.1,0.2))))+rvonmises(100, mu=circular(0), kappa=100)

# circular-linear regression
lm.circular(y=y, x=x, init=c(0.1,0.2), type='c-l', verbose=TRUE)


# plotting y

plot.circular(y, type = "n", cex = 1, bin = 720, stack = T, sep = 0.02, shrink = 1.3, axes = F)
axis.circular(at = circular(seq(1,365, by = 75) * (360/365) * (pi/180), units = "radians", modulo = "2pi", rotation = "counter",  zero = 0), units =  "radians", labels = seq(1,365, by = 75), zero = 0, rotation = "counter",cex = 1, tick = T)
points.circular(y, stack = T, col = "red", pch = 1, cex = 1.2)


