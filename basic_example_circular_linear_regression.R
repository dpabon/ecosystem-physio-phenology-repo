# basic example of circular-linear regression (circular: timing, lineal: climate conditions)
# created by: Daniel E. Pabon-Moreno
# MPI-BGC
# email: dpabon@bgc-jena.mpg.de


library(circular)
library(lubridate)

# function to scale the variables
scal  <- function(x) {
  y  <- (x - mean(x, na.rm = T))/sd(x, na.rm = T)
  return(y)
}

# DOY
DOY <- c(35, 365, 3, 364, 1, 25, 12, 18)


# from DOY to radians
y <- DOY * (360/365) * (pi/180)

# circular object
y <- circular(x, units = "radians", modulo = "2pi", rotation = "clock",  zero = pi/2)

# plotting x

plot.circular(x, type = "n", cex = 1, bin = 720, stack = T, sep = 0.02, shrink = 1.3, axes = F)
axis.circular(at = circular(seq(1,365, by = 75) * (360/365) * (pi/180), units = "radians", modulo = "2pi", rotation = "clock",  zero = pi/2), units =  "radians", labels = seq(1,365, by = 75), zero = pi/2, rotation = "clock",cex = 1, tick = T)
points.circular(x, stack = F, col = "red", pch = 1, cex = 1.2)


# temperature:
set.seed(31)
tair <- rnorm(length(DOY), mean = 27, sd = 5)

# precipitation (mm)
set.seed(31)
precip <- rnorm(length(DOY), mean = 200, sd = 3)


climate_vars <- as.matrix(data.frame(scal(tair), scal(precip)))



# circular-linear regression

regression <- lm.circular(type = "c-l", y = y, x = climate_vars, init = c(0,0), verbose = T)

regression
