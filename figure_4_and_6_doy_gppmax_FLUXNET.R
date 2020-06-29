# libraries ----
library(ncdf4)
library(abind)
library(lubridate)
library(ggplot2)
library(viridis)
library(circular)
library(raster)
library(snow)
library(zoo)
library(rgdal)
library(graticule)
library(MASS)
library(lubridate)
library(maps)
library(proj4)
library(cowplot)
library(gridExtra)
library(grid)
library(plotly)
library(maps)
library(scales)
library(reshape2)
library(rworldmap)
library(mice)
library(Kendall)
library(trend)
library(MASS)
library(glmnet)
library(dimRed)
library(animation)
library(R.utils)
library(gridBase)
library(brew)
library(ggpubr)
library(telegram.bot)

# functions ----
## ggplot2 colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# function to scale the variables
scal  <- function(x) {
  y  <- (x - mean(x, na.rm = T))/sd(x, na.rm = T)
  return(y)
}


# half life function
half.time <- function(init = 1, t, t_half = 30) {
  lambda <- log(2)/t_half
  t <- 1:t
  y <- vector()
  y[1] <- init
  for (i in 2:length(t)) {
    y[i] <- init * exp(1)^(-lambda*t[i])
    #y[i] <- init * ((1/2)^((t[i] - 1)/(1/2)))
  }
  y
  return(y)
}

# function to estimate the mean per day
roll.d <- function(x, by, FUN, na.rm = T) {
  y <- matrix(x, ncol = (length(x)/by), byrow = F)
  y <- apply(y, MARGIN = 2, FUN = FUN, na.rm = na.rm)
  return(y)
}


# circular prediction
circ.recons <- function(x, mu, coef) {
  out <- mu + (2 * atan(sum(coef*x)))
  return(out)
}

# linear prediction
magnitude.recons <- function(x, intersect, coef) {
  y <- intersect + sum(x*coef)
  return(y)
}

## tools
source("bin/gpp/tools.R")


#### Selecting only sites with at least 7 years of duration ----
m5 <- read.csv("data/fluxnet/m5.csv")

View(m5)


## All the data indexed in a single array ----
time.serie <- array(list(list()))

## quality filter value
qc.level <- 2

for (st in 1:nrow(m5)) {
  if (m5$FILTER2[st] == 1) {
    folder <- "data/fluxnet/FLUXNET_sites/F15/"
    nc.file  <- nc_open(paste(folder, m5$SITE_ID[st], ".HH.", m5$DATA_START[st],".",
                              m5$DATA_END[st],".nc", sep = ""))
    # time original for index
    hours.re <- ncvar_get(nc = nc.file, varid = "time", start = 1, count = -1)
    index.days <- as.character(as.Date(hours.re, origin = "1582-10-15 00:00", format = "%Y-%m-%d  %R"))
    # GPP
    # GPP without filter
    gpp.o <- ncvar_get(nc = nc.file, varid = "GPP_NT_VUT_USTAR50", start = 1, count = -1)
    # GPP with filter
    quality.filter.gpp <- ncvar_get(nc = nc.file, varid = "NEE_VUT_USTAR50_QC", start = 1, count = -1)
    gpp.o[which(quality.filter.gpp >= qc.level)] <- NA
    time.serie[st][[1]][[1]] <-  tapply(gpp.o, INDEX = index.days, FUN = quantile, probs = 0.9, na.rm = T)
    
    # time
    hours.re <- ncvar_get(nc = nc.file, varid = "time", start = 1, count = -1)
    time.serie[st][[1]][[2]] <- as.Date(hours.re, origin = "1582-10-15 00:00", format = "%Y-%m-%d  %R")
    unique.days <- unique(as.character(time.serie[st][[1]][[2]]))
    time.serie[st][[1]][[2]] <- as.Date(unique.days)
    length(time.serie[st][[1]][[1]]) == length(time.serie[st][[1]][[2]])
    
    # Tair
    time.serie[st][[1]][[3]] <- tapply(X = ncvar_get(nc = nc.file, varid = "TA_F", start = 1, count = -1), INDEX = index.days, FUN = mean, na.rm = T)
    
    # Precipitation
    time.serie[st][[1]][[4]] <- tapply(X = ncvar_get(nc = nc.file, varid = "P_F", start = 1, count = -1), INDEX = index.days, FUN = sum, na.rm = T)
    
    # Swin
    time.serie[st][[1]][[5]] <- tapply(X = ncvar_get(nc = nc.file, varid = "SW_IN_F", start = 1, count = -1), INDEX = index.days, FUN = mean, na.rm = T)
    
    # VPD
    time.serie[st][[1]][[7]] <- tapply(X = ncvar_get(nc = nc.file, varid = "VPD_F", start = 1, count = -1), INDEX = index.days, FUN = mean, na.rm = T)
    print(paste(st, m5$SITE_ID[st], sep = "-"))
  }
}
m5$mean.temp <- NA
for (i in 1:nrow(m5)) {
  if (m5$FILTER2[i] == 1) {
    y <- year(time.serie[i][[1]][[2]])
    m5$mean.temp[i] <- mean(tapply(X = time.serie[i][[1]][[3]], INDEX = y, FUN = mean, na.rm = T))
  }
}

## detecting incogruencies -----


# depuring site by site

rm.year <- function(x, y.){
  # x Vector with original days.
  # y. years to remove
  # return a vector with the positions correspondig to the year y.
  mid <- lubridate::year(x)
  if (length(y.) == 1) {
    ret <- which(mid == y.)
  }else{
    ret <- match(mid, y.)
    ret <- which(!is.na(ret))
  }
  return(ret)
}

# generating an index list with the necessary modifications
inde <- list()
inde[[1]] <- c(1991)
inde[[2]] <- c(1995)
inde[[3]] <- c(2003, 1996:1998)
inde[[4]] <- c(1996)
inde[[12]] <- c(1998)
inde[[13]] <- c(1998,2000,2001,2004)
inde[[15]] <- c(2007:2010)
inde[[16]] <- c(1997, 2004:2005)
inde[[18]] <- c(2011)
inde[[20]] <- c(2000, 2009,2011:2013)
inde[[22]] <- c(2000)
inde[[23]] <- c(2001)
inde[[25]] <- c(1997:1998)
inde[[29]] <- c(2008:2011)
inde[[31]] <- c(2006)
inde[[32]] <- c(2002,2014)
inde[[33]] <- c(2002)
inde[[37]] <- c(1997:1999, 2009)
inde[[40]] <- c(2013)
inde[[41]] <- c(2013)
inde[[42]] <- c(2013)
inde[[43]] <- c(2014)
inde[[45]] <- c(2014)
inde[[46]] <- c(1995, 1996,2007)
inde[[54]] <- c(2007:2009)
inde[[60]] <- c(2009)
inde[[62]] <- c(2004)
inde[[65]] <- c(2006:2008)
inde[[68]] <- c(2004)
inde[[73]] <- c(1996,1999,2000,2006)
inde[[78]] <- c(2004,2005)
inde[[80]] <- c(2000)
inde[[81]] <- c(1996,1997,2005)
inde[[82]] <- c(1995,2004)

## tier2

inde[[9]] <- 2014 #NL-Loo
inde[[31]] <- c(inde[[31]], 2011,2012,2013) #ZA-Kru
inde[[55]] <- c(2014) #FR-Gri
sites.changed <- c(1,2,3,4,9,12,13,15,16,18,20,22,23,25,29,31,32,33,37,40,41,42,43,45,46,54,55,60,62,65,68,73,78,80,81,82, 9, 31, 55)

for (st in 1:length(inde)) {
  if (is.null(inde[[st]]) == F) {
    if (m5$FILTER2[st] == 1) {
      time.serie[st][[1]][[1]] <- time.serie[st][[1]][[1]][-rm.year(time.serie[st][[1]][[2]], y. = inde[[st]])]
      # index 6 correspond to the time vector modified from vector 2
      time.serie[st][[1]][[6]] <- time.serie[st][[1]][[2]][-rm.year(time.serie[st][[1]][[2]], y. = inde[[st]])]
    }
  }
}

## sites to remove

s.r  <- c(38, 72, 74, 75, 76, 77, 83)

for (i in 1:length(s.r)) {
  time.serie[s.r[i]][[1]][[1]] <- NA
}

# number of sites with information

nrow(m5[-s.r,])
row.names(m5) <- 1:nrow(m5)

length(which(m5$FILTER2 == 1))


## Extracting GPPmax GPPmax timing ----
options(warn = 2)
# template to extract the gppmax
gpp.day.template.f <- list()
for (st in 1:nrow(m5)) {
  if (m5$FILTER2[st] == 1) {
    if (is.element(st, sites.changed)) {
      gpp.nt.vut.ref <- time.serie[st][[1]][[1]]
      gpp.nt.vut.ref.time <- year(time.serie[st][[1]][[6]])
      hours.re <- time.serie[st][[1]][[6]]
    }else{
      gpp.nt.vut.ref <- time.serie[st][[1]][[1]]
      gpp.nt.vut.ref.time <- year(time.serie[st][[1]][[2]])
      hours.re <- time.serie[st][[1]][[2]]
    }
    l.years <- length(unique(gpp.nt.vut.ref.time))
    u.years <- unique(gpp.nt.vut.ref.time)
    all.years  <- matrix(NA, nrow = 365, ncol = l.years)
    
    for (y in 1:l.years) {
      salida <- time.serie[st][[1]][[1]][which(gpp.nt.vut.ref.time == u.years[y])]
      if (leap_year(u.years[y]) == T) {
        salida  <- salida[1:365]
      } else{
        #
      }
      all.years[,y]  <- salida
    }
    # mean seasonal cycle
    msc  <- apply(X = all.years, MARGIN = 1, FUN = mean, na.rm = T)
    if (any(is.nan(msc))) {
      msc[which(is.nan(msc))] <- 0
    }
    if (any(is.na(msc)) == FALSE) {
      # number of cycles
      n.cycles  <- as.integer(which.max(abs(fft(msc))[2:4]))
      # Variance explained of the seasonal patterns
      var_in_cycle <- (2*sum(abs(fft(msc)^2)[2:4]) / sum(abs(fft(msc)[-1]^2)))
      m5[st, "n.cycles"] <- n.cycles
      # timing template
      gpp.day.template  <- vector()
      fy <- fft(msc)
      ny <- length(fy)
      fy[c(1,5:(ny - 3))] <- 0
      yfiltered <- fft(fy,inverse = T)
      yfiltered  <- Re(yfiltered)
      for (d in 1:365) {
        if (is.na(yfiltered[d])) {
          #pass
        }else{
          if (d == 1) {
            if (yfiltered[d] > yfiltered[365] & yfiltered[d] > yfiltered[(d + 1)]) {
              gpp.day.template  <- c(gpp.day.template, d)
            }
          } else if (d == 365) {
            if (yfiltered[d] > yfiltered[(d - 1)] & yfiltered[d] > yfiltered[1]) {
              gpp.day.template  <- c(gpp.day.template, d)
            }
          } else {
            if (yfiltered[d] > yfiltered[(d - 1)] & yfiltered[d] > yfiltered[(d + 1)]) {
              gpp.day.template  <- c(gpp.day.template, d)
            }
          }
        }
      }
      # checking that the number of peaks is equal to the number of cycles (only selecting the highest peaks)
      if (length(gpp.day.template) == n.cycles) {
        gpp.day.template.f[[st]] <- gpp.day.template
      } else {
        tempora  <- yfiltered[gpp.day.template]
        gpp.day.template  <- gpp.day.template[order(tempora, decreasing = T)[1:n.cycles]]
        gpp.day.template <- sort(gpp.day.template)
        gpp.day.template.f[[st]] <- gpp.day.template
      }
    }
  }
}
gpp.day.template.f[[37]]

# extracting gppmax and gppmax timing from original data with quality flags
# sites.info ----
sites.info  <- array(list(list()))
#st <- 23
for (st in 1:nrow(m5)) {
  if (m5$FILTER2[st] == 1) {
    if (is.null(gpp.day.template.f[[st]]) == F) {
      if (is.element(st, sites.changed)) {
        gpp.nt.vut.ref <- time.serie[st][[1]][[1]]
        gpp.nt.vut.ref.time <- year(time.serie[st][[1]][[6]])
        hours.re <- time.serie[st][[1]][[6]]
      }else{
        gpp.nt.vut.ref <- time.serie[st][[1]][[1]]
        gpp.nt.vut.ref.time <- year(time.serie[st][[1]][[2]])
        hours.re <- time.serie[st][[1]][[2]]
      }
      
      l.years <- length(unique(gpp.nt.vut.ref.time))
      u.years <- unique(gpp.nt.vut.ref.time)
      all.years  <- matrix(NA, nrow = 365, ncol = l.years)
      
      for (y in 1:l.years) {
        salida <- gpp.nt.vut.ref[which(gpp.nt.vut.ref.time == u.years[y])]
        if (leap_year(u.years[y]) == T) {
          salida  <- salida[1:365]
        }
        all.years[,y]  <- salida
      }
      gpp.timing  <- vector()
      gpp.max  <- vector()
      gpp.year  <- vector()
      n.cycles <- length(gpp.day.template.f[[st]])
      gpp.day.template <- gpp.day.template.f[[st]]
      for (y in 1:l.years) {
        windows.size  <- 90 / n.cycles
        gpp.temp <- all.years[,y]
        for (d in 1:n.cycles) {
          if ((gpp.day.template[d]  - windows.size) <= 0) {
            back  <- abs(gpp.day.template[d]  - windows.size)
            new.vector  <- c((length(all.years[,y]) - back):length(all.years[,y]),
                             1:(gpp.day.template[d] + windows.size))
            timing <- new.vector[order(gpp.temp[new.vector], decreasing = T)[1:10]]
            gpp.timing  <- c(gpp.timing, timing)
            gp  <- gpp.temp[timing]
            gpp.max  <- c(gpp.max, gp)
            year.t  <- rep(y, times = length(timing))
            gpp.year <- c(gpp.year, year.t)
          } else if ((gpp.day.template[d]  + windows.size) > length(all.years[,y])) {
            advance  <-  (gpp.day.template[d]  + windows.size) - length(all.years[,y])
            new.vector  <- c((gpp.day.template[d] - windows.size):length(all.years[,y]), 1:advance)
            timing <- new.vector[order(gpp.temp[new.vector], decreasing = T)[1:10]]
            gpp.timing  <- c(gpp.timing, timing)
            gp  <- gpp.temp[timing]
            gpp.max  <- c(gpp.max, gp)
            year.t  <- rep(y, times = length(timing))
            gpp.year <- c(gpp.year, year.t)
          } else {
            start  <- gpp.day.template[d] - windows.size
            end  <- gpp.day.template[d] + windows.size
            new.vector  <- start:end
            timing <- new.vector[order(gpp.temp[new.vector], decreasing = T)[1:10]]
            gpp.timing  <- c(gpp.timing, timing)
            gp  <- gpp.temp[timing]
            gpp.max  <- c(gpp.max, gp)
            year.t  <- rep(y, times = length(timing))
            gpp.year <- c(gpp.year, year.t)
          }
        }
      }
      sites.info[st][[1]][[1]]  <- gpp.max
      sites.info[st][[1]][[2]]  <- gpp.timing
      sites.info[st][[1]][[3]]  <- gpp.year
      if (n.cycles == 2) {
        sites.info[st][[1]][[5]]  <- rep(c(1,2), each = 10, times = length(unique(gpp.year)) )
      }
    }
  }
}

gradient <- circular(sites.info[23][[1]][[2]] * (360/365) * (pi/180),
                     units = "radians", modulo = "2pi",
                     rotation = "clock", zero = pi/2)

color <- viridis(length(unique(sites.info[23][[1]][[3]])), direction = -1)

plot(1:130, col = color[sites.info[23][[1]][[3]]])
plot.default(as.numeric(sites.info[23][[1]][[2]]), col = color[sites.info[23][[1]][[3]]])
plot.circular(gradient, stack = T, type = "n", axes = F)
axis.circular(at = circular(seq(1,365, by = 75) * (360/365) * (pi/180), units = "radians", modulo = "2pi", rotation = "clock",  zero = pi/2),
              units =  "radians", labels = seq(1,365, by = 75), zero = pi/2, rotation = "clock",
              cex = 0.7, tick = T)
points.circular(gradient, col = color[sites.info[23][[1]][[3]]], pch = 1, stack = T, cex = 0.3 )



## Extracting index of the year ----

for (st in 1:nrow(m5)) {
  if (m5$FILTER2[st] == 1) {
    if (is.null(sites.info[[st]]) == F) {
      if (is.element(st, sites.changed)) {
        start.ind <- unique(year(time.serie[st][[1]][[6]]))
      }else{
        start.ind <- unique(year(time.serie[st][[1]][[2]]))
      }
      temp.date <- vector()
      for (i in 1:length(sites.info[st][[1]][[3]])) {
        temp.date <- c(temp.date, as.Date(sites.info[st][[1]][[2]][i],
                                          paste((start.ind[sites.info[st][[1]][[3]][i]] - 1), 12, 31, sep = "-")))
      }
      sites.info[st][[1]][[4]] <- temp.date
    }
  }
}

## Summary growing seasons | Figure S6. Supplement 1 ----
# Plot per site (Including mean and sd in a table)

lat.color <- round(sort(m5[which(m5$FILTER2 == 1),]$LAT), digits = 1)
ref <- seq(lat.color[1], lat.color[length(lat.color)], by = 0.1)
lat.color <- sample(lat.color)
lat.color.1 <-  viridis(length(ref), direction = -1, alpha = 0.4)

lat.color.2 <- viridis(length(ref), direction = -1)
mt <- month(seq(as.Date("2005/1/1"), as.Date("2005/12/31"), "days"))

mt <- match(1:12, mt)

plot(1:length(sub.gs.1$LAT), col = lat.color.1)
dev.off()
par(bg = "gray")
pdf("doc/papers/gppmax/images/summary_sites_dist.pdf",width = 5.5, height = 5.5)
gs1 <- sites.info[1][[1]][[2]]
gs1 <- circular(gs1 * (360/365) * (pi/180), units = "radians", modulo = "2pi", rotation = "clock",  zero = pi/2)
plot.circular(gs1, type = "n", cex = 1, bin = 720, stack = T, sep = 0.02, shrink = 1.3, axes = F)
axis.circular(at = circular(seq(1,365, by = 75) * (360/365) * (pi/180), units = "radians", modulo = "2pi", rotation = "clock",  zero = pi/2), units =  "radians", labels = seq(1,365, by = 75), zero = pi/2, rotation = "clock",cex = 1, tick = T)
axis.circular(at = circular(mt * (360/365) * (pi/180), units = "radians", modulo = "2pi", rotation = "clock",  zero = pi/2), units =  "radians", labels = month.abb, zero = pi/2, rotation = "clock",cex = 1, tick = F, tcl.text = 0.35)

gs1.summ <- vector()

for (i in 1:nrow(m5)) {
  if (m5$FILTER2[i] == 1) {
    gs1 <- sites.info[i][[1]][[2]]
    gs1 <- circular(gs1 * (360/365) * (pi/180), units = "radians", modulo = "2pi", rotation = "clock",  zero = pi/2)
    points.circular(gs1, stack = T, col = lat.color.1[match(round(m5$LAT[i], digits = 1),ref)])
    gs1.summ <- c(gs1.summ, gs1)
  }
}

par(mar = c(0,0,0,0), fig = c(0.8,0.83,0.3,0.75), new = T)

image(t(matrix(1:length(ref))), col = lat.color.2, axes = F, ylab = "test")
lat.color <- round(sort(m5[which(m5$FILTER2 == 1),]$LAT), digits = 1)
axis(4, at = seq(0,1, length.out = 15), labels =  round(seq(lat.color[1], lat.color[length(lat.color)], length.out = 15), 1), cex.axis = 1, las = 1)
mtext("Latitude", side = 3, line = 0.5)
dev.off()



## Summary one growing seasons ----

sub.gs.1 <- m5[which(m5$n.cycles == 1 & m5$FILTER2 == 1),]
mean.site <- vector()
sd.site <- vector()

for (i in 1:nrow(sub.gs.1)) {
  index <- as.numeric(rownames(sub.gs.1)[i])
  gs1 <- sites.info[index][[1]][[2]]
  gs1 <- circular(gs1 * (360/365) * (pi/180), units = "radians", modulo = "2pi", rotation = "clock",  zero = pi/2)
  mean.site[i] <- mean.circular(gs1) / (360/365) / (pi/180)
  sd.site[i] <- sd.circular(gs1) / (360/365) / (pi/180)
}


sub.gs.1$mean.site <- mean.site
sub.gs.1$sd.site <- sd.site
#write.csv(sub.gs.1, "data/sub.gs.1.csv", row.names = F)

## Summary Two growing seasons ----

mt <- month(seq(as.Date("2005/1/1"), as.Date("2005/12/31"), "days"))

mt <- match(1:12, mt)
color <- brewer_pal("qual", palette = 2)(2)
dev.off()
sub.gs.2 <- m5[which(m5$n.cycles == 2 & m5$FILTER2 == 1), ]



for (i in 1:nrow(sub.gs.2)) {
  index <- as.numeric(rownames(sub.gs.2)[i])
  gs1 <- sites.info[index][[1]][[2]][which(sites.info[index][[1]][[5]] == 1)]
  gs2 <- sites.info[index][[1]][[2]][which(sites.info[index][[1]][[5]] == 2)]
  
  gs1 <- circular(gs1 * (360/365) * (pi/180), units = "radians", modulo = "2pi", rotation = "clock",  zero = pi/2)
  gs2 <- circular(gs2 * (360/365) * (pi/180), units = "radians", modulo = "2pi", rotation = "clock",  zero = pi/2)
  
  print(paste("site:", m5$SITE_ID[index], "mean1:", mean.circular(gs1) / (360/365) / (pi/180), "sd1:", sd.circular(gs1) / (360/365) / (pi/180),  "mean2:", mean.circular(gs2) / (360/365) / (pi/180), "sd2:", sd.circular(gs2) / (360/365) / (pi/180)))
  
  plot(gs1, axes = F, type = "n", cex = 1, bin = 720, stack = T, sep = 0.02, shrink = 1.3)
  axis.circular(at = circular(seq(1,365, by = 75) * (360/365) * (pi/180), units = "radians", modulo = "2pi", rotation = "clock",  zero = pi/2), units =  "radians", labels = seq(1,365, by = 75), zero = pi/2, rotation = "clock",cex = 0.7, tick = T)
  axis.circular(at = circular(mt * (360/365) * (pi/180), units = "radians", modulo = "2pi", rotation = "clock",  zero = pi/2), units =  "radians", labels = month.abb, zero = pi/2, rotation = "clock",cex = 0.7, tick = F, tcl.text = 0.3)
  points.circular(gs1, stack = T, col = color[1], pch = 1, cex = 1.2)
  arrows.circular(mean.circular(gs1), lwd = 3, length = 0.1, lty = 1, shrink = 0.6, col = color[1])
  points.circular(gs2, stack = T, col = color[2], pch = 1, cex = 1.2)
  arrows.circular(mean.circular(gs2), lwd = 3, length = 0.1, lty = 1, shrink = 0.6, col = color[2])
}

## Extracting Climatic data ----
which(is.element(which(m5$FILTER2 == 1), sites.changed) == T)
## temperature

half_time.g2 <- read.csv("results/gppmax/sensitivity_half_time_g2.csv")

gpp.tair  <- array(list(list()))
for (st in 1:nrow(m5)) {
  if (m5$FILTER2[st] == 1) {
    if (is.null(sites.info[[st]]) == FALSE) {
      gpp.tair.data  <- time.serie[st][[1]][[3]]
      gpp.nt.vut.ref.time <- names(gpp.tair.data)
      n.cycles <- m5$n.cycles[st]
      window.size <- 365
      gpp.tair.tempora.o  <- vector()
      gpp.tair.tempora.w  <- vector()
      years.site.o <- unique(year(time.serie[st][[1]][[2]]))
      for (i in 1:length(sites.info[st][[1]][[4]])) {
        end.t <- which(as.character(gpp.nt.vut.ref.time) ==
                         as.character(as.Date(sites.info[st][[1]][[4]][i])))
        if (year(as.Date(sites.info[st][[1]][[4]][i])) == years.site.o[1]) {
          if (end.t < window.size) {
            start.t <- 1
          }else{
            start.t <- end.t - window.size
          }
        }else{
          start.t <- end.t - window.size
        }
        tempora <- gpp.tair.data[start.t:end.t]
        tair.o  <- mean(gpp.tair.data[start.t:end.t], na.rm = T)
        tempora <- gpp.tair.data[start.t:end.t]
        gpp.tair.tempora.o  <- c(gpp.tair.tempora.o, tair.o)
        
        if (m5$n.cycles[st] == 2) {
          if (sites.info[st][[1]][[5]][i] == 1) {
            v.weight <- half.time(init = 1, t = length(tempora[!is.na(tempora)]), t_half = half_time.g2$half.time[which(as.character(half_time.g2$Site.name) == as.character(m5$SITE_ID[st]) & half_time.g2$n.cycle == 1)])
          }else {
            v.weight <- half.time(init = 1, t = length(tempora[!is.na(tempora)]), t_half = half_time.g2$half.time[which(as.character(half_time.g2$Site.name) == as.character(m5$SITE_ID[st]) & half_time.g2$n.cycle == 2)])
          }
        }else{
          v.weight <- half.time(init = 1, t = length(tempora[!is.na(tempora)]), t_half = m5$half.time[st])
        }
        sum.v.weight <- sum(v.weight)
        tair.w <- (sum(tempora[!is.na(tempora)] * v.weight, na.rm = T)/sum.v.weight)
        gpp.tair.tempora.w <- c(gpp.tair.tempora.w, tair.w)
        
      }
      # original data
      gpp.tair[st][[1]][[1]]  <- gpp.tair.tempora.o
      # weighted data
      gpp.tair[st][[1]][[2]]  <- gpp.tair.tempora.w
    }
  }
}

# Precipitation

gpp.precip  <- array(list(list()))
for (st in 1:nrow(m5)) {
  if (m5$FILTER2[st] == 1) {
    if (is.null(sites.info[[st]]) == FALSE) {
      gpp.precip.data  <- time.serie[st][[1]][[4]]
      gpp.nt.vut.ref.time <- names(gpp.precip.data)
      n.cycles <- m5$n.cycles[st]
      window.size <- 365
      gpp.precip.tempora.o  <- vector()
      gpp.precip.tempora.w  <- vector()
      years.site.o <- unique(year(time.serie[st][[1]][[2]]))
      for (i in 1:length(sites.info[st][[1]][[4]])) {
        end.t <- which(as.character(gpp.nt.vut.ref.time) ==
                         as.character(as.Date(sites.info[st][[1]][[4]][i])))
        if (year(as.Date(sites.info[st][[1]][[4]][i])) == years.site.o[1]) {
          if (end.t < window.size) {
            start.t <- 1
          }else{
            start.t <- end.t - window.size
          }
        }else{
          start.t <- end.t - window.size
        }
        tempora <- gpp.precip.data[start.t:end.t]
        precip.o  <- sum(gpp.precip.data[start.t:end.t], na.rm = T)
        tempora <- gpp.precip.data[start.t:end.t]
        gpp.precip.tempora.o  <- c(gpp.precip.tempora.o, precip.o)
        
        if (m5$n.cycles[st] == 2) {
          if (sites.info[st][[1]][[5]][i] == 1) {
            v.weight <- half.time(init = 1, t = length(tempora[!is.na(tempora)]), t_half = half_time.g2$half.time[which(as.character(half_time.g2$Site.name) == as.character(m5$SITE_ID[st]) & half_time.g2$n.cycle == 1)])
          }else {
            v.weight <- half.time(init = 1, t = length(tempora[!is.na(tempora)]), t_half = half_time.g2$half.time[which(as.character(half_time.g2$Site.name) == as.character(m5$SITE_ID[st]) & half_time.g2$n.cycle == 2)])
          }
        }else{
          v.weight <- half.time(init = 1, t = length(tempora[!is.na(tempora)]), t_half = m5$half.time[st])
        }
        
        sum.v.weight <- sum(v.weight)
        precip.w <- (sum(tempora[!is.na(tempora)] * v.weight, na.rm = T)/sum.v.weight)
        gpp.precip.tempora.w <- c(gpp.precip.tempora.w, precip.w)
      }
      # original data
      gpp.precip[st][[1]][[1]]  <- gpp.precip.tempora.o
      # weighted data
      gpp.precip[st][[1]][[2]]  <- gpp.precip.tempora.w
    }
  }
}

# Shortway incomming radiation
gpp.swrad  <- array(list(list()))
for (st in 1:nrow(m5)) {
  if (m5$FILTER2[st] == 1) {
    if (is.null(sites.info[[st]]) == FALSE) {
      gpp.swrad.data  <- time.serie[st][[1]][[5]]
      gpp.nt.vut.ref.time <- names(gpp.swrad.data)
      n.cycles <- m5$n.cycles[st]
      window.size <- 365
      gpp.swrad.tempora.o  <- vector()
      gpp.swrad.tempora.w  <- vector()
      years.site.o <- unique(year(time.serie[st][[1]][[2]]))
      for (i in 1:length(sites.info[st][[1]][[4]])) {
        end.t <- which(as.character(gpp.nt.vut.ref.time) ==
                         as.character(as.Date(sites.info[st][[1]][[4]][i])))
        if (year(as.Date(sites.info[st][[1]][[4]][i])) == years.site.o[1]) {
          if (end.t < window.size) {
            start.t <- 1
          }else{
            start.t <- end.t - window.size
          }
        }else{
          start.t <- end.t - window.size
        }
        tempora <- gpp.swrad.data[start.t:end.t]
        swrad.o  <- mean(gpp.swrad.data[start.t:end.t], na.rm = T)
        tempora <- gpp.swrad.data[start.t:end.t]
        gpp.swrad.tempora.o  <- c(gpp.swrad.tempora.o, swrad.o)
        if (m5$n.cycles[st] == 2) {
          if (sites.info[st][[1]][[5]][i] == 1) {
            v.weight <- half.time(init = 1, t = length(tempora[!is.na(tempora)]), t_half = half_time.g2$half.time[which(as.character(half_time.g2$Site.name) == as.character(m5$SITE_ID[st]) & half_time.g2$n.cycle == 1)])
          }else {
            v.weight <- half.time(init = 1, t = length(tempora[!is.na(tempora)]), t_half = half_time.g2$half.time[which(as.character(half_time.g2$Site.name) == as.character(m5$SITE_ID[st]) & half_time.g2$n.cycle == 2)])
          }
        }else{
          v.weight <- half.time(init = 1, t = length(tempora[!is.na(tempora)]), t_half = m5$half.time[st])
        }
        
        sum.v.weight <- sum(v.weight)
        swrad.w <- (sum(tempora[!is.na(tempora)] * v.weight, na.rm = T)/sum.v.weight)
        gpp.swrad.tempora.w <- c(gpp.swrad.tempora.w, swrad.w)
        
      }
      # original data
      gpp.swrad[st][[1]][[1]]  <- gpp.swrad.tempora.o
      # weighted data
      gpp.swrad[st][[1]][[2]]  <- gpp.swrad.tempora.w
    }
  }
}

# VPD
gpp.vpd  <- array(list(list()))
for (st in 1:nrow(m5)) {
  if (m5$FILTER2[st] == 1) {
    if (is.null(sites.info[[st]]) == FALSE) {
      gpp.vpd.data  <- time.serie[st][[1]][[7]]
      gpp.nt.vut.ref.time <- names(gpp.vpd.data)
      n.cycles <- m5$n.cycles[st]
      window.size <- 365
      gpp.vpd.tempora.o  <- vector()
      gpp.vpd.tempora.w  <- vector()
      years.site.o <- unique(year(time.serie[st][[1]][[2]]))
      for (i in 1:length(sites.info[st][[1]][[4]])) {
        end.t <- which(as.character(gpp.nt.vut.ref.time) ==
                         as.character(as.Date(sites.info[st][[1]][[4]][i])))
        if (year(as.Date(sites.info[st][[1]][[4]][i])) == years.site.o[1]) {
          if (end.t < window.size) {
            start.t <- 1
          }else{
            start.t <- end.t - window.size
          }
        }else{
          start.t <- end.t - window.size
        }
        tempora <- gpp.vpd.data[start.t:end.t]
        vpd.o  <- mean(gpp.vpd.data[start.t:end.t], na.rm = T)
        tempora <- gpp.vpd.data[start.t:end.t]
        gpp.vpd.tempora.o  <- c(gpp.vpd.tempora.o, vpd.o)
        
        if (m5$n.cycles[st] == 2) {
          if (sites.info[st][[1]][[5]][i] == 1) {
            v.weight <- half.time(init = 1, t = length(tempora[!is.na(tempora)]), t_half = half_time.g2$half.time[which(as.character(half_time.g2$Site.name) == as.character(m5$SITE_ID[st]) & half_time.g2$n.cycle == 1)])
          }else {
            v.weight <- half.time(init = 1, t = length(tempora[!is.na(tempora)]), t_half = half_time.g2$half.time[which(as.character(half_time.g2$Site.name) == as.character(m5$SITE_ID[st]) & half_time.g2$n.cycle == 2)])
          }
        }else{
          v.weight <- half.time(init = 1, t = length(tempora[!is.na(tempora)]), t_half = m5$half.time[st])
        }
        
        sum.v.weight <- sum(v.weight)
        vpd.w <- (sum(tempora[!is.na(tempora)] * v.weight, na.rm = T)/sum.v.weight)
        gpp.vpd.tempora.w <- c(gpp.vpd.tempora.w, vpd.w)
      }
      # original data
      gpp.vpd[st][[1]][[1]]  <- gpp.vpd.tempora.o
      # weighted data
      gpp.vpd[st][[1]][[2]]  <- gpp.vpd.tempora.w
    }
    
  }
}

# reformuling filter 2
m5$FILTER2[which(is.na(m5$n.cycles) == T)] <- 0

m5$n.cycles[which(is.na(m5$n.cycles) == T)] <- 0
## Multiple linear regressions and sd ----
options(warn = 2)
trew <- array(list(list()))
trew.1.1 <- array(list(list()))
trew.1.2 <- array(list(list()))

# loading global climatic space

global.climate.space <- readRDS("global_climatic_space.RDS")


# Circular - Linear regression ----
set.seed(1234)
counter <- 1
for (st in 1:nrow(m5)) {
  if (m5$FILTER2[st] == 1 & m5$n.cycles[st] == 1) {
    if (is.null(sites.info[[st]]) == FALSE) {
      if (is.element(st, sites.changed)) {
        gpp.nt.vut.ref <- time.serie[st][[1]][[1]]
        gpp.nt.vut.ref.time <- year(time.serie[st][[1]][[6]])
        hours.re <- time.serie[st][[1]][[6]]
      }else{
        gpp.nt.vut.ref <- time.serie[st][[1]][[1]]
        gpp.nt.vut.ref.time <- year(time.serie[st][[1]][[2]])
        hours.re <- time.serie[st][[1]][[2]]
      }
      #st <- 1
      print(paste(m5$SITE_ID[st], (counter * 100)/length(which(m5$FILTER2 == 1)), "%"))
      counter <- counter + 1
      timing  <- sites.info[st][[1]][[2]]
      timing.circ  <- circular(timing * (360/365) * (pi/180),
                               units = "radians", modulo = "2pi",
                               rotation = "counter", zero = 0 )

      
      new.climatic <- global.climate.space[which(global.climate.space$site_names == m5$SITE_ID[st]), 1:2]
      
      new.climatic <- as.matrix(new.climatic)

      # bootstraping
      bot <- 1000
      temporal.bot <- matrix(NA, ncol = 6, nrow = bot)
      temporal.bot.l <- matrix(NA, ncol = 5, nrow = bot)
      for (bot.i in 1:bot) {
        index <- tapply(1:length(timing), rep(1:length(unique(year(hours.re))), each = 10), sample, size = 1)
        timing.b <- timing[index]
        climatic.b <- new.climatic[index,]
        timing.circ.b  <- circular(timing.b * (360/365) * (pi/180),units = "radians", modulo = "2pi", rotation = "counter", zero = 0)
        circres <- withTimeout(try(lm.circular(y = timing.circ.b, x = climatic.b,
                                               init = c(0, 0), type = "c-l", verbose = F), silent = T), timeout = 10,
                               onTimeout = "warning")

        if ((class(circres) == "try-error") || is.null(circres)) {
          temporal.bot[bot.i, 1] <- NA
          temporal.bot[bot.i, 2:3] <- NA
          temporal.bot[bot.i, 4] <- NA
          temporal.bot[bot.i, 5:6] <- NA
        } else {
          temporal.bot[bot.i, 1] <- circres$mu
          temporal.bot[bot.i, 2:3] <- circres$coef
          temporal.bot[bot.i, 4] <- circres$kappa
          temporal.bot[bot.i, 5:6] <- circres$p.values
        }

        circres.l <- lm(as.numeric(timing.circ.b) ~ climatic.b)

        temporal.bot.l[bot.i, 1] <- circres.l$coefficients[1]
        temporal.bot.l[bot.i, 2:3] <- circres.l$coefficients[2:3]
        temporal.bot.l[bot.i, 4:5] <-  summary(circres.l)$coefficients[-1,4]

      }

      trew[st][[1]][[1]] <- temporal.bot # raw data

      trew[st][[1]][[2]] <-  apply(temporal.bot[,2:3], MARGIN = 2, mean, na.rm = T) # averages (coefficients)
      trew[st][[1]][[3]] <- apply(temporal.bot[,2:3], MARGIN = 2, FUN = sd, na.rm = T) # sd  coefficients)
      trew[st][[1]][[4]] <- mean.circular(circular(temporal.bot[,1],units = "radians", modulo = "2pi", rotation = "counter", zero = 0), na.rm = T)  # mu mean

      trew[st][[1]][[5]] <- sd.circular(circular(temporal.bot[,1],units = "radians", modulo = "2pi", rotation = "counter", zero = 0), na.rm = T) # mu sd
      trew[st][[1]][[6]] <-  mean(temporal.bot[,4], na.rm = T) # kappa average
      trew[st][[1]][[7]] <- sd(temporal.bot[,4], na.rm = T) # kappa sd
      trew[st][[1]][[8]] <- sd(timing, na.rm = T) # deviation of the number of days (linear)

      trew[st][[1]][[9]] <- apply(temporal.bot[,5:6], MARGIN = 2, median, na.rm = T) # averages (p.values)
      # standard error
      trew[st][[1]][[10]] <- trew[st][[1]][[3]] / length(na.omit(temporal.bot[,2:3]))

      ## Multiple linear regression information

      trew[st][[1]][[11]] <- temporal.bot.l
      # average of coefficients
      trew[st][[1]][[12]] <- apply(temporal.bot.l[,2:3], MARGIN = 2, mean, na.rm = T)
      # sd of coefficients
      trew[st][[1]][[13]] <- apply(temporal.bot.l[,2:3], MARGIN = 2, FUN = sd, na.rm = T)
      # mu mean
      trew[st][[1]][[14]] <- mean(temporal.bot[,1], na.rm = T)
      # mu sd
      trew[st][[1]][[15]] <- sd(temporal.bot[,1], na.rm = T)
      # median p.values
      trew[st][[1]][[16]] <- apply(temporal.bot[,5:6], MARGIN = 2, median, na.rm = T)
      # standard error
      trew[st][[1]][[16]] <- trew[st][[1]][[3]] / length(na.omit(temporal.bot.l[,2:3]))
    }

  } else if (m5$FILTER2[st] == 1 & m5$n.cycles[st] == 2) {
    if (is.null(sites.info[[st]]) == FALSE) {
      for (d in 1:2) {
        if (is.element(st, sites.changed)) {
          gpp.nt.vut.ref <- time.serie[st][[1]][[1]]
          gpp.nt.vut.ref.time <- year(time.serie[st][[1]][[6]])
          hours.re <- time.serie[st][[1]][[6]]
        }else{
          gpp.nt.vut.ref <- time.serie[st][[1]][[1]]
          gpp.nt.vut.ref.time <- year(time.serie[st][[1]][[2]])
          hours.re <- time.serie[st][[1]][[2]]
        }
        #st <- 1
        print(paste(m5$SITE_ID[st], (counter * 100)/length(which(m5$FILTER2 == 1)), "%"))
        counter <- counter + 1
        timing  <- sites.info[st][[1]][[2]][which(sites.info[st][[1]][[5]] == d)]
        timing.circ  <- circular(timing * (360/365) * (pi/180),
                                 units = "radians", modulo = "2pi",
                                 rotation = "counter", zero = 0)
        
        new.climatic <- global.climate.space[which(global.climate.space$site_names == m5$SITE_ID[st]), 1:2]
        
        if (d == 1) {
          new.climatic <- as.matrix(new.climatic[1:(nrow(new.climatic)/2), ])
        }else{
          new.climatic <- as.matrix(new.climatic[((nrow(new.climatic)/2) + 1):nrow(new.climatic), ])
        }

        

        # bootstraping
        bot <- 1000
        temporal.bot <- matrix(NA, ncol = 6, nrow = bot)
        temporal.bot.l <- matrix(NA, ncol = 5, nrow = bot)
        for (bot.i in 1:bot) {
          index <- tapply(1:length(timing), rep(1:length(unique(year(hours.re))), each = 10), sample, size = 1)
          timing.b <- timing[index]
          climatic.b <- new.climatic[index,]
          timing.circ.b  <- circular(timing.b * (360/365) * (pi/180),units = "radians", modulo = "2pi", rotation = "counter", zero = 0)
          circres <- withTimeout(try(lm.circular(y = timing.circ.b, x = climatic.b,
                                                 init = c(0, 0), type = "c-l", verbose = F), silent = T), timeout = 10,
                                 onTimeout = "warning")

          if ((class(circres) == "try-error") || is.null(circres)) {
            temporal.bot[bot.i, 1] <- NA
            temporal.bot[bot.i, 2:5] <- NA
            temporal.bot[bot.i, 6] <- NA
            temporal.bot[bot.i, 7:10] <- NA
          } else {
            temporal.bot[bot.i, 1] <- circres$mu
            temporal.bot[bot.i, 2:3] <- circres$coef
            temporal.bot[bot.i, 4] <- circres$kappa
            temporal.bot[bot.i, 5:6] <- as.numeric(circres$p.values)
          }

          circres.l <- lm(as.numeric(timing.circ.b) ~ climatic.b)

          temporal.bot.l[bot.i, 1] <- circres.l$coefficients[1]
          temporal.bot.l[bot.i, 2:3] <- circres.l$coefficients[2:3]
          temporal.bot.l[bot.i, 4:5] <-  summary(circres.l)$coefficients[-1,4]
        }
        if (d == 1) {
          # raw data
          trew.1.1[st][[1]][[1]] <- temporal.bot
          # averages (coefficients)
          trew.1.1[st][[1]][[2]] <-  apply(temporal.bot[,2:3], MARGIN = 2, mean, na.rm = T)
          # sd  coefficients)
          trew.1.1[st][[1]][[3]] <- apply(temporal.bot[,2:3], MARGIN = 2, FUN = sd, na.rm = T)
          # mu mean
          trew.1.1[st][[1]][[4]] <- mean.circular(circular(temporal.bot[,1],units = "radians", modulo = "2pi", rotation = "counter", zero = 0), na.rm = T)
          # mu sd
          trew.1.1[st][[1]][[5]] <- sd.circular(circular(temporal.bot[,1],units = "radians", modulo = "2pi", rotation = "counter", zero = 0), na.rm = T)
          # kappa average
          trew.1.1[st][[1]][[6]] <-  mean(temporal.bot[,4], na.rm = T)
          # kappa sd
          trew.1.1[st][[1]][[7]] <- sd(temporal.bot[,4], na.rm = T)
          # deviation of the number of days (linear)
          trew.1.1[st][[1]][[8]] <- sd(timing, na.rm = T)
          # averages (p.values)
          trew.1.1[st][[1]][[9]] <- apply(temporal.bot[,5:6], MARGIN = 2, median, na.rm = T)
          trew.1.1[st][[1]][[10]] <- trew.1.1[st][[1]][[3]] / length(na.omit(temporal.bot[,2:3]))

          ## Multiple linear regression information

          trew.1.1[st][[1]][[11]] <- temporal.bot.l
          # average of coefficients
          trew.1.1[st][[1]][[12]] <- apply(temporal.bot.l[,2:3], MARGIN = 2, mean, na.rm = T)
          # sd of coefficients
          trew.1.1[st][[1]][[13]] <- apply(temporal.bot.l[,2:3], MARGIN = 2, FUN = sd, na.rm = T)
          # mu mean
          trew.1.1[st][[1]][[14]] <- mean(temporal.bot[,1], na.rm = T)
          # mu sd
          trew.1.1[st][[1]][[15]] <- sd(temporal.bot[,1], na.rm = T)
          # median p.values
          trew.1.1[st][[1]][[16]] <- apply(temporal.bot[,5:6], MARGIN = 2, median, na.rm = T)
          # standard error
          trew.1.1[st][[1]][[16]] <- trew[st][[1]][[3]] / length(na.omit(temporal.bot.l[,2:3]))

        }else{
          # raw data
          trew.1.2[st][[1]][[1]] <- temporal.bot
          # averages (coefficients)
          trew.1.2[st][[1]][[2]] <-  apply(temporal.bot[,2:3], MARGIN = 2, mean, na.rm = T)
          # sd  coefficients)
          trew.1.2[st][[1]][[3]] <- apply(temporal.bot[,2:3], MARGIN = 2, FUN = sd, na.rm = T)
          # mu mean
          trew.1.2[st][[1]][[4]] <- mean.circular(circular(temporal.bot[,1],units = "radians", modulo = "2pi", rotation = "counter", zero = 0), na.rm = T)
          # mu sd
          trew.1.2[st][[1]][[5]] <- sd.circular(circular(temporal.bot[,1],units = "radians", modulo = "2pi", rotation = "counter", zero = 0), na.rm = T)
          # kappa average
          trew.1.2[st][[1]][[6]] <-  mean(temporal.bot[,4], na.rm = T)
          # kappa sd
          trew.1.2[st][[1]][[7]] <- sd(temporal.bot[,4], na.rm = T)
          # deviation of the number of days (linear)
          trew.1.2[st][[1]][[8]] <- sd(timing, na.rm = T)
          # averages (p.values)
          trew.1.2[st][[1]][[9]] <- apply(temporal.bot[,5:6], MARGIN = 2, median, na.rm = T)
          trew.1.2[st][[1]][[10]] <- trew.1.2[st][[1]][[3]] / length(na.omit(temporal.bot[,2:3]))

          ## Multiple linear regression information

          trew.1.2[st][[1]][[11]] <- temporal.bot.l
          # average of coefficients
          trew.1.2[st][[1]][[12]] <- apply(temporal.bot.l[,2:5], MARGIN = 2, mean, na.rm = T)
          # sd of coefficients
          trew.1.2[st][[1]][[13]] <- apply(temporal.bot.l[,2:5], MARGIN = 2, FUN = sd, na.rm = T)
          # mu mean
          trew.1.2[st][[1]][[14]] <- mean(temporal.bot[,1], na.rm = T)
          # mu sd
          trew.1.2[st][[1]][[15]] <- sd(temporal.bot[,1], na.rm = T)
          # average p.values
          trew.1.2[st][[1]][[16]] <- apply(temporal.bot[,5:6], MARGIN = 2, median, na.rm = T)
          # standard error
          trew.1.2[st][[1]][[16]] <- trew[st][[1]][[3]] / length(na.omit(temporal.bot.l[,2:3]))
        }
      }
    }
  }
}


#
#
# saveRDS(trew, file = "trew_1_global_gs.rds")
# saveRDS(trew.1.1, file = "trew_2_1_global_gs.rds")
# saveRDS(trew.1.2, file = "trew_2_2_global_gs.rds")

trew <- readRDS(file = "trew_1_gs.rds")
trew.1.1 <- readRDS(file = "trew_2_1_gs.rds")
trew.1.2 <- readRDS(file = "trew_2_2_gs.rds")

## summary coefficients 1 growing season ----
color <- brewer_pal("qual", palette = 2)(2)
g1 <- trew[1][[1]][[2]] / sum(abs(trew[1][[1]][[2]]))

coef.comp <- as.data.frame(trew[1][[1]][[1]][,2:3])
# colnames(coef.comp) = c("tair", "swrad", "precip", "vpd")
colnames(coef.comp) = c("PC1", "precip")
coef.comp <- melt(coef.comp)

col <- brewer_pal(type = "qual", palette = 2)(1)
show_col(col)

coef <- ggplot(data = coef.comp, aes(x = variable, y = value, colour = col[1])) +
  geom_boxplot() +
  geom_hline(yintercept = 0) +
  scale_color_brewer(type = "qual", palette = 2, aesthetics = "colour") +
  theme(legend.position = "none")

coef

coef.comp <- as.data.frame(trew[1][[1]][[1]][,5:6])
# colnames(coef.comp) = c("tair", "swrad", "precip", "vpd")
colnames(coef.comp) = c("PC1", "precip")
coef.comp <- melt(coef.comp)


p.values <- ggplot(data = coef.comp, aes(x = variable, y = value, colour = col[1])) +
  geom_boxplot() +
  geom_hline(yintercept = 0.05) +
  scale_color_brewer(type = "qual", palette = 2) +
  scale_color_brewer(type = "qual", palette = 2, aesthetics = "colour") +
  theme(legend.position = "none")

p.values
ggarrange(coef, p.values, ncol = 1, nrow = 2)

# Supplement 3 (EC site plots) ----

dev.off()
sub.gs.1 <- m5[which(m5$n.cycles == 1 & m5$FILTER2 == 1), ]

for (i in 1:nrow(sub.gs.1)) {
  index <- as.numeric(rownames(sub.gs.1)[i])
  # png(paste("doc/papers/doy_gppmax/images/one_gs_comp/", i,"_", m5$SITE_ID[index], ".png", sep = ""), width = 800, height = 550)
  # pdf(paste("doc/papers/doy_gppmax/images/one_gs_comp/", i,"_", m5$SITE_ID[index], ".pdf", sep = ""), paper = "a4r", width = 800, height = 550)
  par(mfrow = c(1,2))
  gs1 <- sites.info[index][[1]][[2]]

  gs1 <- circular(gs1 * (360/365) * (pi/180), units = "radians", modulo = "2pi", rotation = "counter",  zero = 0)

  print(paste("site:", index, m5$SITE_ID[index], "mean1:", mean.circular(gs1) / (360/365) / (pi/180), "sd1:", sd.circular(gs1) / (360/365) / (pi/180)))

  plot.circular(gs1, axes = F, type = "n", cex = 1, bin = 720, stack = T, sep = 0.02, shrink = 1.3, main = m5$SITE_ID[index], rotation = "clock", zero = pi/2)
  axis.circular(at = circular(seq(1,365, by = 75) * (360/365) * (pi/180), units = "radians", modulo = "2pi", rotation = "clock",  zero = pi/2), units =  "radians", labels = seq(1,365, by = 75), zero = pi/2, rotation = "clock",cex = 1, tick = T)
  axis.circular(at = circular(mt * (360/365) * (pi/180), units = "radians", modulo = "2pi", rotation = "clock",  zero = pi/2), units =  "radians", labels = month.abb, zero = pi/2, rotation = "clock",cex = 1, tick = F, tcl.text = 0.4)
  points.circular(gs1, stack = T, col = color[1], pch = 1, cex = 1.2, rotation = "clock", zero = pi/2)
  arrows.circular(mean.circular(gs1), lwd = 3, length = 0.1, lty = 1, shrink = 0.4, col = color[1], rotation = "clock", zero = pi/2)
  mtext("a)", side = 3, adj = 0, cex = 1.5)

  coef.comp <- as.data.frame(trew[index][[1]][[1]][,2:3])
  # colnames(coef.comp) = c("Tair", "SWin", "Precip", "VPD")
  colnames(coef.comp) = c("Tair, SWin, VPD", "Precip")

  coef.comp <- melt(coef.comp)
  coef.comp <- na.omit(coef.comp)

  coef <- ggplot(data = coef.comp, aes(x = variable, y = value, colour = color[1])) +
    geom_boxplot() +
    geom_hline(yintercept = 0) +
    scale_color_brewer(type = "qual", palette = 2, aesthetics = "colour") +
    theme(plot.title = element_text(hjust = 0, size = 17, face = "plain")) +
    labs(color = "Growing season") +
    ylab("Coefficient") +
    xlab("") +
    ggtitle("b)") +
    theme(legend.position = "none")


  # significance value

  coef.comp <- as.data.frame(trew[index][[1]][[1]][,5:6])
  # colnames(coef.comp) = c("Tair", "SWin", "Precip", "VPD")
  colnames(coef.comp) = c("Tair, SWin, VPD", "Precip")

  coef.comp <- melt(coef.comp)
  coef.comp <- na.omit(coef.comp)

  p.values <- ggplot(data = coef.comp, aes(x = variable, y = value, colour = color[1])) +
    geom_boxplot() +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
    scale_color_brewer(type = "qual", palette = 2, aesthetics = "colour") +
    theme(plot.title = element_text(hjust = 0, size = 17, face = "plain")) +
    labs(color = "Growing season") +
    ylab("Significance Value") +
    xlab("") +
    ggtitle("c)") +
    theme(legend.position = "none")

  final <- ggarrange(coef, p.values, ncol = 1, nrow = 2)

  plot.new()              ## suggested by @Josh
  vps <- baseViewports()
  pushViewport(vps$figure) ##   I am in the space of the autocorrelation plot
  vp1 <- plotViewport(c(1,0,0,1)) ## create new vp with margins, you play with this values
  print(final,vp = vp1)
  dev.off()
}



sub.gs.1.original <- m5[which(m5$n.cycles == 1 & m5$FILTER2 == 1), ]
sub.gs.1.original$PC1 <- NA
sub.gs.1.original$PC1.p <- NA
sub.gs.1.original$precip <- NA
sub.gs.1.original$precip.p <- NA

for (i in 1:nrow(sub.gs.1)) {
  index <- which(as.character(m5$SITE_ID) == as.character(sub.gs.1$SITE_ID[i]))
  
  sub.gs.1.original$PC1[i] <- trew[index][[1]][[2]][1]
  sub.gs.1.original$PC1.p[i] <- trew[index][[1]][[9]][1]
  sub.gs.1.original$precip[i] <- trew[index][[1]][[2]][2]
  sub.gs.1.original$precip.p[i] <- trew[index][[1]][[9]][2]
}

## ordering the influence of the parameters

ggplot(data = sub.gs.1.original) +
  geom_point(aes(x = PC1, y = precip, color = IGBP))

intercomp <- sub.gs.1.original

intercomp <- intercomp[,c(17,19)]

influence <- matrix(NA, nrow = nrow(sub.gs.1.original), ncol = 2)

for (i in 1:nrow(influence)) {
  temp.sig <- sub.gs.1.original[i,c(18,20)]
  temp.vec <- sub.gs.1.original[i,c(17,19)]
  temp.vec[which(temp.sig > 0.05)] <- 0
  intercomp[i,which(temp.sig > 0.05)] <- 0
  temp.vec <- abs(temp.vec)
  temp.order <- match(temp.vec, sort(temp.vec, decreasing = T))
  temp.order[which(temp.vec == 0)] <- 0
  influence[i,] <- temp.order
}

intercomp
influence

length(which(intercomp[,1] != 0))

intercomp[which(intercomp[,1] != 0 & intercomp[,2] != 0),]

length(which(intercomp[,1] < 0 & intercomp[,2] == 0))
intercomp[which(intercomp[,1] == 0 & intercomp[,2] == 0),]




## summary coefficients 2 growing seasons ----

g1 <- trew.1.1[68][[1]][[2]] / sum(abs(trew.1.1[68][[1]][[2]]))

g2 <- trew.1.2[68][[1]][[2]] / sum(abs(trew.1.2[68][[1]][[2]]))



coef.comp <- as.data.frame(trew.1.1[68][[1]][[1]][,2:3])
colnames(coef.comp) = c("PC1", "precip")

coef.comp <- melt(coef.comp)

coef.comp$growing.season <-  "g1"

coef.comp.2 <- as.data.frame(trew.1.2[68][[1]][[1]][,2:3])
colnames(coef.comp.2) = c("PC1", "precip")

coef.comp.2 <- melt(coef.comp.2)

coef.comp.2$growing.season <-  "g2"

coef.comp <- rbind(coef.comp, coef.comp.2)

coef <- ggplot(data = coef.comp, aes(x = variable, y = value, color = growing.season)) +
  geom_boxplot() +
  geom_hline(yintercept = 0) +
  scale_color_brewer(type = "qual", palette = 2)

coef

# significance value

coef.comp <- as.data.frame(trew.1.1[68][[1]][[1]][,5:6])
colnames(coef.comp) = c("PC1","precip")

coef.comp <- melt(coef.comp)

coef.comp$growing.season <-  "g1"

coef.comp.2 <- as.data.frame(trew.1.2[68][[1]][[1]][,5:6])
colnames(coef.comp.2) = c("PC1", "precip")

coef.comp.2 <- melt(coef.comp.2)

coef.comp.2$growing.season <-  "g2"

coef.comp <- rbind(coef.comp, coef.comp.2)

p.values <- ggplot(data = coef.comp, aes(x = variable, y = value, color = growing.season)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.05) +
  scale_color_brewer(type = "qual", palette = 2)

ggarrange(coef, p.values, ncol = 1, nrow = 2, common.legend = T, legend = "right")


color <- brewer_pal("qual", palette = 2)(2)


dev.off()
sub.gs.2 <- m5[which(m5$n.cycles == 2 & m5$FILTER2 == 1), ]

# 
# 
# for (i in 1:nrow(sub.gs.2)) {
#   index <- as.numeric(rownames(sub.gs.2)[i])
#   # png(paste("doc/papers/doy_gppmax/images/two_gs_comp/", i,"_", m5$SITE_ID[index], ".png", sep = ""), width = 800, height = 550)
#   pdf(paste("doc/papers/doy_gppmax/images/two_gs_comp/", i,"_", m5$SITE_ID[index], ".pdf", sep = ""), paper = "a4r", width = 800, height = 550)
#   par(mfrow = c(1,2))
#   gs1 <- sites.info[index][[1]][[2]][which(sites.info[index][[1]][[5]] == 1)]
#   gs2 <- sites.info[index][[1]][[2]][which(sites.info[index][[1]][[5]] == 2)]
#   
#   gs1 <- circular(gs1 * (360/365) * (pi/180), units = "radians", modulo = "2pi", rotation = "clock",  zero = pi/2)
#   gs2 <- circular(gs2 * (360/365) * (pi/180), units = "radians", modulo = "2pi", rotation = "clock",  zero = pi/2)
#   
#   print(paste("site:", m5$SITE_ID[index], "mean1:", mean.circular(gs1) / (360/365) / (pi/180), "sd1:", sd.circular(gs1) / (360/365) / (pi/180),  "mean2:", mean.circular(gs2) / (360/365) / (pi/180), "sd2:", sd.circular(gs2) / (360/365) / (pi/180)))
#   
#   plot(gs1, axes = F, type = "n", cex = 1, bin = 720, stack = T, sep = 0.02, shrink = 1.3, main = m5$SITE_ID[index])
#   axis.circular(at = circular(seq(1,365, by = 75) * (360/365) * (pi/180), units = "radians", modulo = "2pi", rotation = "clock",  zero = pi/2), units =  "radians", labels = seq(1,365, by = 75), zero = pi/2, rotation = "clock",cex = 1, tick = T)
#   axis.circular(at = circular(mt * (360/365) * (pi/180), units = "radians", modulo = "2pi", rotation = "clock",  zero = pi/2), units =  "radians", labels = month.abb, zero = pi/2, rotation = "clock",cex = 1, tick = F, tcl.text = 0.4)
#   points.circular(gs1, stack = T, col = color[1], pch = 1, cex = 1.2)
#   arrows.circular(mean.circular(gs1), lwd = 3, length = 0.1, lty = 1, shrink = 0.4, col = color[1])
#   points.circular(gs2, stack = T, col = color[2], pch = 1, cex = 1.2)
#   arrows.circular(mean.circular(gs2), lwd = 3, length = 0.1, lty = 1, shrink = 0.4, col = color[2])
#   mtext("a)", side = 3, adj = 0, cex = 1.5)
#   
#   coef.comp <- as.data.frame(trew.1.1[index][[1]][[1]][,2:3])
#   colnames(coef.comp) = c("Tair, SWin, VPD", "Precip")
#   
#   coef.comp <- melt(coef.comp)
#   
#   coef.comp$growing.season <-  "g1"
#   
#   coef.comp.2 <- as.data.frame(trew.1.2[index][[1]][[1]][,2:3])
#   colnames(coef.comp.2) = c("Tair, SWin, VPD", "Precip")
#   
#   coef.comp.2 <- melt(coef.comp.2)
#   
#   coef.comp.2$growing.season <-  "g2"
#   
#   coef.comp <- rbind(coef.comp, coef.comp.2)
#   
#   coef <- ggplot(data = coef.comp, aes(x = variable, y = value, color = growing.season)) +
#     geom_boxplot() +
#     geom_hline(yintercept = 0) +
#     labs(color = "Growing season") +
#     scale_color_brewer(type = "qual", palette = 2) +
#     theme(plot.title = element_text(hjust = 0, size = 17, face = "plain")) +
#     ylab("Coefficient") +
#     xlab("") + ggtitle("b)")
#   
#   
#   # significance value
#   
#   coef.comp <- as.data.frame(trew.1.1[index][[1]][[1]][,5:6])
#   colnames(coef.comp) = c("Tair, SWin, VPD", "Precip")
#   
#   coef.comp <- melt(coef.comp)
#   
#   coef.comp$growing.season <-  "g1"
#   
#   coef.comp.2 <- as.data.frame(trew.1.2[index][[1]][[1]][,5:6])
#   colnames(coef.comp.2) = c("Tair, SWin, VPD", "Precip")
#   
#   coef.comp.2 <- melt(coef.comp.2)
#   
#   coef.comp.2$growing.season <-  "g2"
#   
#   coef.comp <- rbind(coef.comp, coef.comp.2)
#   
#   p.values <- ggplot(data = coef.comp, aes(x = variable, y = value, color = growing.season)) +
#     geom_boxplot() +
#     geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
#     labs(color = "Growing season") +
#     scale_color_brewer(type = "qual", palette = 2) +
#     theme(plot.title = element_text(hjust = 0, size = 17, face = "plain")) +
#     ylab("Significance Value") +
#     xlab("") + ggtitle("c)")
#   
#   
#   final <- ggarrange(coef, p.values, ncol = 1, nrow = 2, common.legend = T, legend = "right")
#   
#   plot.new()              ## suggested by @Josh
#   vps <- baseViewports()
#   pushViewport(vps$figure) ##   I am in the space of the autocorrelation plot
#   vp1 <- plotViewport(c(1,0,0,1)) ## create new vp with margins, you play with this values
#   print(final,vp = vp1)
#   dev.off()
# }

# fiting -----
# timing
rsqrt <- vector()
name <- vector()
for (st in 1:nrow(m5)) {
  if (m5$FILTER2[st] == 1 & m5$n.cycles[st] == 1) {
    if (is.null(sites.info[[st]]) == F) {
      timing  <- sites.info[st][[1]][[2]]
      timing.circ  <- circular(timing * (360/365) * (pi/180),
                               units = "radians", modulo = "2pi",
                               rotation = "clock", zero = pi/2 )
      climatic  <- matrix(NA, nrow = length(sites.info[st][[1]][[1]]), ncol = 4)
      climatic[,1]  <- scal(gpp.tair[st][[1]][[2]])
      climatic[,2]  <- scal(gpp.swrad[st][[1]][[2]])
      climatic[,3]  <- scal(gpp.precip[st][[1]][[2]])
      climatic[,4] <- scal(gpp.vpd[st][[1]][[2]])
      sim <- vector()
      for (i in 1:nrow(climatic)) {
        x <- climatic[i,]
        sim <- c(sim, as.numeric(circ.recons(x, mu = as.numeric(trew[st][[1]][[4]]), coef = trew[st][[1]][[2]])))
      }
      #sim <- sim / (360/365) / (pi/180)
      test <- circular(sim, units = "radians", modulo = "2pi",
                       rotation = "clock", zero = pi/2 )
      ex <- cor.circular(timing.circ, test, test = T)$cor
      plot(as.numeric(timing.circ))
      points(as.numeric(sim), col = "red")
      rsqrt <- c(rsqrt, cor.circular(timing.circ, test, test = T)$cor)
      name <- c(name, st)
      # png(paste("data/fluxnet/fluxnet_circular_orig_vs_predicted/gppsat/", st, ".png", sep = ""),
      # width = 1056, height = 668)
      # plot.circular(timing.circ, stack = T, shrink = 2, main = paste(st, m5$SITE_ID[st], "\n",
      #                                                                "cor.circular = ", ex), axes = F)
      # axis.circular(at = circular(seq(1,365, by = 75) * (360/365) * (pi/180), units = "radians", modulo = "2pi", rotation = "clock",  zero = pi/2),
      #               units =  "radians", labels = seq(1,365, by = 75), zero = pi/2, rotation = "clock",
      #               cex = 0.7, tick = T)
      # points(circular(sim, units = "radians", modulo = "2pi",
      #                 rotation = "clock", zero = pi/2), col = "red", pch = 21, stack = T )
      # arrows.circular(trew[st][[1]][[4]],
      #                 rho.circular(timing.circ) - 0.1, lwd = 1, length = 0.1, col = 2, lty = 1)
      # dev.off()
    }
  }
}



for (st in 1:nrow(m5)) {
  if (m5$FILTER2[st] == 1) {
    magnitude  <- scal(sites.info[st][[1]][[1]])
    climatic  <- matrix(NA, nrow = length(sites.info[st][[1]][[1]]), ncol = 4)
    climatic[,1]  <- scal(gpp.tair[st][[1]][[2]])
    climatic[,2]  <- scal(gpp.swrad[st][[1]][[2]])
    climatic[,3]  <- scal(gpp.precip[st][[1]][[2]])
    climatic[,4] <- scal(gpp.vpd[st][[1]][[2]])
    sim <- vector()
    for (i in 1:nrow(climatic)) {
      x <- climatic[i,]
      sim <- c(sim, magnitude.recons(x, intersect =  coef(trew[st][[1]][[2]], s = "lambda.min")[1],
                                     coef = coef(trew[st][[1]][[2]], s = "lambda.min")[-1]))
    }
    #sim <- sim / (360/365) / (pi/180)
    rsqrt.2 <- rsquare(true = magnitude, predicted = sim)
    rsqrt <- c(rsqrt, rsquare(true = magnitude, predicted = sim))
    name <- c(name, st)
    # png(paste("data/fluxnet/fluxnet_magni_orig_vs_predicted/", st, ".png", sep = ""),
    #     width = 1056, height = 668)
    # plot(magnitude, main = paste(st, m5$SITE_ID[st], "\n", "R2=", rsqrt.2))
    # points(sim, col = "red", pch = 21)
    # dev.off()
    
  }
}


# recovering timing (Empirical), circular vs linear (Figure 4) ----
#

color <- brewer_pal("qual", palette = 2)(3)


sites.eg <- c(1, 23)

pdf("circular_vs_linear_empirical.pdf", width = 9, height = 6.5)
# png("circular_vs_linear_empirical.png", width = 900, height = 650)
par(mfrow = c(1,2), mar = c(1, 1, 1, 1))
for (st in sites.eg) {
  if (m5$FILTER2[st] == 1) {
    # circular
    timing  <- sites.info[st][[1]][[2]]
    timing.circ  <- circular(timing * (360/365) * (pi/180),
                             units = "radians", modulo = "2pi",
                             rotation = "counter", zero = 0)
    climatic  <- matrix(NA, nrow = length(sites.info[st][[1]][[1]]), ncol = 4)
    climatic[,1]  <- scal(gpp.tair[st][[1]][[2]])
    climatic[,2]  <- scal(gpp.swrad[st][[1]][[2]])
    climatic[,3]  <- scal(gpp.precip[st][[1]][[2]])
    climatic[,4] <- scal(gpp.vpd[st][[1]][[2]])
    # pca
    #
    climatic_pca <- prcomp(climatic[,c(1,2,4)], center = T, scale. = T)
    new.climatic <- matrix(NA, nrow = length(sites.info[st][[1]][[1]]), ncol = 2)
    new.climatic[,1] <- scale(climatic_pca$x[,1], center = T, scale = T)
    new.climatic[,2] <- climatic[,3]
    
    
    sim <- vector()
    for (i in 1:nrow(climatic)) {
      x <- new.climatic[i,]
      sim <- c(sim, as.numeric(circ.recons(x, mu = as.numeric(trew[st][[1]][[4]]), coef = trew[st][[1]][[2]])))
    }
    
    test <- circular(sim, units = "radians", modulo = "2pi",
                     rotation = "counter", zero = 0)
    ex <- cor.circular(timing.circ, test, test = T)$cor
    
    rsqrt <- c(ex, cor.circular(timing.circ, test, test = T)$cor)
    name <- c(name, st)
    
    # linear
    magnitude  <- as.numeric(timing.circ)
    climatic  <- matrix(NA, nrow = length(sites.info[st][[1]][[1]]), ncol = 4)
    climatic[,1]  <- scal(gpp.tair[st][[1]][[2]])
    climatic[,2]  <- scal(gpp.swrad[st][[1]][[2]])
    climatic[,3]  <- scal(gpp.precip[st][[1]][[2]])
    climatic[,4] <- scal(gpp.vpd[st][[1]][[2]])
    
    # pca
    #
    climatic_pca <- prcomp(climatic[,c(1,2,4)], center = T, scale. = T)
    new.climatic <- matrix(NA, nrow = length(sites.info[st][[1]][[1]]), ncol = 2)
    new.climatic[,1] <- scale(climatic_pca$x[,1], center = T, scale = T)
    new.climatic[,2] <- climatic[,3]
    
    sim.2 <- vector()
    for (i in 1:nrow(climatic)) {
      x <- new.climatic[i,]
      
      sim.2 <- c(sim.2, magnitude.recons(x, intersect =  trew[st][[1]][[14]],
                                         coef = trew[st][[1]][[12]]))
    }
    test <- circular(sim.2, units = "radians", modulo = "2pi",
                     rotation = "counter", zero = 0)
    rsqrt.2 <- cor.circular(timing.circ, test, test = T)$cor
    l.cor <- cor(as.numeric(timing.circ), sim.2)
    name <- c(name, st)
    
    plot.circular(circular(timing.circ, units = "radians", modulo = "2pi",
                           rotation = "clock", zero = pi/2), stack = T, shrink = 1.3, axes = F,
                  col = color[1], cex = 1.2, pch = 1)
    axis.circular(at = circular(seq(1,365, by = 75) * (360/365) * (pi/180), units = "radians", modulo = "2pi", rotation = "clock",  zero = pi/2), units =  "radians", labels = seq(1,365, by = 75), zero = pi/2, rotation = "clock",cex = 1, tick = T)
    axis.circular(at = circular(mt * (360/365) * (pi/180), units = "radians", modulo = "2pi", rotation = "clock",  zero = pi/2), units =  "radians", labels = month.abb, zero = pi/2, rotation = "clock",cex = 1, tick = F, tcl.text = 0.35)
    
    points(circular(sim, units = "radians", modulo = "2pi",
                    rotation = "clock", zero = pi/2), col = color[2], pch = 1, stack = T, cex = 1.2 )
    
    points(circular(sim.2, units = "radians", modulo = "2pi",
                    rotation = "clock", zero = pi/2), col = color[3], pch = 1, stack = T, cex = 1.2 )
    arrows.circular(mean.circular(timing.circ), lwd = 3, length = 0.1, lty = 1, shrink = 0.6, col = color[1],
                    rotation = "clock", zero = pi/2)
    
    
    
    
    mtext(paste("cor.circular (JS) = ", round(ex, digits = 3), "\n", "cor.linear (Pearson) =", round(l.cor, digits = 3)), side = 1, line = -3.5)
    if (st == 1) {
      let <- "a)"
      legend(x =  0.8, y = -0.6, legend = c("Original", "Circular", "Linear"), col = color, pch = 19, bty = "n")
    }else{
      legend(x =  0.8, y = -0.6, legend = c("Original", "Circular", "Linear"), col = color, pch = 19, bty = "n")
      let <- "b)"
    }
    mtext(let, side = 3, adj = 0, line = -9.2, cex = 2, font = 2)
  }
}
dev.off()
