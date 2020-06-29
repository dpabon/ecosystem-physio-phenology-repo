# libraries ----
setwd("~/phd/")
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
library(scales)
library(reshape2)
library(rworldmap)
library(mice)
library(Kendall)
library(trend)
library(glmnet)
library(dimRed)
library(animation)
library(R.utils)
library(gridBase)
library(brew)
library(ggpubr)
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

# r^2
rsquare <- function(true, predicted) {
  sse <- sum((predicted - true)^2)
  sst <- sum((true - mean(true))^2)
  rsq <- 1 - (sse / sst)
  
  # For this post, impose floor...
  #if (rsq < 0) rsq <- 0
  
  return(rsq)
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


## tools
source("bin/gpp/tools.R")

# base map
wmap <- readOGR(dsn = "data/tools/shapefiles/ne_110m_land")
# convert to dataframe
wmap_robin <- spTransform(wmap, CRS("+proj=robin"))
wmap_df_robin <- fortify(wmap_robin)

theme_opts <- list(theme(panel.grid.minor = element_blank(),
                         panel.grid.major = element_blank(),
                         panel.background = element_blank(),
                         plot.background = element_rect(fill = "white"),
                         panel.border = element_blank(),
                         axis.line = element_blank(),
                         axis.text.x = element_blank(),
                         axis.text.y = element_blank(),
                         axis.ticks = element_blank(),
                         axis.title.x = element_blank(),
                         axis.title.y = element_blank(),
                         plot.title = element_blank()))



#### Selecting only sites with at least 7 years of duration ----
m5 <- read.csv("m5.csv")

View(m5)


## All the data indexed in a single array ----
time.serie <- array(list(list()))

## quality filter value 
qc.level <- 2


#st <- 23
for (st in 1:nrow(m5)) {
  if (m5$FILTER2[st] == 1) {
    folder <- "FLUXNET_sites/F15/"
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

## Sites distribution -----

world <- map_data("world")
map <- ggplot() + geom_polygon(data = world, aes(x = long, y = lat, group = group)) + theme_map() +
  coord_fixed(1.3) +
  geom_point(data = m5[which(m5$FILTER2 == 1),], aes(x = LON, y = LAT, color = IGBP), size = 1) 
map

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
      # Variance explained of the seassonal patterns
      # Measure of the confidence that really exist a growing season. (TO INCORPORATE IN THE RESULTS)
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

#write.csv(m5, "data/fluxnet/m5.csv", row.names = F)

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
        #y <- 1 # to remove
        windows.size  <- 90 / n.cycles
        gpp.temp <- all.years[,y]
        
        for (d in 1:n.cycles) {
          if ((gpp.day.template[d]  - windows.size) <= 0) {
            back  <- abs(gpp.day.template[d]  - windows.size)
            new.vector  <- c((length(all.years[,y]) - back):length(all.years[,y]),
                             1:(gpp.day.template[d] + windows.size))
            timing <- new.vector[order(gpp.temp[new.vector], decreasing = T)[1:10]] #new.vector[which.max(gpp.temp[new.vector])] # bug fixed
            gpp.timing  <- c(gpp.timing, timing)
            gp  <- gpp.temp[timing]
            gpp.max  <- c(gpp.max, gp)
            year.t  <- rep(y, times = length(timing))
            gpp.year <- c(gpp.year, year.t)
          } else if ((gpp.day.template[d]  + windows.size) > length(all.years[,y])) {
            advance  <-  (gpp.day.template[d]  + windows.size) - length(all.years[,y])
            new.vector  <- c((gpp.day.template[d] - windows.size):length(all.years[,y]), 1:advance)
            timing <- new.vector[order(gpp.temp[new.vector], decreasing = T)[1:10]] #new.vector[which.max(gpp.temp[new.vector])]
            gpp.timing  <- c(gpp.timing, timing)
            gp  <- gpp.temp[timing]
            gpp.max  <- c(gpp.max, gp)
            year.t  <- rep(y, times = length(timing))
            gpp.year <- c(gpp.year, year.t)
          } else {
            start  <- gpp.day.template[d] - windows.size
            end  <- gpp.day.template[d] + windows.size
            new.vector  <- start:end
            timing <- new.vector[order(gpp.temp[new.vector], decreasing = T)[1:10]] #new.vector[which.max(gpp.temp[new.vector])] 
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

?points.circular


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


## CROSS Validation ----

###CROSS DATA ----
m5 <- read.csv("data/fluxnet/m5.csv")
vegetation <- "MF"

summary(m5$IGBP[which(m5$FILTER2 == 1)])

m5[-which(m5$FILTER2 == 1 & m5$IGBP == vegetation),"FILTER2"] <- 0

half.time.DBF <- round(mean(m5$half.time[which(m5$FILTER2 == 1)]), digits = 0)

# all are DBF half.time
m5$half.time <- half.time.DBF


## Extracting Climatic data ----
half_time.g2 <- read.csv("results/gppmax/sensitivity_half_time_g2.csv")

gpp.tair  <- array(list(list()))
for (st in 1:nrow(m5)) {
  #st <- 1
  if (m5$FILTER2[st] == 1) {
    if (is.null(sites.info[[st]]) == FALSE) {
      gpp.tair.data  <- time.serie[st][[1]][[3]]
      gpp.nt.vut.ref.time <- names(gpp.tair.data)
      #hours.re <- time.serie[st][[1]][[2]][seq(1, length(time.serie[st][[1]][[2]]),
      #by = 48)]
      n.cycles <- m5$n.cycles[st]
      window.size <- 365
      gpp.tair.tempora.o  <- vector()
      gpp.tair.tempora.w  <- vector()
      years.site.o <- unique(year(time.serie[st][[1]][[2]]))
      for (i in 1:length(sites.info[st][[1]][[4]])) {
        #i <- 108
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
      #hours.re <- time.serie[st][[1]][[2]][seq(1, length(time.serie[st][[1]][[2]]),
      #by = 48)]
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
      #hours.re <- time.serie[st][[1]][[2]][seq(1, length(time.serie[st][[1]][[2]]),
      #by = 48)]
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
      #hours.re <- time.serie[st][[1]][[2]][seq(1, length(time.serie[st][[1]][[2]]),
      #by = 48)]
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
# 
# 
# 
# ## Cross validation per vegetation type ----

dbf.data <- m5[which(m5$FILTER2 == 1), ]

# preparing the data

site.id <- rep(m5$SITE_ID[which(m5$FILTER2 == 1)], m5$RANGE_ANALYSIS[which(m5$FILTER2 == 1)] * 10)

year <- vector()

for (i in 1:nrow(dbf.data)) {
  temp <- rep(1:dbf.data$RANGE_ANALYSIS[i], each = 10)
  year <- c(year, temp)
}

year

# tair <- vector()
PC1 <- vector()
precip <- vector()
# swir <- vector()
# vpd <- vector()
DOYmax <- vector()
options(warn = 2)

for (i in 1:nrow(m5)) {
  if (m5$FILTER2[i] == 1) {
    
    to_reduce <- matrix(data = NA, nrow = length(gpp.tair[i][[1]][[2]]), ncol = 3)
    to_reduce[,1] <- gpp.tair[i][[1]][[2]]
    to_reduce[,2] <- gpp.swrad[i][[1]][[2]]
    to_reduce[,3] <- gpp.vpd[i][[1]][[2]]
    climatic_pca <- prcomp(to_reduce, center = T, scale. = T)
    
    PC1 <- c(PC1, climatic_pca$x[,1])
    precip <- c(precip, gpp.precip[i][[1]][[2]])
    DOYmax  <- c(DOYmax, sites.info[i][[1]][[2]])
  }
}


cross.dbf <- data.frame(site.id, year, DOYmax, PC1, precip)
cross.dbf$index.id <- paste0(cross.dbf$site.id, cross.dbf$year)
cross.dbf$index.out <- 1:nrow(cross.dbf)

# saveRDS(cross.dbf, file = paste("data/cross_validation_",vegetation ,".RDS", sep = ""))


## DBF-CROSS ----
m5 <- read.csv("data/fluxnet/m5.csv")
dbf.lat <- m5$LAT[which(m5$FILTER2 == 1 & m5$IGBP == "DBF")]

dbf.lat <- rep(dbf.lat, m5$RANGE_ANALYSIS[which(m5$FILTER2 == 1 & m5$IGBP == "DBF") ] * 10)

cross.dbf <- readRDS("data/cross_validation_DBF.RDS")

cross.dbf.predict <- readRDS("data/cross_dbf_predicted.RDS")

observed <- circular(cross.dbf$DOYmax * (360/365) * (pi/180),units = "radians", modulo = "2pi", rotation = "counter", zero = 0)

predicted <- circular(cross.dbf.predict,units = "radians", modulo = "2pi", rotation = "counter", zero = 0)

cor.dbf <- cor.circular(predicted, observed, test = T)

observed <- as.numeric(observed) / (360/365) / (pi/180)
predicted <- as.numeric(predicted) / (360/365) / (pi/180)
input.gg <- data.frame(observed, predicted, dbf.lat)

#lab <- paste("italic(R) ^ 2 ==", round(cor.dbf$cor, digits = 2))
lab <- paste("JS.Cor(cv) =", round(cor.dbf$cor, digits = 2))

scatter.dbf <- ggplot(data = input.gg, aes(x = predicted, y = observed)) +
  geom_abline(slope = 1, intercept = 0, col = "red", size = 1.5) +
  geom_point(aes(colour = dbf.lat), alpha = 0.6) +
  theme_gray() + 
  #annotate("text", x = (120 + 20), y = 259, label = lab, size = 5 ) +
  # geom_smooth() +
  scale_color_viridis(name = "Latitude") +
  ylab(expression("Observed "*"(DOY"["GPPmax"]*")")) + 
  xlab(expression("Predicted "*"(DOY"["GPPmax"]*")")) +
  xlim(120,260) +
  ylim(120, 260) +
  theme_bw() +
  ggtitle(paste("DBF", lab, sep = "  ")) 

scatter.dbf

VT <- m5[which(m5$IGBP == "DBF" & m5$FILTER2 == 1),]

# Robinson transformation
places_robin_df <- project(cbind(VT$LON, VT$LAT), proj="+init=esri:54030") 

places_robin_df <- as.data.frame(places_robin_df)
names(places_robin_df) <- c("LONGITUDE", "LATITUDE")



map.dbf <- ggplot(data = wmap_df_robin, aes(long,lat, group = group, fill = hole)) +
  geom_polygon() + 
  geom_point(data = places_robin_df, aes(LONGITUDE, LATITUDE, group = NULL, fill = NULL), color = "#F4511E", size = 1.5) +
  coord_equal() + 
  theme_opts +
  scale_fill_manual(values = c("#262626", "white"), guide = "none") # change colors & remove legend

ggarrange(scatter.dbf, map.dbf, widths = c(1,1.5))

# plot.circular(observed, stack = T, shrink = 1.8)
# points.circular(predicted, stack = T, shrink = 1.8, pch = 21, col = "red")



## EBF-CROSS  ----
m5 <- read.csv("data/m5.csv")

VT <- m5[which(m5$IGBP == "EBF" & m5$FILTER2 == 1),]
ebf.lat <- m5$LAT[which(m5$FILTER2 == 1 & m5$IGBP == "EBF")]
ebf.lat <- rep(ebf.lat, m5$RANGE_ANALYSIS[which(m5$FILTER2 == 1 & m5$IGBP == "EBF") ] * 10)

observed.ebf <- readRDS("data/cross_validation_EBF.RDS")
observed.ebf <- circular(observed.ebf$DOYmax * (360/365) * (pi/180),units = "radians", modulo = "2pi", rotation = "counter", zero = 0)

predicted.ebf <- readRDS("data/cross_ebf_predicted.RDS")
predicted.ebf <- circular(predicted.ebf, units = "radians", modulo = "2pi", rotation = "counter", zero = 0)

cor.ebf <- cor.circular(predicted.ebf, observed.ebf, test = T)

plot.circular(observed.ebf, stack = T)
points.circular(predicted.ebf, stack = T, pch = 21, col = "red")

#predicted.ebf[which(predicted.ebf > 4)] <- predicted.ebf[which(predicted.ebf > 4)] - (2*pi)
#observed.ebf[observed.ebf > 4] <- observed.ebf[observed.ebf > 4] - (2*pi) 

predicted.ebf <- predicted.ebf / (360/365) / (pi/180)
observed.ebf <- observed.ebf / (360/365) / (pi/180)

input.gg <- data.frame(predicted.ebf, observed.ebf, ebf.lat)

lab <- paste("JS.Cor(cv) =", round(cor.ebf$cor, digits = 2))

scatter.ebf <- ggplot(data = input.gg, aes(x = predicted.ebf, y = observed.ebf)) +
  geom_abline(slope = 1, intercept = 0, col = "red", size = 1.5) +
  geom_point(aes(colour = ebf.lat), alpha = 0.6) +
  theme_gray() + 
  #annotate("text", x = (0 + 50), y = 364, label = lab, size = 5) +
  # geom_smooth() +
  scale_color_viridis(name = "Latitude") +
  ylab(expression("Observed "*"(DOY"["GPPmax"]*")")) + 
  xlab(expression("Predicted "*"(DOY"["GPPmax"]*")")) +
  xlim(0,365) +
  ylim(0, 365) +
  theme_bw() +
  ggtitle(paste("EBF", lab, sep = "  "))

scatter.ebf

# Robinson transformation
places_robin_df <- project(cbind(VT$LON, VT$LAT), proj="+init=esri:54030") 

places_robin_df <- as.data.frame(places_robin_df)
names(places_robin_df) <- c("LONGITUDE", "LATITUDE")

wmap_robin <- spTransform(wmap, CRS("+proj=robin"))
wmap_df_robin <- fortify(wmap_robin)

map.ebf <- ggplot(data = wmap_df_robin, aes(long,lat, group=group, fill=hole)) + 
  geom_polygon() + 
  geom_point(data = places_robin_df, aes(LONGITUDE, LATITUDE, group = NULL, fill = NULL), color="#F4511E", size = 1) +
  coord_equal() + 
  theme_opts +
  scale_fill_manual(values = c("#262626", "#e6e8ed"), guide="none") # change colors & remove legend


ggarrange(scatter.ebf, map.ebf)


## ENF-CROSS ----
m5 <- read.csv("data/m5.csv")

VT <- m5[which(m5$IGBP == "ENF" & m5$FILTER2 == 1),]
enf.lat <- m5$LAT[which(m5$FILTER2 == 1 & m5$IGBP == "ENF")]
enf.lat <- rep(enf.lat, m5$RANGE_ANALYSIS[which(m5$FILTER2 == 1 & m5$IGBP == "ENF") ] * 10)

observed.enf <- readRDS("data/cross_validation_ENF.RDS")
observed.enf <- circular(observed.enf$DOYmax * (360/365) * (pi/180),units = "radians", modulo = "2pi", rotation = "counter", zero = 0)

predicted.enf <- readRDS("data/cross_enf_predicted.RDS")
predicted.enf <- circular(predicted.enf, units = "radians", modulo = "2pi", rotation = "counter", zero = 0)

cor.enf <- cor.circular(predicted.enf, observed.enf, test = T)

plot.circular(observed.enf, stack = T)
points.circular(predicted.enf, stack = T, pch = 21, col = "red")

predicted.enf <- predicted.enf / (360/365) / (pi/180)
observed.enf <- observed.enf / (360/365) / (pi/180)

input.gg <- data.frame(predicted.enf, observed.enf, enf.lat)

lab <- paste("JS.Cor(cv) =", round(cor.enf$cor, digits = 2))

scatter.enf <- ggplot(data = input.gg, aes(x = predicted.enf, y = observed.enf)) +
  geom_abline(slope = 1, intercept = 0, col = "red", size = 1.5) +
  geom_point(aes(colour = enf.lat), alpha = 0.6) +
  theme_gray() + 
  #annotate("text", x = (80 + 30), y = 274, label = lab, size = 5) +
  # geom_smooth() +
  scale_color_viridis(name = "Latitude") +
  ylab(expression("Observed "*"(DOY"["GPPmax"]*")")) + 
  xlab(expression("Predicted "*"(DOY"["GPPmax"]*")")) +
  xlim(80,275) +
  ylim(80, 275) +
  theme_bw() +
  ggtitle(paste("ENF", lab, sep = "  "))

scatter.enf

# Robinson transformation
places_robin_df <- project(cbind(VT$LON, VT$LAT), proj="+init=esri:54030") 

places_robin_df <- as.data.frame(places_robin_df)
names(places_robin_df) <- c("LONGITUDE", "LATITUDE")

wmap_robin <- spTransform(wmap, CRS("+proj=robin"))
wmap_df_robin <- fortify(wmap_robin)

map.enf <- ggplot(data = wmap_df_robin, aes(long,lat, group = group, fill = hole)) + 
  geom_polygon() + 
  geom_point(data = places_robin_df, aes(LONGITUDE, LATITUDE, group = NULL, fill = NULL), color = "#F4511E", size = 1) +
  coord_equal() + 
  theme_opts +
  scale_fill_manual(values = c("#262626", "#e6e8ed"), guide="none") # change colors & remove legend


ggarrange(scatter.enf, map.enf)

## GRAS-CROSS ----
m5 <- read.csv("data/m5.csv")

VT <- m5[which(m5$IGBP == "GRA" & m5$FILTER2 == 1),]
gras.lat <- m5$LAT[which(m5$FILTER2 == 1 & m5$IGBP == "GRA")]
gras.lat <- rep(gras.lat, m5$RANGE_ANALYSIS[which(m5$FILTER2 == 1 & m5$IGBP == "GRA") ] * 10)

observed.gras <- readRDS("data/cross_validation_GRA.RDS")
observed.gras <- circular(observed.gras$DOYmax * (360/365) * (pi/180),units = "radians", modulo = "2pi", rotation = "counter", zero = 0)

predicted.gras <- readRDS("data/cross_gra_predicted.RDS")
predicted.gras <- circular(predicted.gras, units = "radians", modulo = "2pi", rotation = "counter", zero = 0)

cor.gras <- cor.circular(predicted.gras, observed.gras, test = T)

plot.circular(observed.gras, stack = T, shrink = 1.2)
points.circular(predicted.gras, stack = T, pch = 21, col = "red")

predicted.gras <- predicted.gras / (360/365) / (pi/180)
observed.gras <- observed.gras / (360/365) / (pi/180)

input.gg <- data.frame(predicted.gras, observed.gras, gras.lat)

lab <- paste("JS.Cor(cv) =", round(cor.gras$cor, digits = 2))

scatter.gras <- ggplot(data = input.gg, aes(x = predicted.gras, y = observed.gras)) +
  geom_abline(slope = 1, intercept = 0, col = "red", size = 1.5) +
  geom_point(aes(colour = gras.lat), alpha = 0.6) +
  theme_gray() + 
  #annotate("text", x = (20 + 30), y = (307 - 1), label = lab, size = 5) +
  # geom_smooth() +
  scale_color_viridis(name = "Latitude") +
  ylab(expression("Observed "*"(DOY"["GPPmax"]*")")) + 
  xlab(expression("Predicted "*"(DOY"["GPPmax"]*")")) +
  xlim(20,307) +
  ylim(20,307) +
  theme_bw() +
  ggtitle(paste("GRAS", lab, sep = "  "))

scatter.gras

# Robinson transformation
places_robin_df <- project(cbind(VT$LON, VT$LAT), proj="+init=esri:54030") 

places_robin_df <- as.data.frame(places_robin_df)
names(places_robin_df) <- c("LONGITUDE", "LATITUDE")

wmap_robin <- spTransform(wmap, CRS("+proj=robin"))
wmap_df_robin <- fortify(wmap_robin)

map.gras <- ggplot(data = wmap_df_robin, aes(long,lat, group = group, fill = hole)) + 
  geom_polygon() + 
  geom_point(data = places_robin_df, aes(LONGITUDE, LATITUDE, group = NULL, fill = NULL), color="#F4511E", size = 1) +
  coord_equal() + 
  theme_opts +
  scale_fill_manual(values = c("#262626", "#e6e8ed"), guide="none") # change colors & remove legend


ggarrange(scatter.gras, map.gras)


## MF-CROSS ----
m5 <- read.csv("data/fluxnet/m5.csv")

VT <- m5[which(m5$IGBP == "MF" & m5$FILTER2 == 1),]
mf.lat <- m5$LAT[which(m5$FILTER2 == 1 & m5$IGBP == "MF")]
mf.lat <- rep(mf.lat, m5$RANGE_ANALYSIS[which(m5$FILTER2 == 1 & m5$IGBP == "MF") ] * 10)

observed.mf <- readRDS("data/cross_validation_MF.RDS")
observed.mf <- circular(observed.mf$DOYmax * (360/365) * (pi/180),units = "radians", modulo = "2pi", rotation = "counter", zero = 0)

predicted.mf <- readRDS("data/cross_mf_predicted.RDS")
predicted.mf <- circular(predicted.mf, units = "radians", modulo = "2pi", rotation = "counter", zero = 0)

cor.mf <- cor.circular(predicted.mf, observed.mf, test = T)

plot.circular(observed.mf, stack = T, shrink = 1.2)
points.circular(predicted.mf, stack = T, pch = 21, col = "red")

predicted.mf <- predicted.mf / (360/365) / (pi/180)
observed.mf <- observed.mf / (360/365) / (pi/180)

input.gg <- data.frame(predicted.mf, observed.mf, mf.lat)

lab <- paste("JS.Cor(cv) =", round(cor.mf$cor, digits = 2))

scatter.mf <- ggplot(data = input.gg, aes(x = predicted.mf, y = observed.mf)) +
  geom_abline(slope = 1, intercept = 0, col = "red", size = 1.5) +
  geom_point(aes(colour = mf.lat), alpha = 0.6) +
  theme_gray() + 
  #annotate("text", x = (120 + 20), y = 267 - 1, label = lab, size = 5) +
  # geom_smooth() +
  scale_color_viridis(name = "Latitude") +
  ylab(expression("Observed "*"(DOY"["GPPmax"]*")")) + 
  xlab(expression("Predicted "*"(DOY"["GPPmax"]*")")) +
  xlim(120,267) +
  ylim(120,267) +
  theme_bw() +
  ggtitle(paste("MF", lab, sep = "  "))
  

scatter.mf

# Robinson transformation
places_robin_df <- project(cbind(VT$LON, VT$LAT), proj="+init=esri:54030") 

places_robin_df <- as.data.frame(places_robin_df)
names(places_robin_df) <- c("LONGITUDE", "LATITUDE")

wmap_robin <- spTransform(wmap, CRS("+proj=robin"))
wmap_df_robin <- fortify(wmap_robin)

map.mf <- ggplot(data = wmap_df_robin, aes(long,lat, group = group, fill = hole)) + 
  geom_polygon() + 
  geom_point(data = places_robin_df, aes(LONGITUDE, LATITUDE, group = NULL, fill = NULL), color= "#F4511E", size = 1) +
  coord_equal() + 
  theme_opts +
  scale_fill_manual(values = c("#262626", "#e6e8ed"), guide = "none") # change colors & remove legend


ggarrange(scatter.mf, map.mf)


## ALL plots ----
ggarrange(scatter.dbf, map.dbf, scatter.ebf, map.ebf, scatter.gras, map.gras, scatter.mf, map.mf, scatter.enf, map.enf, ncol = 2, nrow = 5)

phenocross <- ggarrange(scatter.dbf, scatter.ebf, scatter.gras, scatter.mf, scatter.enf,  ncol = 2, nrow = 3, align = "hv")


phenocross
ggsave(phenocross, filename = "figure_7.png", width = 8, height = 10)
ggsave(phenocross, filename = "figure_7.pdf", width = 8, height = 10)

