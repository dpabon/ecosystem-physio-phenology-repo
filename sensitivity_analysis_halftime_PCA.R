# Upscaling the half time sensitivity analysis Defining cluster sections ----

# name.R n.proc = 15, current job = $i

args <- commandArgs(trailingOnly = TRUE)

chunk2 <- function(x, n) split(x, cut(seq_along(x), n, labels = FALSE))

n.proc <- as.numeric(args[1])
# n.proc <- 16
c.job <- as.numeric(args[2])
# c.job <- 1

u <- chunk2(2:365, n.proc)
section <- as.numeric(u[[c.job]])

print(paste("Partition =", c.job, "my section is from:", section[1], "to:", section[length(section)]))

# Libraries and functions ----
setwd("/Net/Groups/BGI/scratch/dpabon/")
library(methods)
library(lubridate)
library(circular)
library(MASS)
library(glmnet)
library(ncdf4)
library(zoo)
library(R.utils)
#library(bigleaf)

## ggplot2 colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# function to scale the variables
scal <- function(x) {
  y <- (x - mean(x, na.rm = T))/sd(x, na.rm = T)
  return(y)
}

# r^2
rsquare <- function(true, predicted) {
  sse <- sum((predicted - true)^2)
  sst <- sum((true - mean(true))^2)
  rsq <- 1 - (sse/sst)
  
  # For this post, impose floor... if (rsq < 0) rsq <- 0
  
  return(rsq)
}

# half life function
half.time <- function(init = 1, t, t_half = 30) {
  lambda <- log(2)/t_half
  t <- 1:t
  y <- vector()
  y[1] <- init
  for (i in 2:length(t)) {
    y[i] <- init * exp(1)^(-lambda * t[i])
    # y[i] <- init * ((1/2)^((t[i] - 1)/(1/2)))
  }
  y
  return(y)
}

simpleWAI <- function(precip, pet, theta=0.05, awc=100) {
  
  n <- length(precip)
  
  WAI <- rep(NA, n)
  WAI[1] <- awc
  
  ETmodWAI <- rep(NA, n)
  input <- rep(NA, n)
  
  for (i in 1:(n-1)) {
    input[i] <- min(precip[i], awc - WAI[i]) #Recharge = min (precip, water deficit previous state)
    ETsupply <- theta * WAI[i] #Evapotranspiration supply, theta is a parameter from literature
    ETmodWAI[i] <- min(pet[i], ETsupply)
    WAI[i+1] <- WAI[i] + input[i] - ETmodWAI[i] 
    
    if (is.na(WAI[i+1])) 
    {
      WAI[i+1] <- WAI[i]
    }
    
  } 
  
  data.frame(WAI=WAI, ETmodWAI=ETmodWAI, input=input)
}

# prediction using circular model
circ.recons <- function(x, mu, coef) {
  out <- mu + (2 * atan(sum(coef * x)))
  return(out)
}

# Reading and preparing data ---- Selecting only sites with at least nine years
# of duration -
m5 <- read.csv("data/fluxnet/m5.csv")
m5$id <- 1:nrow(m5)
# m5 <- m5[-which(m5$FILTER2 == 0),] View(m5)


# only for B15

## All the data indexed in a single array --
time.serie <- array(list(list()))
## quality filter value 
qc.level <- 2


# st <- 23
for (st in 1:nrow(m5)) {
  if (m5$FILTER2[st] == 1) {
    folder <- "data/fluxnet/FLUXNET_sites/F15/"
    nc.file <- nc_open(paste(folder, m5$SITE_ID[st], ".HH.", m5$DATA_START[st],
                             ".", m5$DATA_END[st], ".nc", sep = ""))
    # GPP
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
    time.serie[st][[1]][[2]] <- as.Date(hours.re, origin = "1582-10-15 00:00",
                                        format = "%Y-%m-%d  %R")
    unique.days <- unique(as.character(time.serie[st][[1]][[2]]))
    time.serie[st][[1]][[2]] <- as.Date(unique.days)
    length(time.serie[st][[1]][[1]]) == length(time.serie[st][[1]][[2]])
    # index days used to estimate the mean and the sum of each day
    index.days <- as.character(as.Date(hours.re, origin = "1582-10-15 00:00",
                                       format = "%Y-%m-%d  %R"))
    
    # Tair
    time.serie[st][[1]][[3]] <- tapply(X = ncvar_get(nc = nc.file, varid = "TA_F",
                                                     start = 1, count = -1), INDEX = index.days, FUN = mean, na.rm = T)
    
    # Precipitation
    time.serie[st][[1]][[4]] <- tapply(X = ncvar_get(nc = nc.file, varid = "P_F",
                                                     start = 1, count = -1), INDEX = index.days, FUN = sum, na.rm = T)
    
    # Swin
    time.serie[st][[1]][[5]] <- tapply(X = ncvar_get(nc = nc.file, varid = "SW_IN_F",
                                                     start = 1, count = -1), INDEX = index.days, FUN = mean, na.rm = T)
    
    # VPD
    time.serie[st][[1]][[7]] <- tapply(X = ncvar_get(nc = nc.file, varid = "VPD_F",
                                                     start = 1, count = -1), INDEX = index.days, FUN = mean, na.rm = T)
    print(paste(st, m5$SITE_ID[st], sep = "-"))
    
  }
}

## detecting incogruencies -


# depuring site by site

rm.year <- function(x, y.) {
  # x Vector with original days.  y. years to remove return a vector with the
  # positions correspondig to the year y.
  mid <- lubridate::year(x)
  if (length(y.) == 1) {
    ret <- which(mid == y.)
  } else {
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
inde[[13]] <- c(1998, 2000, 2001, 2004)
inde[[15]] <- c(2007:2010)
inde[[16]] <- c(1997, 2004:2005)
inde[[18]] <- c(2011)
inde[[20]] <- c(2000, 2009, 2011:2013)
inde[[22]] <- c(2000)
inde[[23]] <- c(2001)
inde[[25]] <- c(1997:1998)
inde[[29]] <- c(2008:2011)
inde[[31]] <- c(2006)
inde[[32]] <- c(2002, 2014)
inde[[33]] <- c(2002)
inde[[37]] <- c(1997:1999, 2009)
inde[[40]] <- c(2013)
inde[[41]] <- c(2013)
inde[[42]] <- c(2013)
inde[[43]] <- c(2014)
inde[[45]] <- c(2014)
inde[[46]] <- c(1995, 1996, 2007)
inde[[54]] <- c(2007:2009)
inde[[60]] <- c(2009)
inde[[62]] <- c(2004)
inde[[65]] <- c(2006:2008)
inde[[68]] <- c(2004)
inde[[73]] <- c(1996, 1999, 2000, 2006)
inde[[78]] <- c(2004, 2005)
inde[[80]] <- c(2000)
inde[[81]] <- c(1996, 1997, 2005)
inde[[82]] <- c(1995, 2004)

## tier2

inde[[9]] <- 2014  #NL-Loo
inde[[31]] <- c(inde[[31]], 2011, 2012, 2013)  #ZA-Kru
inde[[55]] <- c(2014)  #FR-Gri
sites.changed <- c(1, 2, 3, 4, 9, 12, 13, 15, 16, 18, 20, 22, 23, 25, 29, 31, 32,
                   33, 37, 40, 41, 42, 43, 45, 46, 54, 55, 60, 62, 65, 68, 73, 78, 80, 81, 82, 9,
                   31, 55)

for (st in 1:length(inde)) {
  if (is.null(inde[[st]]) == F) {
    if (m5$FILTER2[st] == 1) {
      time.serie[st][[1]][[1]] <- time.serie[st][[1]][[1]][-rm.year(time.serie[st][[1]][[2]],
                                                                    y. = inde[[st]])]
      # index 6 correspond to the time vector modified from vector 2
      time.serie[st][[1]][[6]] <- time.serie[st][[1]][[2]][-rm.year(time.serie[st][[1]][[2]],
                                                                    y. = inde[[st]])]
    }
  }
}

## sites to remove

s.r <- c(38, 72, 74, 75, 76, 77, 83)

for (i in 1:length(s.r)) {
  time.serie[s.r[i]][[1]][[1]] <- NA
}

# number of sites with information


row.names(m5) <- 1:nrow(m5)



# Extracting GPPmax, DOYmax and index of the year ----

options(warn = 2)
# template to extract the gppmax
gpp.day.template.f <- list()
for (st in 1:nrow(m5)) {
  if (m5$FILTER2[st] == 1) {
    if (is.element(st, sites.changed)) {
      gpp.nt.vut.ref <- time.serie[st][[1]][[1]]
      gpp.nt.vut.ref.time <- year(time.serie[st][[1]][[6]])
      hours.re <- time.serie[st][[1]][[6]]
    } else {
      gpp.nt.vut.ref <- time.serie[st][[1]][[1]]
      gpp.nt.vut.ref.time <- year(time.serie[st][[1]][[2]])
      hours.re <- time.serie[st][[1]][[2]]
    }
    l.years <- length(unique(gpp.nt.vut.ref.time))
    u.years <- unique(gpp.nt.vut.ref.time)
    all.years <- matrix(NA, nrow = 365, ncol = l.years)
    
    for (y in 1:l.years) {
      # salida <- tapply(time.serie[st][[1]][[1]][which(gpp.nt.vut.ref.time ==
      # u.years[y])], INDEX = hours.re[which(gpp.nt.vut.ref.time == u.years[y])], FUN =
      # quantile, na.rm = T, probs = 0.9)
      salida <- time.serie[st][[1]][[1]][which(gpp.nt.vut.ref.time == u.years[y])]
      if (leap_year(u.years[y]) == T) {
        salida <- salida[1:365]
      } else {
        #
      }
      all.years[, y] <- salida
    }
    # mean seasonal cycle
    msc <- apply(X = all.years, MARGIN = 1, FUN = mean, na.rm = T)
    if (any(is.nan(msc)) == T) {
      msc[which(is.nan(msc))] <- 0
    }
    if (any(is.na(msc)) == FALSE) {
      # number of cycles
      n.cycles <- as.integer(which.max(abs(fft(msc))[2:4]))
      # Variance explained of the seassonal patterns Measure of the confidence that
      # really exist a growing season. (TO INCORPORATE IN THE RESULTS)
      var_in_cycle <- (2 * sum(abs(fft(msc)^2)[2:4])/sum(abs(fft(msc)[-1]^2)))
      m5[st, "n.cycles"] <- n.cycles
      # timing template
      gpp.day.template <- vector()
      fy <- fft(msc)
      ny <- length(fy)
      fy[c(1, 5:(ny - 3))] <- 0
      yfiltered <- fft(fy, inverse = T)
      yfiltered <- Re(yfiltered)
      
      for (d in 1:365) {
        if (is.na(yfiltered[d])) {
          # pass
        } else {
          if (d == 1) {
            if (yfiltered[d] > yfiltered[365] & yfiltered[d] > yfiltered[(d +
                                                                          1)]) {
              gpp.day.template <- c(gpp.day.template, d)
            }
          } else if (d == 365) {
            if (yfiltered[d] > yfiltered[(d - 1)] & yfiltered[d] > yfiltered[1]) {
              gpp.day.template <- c(gpp.day.template, d)
            }
          } else {
            if (yfiltered[d] > yfiltered[(d - 1)] & yfiltered[d] > yfiltered[(d +
                                                                              1)]) {
              gpp.day.template <- c(gpp.day.template, d)
            }
          }
        }
      }
      # checking that the number of peaks is equal to the number of cycles (only
      # selecting the highest peaks)
      if (length(gpp.day.template) == n.cycles) {
        gpp.day.template.f[[st]] <- gpp.day.template
      } else {
        tempora <- yfiltered[gpp.day.template]
        gpp.day.template <- gpp.day.template[order(tempora, decreasing = T)[1:n.cycles]]
        gpp.day.template.f[[st]] <- gpp.day.template
      }
    }
  }
}
# gpp.day.template.f[[24]] length(gpp.day.template.f)

# extracting gppmax and gppmax timing from original data with quality flags

sites.info <- array(list(list()))
# st <- 23
for (st in 1:nrow(m5)) {
  if (m5$FILTER2[st] == 1) {
    if (is.null(gpp.day.template.f[[st]]) == F) {
      if (is.element(st, sites.changed)) {
        gpp.nt.vut.ref <- time.serie[st][[1]][[1]]
        gpp.nt.vut.ref.time <- year(time.serie[st][[1]][[6]])
        hours.re <- time.serie[st][[1]][[6]]
      } else {
        gpp.nt.vut.ref <- time.serie[st][[1]][[1]]
        gpp.nt.vut.ref.time <- year(time.serie[st][[1]][[2]])
        hours.re <- time.serie[st][[1]][[2]]
      }
      
      l.years <- length(unique(gpp.nt.vut.ref.time))
      u.years <- unique(gpp.nt.vut.ref.time)
      all.years <- matrix(NA, nrow = 365, ncol = l.years)
      
      for (y in 1:l.years) {
        # salida <- tapply(gpp.nt.vut.ref[which(gpp.nt.vut.ref.time == u.years[y])],
        # INDEX = hours.re[which(gpp.nt.vut.ref.time == u.years[y])], FUN = quantile,
        # na.rm = T, probs = 0.9)
        salida <- gpp.nt.vut.ref[which(gpp.nt.vut.ref.time == u.years[y])]
        if (leap_year(u.years[y]) == T) {
          salida <- salida[1:365]
        } else {
          #
        }
        all.years[, y] <- salida
      }
      gpp.timing <- vector()
      gpp.max <- vector()
      gpp.year <- vector()
      n.cycles <- length(gpp.day.template.f[[st]])
      gpp.day.template <- gpp.day.template.f[[st]]
      for (y in 1:l.years) {
        # y <- 1 # to removed
        windows.size <- 90/n.cycles
        gpp.temp <- all.years[, y]
        
        for (d in 1:n.cycles) {
          if ((gpp.day.template[d] - windows.size) <= 0) {
            back <- abs(gpp.day.template[d] - windows.size)
            new.vector <- c((length(all.years[, y]) - back):length(all.years[,
                                                                             y]), 1:(gpp.day.template[d] + windows.size))
            timing <- new.vector[order(gpp.temp[new.vector], decreasing = T)[1:10]]  #new.vector[which.max(gpp.temp[new.vector])] # bug fixed
            gpp.timing <- c(gpp.timing, timing)
            gp <- gpp.temp[timing]
            gpp.max <- c(gpp.max, gp)
            year.t <- rep(y, times = length(timing))
            gpp.year <- c(gpp.year, year.t)
          } else if ((gpp.day.template[d] + windows.size) > length(all.years[,
                                                                             y])) {
            advance <- (gpp.day.template[d] + windows.size) - length(all.years[,
                                                                               y])
            new.vector <- c((gpp.day.template[d] - windows.size):length(all.years[,
                                                                                  y]), 1:advance)
            timing <- new.vector[order(gpp.temp[new.vector], decreasing = T)[1:10]]  #new.vector[which.max(gpp.temp[new.vector])]
            gpp.timing <- c(gpp.timing, timing)
            gp <- gpp.temp[timing]
            gpp.max <- c(gpp.max, gp)
            year.t <- rep(y, times = length(timing))
            gpp.year <- c(gpp.year, year.t)
          } else {
            start <- gpp.day.template[d] - windows.size
            end <- gpp.day.template[d] + windows.size
            new.vector <- start:end
            timing <- new.vector[order(gpp.temp[new.vector], decreasing = T)[1:10]]  #new.vector[which.max(gpp.temp[new.vector])]
            gpp.timing <- c(gpp.timing, timing)
            gp <- gpp.temp[timing]
            gpp.max <- c(gpp.max, gp)
            year.t <- rep(y, times = length(timing))
            gpp.year <- c(gpp.year, year.t)
          }
        }
      }
      sites.info[st][[1]][[1]] <- gpp.max
      sites.info[st][[1]][[2]] <- gpp.timing
      sites.info[st][[1]][[3]] <- gpp.year
      if (n.cycles == 2) {
        sites.info[st][[1]][[5]]  <- rep(c(1,2), each = 10, times = length(unique(gpp.year)) )
      }
    }
  }
}

## Extracting index of the year --

for (st in 1:nrow(m5)) {
  if (m5$FILTER2[st] == 1) {
    if (is.null(sites.info[[st]]) == F) {
      if (is.element(st, sites.changed)) {
        start.ind <- unique(year(time.serie[st][[1]][[6]]))
      } else {
        start.ind <- unique(year(time.serie[st][[1]][[2]]))
      }
      temp.date <- vector()
      for (i in 1:length(sites.info[st][[1]][[3]])) {
        temp.date <- c(temp.date, as.Date(sites.info[st][[1]][[2]][i], paste((start.ind[sites.info[st][[1]][[3]][i]] -
                                                                                1), 12, 31, sep = "-")))
      }
      sites.info[st][[1]][[4]] <- temp.date
    }
  }
}

# Half-time sensitivity (Deus ex-machina) ----

site <- vector()
length.days <- vector()
rsqrt.2 <- vector()
# deus.ex <- round(exp(seq(log(2), log(365), length.out = 5)))
deus.ex <- section
## deus ex machina ----
for (deus in 1:length(deus.ex)) {
  # deus <- 30
  
  ## Extracting Climatic data ----
  gpp.tair <- array(list(list()))
  for (st in 1:nrow(m5)) {
    # st <- 1
    if (m5$FILTER2[st] == 1) {
      if (is.null(sites.info[[st]]) == FALSE) {
        gpp.tair.data <- time.serie[st][[1]][[3]]
        gpp.nt.vut.ref.time <- names(gpp.tair.data)
        # hours.re <- time.serie[st][[1]][[2]][seq(1, length(time.serie[st][[1]][[2]]),
        # by = 48)]
        n.cycles <- m5$n.cycles[st]
        window.size <- 365
        gpp.tair.tempora.o <- vector()
        gpp.tair.tempora.w <- vector()
        years.site.o <- unique(year(time.serie[st][[1]][[2]]))
        for (i in 1:length(sites.info[st][[1]][[4]])) {
          # i <- 108
          end.t <- which(as.character(gpp.nt.vut.ref.time) == as.character(as.Date(sites.info[st][[1]][[4]][i])))
          if (year(as.Date(sites.info[st][[1]][[4]][i])) == years.site.o[1]) {
            if (end.t < window.size) {
              start.t <- 1
            } else {
              start.t <- end.t - window.size
            }
          } else {
            start.t <- end.t - window.size
          }
          tempora <- gpp.tair.data[start.t:end.t]
          tair.o <- mean(gpp.tair.data[start.t:end.t], na.rm = T)
          tempora <- gpp.tair.data[start.t:end.t]
          gpp.tair.tempora.o <- c(gpp.tair.tempora.o, tair.o)
          v.weight <- half.time(init = 1, t = length(tempora[!is.na(tempora)]),
                                t_half = deus.ex[deus])
          sum.v.weight <- sum(v.weight)
          tair.w <- (sum(tempora[!is.na(tempora)] * v.weight, na.rm = T)/sum.v.weight)
          gpp.tair.tempora.w <- c(gpp.tair.tempora.w, tair.w)
        }
        # original data
        gpp.tair[st][[1]][[1]] <- gpp.tair.tempora.o
        # weighted data
        gpp.tair[st][[1]][[2]] <- gpp.tair.tempora.w
      }
    }
  }
  
  # Precipitation
  
  gpp.precip <- array(list(list()))
  for (st in 1:nrow(m5)) {
    if (m5$FILTER2[st] == 1) {
      if (is.null(sites.info[[st]]) == FALSE) {
        gpp.precip.data <- time.serie[st][[1]][[4]]
        gpp.nt.vut.ref.time <- names(gpp.precip.data)
        # hours.re <- time.serie[st][[1]][[2]][seq(1, length(time.serie[st][[1]][[2]]),
        # by = 48)]
        n.cycles <- m5$n.cycles[st]
        window.size <- 365
        gpp.precip.tempora.o <- vector()
        gpp.precip.tempora.w <- vector()
        years.site.o <- unique(year(time.serie[st][[1]][[2]]))
        for (i in 1:length(sites.info[st][[1]][[4]])) {
          end.t <- which(as.character(gpp.nt.vut.ref.time) == as.character(as.Date(sites.info[st][[1]][[4]][i])))
          if (year(as.Date(sites.info[st][[1]][[4]][i])) == years.site.o[1]) {
            if (end.t < window.size) {
              start.t <- 1
            } else {
              start.t <- end.t - window.size
            }
          } else {
            start.t <- end.t - window.size
          }
          tempora <- gpp.precip.data[start.t:end.t]
          precip.o <- sum(gpp.precip.data[start.t:end.t], na.rm = T)
          tempora <- gpp.precip.data[start.t:end.t]
          gpp.precip.tempora.o <- c(gpp.precip.tempora.o, precip.o)
          # linear weighted v.weight <- seq(from = 0.01 to = 1, length.out = (window.size +
          # 1)) exponential weighted
          v.weight <- half.time(init = 1, t = length(tempora[!is.na(tempora)]),
                                t_half = deus.ex[deus])
          sum.v.weight <- sum(v.weight)
          precip.w <- (sum(tempora[!is.na(tempora)] * v.weight, na.rm = T)/sum.v.weight)
          gpp.precip.tempora.w <- c(gpp.precip.tempora.w, precip.w)
        }
        # original data
        gpp.precip[st][[1]][[1]] <- gpp.precip.tempora.o
        # weighted data
        gpp.precip[st][[1]][[2]] <- gpp.precip.tempora.w
      }
    }
  }
  
  # Shortway incomming radiation
  gpp.swrad <- array(list(list()))
  for (st in 1:nrow(m5)) {
    if (m5$FILTER2[st] == 1) {
      if (is.null(sites.info[[st]]) == FALSE) {
        gpp.swrad.data <- time.serie[st][[1]][[5]]
        gpp.nt.vut.ref.time <- names(gpp.swrad.data)
        # hours.re <- time.serie[st][[1]][[2]][seq(1, length(time.serie[st][[1]][[2]]),
        # by = 48)]
        n.cycles <- m5$n.cycles[st]
        window.size <- 365
        gpp.swrad.tempora.o <- vector()
        gpp.swrad.tempora.w <- vector()
        years.site.o <- unique(year(time.serie[st][[1]][[2]]))
        for (i in 1:length(sites.info[st][[1]][[4]])) {
          end.t <- which(as.character(gpp.nt.vut.ref.time) == as.character(as.Date(sites.info[st][[1]][[4]][i])))
          if (year(as.Date(sites.info[st][[1]][[4]][i])) == years.site.o[1]) {
            if (end.t < window.size) {
              start.t <- 1
            } else {
              start.t <- end.t - window.size
            }
          } else {
            start.t <- end.t - window.size
          }
          tempora <- gpp.swrad.data[start.t:end.t]
          swrad.o <- mean(gpp.swrad.data[start.t:end.t], na.rm = T)
          tempora <- gpp.swrad.data[start.t:end.t]
          gpp.swrad.tempora.o <- c(gpp.swrad.tempora.o, swrad.o)
          # linear weighted v.weight <- seq(from = 0.01 to = 1, length.out = (window.size +
          # 1)) exponential weighted
          v.weight <- half.time(init = 1, t = length(tempora[!is.na(tempora)]),
                                t_half = deus.ex[deus])
          sum.v.weight <- sum(v.weight)
          swrad.w <- (sum(tempora[!is.na(tempora)] * v.weight, na.rm = T)/sum.v.weight)
          gpp.swrad.tempora.w <- c(gpp.swrad.tempora.w, swrad.w)
          
        }
        # original data
        gpp.swrad[st][[1]][[1]] <- gpp.swrad.tempora.o
        # weighted data
        gpp.swrad[st][[1]][[2]] <- gpp.swrad.tempora.w
      }
    }
  }
  
  # VPD
  gpp.vpd <- array(list(list()))
  for (st in 1:nrow(m5)) {
    if (m5$FILTER2[st] == 1) {
      if (is.null(sites.info[[st]]) == FALSE) {
        gpp.vpd.data <- time.serie[st][[1]][[7]]
        gpp.nt.vut.ref.time <- names(gpp.vpd.data)
        # hours.re <- time.serie[st][[1]][[2]][seq(1, length(time.serie[st][[1]][[2]]),
        # by = 48)]
        n.cycles <- m5$n.cycles[st]
        window.size <- 365
        gpp.vpd.tempora.o <- vector()
        gpp.vpd.tempora.w <- vector()
        years.site.o <- unique(year(time.serie[st][[1]][[2]]))
        for (i in 1:length(sites.info[st][[1]][[4]])) {
          end.t <- which(as.character(gpp.nt.vut.ref.time) == as.character(as.Date(sites.info[st][[1]][[4]][i])))
          if (year(as.Date(sites.info[st][[1]][[4]][i])) == years.site.o[1]) {
            if (end.t < window.size) {
              start.t <- 1
            } else {
              start.t <- end.t - window.size
            }
          } else {
            start.t <- end.t - window.size
          }
          tempora <- gpp.vpd.data[start.t:end.t]
          vpd.o <- mean(gpp.vpd.data[start.t:end.t], na.rm = T)
          tempora <- gpp.vpd.data[start.t:end.t]
          gpp.vpd.tempora.o <- c(gpp.vpd.tempora.o, vpd.o)
          # linear weighted v.weight <- seq(from = 0.01 to = 1, length.out = (window.size +
          # 1)) exponential weighted
          v.weight <- half.time(init = 1, t = length(tempora[!is.na(tempora)]),
                                t_half = deus.ex[deus])
          sum.v.weight <- sum(v.weight)
          vpd.w <- (sum(tempora[!is.na(tempora)] * v.weight, na.rm = T)/sum.v.weight)
          gpp.vpd.tempora.w <- c(gpp.vpd.tempora.w, vpd.w)
        }
        # original data
        gpp.vpd[st][[1]][[1]] <- gpp.vpd.tempora.o
        # weighted data
        gpp.vpd[st][[1]][[2]] <- gpp.vpd.tempora.w
      }
    }
  }
  
  # reformuling filter 2
  m5$FILTER2[which(is.na(m5$n.cycles) == T)] <- 0
  
  m5$n.cycles[which(is.na(m5$n.cycles) == T)] <- 0
  
  
  ## Removing LT sites (temporal) ---- problems with the variable radiation ----
  m5$FILTER2[which(m5$DATASET == "LT")] <- 0
  # two growing seassons
  m5$FILTER2[which(m5$n.cycles == 2)] <- 0
  # temporal filter for sites with other problems m5$FILTER2[16] <- 0 # without
  # observations in vpd m5$FILTER2[49] <- 0 # without observations in precip.
  # m5$FILTER2[65] <- 0 # without observations in precip. m5$FILTER2[45] <- 0 # too
  # many NAs in vpd only 20 values after removed the NAs Multiple linear
  # regressions and sd ----
  options(warn = 2)
  trew <- array(list(list()))
  # st <- 3
  for (st in 1:nrow(m5)) {
    if (m5$FILTER2[st] == 1) {
      if (m5$FILTER2[st] == 1 & m5$n.cycles[st] == 1) {
        if (is.null(sites.info[[st]]) == FALSE) {
          if (is.element(st, sites.changed)) {
            print(paste("deus=", deus.ex[deus], "site:", st, m5$SITE_ID[st],
                        sep = " "))
            gpp.nt.vut.ref <- time.serie[st][[1]][[1]]
            gpp.nt.vut.ref.time <- year(time.serie[st][[1]][[6]])
            hours.re <- time.serie[st][[1]][[6]]
          } else {
            print(paste("deus=", deus.ex[deus], "site:", st, m5$SITE_ID[st],
                        sep = " "))
            gpp.nt.vut.ref <- time.serie[st][[1]][[1]]
            gpp.nt.vut.ref.time <- year(time.serie[st][[1]][[2]])
            hours.re <- time.serie[st][[1]][[2]]
          }
          # print(st)
          timing <- sites.info[st][[1]][[2]]
          timing.circ <- circular(timing * (360/365) * (pi/180), units = "radians",
                                  modulo = "2pi", rotation = "clock", zero = pi/2)
          
          climatic <- matrix(NA, nrow = length(sites.info[st][[1]][[1]]),
                             ncol = 4)
          climatic[, 1] <- scal(gpp.tair[st][[1]][[2]])
          climatic[, 2] <- scal(gpp.swrad[st][[1]][[2]])
          climatic[, 3] <- scal(gpp.precip[st][[1]][[2]])
          climatic[, 4] <- scal(gpp.vpd[st][[1]][[2]])
          
          if (any(is.na(climatic))) {
            timing <- sites.info[st][[1]][[2]][-which(is.na(climatic), arr.ind = T)[,1]]
            climatic <- na.omit(climatic)
          }
          # pca
          
          climatic_pca <- prcomp(climatic[,c(1,2,4)], center = T, scale. = T)
          
          new.climatic <- matrix(NA, nrow = length(sites.info[st][[1]][[1]]),
                                 ncol = 2)
          
          new.climatic[,1] <- scale(climatic_pca$x[,1], center = T, scale = T)
          new.climatic[,2] <- climatic[,3]
          
          
          bot <- 1000
          # number of columns is equal to the number of variables in the climatic matrix + 1
          temporal.bot <- matrix(NA, ncol = 3, nrow = bot)
          for (bot.i in 1:bot) {
            # bot.i <- 1
            index <- tapply(1:length(timing), rep(1:length(unique(year(hours.re))),
                                                  each = 10), sample, size = 1)
            timing.b <- timing[index]
            climatic.b <- new.climatic[index, ]
            timing.circ.b <- circular(timing.b * (360/365) * (pi/180), units = "radians",
                                      modulo = "2pi", rotation = "counter", zero = 0)
            
            circres <- withTimeout(try(lm.circular(y = timing.circ.b, x = climatic.b,
                                                   init = c(0, 0), type = "c-l", verbose = F), silent = T), timeout = 10,
                                   onTimeout = "warning")
            
            if ((class(circres) == "try-error") || is.null(circres)) {
              temporal.bot[bot.i, 1] <- NA
              temporal.bot[bot.i, 2:3] <- NA
            } else {
              temporal.bot[bot.i, 1] <- circres$mu
              temporal.bot[bot.i, 2:3] <- circres$coef
            }
          }
          trew[st][[1]][[1]] <- temporal.bot  # raw data
          trew[st][[1]][[2]] <- apply(temporal.bot[, 2:3], MARGIN = 2, mean,
                                      na.rm = T)  # averages (coefficients)
          trew[st][[1]][[3]] <- apply(temporal.bot[, 2:3], MARGIN = 2, FUN = sd,
                                      na.rm = T)  # sd  coefficients)
          trew[st][[1]][[4]] <- mean.circular(circular(temporal.bot[, 1],
                                                       units = "radians", modulo = "2pi", rotation = "counter", zero = 0),
                                              na.rm = T)  # mu mean
          #trew[st][[1]][[5]] <- sd.circular(circular(temporal.bot[, 1], units = "radians",
          #                                           modulo = "2pi", rotation = "counter", zero = 0), na.rm = T)  # mu sd
          #trew[st][[1]][[6]] <- mean(temporal.bot[, 6], na.rm = T)  # kappa average
          #trew[st][[1]][[7]] <- sd(temporal.bot[, 6], na.rm = T)  # kappa sd
          trew[st][[1]][[8]] <- sd(timing, na.rm = T)  # deviation of the number of days (linear)
        }
      }
    }
    # testing angular direction ----
    # fiting ----- timing
    rsqrt <- vector()
    name <- vector()
    for (st in 1:nrow(m5)) {
      if (m5$FILTER2[st] == 1 & m5$n.cycles[st] == 1) {
        if (is.null(unlist(trew[st])) == F) {
          if (all(is.nan(trew[st][[1]][[2]])) == FALSE) {
            timing <- sites.info[st][[1]][[2]]
            timing.circ <- circular((timing * (360/365) * (pi/180)), units = "radians",
                                    modulo = "2pi", rotation = "counter", zero = 0)
            climatic <- matrix(NA, nrow = length(sites.info[st][[1]][[1]]),
                               ncol = 4)
            climatic[, 1] <- scal(gpp.tair[st][[1]][[2]])
            climatic[, 2] <- scal(gpp.swrad[st][[1]][[2]])
            climatic[, 3] <- scal(gpp.precip[st][[1]][[2]])
            climatic[, 4] <- scal(gpp.vpd[st][[1]][[2]])
            
            # pca
            climatic_pca <- prcomp(climatic[,c(1,2,4)], center = T, scale. = T)
            
            new.climatic <- matrix(NA, nrow = length(sites.info[st][[1]][[1]]),
                                   ncol = 2)
            
            new.climatic[,1] <- scale(climatic_pca$x[,1], center = T, scale = T)
            new.climatic[,2] <- climatic[,3]
            
            
            sim <- vector()
            for (i in 1:nrow(climatic)) {
              # i <- 1
              x <- new.climatic[i, ]
              sim <- c(sim, as.numeric(circ.recons(x, mu = as.numeric(trew[st][[1]][[4]]),
                                                   coef = trew[st][[1]][[2]])))
            }
            # sim <- sim / (360/365) / (pi/180)
            # test <- circular(sim, units = "radians",modulo = "2pi", rotation = "clock", zero = pi/2)
            test <- circular(sim, units = "radians", modulo = "2pi", zero = 0, rotation = "counter")
            #timing.circ <- conversion.circular(timing.circ, units = "radians",
            #                        modulo = "2pi", rotation = "counter", zero = 0)
            ex <- cor.circular(timing.circ, test, test = T)$cor
            # plot(as.numeric(timing.circ), points(as.numeric(sim), col = 'red'))
            rsqrt <- c(rsqrt, cor.circular(timing.circ, test, test = T)$cor)
            name <- c(name, st)
          } else {
            rsqrt <- c(rsqrt, NA)
            name <- c(name, st)
            print(paste("deus:", deus.ex[deus], st, m5$SITE_ID[st], "NA",
                        sep = " "))
          }
        }
      }
    }
  }
  length.days <- c(length.days, rep(deus.ex[deus], length(rsqrt)))
  site <- c(site, name)
  rsqrt.2 <- c(rsqrt.2, rsqrt)
  print(paste("deus:", deus.ex[deus], "finished", sep = " "))
}


test.data.frame <- data.frame(sites = m5$SITE_ID[site], half.time = length.days,
                              jama = rsqrt.2)
print(head(test.data.frame))
write.csv(test.data.frame, paste("results/sensitivity_half_time_", c.job, ".csv", sep = ""))
