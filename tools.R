# extract specific years from an array and a vectors days.
year.extract <- function(x, date, extr.years, time = 3){
  # time dimension should be in third position
  if (is.Date(date)) {
    require(lubridate)
    years <- year(date)
    years <- as.numeric(years)
  }else{
    years <- date
  }
  posit <- match(years, extr.years)
  y <- x[,,which(is.na(posit) == FALSE)]
  return(y)
  
}
# Mean seasonal cycle ---------
msc <- function(x, date, u.years, na.rm = T){
  # x = could be a vector or an multidimmensional array
  # u.years = vector with specific years
  # a priori all the years have 365 days
  if (is.Date(date)) {
    require(lubridate)
    years <- year(date)
    years <- as.numeric(years)
  }else{
    years <- date
  }
  if (missing(u.years) == F) {
    posit <- match(years, u.years)
    x <- x[,,which(is.na(posit) == FALSE)]
  }else{
    x <- x
  }
  lenght.year <- length(which(years == years[1]))
  
  if (is.array(x) == T) {
    d.msc <- array(data = NA, dim = c(dim(x)[1], dim(x)[2], lenght.year))
    for (i in 1:dim(x)[2]) {
      for (j in 1:dim(x)[1]) {
        for (d in 1:lenght.year) {
          d.msc[j,i,d] <- mean(x[j,i, c(seq(d, dim(x)[3], by = lenght.year))], na.rm = na.rm)
          
        }
      }
    }
    return(d.msc)
  }else{
    for (d in 1:lenght.year) {
      d.msc <- rep(NA, lenght.year)
      d.msc[d] <- mean(x[seq(d, length(x), by = lenght.year)])
    }
  }
  
}

# Number of cycles per year -------------

n.cycles <- function(x, max = 3){
  max <- max + 1
  n.cycles  <- array(data = NA, dim = c(dim(x)[1], dim(x)[2]))
  for (i in 1:dim(x)[2]) {
    for (j in 1:dim(x)[1]) {
      if (anyNA(x[j,i,])) {
        
      }else {
        n.cycles[j,i]  <- as.integer(which.max(abs(fft(x[j,i,]))[2:max]))
      }
      }
  }
  return(n.cycles)
}


# day template -------------

template <- function(msc, cycles, date){
  # msc = Mean seasonal cycle object
  # cycles = number of cycle object
  if (is.Date(date)) {
    require(lubridate)
    years <- year(date)
    years <- as.numeric(years)
  }else{
    years <- date
  }
  year.l <- length(which(years == years[1]))
  template <- array(data = list(), dim = c(dim(msc)[1], dim(msc)[2]))
  
  for (i in 1:dim(msc)[2]) {
    for (j in 1:dim(msc)[1]) {
      if (anyNA(msc[j,i,])) {
        
      }else {
        fy <- fft(msc[j,i,])
        ny <- length(fy)
        fy[c(1,5:(ny - 3))] <- 0
        yfiltered <- fft(fy,inverse = T)
        yfiltered  <- Re(yfiltered)
        count  <- 1
        for (d in 1:year.l) {
          
          if (d == 1) {
            if (yfiltered[d] > yfiltered[year.l] & yfiltered[d] > yfiltered[(d + 1)]) {
              template[j,i][[1]][count]  <- d
              count  <- count + 1
            }
          } else if (d == year.l) {
            if (yfiltered[d] > yfiltered[(d - 1)] & yfiltered[d] > yfiltered[1]) {
              template[j,i][[1]][count]  <- d
              count  <- count + 1
            }
          } else {
            if (yfiltered[d] > yfiltered[(d - 1)] & yfiltered[d] > yfiltered[(d + 1)]) {
              template[j,i][[1]][count]  <- d
              count  <- count + 1
            }
          }
          
        }
        
        # checking that the number of peaks is equal to the number of cycles (only selecting the highest picks)
        if (length(template[j,i][[1]]) == cycles[j,i]) {
          #
        }else {
          tempora  <- yfiltered[template[j,i][[1]]]
          template[j,i][[1]]  <- template[j,i][[1]][order(tempora, decreasing = T)[1:cycles[j,i]]]
        } 
      }
    }
  }
  return(template)
}

# timing extraction -------------

timing <- function(original, template, n.cycles, date){
  # original data
  # template object
  # mean seasonal cycle object
  # date = vector of classe date or numberic
  # [position][[1]][[1]][[1]]
  # position - 1 -  cycle - year
  if (is.Date(date)) {
    require(lubridate)
    years <- year(date)
    years <- as.numeric(years)
  }else{
    years <- date
  }
  gpp.timing  <- array(data = list(list()), dim =  c(dim(original)[1],dim(original)[2], length(unique(years))))
  # DEFININTION: base number of days giving the number of measures per year
  base.window <- round(length(which(years == years[1])) / 4) - 1
  # if it's an array
  for (y in 1:length(unique(years))) {
    # current year
    temp.year <- unique(years)[y]
    
    for (i in 1:dim(original)[2]) {
      for (j in 1:dim(original)[1]) {
        if (all(is.na(n.cycles[j,i]))) {
          gpp.timing[j,i,y][[1]] <- NA
        }else{
          
          windows.size <- round(base.window / n.cycles[j,i])
          gpp.temp <- original[j,i,which(years == temp.year)]
          date.temp <- date[which(years == temp.year)]
          
          for (d in 1:n.cycles[j,i]) {
            if ((template[j,i][[1]][d] - windows.size) <= 0 ) {
              back <- abs(template[j,i][[1]][d] - windows.size)
              new.vector <- c((length(gpp.temp) - back):length(gpp.temp),
              1:(template[j,i][[1]][d] + windows.size))
              index.day <- new.vector[which.max(gpp.temp[new.vector])]
              gpp.timing[j,i,y][[1]][[d]] <- yday(date.temp[index.day])
              
            }else if ((template[j,i][[1]][d] + windows.size) > length(gpp.temp)) {
              advance <- (template[j,i][[1]][d] + windows.size) - length(gpp.temp)
              new.vector <- c((template[j,i][[1]][d] - windows.size):length(gpp.temp),
                              1:advance)
              index.day <- new.vector[which.max(gpp.temp[new.vector])]
              gpp.timing[j,i,y][[1]][[d]] <- yday(date.temp[index.day])
              
            }else {
              start <- template[j,i][[1]][d] - windows.size
              end <- template[j,i][[1]][d] + windows.size
              new.vector <- start:end
              index.day <- new.vector[which.max(gpp.temp[new.vector])]
              gpp.timing[j,i,y][[1]][[d]] <- yday(date.temp[index.day])
            }
          }
        }
      }
    }
  }
  return(gpp.timing)
}



# maximum extraction -----

maximum <- function(original, template, n.cycles, date){
  # original data
  # template object
  # mean seasonal cycle object
  # date = vector of classe date or numberic
  # [position][[1]][[1]][[1]]
  # position - 1 -  cycle - year
  
  if (is.Date(date)) {
    require(lubridate)
    years <- year(date)
    years <- as.numeric(years)
  }else{
    years <- date
  }
  gpp.max  <- array(data = list(list()), dim =  c(720,360, length(unique(years))))
  
  # DEFININTION: base number of days giving the number of measures per year
  base.window <- round(length(which(years == years[1])) / 4) - 1
  # if it's an array
  for (y in 1:length(unique(years))) {
    # current year
    temp.year <- unique(years)[y]
    
    for (i in 1:dim(original)[2]) {
      for (j in 1:dim(original)[1]) {
        if (all(is.na(n.cycles[j,i]))) {
          gpp.max[j,i,y][[1]] <- NA
        }else{
          
          windows.size <- round(base.window / n.cycles[j,i])
          gpp.temp <- original[j,i,which(years == temp.year)]
          
          for (d in 1:n.cycles[j,i]) {
            if ((template[j,i][[1]][d] - windows.size) <= 0 ) {
              back <- abs(template[j,i][[1]][d] - windows.size)
              new.vector <- c((length(gpp.temp) - back):length(gpp.temp),
                              1:(template[j,i][[1]][d] + windows.size))
              gpp.timing <- new.vector[which.max(gpp.temp[new.vector])]
              gpp.max[j,i,y][[1]][[d]] <- gpp.temp[gpp.timing]
              
            }else if ((template[j,i][[1]][d] + windows.size) > length(gpp.temp)) {
              advance <- (template[j,i][[1]][d] + windows.size) - length(gpp.temp)
              new.vector <- c((template[j,i][[1]][d] - windows.size):length(gpp.temp),
                              1:advance)
              gpp.timing <- new.vector[which.max(gpp.temp[new.vector])]
              gpp.max[j,i,y][[1]][[d]] <- gpp.temp[gpp.timing]
              
            }else {
              start <- template[j,i][[1]][d] - windows.size
              end <- template[j,i][[1]][d] + windows.size
              new.vector <- start:end
              gpp.timing <- new.vector[which.max(gpp.temp[new.vector])]
              gpp.max[j,i,y][[1]][[d]] <- gpp.temp[gpp.timing]
            }
          }
        }
      }
    }
  }
  return(gpp.max)
}




