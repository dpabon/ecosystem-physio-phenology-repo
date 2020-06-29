library(GGally)

climatic.allsites <- matrix(data = NA,nrow = 1, ncol = 5)
climatic.allsites <- as.data.frame(climatic.allsites)
colnames(climatic.allsites) <- c("Tair", "SWin", "Precip", "VPD", "site_name")

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
      print(st)
      timing  <- sites.info[st][[1]][[2]]
      timing.circ  <- circular(timing * (360/365) * (pi/180),
                               units = "radians", modulo = "2pi",
                               rotation = "counter", zero = 0 )
      
      climatic  <- matrix(NA, nrow = length(sites.info[st][[1]][[1]]), ncol = 4)
      climatic[,1]  <- gpp.tair[st][[1]][[2]]
      climatic[,2]  <- gpp.swrad[st][[1]][[2]]
      climatic[,3]  <- gpp.precip[st][[1]][[2]]
      climatic[,4] <- gpp.vpd[st][[1]][[2]]
      if (any(is.na(climatic))) {
        timing <- sites.info[st][[1]][[2]][-which(is.na(climatic), arr.ind = T)[,1]]
        climatic <- na.omit(climatic)
      }
      
      climatic <- as.data.frame(climatic)
      colnames(climatic) <- c("Tair", "SWin", "Precip", "VPD")
      
      climatic$site_name <- rep(m5$SITE_ID[st], nrow(climatic))
      
      climatic.allsites <- rbind(climatic.allsites, climatic)

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
        print(st)
        timing  <- sites.info[st][[1]][[2]][which(sites.info[st][[1]][[5]] == d)]
        timing.circ  <- circular(timing * (360/365) * (pi/180),
                                 units = "radians", modulo = "2pi",
                                 rotation = "counter", zero = 0)
        
        climatic  <- matrix(NA, nrow = length(sites.info[st][[1]][[1]][which(sites.info[st][[1]][[5]] == d)]), ncol = 4)
        climatic[,1]  <- gpp.tair[st][[1]][[2]][which(sites.info[st][[1]][[5]] == d)]
        climatic[,2]  <- gpp.swrad[st][[1]][[2]][which(sites.info[st][[1]][[5]] == d)]
        climatic[,3]  <- gpp.precip[st][[1]][[2]][which(sites.info[st][[1]][[5]] == d)]
        climatic[,4] <- gpp.vpd[st][[1]][[2]][which(sites.info[st][[1]][[5]] == d)]
      
        if (any(is.na(climatic))) {
          timing <- sites.info[st][[1]][[2]][-which(is.na(climatic), arr.ind = T)[,1]]
          climatic <- na.omit(climatic)
        }
        
        climatic <- as.data.frame(climatic)
        colnames(climatic) <- c("Tair", "SWin", "Precip", "VPD")
        climatic$site_name <- rep(m5$SITE_ID[st], nrow(climatic))
        
        
        climatic.allsites <- rbind(climatic.allsites, climatic)
      }
      
    }
  }
}

nrow(climatic.allsites)
climatic.allsites <- climatic.allsites[-1,]
climatic.allsites[,1] <- scale(climatic.allsites[,1], center = T, scale = T)
climatic.allsites[,2] <- scale(climatic.allsites[,2], center = T, scale = T)
climatic.allsites[,3] <- scale(climatic.allsites[,3], center = T, scale = T)
climatic.allsites[,4] <- scale(climatic.allsites[,4], center = T, scale = T)

plot(climatic.allsites[,1:4])

climatic_pca <- prcomp(climatic.allsites[,1:4], center = T, scale. = T)

climatic_pca

new.climatic.space <- scale(climatic_pca$x[,1:2], center = T, scale = T)

new.climatic.space <- as.data.frame(new.climatic.space)

new.climatic.space

new.climatic.space$site_names <- climatic.allsites$site_name


saveRDS(new.climatic.space, "data/fluxnet/global_climatic_space.RDS")
