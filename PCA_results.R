library(ggplto2)
library(reshape2)


climatic.pca.contrib <- matrix(data = NA, nrow = 1, ncol = 4)
climatic.pca.contrib <- as.data.frame(climatic.pca.contrib)
colnames(climatic.pca.contrib) <- c("PC1", "PC2", "PC3", "site_name")

pc1.climatic.contrib <- matrix(data = NA, nrow = 1, ncol = 4)
pc1.climatic.contrib <- as.data.frame(pc1.climatic.contrib)
colnames(pc1.climatic.contrib) <- c("Tair", "SWin", "VPD", "site_name")

# st <- 65
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
      climatic[,1]  <- scal(gpp.tair[st][[1]][[2]])
      climatic[,2]  <- scal(gpp.swrad[st][[1]][[2]])
      climatic[,3]  <- scal(gpp.precip[st][[1]][[2]])
      climatic[,4] <- scal(gpp.vpd[st][[1]][[2]])

      if (any(is.na(climatic))) {
        timing <- sites.info[st][[1]][[2]][-which(is.na(climatic), arr.ind = T)[,1]]
        climatic <- na.omit(climatic)
      }

      # pca

      climatic_pca <- prcomp(climatic[,c(1,2,4)], center = T, scale. = T)
      
      temp.pc1 <- data.frame(Tair = climatic_pca$rotation[1,1], SWin = climatic_pca$rotation[2,1], VPD = climatic_pca$rotation[3,1], site_name = as.character(m5$SITE_ID[st]))
      
      pc1.climatic.contrib <- rbind(pc1.climatic.contrib, temp.pc1)
      
      temp.pcas <- summary(climatic_pca)
      
      temp.pcas <- data.frame(PC1 = temp.pcas$importance[2,1], PC2 = temp.pcas$importance[2,2], PC3 = temp.pcas$importance[2,3], site_name = m5$SITE_ID[st])
      
      climatic.pca.contrib <- rbind(climatic.pca.contrib, temp.pcas)

      # bootstraping
      # bot <- 1000
      # temporal.bot <- matrix(NA, ncol = 6, nrow = bot)
      # temporal.bot.l <- matrix(NA, ncol = 5, nrow = bot)
      # for (bot.i in 1:bot) {
      #   index <- tapply(1:length(timing), rep(1:length(unique(year(hours.re))), each = 10), sample, size = 1)
      #   timing.b <- timing[index]
      #   climatic.b <- new.climatic[index,]
      #   timing.circ.b  <- circular(timing.b * (360/365) * (pi/180),units = "radians", modulo = "2pi", rotation = "counter", zero = 0)
      #   circres <- withTimeout(try(lm.circular(y = timing.circ.b, x = climatic.b,
      #                                          init = c(0, 0), type = "c-l", verbose = F), silent = T), timeout = 10,
      #                          onTimeout = "warning")
      # 
      #   if ((class(circres) == "try-error") || is.null(circres)) {
      #     temporal.bot[bot.i, 1] <- NA
      #     temporal.bot[bot.i, 2:3] <- NA
      #     temporal.bot[bot.i, 4] <- NA
      #     temporal.bot[bot.i, 5:6] <- NA
      #   } else {
      #     temporal.bot[bot.i, 1] <- circres$mu
      #     temporal.bot[bot.i, 2:3] <- circres$coef
      #     temporal.bot[bot.i, 4] <- circres$kappa
      #     temporal.bot[bot.i, 5:6] <- circres$p.values
      #   }
      # 
      #   circres.l <- lm(as.numeric(timing.circ.b) ~ climatic.b)
      # 
      #   temporal.bot.l[bot.i, 1] <- circres.l$coefficients[1]
      #   temporal.bot.l[bot.i, 2:3] <- circres.l$coefficients[2:3]
      #   temporal.bot.l[bot.i, 4:5] <-  summary(circres.l)$coefficients[-1,4]
      # 
      # }
      # 
      # trew[st][[1]][[1]] <- temporal.bot # raw data
      # 
      # trew[st][[1]][[2]] <-  apply(temporal.bot[,2:3], MARGIN = 2, mean, na.rm = T) # averages (coefficients)
      # trew[st][[1]][[3]] <- apply(temporal.bot[,2:3], MARGIN = 2, FUN = sd, na.rm = T) # sd  coefficients)
      # trew[st][[1]][[4]] <- mean.circular(circular(temporal.bot[,1],units = "radians", modulo = "2pi", rotation = "counter", zero = 0), na.rm = T)  # mu mean
      # 
      # trew[st][[1]][[5]] <- sd.circular(circular(temporal.bot[,1],units = "radians", modulo = "2pi", rotation = "counter", zero = 0), na.rm = T) # mu sd
      # trew[st][[1]][[6]] <-  mean(temporal.bot[,4], na.rm = T) # kappa average
      # trew[st][[1]][[7]] <- sd(temporal.bot[,4], na.rm = T) # kappa sd
      # trew[st][[1]][[8]] <- sd(timing, na.rm = T) # deviation of the number of days (linear)
      # 
      # trew[st][[1]][[9]] <- apply(temporal.bot[,5:6], MARGIN = 2, median, na.rm = T) # averages (p.values)
      # # standard error
      # trew[st][[1]][[10]] <- trew[st][[1]][[3]] / length(na.omit(temporal.bot[,2:3]))
      # 
      # ## Multiple linear regression information
      # 
      # trew[st][[1]][[11]] <- temporal.bot.l
      # # average of coefficients
      # trew[st][[1]][[12]] <- apply(temporal.bot.l[,2:3], MARGIN = 2, mean, na.rm = T)
      # # sd of coefficients
      # trew[st][[1]][[13]] <- apply(temporal.bot.l[,2:3], MARGIN = 2, FUN = sd, na.rm = T)
      # # mu mean
      # trew[st][[1]][[14]] <- mean(temporal.bot[,1], na.rm = T)
      # # mu sd
      # trew[st][[1]][[15]] <- sd(temporal.bot[,1], na.rm = T)
      # # median p.values
      # trew[st][[1]][[16]] <- apply(temporal.bot[,5:6], MARGIN = 2, median, na.rm = T)
      # # standard error
      # trew[st][[1]][[16]] <- trew[st][[1]][[3]] / length(na.omit(temporal.bot.l[,2:3]))
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
        climatic[,1]  <- scal(gpp.tair[st][[1]][[2]][which(sites.info[st][[1]][[5]] == d)])
        climatic[,2]  <- scal(gpp.swrad[st][[1]][[2]][which(sites.info[st][[1]][[5]] == d)])
        climatic[,3]  <- scal(gpp.precip[st][[1]][[2]][which(sites.info[st][[1]][[5]] == d)])
        climatic[,4] <- scal(gpp.vpd[st][[1]][[2]][which(sites.info[st][[1]][[5]] == d)])
        if (any(is.na(climatic))) {
          timing <- sites.info[st][[1]][[2]][-which(is.na(climatic), arr.ind = T)[,1]]
          climatic <- na.omit(climatic)
        }

        # pca

        climatic_pca <- prcomp(climatic[,c(1,2,4)], center = T, scale. = T)
        
        climatic_pca <- prcomp(climatic[,c(1,2,4)], center = T, scale. = T)
        
        temp.pc1 <- data.frame(Tair = climatic_pca$rotation[1,1], SWin = climatic_pca$rotation[2,1], VPD = climatic_pca$rotation[3,1], site_name = as.character(m5$SITE_ID[st]))
        
        pc1.climatic.contrib <- rbind(pc1.climatic.contrib, temp.pc1)
        
        temp.pcas <- summary(climatic_pca)
        
        temp.pcas <- data.frame(PC1 = temp.pcas$importance[2,1], PC2 = temp.pcas$importance[2,2], PC3 = temp.pcas$importance[2,3], site_name = m5$SITE_ID[st])
        
        climatic.pca.contrib <- rbind(climatic.pca.contrib, temp.pcas)

        # bootstraping
        # bot <- 1000
        # temporal.bot <- matrix(NA, ncol = 6, nrow = bot)
        # temporal.bot.l <- matrix(NA, ncol = 5, nrow = bot)
        # for (bot.i in 1:bot) {
        #   index <- tapply(1:length(timing), rep(1:length(unique(year(hours.re))), each = 10), sample, size = 1)
        #   timing.b <- timing[index]
        #   climatic.b <- new.climatic[index,]
        #   timing.circ.b  <- circular(timing.b * (360/365) * (pi/180),units = "radians", modulo = "2pi", rotation = "counter", zero = 0)
        #   circres <- withTimeout(try(lm.circular(y = timing.circ.b, x = climatic.b,
        #                                          init = c(0, 0), type = "c-l", verbose = F), silent = T), timeout = 10,
        #                          onTimeout = "warning")
        # 
        #   if ((class(circres) == "try-error") || is.null(circres)) {
        #     temporal.bot[bot.i, 1] <- NA
        #     temporal.bot[bot.i, 2:5] <- NA
        #     temporal.bot[bot.i, 6] <- NA
        #     temporal.bot[bot.i, 7:10] <- NA
        #   } else {
        #     temporal.bot[bot.i, 1] <- circres$mu
        #     temporal.bot[bot.i, 2:3] <- circres$coef
        #     temporal.bot[bot.i, 4] <- circres$kappa
        #     temporal.bot[bot.i, 5:6] <- as.numeric(circres$p.values)
        #   }
        # 
        #   circres.l <- lm(as.numeric(timing.circ.b) ~ climatic.b)
        # 
        #   temporal.bot.l[bot.i, 1] <- circres.l$coefficients[1]
        #   temporal.bot.l[bot.i, 2:3] <- circres.l$coefficients[2:3]
        #   temporal.bot.l[bot.i, 4:5] <-  summary(circres.l)$coefficients[-1,4]
        # }
        # if (d == 1) {
        #   # raw data
        #   trew.1.1[st][[1]][[1]] <- temporal.bot
        #   # averages (coefficients)
        #   trew.1.1[st][[1]][[2]] <-  apply(temporal.bot[,2:3], MARGIN = 2, mean, na.rm = T)
        #   # sd  coefficients)
        #   trew.1.1[st][[1]][[3]] <- apply(temporal.bot[,2:3], MARGIN = 2, FUN = sd, na.rm = T)
        #   # mu mean
        #   trew.1.1[st][[1]][[4]] <- mean.circular(circular(temporal.bot[,1],units = "radians", modulo = "2pi", rotation = "counter", zero = 0), na.rm = T)
        #   # mu sd
        #   trew.1.1[st][[1]][[5]] <- sd.circular(circular(temporal.bot[,1],units = "radians", modulo = "2pi", rotation = "counter", zero = 0), na.rm = T)
        #   # kappa average
        #   trew.1.1[st][[1]][[6]] <-  mean(temporal.bot[,4], na.rm = T)
        #   # kappa sd
        #   trew.1.1[st][[1]][[7]] <- sd(temporal.bot[,4], na.rm = T)
        #   # deviation of the number of days (linear)
        #   trew.1.1[st][[1]][[8]] <- sd(timing, na.rm = T)
        #   # averages (p.values)
        #   trew.1.1[st][[1]][[9]] <- apply(temporal.bot[,5:6], MARGIN = 2, median, na.rm = T)
        #   trew.1.1[st][[1]][[10]] <- trew.1.1[st][[1]][[3]] / length(na.omit(temporal.bot[,2:3]))
        # 
        #   ## Multiple linear regression information
        # 
        #   trew.1.1[st][[1]][[11]] <- temporal.bot.l
        #   # average of coefficients
        #   trew.1.1[st][[1]][[12]] <- apply(temporal.bot.l[,2:3], MARGIN = 2, mean, na.rm = T)
        #   # sd of coefficients
        #   trew.1.1[st][[1]][[13]] <- apply(temporal.bot.l[,2:3], MARGIN = 2, FUN = sd, na.rm = T)
        #   # mu mean
        #   trew.1.1[st][[1]][[14]] <- mean(temporal.bot[,1], na.rm = T)
        #   # mu sd
        #   trew.1.1[st][[1]][[15]] <- sd(temporal.bot[,1], na.rm = T)
        #   # median p.values
        #   trew.1.1[st][[1]][[16]] <- apply(temporal.bot[,5:6], MARGIN = 2, median, na.rm = T)
        #   # standard error
        #   trew.1.1[st][[1]][[16]] <- trew[st][[1]][[3]] / length(na.omit(temporal.bot.l[,2:3]))
        # 
        # }else{
        #   # raw data
        #   trew.1.2[st][[1]][[1]] <- temporal.bot
        #   # averages (coefficients)
        #   trew.1.2[st][[1]][[2]] <-  apply(temporal.bot[,2:3], MARGIN = 2, mean, na.rm = T)
        #   # sd  coefficients)
        #   trew.1.2[st][[1]][[3]] <- apply(temporal.bot[,2:3], MARGIN = 2, FUN = sd, na.rm = T)
        #   # mu mean
        #   trew.1.2[st][[1]][[4]] <- mean.circular(circular(temporal.bot[,1],units = "radians", modulo = "2pi", rotation = "counter", zero = 0), na.rm = T)
        #   # mu sd
        #   trew.1.2[st][[1]][[5]] <- sd.circular(circular(temporal.bot[,1],units = "radians", modulo = "2pi", rotation = "counter", zero = 0), na.rm = T)
        #   # kappa average
        #   trew.1.2[st][[1]][[6]] <-  mean(temporal.bot[,4], na.rm = T)
        #   # kappa sd
        #   trew.1.2[st][[1]][[7]] <- sd(temporal.bot[,4], na.rm = T)
        #   # deviation of the number of days (linear)
        #   trew.1.2[st][[1]][[8]] <- sd(timing, na.rm = T)
        #   # averages (p.values)
        #   trew.1.2[st][[1]][[9]] <- apply(temporal.bot[,5:6], MARGIN = 2, median, na.rm = T)
        #   trew.1.2[st][[1]][[10]] <- trew.1.2[st][[1]][[3]] / length(na.omit(temporal.bot[,2:3]))
        # 
        #   ## Multiple linear regression information
        # 
        #   trew.1.2[st][[1]][[11]] <- temporal.bot.l
        #   # average of coefficients
        #   trew.1.2[st][[1]][[12]] <- apply(temporal.bot.l[,2:5], MARGIN = 2, mean, na.rm = T)
        #   # sd of coefficients
        #   trew.1.2[st][[1]][[13]] <- apply(temporal.bot.l[,2:5], MARGIN = 2, FUN = sd, na.rm = T)
        #   # mu mean
        #   trew.1.2[st][[1]][[14]] <- mean(temporal.bot[,1], na.rm = T)
        #   # mu sd
        #   trew.1.2[st][[1]][[15]] <- sd(temporal.bot[,1], na.rm = T)
        #   # average p.values
        #   trew.1.2[st][[1]][[16]] <- apply(temporal.bot[,5:6], MARGIN = 2, median, na.rm = T)
        #   # standard error
        #   trew.1.2[st][[1]][[16]] <- trew[st][[1]][[3]] / length(na.omit(temporal.bot.l[,2:3]))
        # }
      }
    }
  }
}

climatic.pca.contrib <- climatic.pca.contrib[-1,]
pc1.climatic.contrib <- pc1.climatic.contrib[-1,]

# write.csv(climatic.pca.contrib, file = "results/gppmax/climatic_pca_contribution.csv", row.names = F)
# 
# write.csv(pc1.climatic.contrib, file = "results/gppmax/pc1_variables_contribution.csv", row.names = F)

climatic.pca.contrib.plot <- melt(climatic.pca.contrib[,1:3])

ggplot(data = climatic.pca.contrib.plot) +
  geom_boxplot(aes(x = variable, y = value)) +
  xlab("Principal Component") +
  ylab("Proportion of variance explained") + 
  theme_bw(base_size = 18)

# ggsave(filename = "doc/papers/doy_gppmax/images/summary_pcas_contribution.png", width = 12, height = 8)

for (i in 1:nrow(pc1.climatic.contrib)) {
  if (pc1.climatic.contrib[i,1] < 1 & pc1.climatic.contrib[i,2] < 1 & pc1.climatic.contrib[i,3] < 1) {
    pc1.climatic.contrib[i,1] <- abs(pc1.climatic.contrib[i,1])
    pc1.climatic.contrib[i,2] <- abs(pc1.climatic.contrib[i,2])
    pc1.climatic.contrib[i,3] <- abs(pc1.climatic.contrib[i,3])
  }
}

pc1.climatic.contrib[,1:3]


pc1.climatic.contrib.plot <- melt(pc1.climatic.contrib[,1:3])

ggplot(pc1.climatic.contrib.plot) +
  geom_boxplot(aes(x = variable, y = value)) +
  xlab("Climatic variable") +
  ylab("Contribution to PC1") + 
  theme_bw(base_size = 18)

# ggsave(filename = "doc/papers/doy_gppmax/images/summary_pc1_contribution.png", width = 12, height = 8)


