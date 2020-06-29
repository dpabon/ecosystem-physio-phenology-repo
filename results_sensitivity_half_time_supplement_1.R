library(ggplot2)
library(plotly)
library(gtools)
library(cowplot)
library(ggpubr)
setwd("~/phd/")
m5 <- read.csv("m5.csv")
files <- list.files(path = "data/half_time_sensitivity_365/", full.names = T, pattern = "*.csv")
files <- mixedsort(sort(files))
test.data.frame <- read.csv(files[1], row.names = 1)
for (i in 2:length(files)) {
  test.data.frame.2 <- read.csv(files[i], row.names = 1)
  test.data.frame <- rbind(test.data.frame, test.data.frame.2)
}




#View(test.data.frame)
sites <- unique(test.data.frame$sites)

p <- ggplot(data = test.data.frame, aes(x = half.time, y = jama, group = sites, colour = sites)) +
  geom_line() +
  xlab("Half-time") +
  ylab("Correlation Coefficient (JS)")

ggplotly(p)

ggsave(plot = p, filename = "doc/papers/doy_gppmax/images/append3_halftime.pdf",  width = 11, height = 6.5)
ggsave(plot = p, filename = "doc/papers/doy_gppmax/images/append3_halftime.png",  width = 11, height = 6.5)
# Summary per site ----


# y <- match(m5$SITE_ID, unique(test.data.frame$sites))
# 
# m5$SITE_ID[which(is.na(y))]
# 
# for (i in 1:nrow(m5)) {
#   if (m5$FILTER2[i] == 1 & m5$n.cycles[i] == 1) {
#     subsection <- test.data.frame[which(test.data.frame$sites == as.character(m5$SITE_ID[i])),]
#     m5[i, "half.time"] <- subsection[which.max(subsection[,"jama"]), "half.time"]
#     m5[i, "j.s"] <- max(subsection[,"jama"])
#   }
# 
# }
# write.csv(m5, "~/phd/data/fluxnet/m5.csv", row.names = F)
View(m5)



# Summary per vegetation type and climate clase ----

maximum <- vector()
igbp.2 <- vector()
koppen.2 <- vector()


igbp <- match(sites, m5$SITE_ID)
koppen.2 <- as.factor(m5$kop[igbp])
igbp.2 <- as.factor(m5$IGBP[igbp])
lat <- m5$LAT[igbp]
jama.2 <- vector()

plot(test.data.frame[which(test.data.frame$sites == "US-Ha1"), "jama"], type = "n")
for (i in 1:length(sites)) {
  
  abline(v = which.max(test.data.frame[which(test.data.frame$sites == sites[i]), "jama"]), lwd = 2)
  subsection <- test.data.frame[which(test.data.frame$sites == sites[i]),]
  maximum[i] <- subsection[which.max(subsection[,"jama"]), "half.time"]
  jama.2[i] <- max(subsection[,"jama"])
}

## How many sites have a half-time value greather than 100 

m5[which(m5$half.time > 100),]

## How many sites have a half-time value less than 10

m5[which(m5$half.time < 10),]

hist(m5$half.time[which(m5$IGBP == "EBF")])

## Distribution of the maximum JS correlation coefficient

hist(m5$j.s)

max(m5$j.s, na.rm = T)
quantile(m5$j.s, probs = 0.5, na.rm = T)

# Which sites with a JS less than 0.8

m5[which(m5$j.s < 0.8),]

plot(density(maximum))
hist(maximum)

tyu <- data.frame(maximum, jama.2, koppen.2, igbp.2, lat)

m5.temp <- m5[igbp,]

m5.temp$maximum <- maximum

m5.temp$jama <- jama.2

m5.temp[which(m5.temp$jama < 0.4),]

m5.temp[which(m5.temp$maximum > 200),]

half.time.js.hex <- ggplot(tyu, aes(maximum, jama.2, colour = lat)) +
  geom_hex() + theme_bw() + theme(axis.text = element_text(size = 14),
                                    axis.title = element_text(size = 14,face = "bold")) +
  xlab("Half-time") + ylab("Correlation Coefficient (JS)") +
  scale_fill_viridis_c()

# ggsave(plot = half.time.js.hex, filename = "doc/papers/doy_gppmax/images/half_time_js_hex.pdf",  width = 9, height = 6.5)
# ggsave(plot = half.time.js.hex, filename = "doc/papers/doy_gppmax/images/half_time_js_hex.png",  width = 9, height = 6.5)

half.time.js <- ggplot(tyu, aes(jama.2)) +
  geom_histogram() + theme_bw() + theme(axis.text = element_text(size = 14),
                           axis.title = element_text(size = 14,face = "bold"), plot.title = element_text(size = 15, face = "bold")) +
  xlab("Correlation Coefficient (JS)") + ylab("Number of sites") +
  ggtitle("a)")

ggsave(plot = half.time.js, filename = "doc/papers/doy_gppmax/images/append3_halftime_a.pdf",  width = 11, height = 6.5)
ggsave(plot = half.time.js, filename = "doc/papers/doy_gppmax/images/append3_halftime_a.png",  width = 11, height = 6.5)

second.data <- data.frame(maxi = maximum, kop = koppen.2, igbp.2 = igbp.2)

half.time.days <- ggplot(data = second.data, aes(maxi)) +
  geom_histogram() + theme_bw() + theme(axis.text = element_text(size = 14),
                                    axis.title = element_text(size = 14,face = "bold"), plot.title = element_text(size = 15, face = "bold")) +
  xlab("Half-time") + ylab("Number of sites") +
  ggtitle("b)")

comp.half <- ggarrange(half.time.js, half.time.days, ncol = 1, nrow = 2)

ggsave(plot = comp.half, filename = "doc/papers/doy_gppmax/images/figure6_half_time.pdf", width = 11, height = 6.5)
ggsave(plot = comp.half, filename = "doc/papers/doy_gppmax/images/figure6_half_time.png", width = 11, height = 6.5)

ggsave(plot = half.time.days, filename = "doc/papers/doy_gppmax/images/append3_halftime_b.pdf",  width = 11, height = 6.5)
ggsave(plot = half.time.days, filename = "doc/papers/doy_gppmax/images/append3_halftime_b.png",  width = 11, height = 6.5)



#View(second.data)
#second.data <- second.data[is.element(second.data$igbp.2, c("DBF", "EBF", "ENF", "GRA", "MF")),]

  half.time.kop <- ggplot(second.data, aes(x = kop, y = maximum)) +
  geom_violin() +
  geom_jitter(shape = 16, position=position_jitter(0.2)) +
  geom_boxplot(width = 0.1) + theme_bw() + theme(axis.text = element_text(size = 14),
                                                 axis.title = element_text(size = 14,face = "bold")) +
  xlab("Koppengeiger classes") + ylab("Half-time")

ggsave(plot = half.time.kop, filename = "doc/papers/doy_gppmax/images/append3_halftime_c.pdf",  width = 11, height = 6.5)
ggsave(plot = half.time.kop, filename = "doc/papers/doy_gppmax/images/append3_halftime_c.png",  width = 11, height = 6.5)


half.time.vt <- ggplot(second.data, aes(x = igbp.2, y = maxi)) +
  geom_violin() +
  geom_jitter(shape = 16, position = position_jitter(0.2)) +
  geom_boxplot(width = 0.1) + theme_bw() + theme(axis.text = element_text(size = 14),
                                                 axis.title = element_text(size = 14,face = "bold")) +
  xlab("Vegetation type") + ylab("Half-time")

ggsave(plot = half.time.vt, filename = "doc/papers/doy_gppmax/images/append3_halftime_d.pdf",  width = 11, height = 6.5)
ggsave(plot = half.time.vt, filename = "doc/papers/doy_gppmax/images/append3_halftime_d.png",  width = 11, height = 6.5)


subm5 <- m5[which(m5$FILTER2 == 1 & m5$n.cycles == 1),]

write.csv(subm5, "doc/papers/doy_gppmax/images/append3_halftime.csv",row.names = F)

## Two growing seasons ----


m5 <- read.csv("~/phd/data/fluxnet/m5.csv")
files <- list.files(path = "~/phd/data/fluxnet/half_time_sensitivity_365/2_seasons_180/1_season/", full.names = T)
files <- mixedsort(sort(files))
test.data.frame.1.1 <- read.csv(files[1], row.names = 1)

for (i in 2:length(files)) {
  test.data.frame.2 <- read.csv(files[i], row.names = 1)
  test.data.frame.1.1 <- rbind(test.data.frame.1.1, test.data.frame.2)
}



test.data.frame.1.1[which.max(test.data.frame.1.1$jama),]
#View(test.data.frame)
sites <- unique(test.data.frame.1.1$sites)

#test.data.frame <- test.data.frame[which(test.data.frame$length.days <= 200),]
p.1.1 <- ggplot(data = test.data.frame.1.1, aes(x = half.time, y = jama, group = sites, colour = sites)) +
  geom_line() +
  theme(plot.title = element_text(hjust = 0, size = 17, face = "plain")) +
  xlab("Half-time") +
  ylab("Correlation Coefficient (JS)") +
  ggtitle("a)")


p.1.1


files <- list.files(path = "~/phd/data/fluxnet/half_time_sensitivity_365/2_seasons_180/2_season/", full.names = T)
files <- mixedsort(sort(files))
test.data.frame.1.2 <- read.csv(files[1], row.names = 1)
for (i in 2:length(files)) {
  test.data.frame.2 <- read.csv(files[i], row.names = 1)
  test.data.frame.1.2 <- rbind(test.data.frame.1.2, test.data.frame.2)
}

test.data.frame.1.2[which.max(test.data.frame.1.2$jama),]

## plots for only one sites ----
sites <- unique(test.data.frame.1.2$sites)

test.data.frame.1.1$GS <- "GS1"
test.data.frame.1.2$GS <- "GS2"

test.data.frame.1.site <- rbind(test.data.frame.1.1, test.data.frame.1.2)


half.time.1.site <- ggplot(data = test.data.frame.1.site, aes(x = half.time, y = jama, colour = GS)) +
  geom_line() +
  theme(plot.title = element_text(hjust = 0, size = 17, face = "plain")) +
  xlab("Half-time") +
  ylab("Correlation Coefficient (JS)") +
  scale_colour_discrete(name  ="",
                        breaks=c("GS1", "GS2"),
                        labels=c("First growing season", "Second growing season"))

ggsave(plot = half.time.1.site, filename = "doc/papers/doy_gppmax/images/append4.pdf",  width = 11, height = 6.5)
ggsave(plot = half.time.1.site, filename = "doc/papers/doy_gppmax/images/append4.png",  width = 11, height = 6.5)


## plots for two sites ----
#test.data.frame <- test.data.frame[which(test.data.frame$length.days <= 200),]
p.1.2 <- ggplot(data = test.data.frame.1.2, aes(x = half.time, y = jama, group = sites, colour = sites)) +
  geom_line() +
  theme(plot.title = element_text(hjust = 0, size = 17, face = "plain")) +
  xlab("Half-time") +
  ylab("Correlation Coefficient (JS)") +
  ggtitle("b)")


g.s.2 <- ggarrange(p.1.1, p.1.2, ncol = 1, nrow = 2, common.legend = T, legend = "bottom")

# ggsave(plot = g.s.2, filename = "doc/papers/gppmax/images/append4.pdf",  width = 11, height = 6.5)
# ggsave(plot = g.s.2, filename = "doc/papers/gppmax/images/append4.png",  width = 11, height = 6.5)

maximum <- vector()
igbp.2 <- vector()
koppen.2 <- vector()
sites <- unique(test.data.frame.1.1$sites)

igbp <- match(sites, m5$SITE_ID)
koppen.2 <- as.factor(m5$kop[igbp])
igbp.2 <- as.factor(m5$IGBP[igbp])
jama.2 <- vector()

plot(test.data.frame.1.1[which(test.data.frame.1.1$sites == "FR-Pue"), "jama"], type = "n")
for (i in 1:length(sites)) {
  
  abline(v = which.max(test.data.frame.1.1[which(test.data.frame.1.1$sites == sites[i]), "jama"]), lwd = 2)
  subsection <- test.data.frame.1.1[which(test.data.frame.1.1$sites == sites[i]),]
  maximum[i] <- subsection[which.max(subsection[,"jama"]), "half.time"]
  jama.2[i] <- max(subsection[,"jama"])
}

tyu.1 <- data.frame(maximum, jama.2)

sites
####

maximum <- vector()
igbp.2 <- vector()
koppen.2 <- vector()
sites <- unique(test.data.frame.1.2$sites)

igbp <- match(sites, m5$SITE_ID)
koppen.2 <- as.factor(m5$kop[igbp])
igbp.2 <- as.factor(m5$IGBP[igbp])
jama.2 <- vector()

plot(test.data.frame.1.2[which(test.data.frame.1.2$sites == "FR-Pue"), "jama"], type = "n")
for (i in 1:length(sites)) {
  
  abline(v = which.max(test.data.frame.1.2[which(test.data.frame.1.2$sites == sites[i]), "jama"]), lwd = 2)
  subsection <- test.data.frame.1.2[which(test.data.frame.1.2$sites == sites[i]),]
  maximum[i] <- subsection[which.max(subsection[,"jama"]), "half.time"]
  jama.2[i] <- max(subsection[,"jama"])
}

tyu.2 <- data.frame(maximum, jama.2)

# Summary

summary.g.2 <- data.frame(Site.name = sites, half.time = c(tyu.1$maximum, tyu.2$maximum), js = c(tyu.1$jama.2, tyu.2$jama.2), n.cycle = c(1,2))

write.csv(summary.g.2, "results/gppmax/sensitivity_half_time_g2.csv", row.names = F)
