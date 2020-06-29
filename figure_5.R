#library(BiocManager)
#library("Rgraphviz")
influence <- read.csv("~/phd/data/prob_comparison/influence.csv")
intercomp <- read.csv("~/phd/data/prob_comparison/intercomp.cvs")
sub.gs.1 <- read.csv("~/phd/data/prob_comparison/sub.gs.1.origina.csv")



# influence of the climatic variables ----

influence
# PC1
length(which(influence[,1] == 1))
# # radiation
# length(which(influence[,2] == 1))
# precip
length(which(influence[,2] == 1))
# vpd
# length(which(influence[,4] == 1))
# summary barplot ----


library(RColorBrewer)
library(ggplot2)
colores <- brewer.pal(3, "BrBG")



# tair.sites <- length(which(intercomp$tair != 0))
# swin.sites <- length(which(intercomp$swrad != 0))
# precip.sites <- length(which(intercomp$precip != 0))
# vpd.sites <- length(which(intercomp$vpd != 0))

PC1.sites <- length(which(intercomp$PC1 != 0))
precip.sites <- length(which(intercomp$precip != 0))

input <- data.frame(climate = c("PC1", "Precip"), value = c(PC1.sites,precip.sites))
# input$climate <- factor(input$climate, levels = c(c("SWin", "Tair", "VPD", "Precip")))

ggplot(data = input, aes(x = climate, y = value, fill = climate)) +
  geom_bar(stat = "identity") +
  xlab(expression("DOY"["GPPmax"]*" driver")) +
  ylab("No. of sites \n (where climate variable is significant)") +
  scale_fill_manual(values = c(colores[3], colores[3], colores[3], colores[3])) +
  theme_bw() +
  theme(legend.position = "none")

ggsave("~/phd/doc/papers/doy_gppmax/images/summary_bar_plot.png", width = 6.5, height = 7)
ggsave("~/phd/doc/papers/doy_gppmax/images/summary_bar_plot.pdf", width = 6.5, height = 7)


igbp <- unique(sub.gs.1$IGBP)

colores <- brewer.pal(3, "BrBG")

values <- vector()

for (i in 1:length(igbp)) {
  # tair.sites <- length(which(intercomp$tair != 0 & sub.gs.1$IGBP == igbp[i]))
  # swin.sites <- length(which(intercomp$swrad != 0 & sub.gs.1$IGBP == igbp[i]))
  # precip.sites <- length(which(intercomp$precip != 0 & sub.gs.1$IGBP == igbp[i]))
  # vpd.sites <- length(which(intercomp$vpd != 0 & sub.gs.1$IGBP == igbp[i]))
  
  PC1.sites <- length(which(intercomp$PC1 != 0 & sub.gs.1$IGBP == igbp[i]))
  precip.sites <- length(which(intercomp$precip != 0 & sub.gs.1$IGBP == igbp[i]))
  values <- c(values, PC1.sites, precip.sites)
}

climate <- rep(c("PC1", "Precip"), time = length(igbp))

igbp <- rep(igbp, each = 2)
number_ticks <- function(n) {function(limits) pretty(limits, n)}

input <- data.frame(igbp, climate, values)

out <- ggplot(data = input, aes(x = igbp, y = values, fill = climate)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  theme_bw() +
  xlab("Vegetation type") +
  ylab("Cumulative number of sites") +
  labs(fill = expression("DOY"["GPPmax"] * "\n driver")) +
  scale_fill_discrete(labels = c("Tair, SWin, VPD", "Precip")) +
  scale_y_continuous(breaks=number_ticks(10))

out
ggsave("figure_5.png", plot = out, width = 6.5, height = 7)
ggsave("figure_5.pdf", plot = out, width = 6.5, height = 7)


input2 <- input

for (i in 1:nrow(input2)) {
  input2$values[i] <- input$values[i] / sum(input$values[which(input$igbp == input$igbp[i])])
}

ggplot(data = input2, aes(x = igbp, y = values, fill = climate)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  xlab("Vegetation type") +
  ylab("Cumulative contribution") +
  labs(fill = "DOYmax \n driver")



# summary statistics (sense) ----

head(intercomp)

length(which(intercomp$PC2 > 0))

length(which(intercomp$PC1 < 0))

length(which(intercomp$precip < 0))


