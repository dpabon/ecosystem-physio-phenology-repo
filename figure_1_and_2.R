library(circular)
library(lubridate)

# REFERENCE DISTRIBUTION
set.seed(1234)
dist <- rvonmises(n = 20, mu = 0, kappa = 20)

mu <- 0
b1 <- -0.3
x1 <- 0.6
b2 <- 0.3
x2 <- 0.6
y <- circular( mu + (2 * atan(sum(c(b1*x1), (b2*x2)))))
pdf("raw_conceptual.pdf")
plot.circular(y, zero = pi/2, rotation = "clock", shrink = 2, type = "n", cex = 1.3, tcl.text = 0.16)
points.circular(y, zero = pi/2, cex = 2, shrink = 1.5, rotation = "clock")
lines(density.circular(dist, bw = 20), rotation = "clock", zero = pi/2)
#text(1.4,-1, labels = expression(paste(2.94, "=", pi, - 2*atan("0.4 x -0.2 + 0.6 x 0.3"))) )

# increasing x1
x1 <- 1.3
x2 <- 0.6
y <- circular( mu + (2 * atan(sum(b1*x1 + b2*x2))))
points.circular(y, zero = pi/2, rotation = "clock", cex = 2, bg = "#5ab4ac", pch = 21)
#text(-1.6,-0.7, labels = "/mu = /pi- 2 * atan(5.8 * -0.2 + 0.6 * 0.3)")
# increasing x2
x1 <- 0.6
x2 <- 1.3
y <- circular( mu + (2 * atan(sum(b1*x1 + b2*x2))))
points(y, zero = pi/2, rotation = "clock", cex = 2, bg = "#d8b365", pch = 21)
dev.off()
#text(1.7,-0.6, labels = "/mu = /pi - 2 * atan(0.4 * -0.2 + 4.0 * 0.3)")

#legend("toplef", legend = c("reference", "increase of x1", "increase of x2"), col = c("black", "red", "blue"), pch = 19,bty = "n")

## Figure 1

mt <- month(seq(as.Date("2005/1/1"), as.Date("2005/12/31"), "days"))

mt <- match(1:12, mt)



pdf("conceptual_figure1.pdf")

plot.circular(y, zero = pi/2, rotation = "clock", type = "n", cex = 0.8, shrink = 2.35)

#axis.circular(at = circular(seq(1,365, by = 75) * (360/365) * (pi/180), units = "radians", modulo = "2pi", rotation = "clock",  zero = pi/2), units =  "radians", labels = seq(1,365, by = 75), zero = pi/2, rotation = "clock",cex = 1, tick = F, tcl.text = 0.46)

axis.circular(at = circular(mt * (360/365) * (pi/180), units = "radians", modulo = "2pi", rotation = "clock",  zero = pi/2), units =  "radians", labels = month.abb, zero = pi/2, rotation = "clock", cex = 0.8, tick = F, tcl.text = 0.3)

set.seed(45454)
northern <- rvonmises(n = 300, mu = 2.9, kappa = 10)
lines(density.circular(northern, bw = 3, adjust = 10 ), zero = pi/2, rotation = "clock", col = "#67a9cf", lwd = 3)


southern <- rvonmises(n = 300, mu = 6.2, kappa = 12)
#plot.circular(southern, stack = T)
lines(density.circular(southern, bw = 2, adjust = 12 ), zero = pi/2, rotation = "clock", col = "#ef8a62", lwd = 3)


dev.off()
