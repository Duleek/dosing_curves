library(drc)

data <- read.csv("Drug dosing - cyclophosphamide.csv", head = TRUE)
head(data)

x <- log10(data$ï..Drug.Concentration)
y1 <- c(data$BT20_Day1)
y2 <- c(data$CAL51_Day1)
y3 <- c(data$HCC1937_Day1)
y4 <- c(data$HDQP1_Day1)

plot(y1 ~ x, main = "Cyclophosphamide day 1 exposure efficacy with CellTiter-Glo",
	ylim = c(0,150), col = "blue",
	xlab="Dosing log10(nM)", ylab="Cell Viability (%)", pch=19)

points(x, y2, col="red")
points(x, y3, col="green") 
points(x, y4, col="orange")

fit1 <- nls(y1 ~ SSlogis(x, Asym, xmid, scal), data = data.frame(x, y1))
summary(fit1)

fit2 <- nls(y2 ~ SSlogis(x, Asym, xmid, scal), data = data.frame(x, y2))
summary(fit2)

fit3 <- nls(y3 ~ SSlogis(x, Asym, xmid, scal), data = data.frame(x, y3))
summary(fit3)

fit4 <- nls(y4 ~ SSlogis(x, Asym, xmid, scal), data = data.frame(x, y4))
summary(fit4)

#if 4-paramater curve does not work due to too few data points for proper fit
#fit3 <- drm(y3 ~ x, data = data.frame(x, y3), fct = L.4())


fit1 <- drm(y1 ~ x, data = data.frame(x, y1), fct = L.4())
fit2 <- drm(y2 ~ x, data = data.frame(x, y2), fct = L.4())
fit3 <- drm(y3 ~ x, data = data.frame(x, y3), fct = L.4())
fit4 <- drm(y4 ~ x, data = data.frame(x, y4), fct = L.4())

lines(seq(0, 7, length.out = 100), col = "blue", 
	predict(fit1, newdata = data.frame(x = seq(0, 7, length.out = 100))))

lines(seq(0, 7, length.out = 100), col = "red", 
	predict(fit2, newdata = data.frame(x = seq(0, 7, length.out = 100))))

lines(seq(0, 7, length.out = 100), col = "green", 
	predict(fit3, newdata = data.frame(x = seq(0, 7, length.out = 100))))

lines(seq(0, 7, length.out = 100), col = "orange", 
	predict(fit4, newdata = data.frame(x = seq(0, 7, length.out = 100))))


legend(4.5, 150, legend=c("BT20 viability", "CAL51 viability", "HCC1937 viability", "HDQP1 viability"),
       col=c("blue", "red", "green", "orange"), lty=1:2, cex=0.8)
