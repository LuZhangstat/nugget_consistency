#setwd("/Users/luzhang/Documents/github/nugget_consistency")
rm(list = ls())
load("./data/simdata.RData")
source("./projects/utils.R")
# need to compare the MLE of parameter sigma^2 tau^2 phi

# First obtain the MLE of phi and delta = tau^2 / sigma^2 through gradient decent
# of the profile likelihood, then calculate the sigma^2
library("optimx")
# tt <- optimx(c(phi0.5_1, tau22), neg_loglik_0.5_profile,
#              X = Sim_data_0.5["0.2", "0.5_0.15", 1, 1:n1], D = D[1:n1, 1:n1],
#              n = n1, method = c('L-BFGS-B'), lower = c(0, 0))
# # 0.912 seconds
# 
# tt <- optimx(c(phi0.5_1, tau22), neg_loglik_0.5_profile,
#              X = Sim_data_0.5["0.2", "0.5_0.15", 1, 1:n2], D = D[1:n2, 1:n2],
#              n = n2, method = c('L-BFGS-B'), lower = c(0, 0))
# # 7.784 seconds
# 
# tt <- optimx(c(phi0.5_1, tau22), neg_loglik_0.5_profile,
#              X = Sim_data_0.5["0.2", "0.5_0.15", 1, 1:n3], D = D[1:n3, 1:n3],
#              n = n3, method = c('L-BFGS-B'), lower = c(0, 0))
# 54.01 seconds

load("./data/sim2_0.5_0_0.4results.RData")
load("./data/sim2_0.5_0.2_0.4results.RData")
load("./data/sim2_0.5_0.8_0.4results.RData")

library(ggplot2)

#-------------- check plots of each setting -----------------------#
par(mfrow = c(1, 3))
hist(MLE_0.5_0_0.4_n1[, "sigma2"], xlim = c(0.2, 2.3))
hist(MLE_0.5_0_0.4_n2[, "sigma2"], xlim = c(0.2, 2.3))
hist(MLE_0.5_0_0.4_n3[, "sigma2"], xlim = c(0.2, 2.3))
par(mfrow = c(1, 3))
hist(MLE_0.5_0_0.4_n1[, "phi"], xlim = c(1, 20))
hist(MLE_0.5_0_0.4_n2[, "phi"], xlim = c(1, 20))
hist(MLE_0.5_0_0.4_n3[, "phi"], xlim = c(1, 20))
par(mfrow = c(1, 3))
hist(MLE_0.5_0_0.4_n1[, "phi"] * MLE_0.5_0_0.4_n1[, "sigma2"], 
     xlim = c(1, 20))
hist(MLE_0.5_0_0.4_n2[, "phi"] * MLE_0.5_0_0.4_n2[, "sigma2"], 
     xlim = c(1, 20))
hist(MLE_0.5_0_0.4_n3[, "phi"] * MLE_0.5_0_0.4_n3[, "sigma2"], 
     xlim = c(1, 20))




par(mfrow = c(1, 3))
hist(MLE_0.5_0.2_0.4_n1[, "tau2"], xlim = c(0, 0.4))
hist(MLE_0.5_0.2_0.4_n2[, "tau2"], xlim = c(0, 0.4))
hist(MLE_0.5_0.2_0.4_n3[, "tau2"], xlim = c(0, 0.4))
par(mfrow = c(1, 3))
hist(MLE_0.5_0.2_0.4_n1[, "sigma2"], xlim = c(0.2, 2.3))
hist(MLE_0.5_0.2_0.4_n2[, "sigma2"], xlim = c(0.2, 2.3))
hist(MLE_0.5_0.2_0.4_n3[, "sigma2"], xlim = c(0.2, 2.3))
par(mfrow = c(1, 3))
hist(MLE_0.5_0.2_0.4_n1[, "phi"], xlim = c(1, 20))
hist(MLE_0.5_0.2_0.4_n2[, "phi"], xlim = c(1, 20))
hist(MLE_0.5_0.2_0.4_n3[, "phi"], xlim = c(1, 20))
par(mfrow = c(1, 3))
hist(MLE_0.5_0.2_0.4_n1[, "phi"] * MLE_0.5_0.2_0.4_n1[, "sigma2"], 
     xlim = c(1, 20))
hist(MLE_0.5_0.2_0.4_n2[, "phi"] * MLE_0.5_0.2_0.4_n2[, "sigma2"], 
     xlim = c(1, 20))
hist(MLE_0.5_0.2_0.4_n3[, "phi"] * MLE_0.5_0.2_0.4_n3[, "sigma2"], 
     xlim = c(1, 20))


par(mfrow = c(1, 3))
hist(MLE_0.5_0.8_0.4_n1[, "tau2"], xlim = c(0.4, 1.2))
hist(MLE_0.5_0.8_0.4_n2[, "tau2"], xlim = c(0.4, 1.2))
hist(MLE_0.5_0.8_0.4_n3[, "tau2"], xlim = c(0.4, 1.2))
par(mfrow = c(1, 3))
hist(MLE_0.5_0.8_0.4_n1[, "sigma2"], xlim = c(0.2, 2.5))
hist(MLE_0.5_0.8_0.4_n2[, "sigma2"], xlim = c(0.2, 2.5))
hist(MLE_0.5_0.8_0.4_n3[, "sigma2"], xlim = c(0.2, 2.5))
par(mfrow = c(1, 3))
hist(MLE_0.5_0.8_0.4_n1[, "phi"], xlim = c(1, 21))
hist(MLE_0.5_0.8_0.4_n2[, "phi"], xlim = c(1, 21))
hist(MLE_0.5_0.8_0.4_n3[, "phi"], xlim = c(1, 21))
par(mfrow = c(1, 3))
hist(MLE_0.5_0.8_0.4_n1[, "phi"] * MLE_0.5_0.8_0.4_n1[, "sigma2"], 
     xlim = c(1, 20))
hist(MLE_0.5_0.8_0.4_n2[, "phi"] * MLE_0.5_0.8_0.4_n2[, "sigma2"], 
     xlim = c(1, 20))
hist(MLE_0.5_0.8_0.4_n3[, "phi"] * MLE_0.5_0.8_0.4_n3[, "sigma2"], 
     xlim = c(1, 20))


#-------------- check plots of phi, sigma2 and phi*sigma2 -------------------#
width <- 1000
height <- 1000
pointsize <- 16

png(paste("./pics/hist_tau2.png", sep = ""), 
    width = width, height = height * 2 / 3, 
    pointsize = pointsize, family = "Courier")
par(mfrow = c(2, 3))
hist(MLE_0.5_0.2_0.4_n1[, "tau2"], xlim = c(0, 0.4), breaks = 15,
     main = expression(tau^2 == 0.2 ~ "," ~ n == 400), xlab = expression(tau^2))
hist(MLE_0.5_0.2_0.4_n2[, "tau2"], xlim = c(0, 0.4), breaks = 15,
     main = expression(tau^2 == 0.2 ~ "," ~ n == 900), xlab = expression(tau^2))
hist(MLE_0.5_0.2_0.4_n3[, "tau2"], xlim = c(0, 0.4), breaks = 15,
     main = expression(tau^2 == 0.2 ~ "," ~ n == 1600), xlab = expression(tau^2))
hist(MLE_0.5_0.8_0.4_n1[, "tau2"], xlim = c(0.4, 1.2), breaks = 15,
     main = expression(tau^2 == 0.8 ~ "," ~ n == 400), xlab = expression(tau^2))
hist(MLE_0.5_0.8_0.4_n2[, "tau2"], xlim = c(0.4, 1.2), breaks = 15,
     main = expression(tau^2 == 0.8 ~ "," ~ n == 900), xlab = expression(tau^2))
hist(MLE_0.5_0.8_0.4_n3[, "tau2"], xlim = c(0.4, 1.2), breaks = 15,
     main = expression(tau^2 == 0.8 ~ "," ~ n == 1600), xlab = expression(tau^2))
dev.off()


png(paste("./pics/hist_sigma2.png", sep = ""), 
    width = width, height = height, pointsize = pointsize, family = "Courier")
par(mfrow = c(3, 3))
hist(MLE_0.5_0_0.4_n1[, "sigma2"], xlim = c(0.2, 2.5), breaks = 15,
     main = expression(tau^2 == 0 ~ "," ~ n == 400), xlab = expression(sigma^2))
hist(MLE_0.5_0_0.4_n2[, "sigma2"], xlim = c(0.2, 2.5), breaks = 15,
     main = expression(tau^2 == 0 ~ "," ~ n == 900), xlab = expression(sigma^2))
hist(MLE_0.5_0_0.4_n3[, "sigma2"], xlim = c(0.2, 2.5), breaks = 15,
     main = expression(tau^2 == 0 ~ "," ~ n == 1600), xlab = expression(sigma^2))
hist(MLE_0.5_0.2_0.4_n1[, "sigma2"], xlim = c(0.2, 2.5), breaks = 15,
     main = expression(tau^2 == 0.2 ~ "," ~ n == 400), xlab = expression(sigma^2))
hist(MLE_0.5_0.2_0.4_n2[, "sigma2"], xlim = c(0.2, 2.5), breaks = 15,
     main = expression(tau^2 == 0.2 ~ "," ~ n == 900), xlab = expression(sigma^2))
hist(MLE_0.5_0.2_0.4_n3[, "sigma2"], xlim = c(0.2, 2.5), breaks = 15,
     main = expression(tau^2 == 0.2 ~ "," ~ n == 1600), xlab = expression(sigma^2))
hist(MLE_0.5_0.8_0.4_n1[, "sigma2"], xlim = c(0.2, 2.5), breaks = 15,
     main = expression(tau^2 == 0.8 ~ "," ~ n == 400), xlab = expression(sigma^2))
hist(MLE_0.5_0.8_0.4_n2[, "sigma2"], xlim = c(0.2, 2.5), breaks = 15,
     main = expression(tau^2 == 0.8 ~ "," ~ n == 900), xlab = expression(sigma^2))
hist(MLE_0.5_0.8_0.4_n3[, "sigma2"], xlim = c(0.2, 2.5), breaks = 15,
     main = expression(tau^2 == 0.8 ~ "," ~ n == 1600), xlab = expression(sigma^2))
dev.off()

png(paste("./pics/hist_phi.png", sep = ""), 
    width = width, height = height, pointsize = pointsize, family = "Courier")
par(mfrow = c(3, 3))
hist(MLE_0.5_0_0.4_n1[, "phi"], xlim = c(1, 21), breaks = 15,
     main = expression(tau^2 == 0 ~ "," ~ n == 400), xlab = expression(phi))
hist(MLE_0.5_0_0.4_n2[, "phi"], xlim = c(1, 21), breaks = 15,
     main = expression(tau^2 == 0 ~ "," ~ n == 900), xlab = expression(phi))
hist(MLE_0.5_0_0.4_n3[, "phi"], xlim = c(1, 21), breaks = 15,
     main = expression(tau^2 == 0 ~ "," ~ n == 1600), xlab = expression(phi))
hist(MLE_0.5_0.2_0.4_n1[, "phi"], xlim = c(1, 21), breaks = 15,
     main = expression(tau^2 == 0.2 ~ "," ~ n == 400), xlab = expression(phi))
hist(MLE_0.5_0.2_0.4_n2[, "phi"], xlim = c(1, 21), breaks = 15,
     main = expression(tau^2 == 0.2 ~ "," ~ n == 900), xlab = expression(phi))
hist(MLE_0.5_0.2_0.4_n3[, "phi"], xlim = c(1, 21), breaks = 15,
     main = expression(tau^2 == 0.2 ~ "," ~ n == 1600), xlab = expression(phi))
hist(MLE_0.5_0.8_0.4_n1[, "phi"], xlim = c(1, 21), breaks = 15,
     main = expression(tau^2 == 0.8 ~ "," ~ n == 400), xlab = expression(phi))
hist(MLE_0.5_0.8_0.4_n2[, "phi"], xlim = c(1, 21), breaks = 15,
     main = expression(tau^2 == 0.8 ~ "," ~ n == 900), xlab = expression(phi))
hist(MLE_0.5_0.8_0.4_n3[, "phi"], xlim = c(1, 21), breaks = 15,
     main = expression(tau^2 == 0.8 ~ "," ~ n == 1600), xlab = expression(phi))
dev.off()

png(paste("./pics/hist_phisigma2.png", sep = ""), 
    width = width, height = height, pointsize = pointsize, family = "Courier")
par(mfrow = c(3, 3))
hist(MLE_0.5_0_0.4_n1[, "phi"] * MLE_0.5_0_0.4_n1[, "sigma2"], 
     xlim = c(1, 20), breaks = 15 , xlab = expression(sigma^2 * phi),
     main = expression(tau^2 == 0 ~ "," ~ n == 400))
hist(MLE_0.5_0_0.4_n2[, "phi"] * MLE_0.5_0_0.4_n2[, "sigma2"], 
     xlim = c(1, 20), breaks = 15, xlab = expression(sigma^2 * phi),
     main = expression(tau^2 == 0 ~ "," ~ n == 900))
hist(MLE_0.5_0_0.4_n3[, "phi"] * MLE_0.5_0_0.4_n3[, "sigma2"], 
     xlim = c(1, 20), breaks = 15, xlab = expression(sigma^2 * phi),
     main = expression(tau^2 == 0 ~ "," ~ n == 1600))
hist(MLE_0.5_0.2_0.4_n1[, "phi"] * MLE_0.5_0.2_0.4_n1[, "sigma2"], 
     xlim = c(1, 20), breaks = 15, xlab = expression(sigma^2 * phi),
     main = expression(tau^2 == 0.2 ~ "," ~ n == 400))
hist(MLE_0.5_0.2_0.4_n2[, "phi"] * MLE_0.5_0.2_0.4_n2[, "sigma2"], 
     xlim = c(1, 20), breaks = 15, xlab = expression(sigma^2 * phi),
     main = expression(tau^2 == 0.2 ~ "," ~ n == 900))
hist(MLE_0.5_0.2_0.4_n3[, "phi"] * MLE_0.5_0.2_0.4_n3[, "sigma2"], 
     xlim = c(1, 20), breaks = 15, xlab = expression(sigma^2 * phi),
     main = expression(tau^2 == 0.2 ~ "," ~ n == 1600))
hist(MLE_0.5_0.8_0.4_n1[, "phi"] * MLE_0.5_0.8_0.4_n1[, "sigma2"], 
     xlim = c(1, 20), breaks = 15, xlab = expression(sigma^2 * phi),
     main = expression(tau^2 == 0.8 ~ "," ~ n == 400))
hist(MLE_0.5_0.8_0.4_n2[, "phi"] * MLE_0.5_0.8_0.4_n2[, "sigma2"], 
     xlim = c(1, 20), breaks = 15, xlab = expression(sigma^2 * phi),
     main = expression(tau^2 == 0.8 ~ "," ~ n == 900))
hist(MLE_0.5_0.8_0.4_n3[, "phi"] * MLE_0.5_0.8_0.4_n3[, "sigma2"], 
     xlim = c(1, 20), breaks = 15, xlab = expression(sigma^2 * phi),
     main = expression(tau^2 == 0.8 ~ "," ~ n == 1600))
dev.off()

width <- 5.0
height <- 5.0
pointsize <- 16

setEPS()
postscript("./pics/hist-L_05_02_04_tau400.eps", width = width, height = height)
par(mfrow = c(1, 1))
hist(MLE_0.5_0.2_0.4_n1[, "tau2"], xlim = c(0, 0.4), breaks = 15,
     main = expression({tau[0]}^2 == 0.2 ~ "," ~ n == 400), xlab = expression(tau^2))
dev.off()

postscript("./pics/hist-L_05_02_04_tau900.eps", width = width, height = height)
par(mfrow = c(1, 1))
hist(MLE_0.5_0.2_0.4_n2[, "tau2"], xlim = c(0, 0.4), breaks = 15,
     main = expression(tau[0]^2 == 0.2 ~ "," ~ n == 900), xlab = expression(tau^2))
dev.off()

postscript("./pics/hist-L_05_02_04_tau1600.eps", width = width, height = height)
par(mfrow = c(1, 1))
hist(MLE_0.5_0.2_0.4_n3[, "tau2"], xlim = c(0, 0.4), breaks = 15,
     main = expression(tau[0]^2 == 0.2 ~ "," ~ n == 1600), xlab = expression(tau^2))
dev.off()

postscript("./pics/hist-L_05_02_04_sigma400.eps", width = width, height = height)
par(mfrow = c(1, 1))
hist(MLE_0.5_0.2_0.4_n1[, "sigma2"], xlim = c(0.2, 2.5), breaks = 15,
     main = expression(tau[0]^2 == 0.2 ~ "," ~ n == 400), xlab = expression(sigma^2))
dev.off()

postscript("./pics/hist-L_05_02_04_sigma900.eps", width = width, height = height)
par(mfrow = c(1, 1))
hist(MLE_0.5_0.2_0.4_n2[, "sigma2"], xlim = c(0.2, 2.5), breaks = 15,
     main = expression(tau[0]^2 == 0.2 ~ "," ~ n == 900), xlab = expression(sigma^2))
dev.off()

postscript("./pics/hist-L_05_02_04_sigma1600.eps", width = width, height = height)
par(mfrow = c(1, 1))
hist(MLE_0.5_0.2_0.4_n3[, "sigma2"], xlim = c(0.2, 2.5), breaks = 15,
     main = expression(tau[0]^2 == 0.2 ~ "," ~ n == 1600), xlab = expression(sigma^2))
dev.off()

postscript("./pics/hist-L_05_02_04_phi400.eps", width = width, height = height)
par(mfrow = c(1, 1))
hist(MLE_0.5_0.2_0.4_n1[, "phi"], xlim = c(1, 21), breaks = 15,
     main = expression(tau[0]^2 == 0.2 ~ "," ~ n == 400), xlab = expression(phi))
dev.off()

postscript("./pics/hist-L_05_02_04_phi900.eps", width = width, height = height)
par(mfrow = c(1, 1))
hist(MLE_0.5_0.2_0.4_n2[, "phi"], xlim = c(1, 21), breaks = 15,
     main = expression(tau[0]^2 == 0.2 ~ "," ~ n == 900), xlab = expression(phi))
dev.off()

postscript("./pics/hist-L_05_02_04_phi1600.eps", width = width, height = height)
par(mfrow = c(1, 1))
hist(MLE_0.5_0.2_0.4_n3[, "phi"], xlim = c(1, 21), breaks = 15,
     main = expression(tau[0]^2 == 0.2 ~ "," ~ n == 1600), xlab = expression(phi))
dev.off()

postscript("./pics/hist-L_05_02_04_c400.eps", width = width, height = height)
par(mfrow = c(1, 1))
hist(MLE_0.5_0.2_0.4_n1[, "phi"] * MLE_0.5_0.2_0.4_n1[, "sigma2"], 
     xlim = c(1, 20), breaks = 15, xlab = expression(sigma^2 * phi),
     main = expression(tau[0]^2 == 0.2 ~ "," ~ n == 400))
dev.off()

postscript("./pics/hist-L_05_02_04_c900.eps", width = width, height = height)
par(mfrow = c(1, 1))
hist(MLE_0.5_0.2_0.4_n2[, "phi"] * MLE_0.5_0.2_0.4_n2[, "sigma2"], 
     xlim = c(1, 20), breaks = 15, xlab = expression(sigma^2 * phi),
     main = expression(tau[0]^2 == 0.2 ~ "," ~ n == 900))
dev.off()

postscript("./pics/hist-L_05_02_04_c1600.eps", width = width, height = height)
par(mfrow = c(1, 1))
hist(MLE_0.5_0.2_0.4_n3[, "phi"] * MLE_0.5_0.2_0.4_n3[, "sigma2"], 
     xlim = c(1, 20), breaks = 15, xlab = expression(sigma^2 * phi),
     main = expression(tau[0]^2 == 0.2 ~ "," ~ n == 1600))
dev.off()


####### update plots for jrssb #######
width <- 5.0
height <- 5.0

png("./pics/fig4/a.png",
    width = width, height = height, units = "in", res = 1000)
par(mfrow = c(1, 1))
hist(MLE_0.5_0.2_0.4_n1[, "tau2"], xlim = c(0, 0.4), breaks = 15,
     main = expression({tau[0]}^2 == 0.2 ~ "," ~ n == 400), xlab = expression(tau^2))
dev.off()

png("./pics/fig4/b.png",
    width = width, height = height, units = "in", res = 1000)
par(mfrow = c(1, 1))
hist(MLE_0.5_0.2_0.4_n2[, "tau2"], xlim = c(0, 0.4), breaks = 15,
     main = expression(tau[0]^2 == 0.2 ~ "," ~ n == 900), xlab = expression(tau^2))
dev.off()

png("./pics/fig4/c.png",
    width = width, height = height, units = "in", res = 1000)
par(mfrow = c(1, 1))
hist(MLE_0.5_0.2_0.4_n3[, "tau2"], xlim = c(0, 0.4), breaks = 15,
     main = expression(tau[0]^2 == 0.2 ~ "," ~ n == 1600), xlab = expression(tau^2))
dev.off()

png("./pics/fig4/d.png",
    width = width, height = height, units = "in", res = 1000)
par(mfrow = c(1, 1))
hist(MLE_0.5_0.2_0.4_n1[, "sigma2"], xlim = c(0.2, 2.5), breaks = 15,
     main = expression(tau[0]^2 == 0.2 ~ "," ~ n == 400), xlab = expression(sigma^2))
dev.off()

png("./pics/fig4/e.png",
    width = width, height = height, units = "in", res = 1000)
par(mfrow = c(1, 1))
hist(MLE_0.5_0.2_0.4_n2[, "sigma2"], xlim = c(0.2, 2.5), breaks = 15,
     main = expression(tau[0]^2 == 0.2 ~ "," ~ n == 900), xlab = expression(sigma^2))
dev.off()

png("./pics/fig4/f.png",
    width = width, height = height, units = "in", res = 1000)
par(mfrow = c(1, 1))
hist(MLE_0.5_0.2_0.4_n3[, "sigma2"], xlim = c(0.2, 2.5), breaks = 15,
     main = expression(tau[0]^2 == 0.2 ~ "," ~ n == 1600), xlab = expression(sigma^2))
dev.off()

png("./pics/fig4/g.png",
    width = width, height = height, units = "in", res = 1000)
par(mfrow = c(1, 1))
hist(MLE_0.5_0.2_0.4_n1[, "phi"], xlim = c(1, 21), breaks = 15,
     main = expression(tau[0]^2 == 0.2 ~ "," ~ n == 400), xlab = expression(phi))
dev.off()

png("./pics/fig4/h.png",
    width = width, height = height, units = "in", res = 1000)
par(mfrow = c(1, 1))
hist(MLE_0.5_0.2_0.4_n2[, "phi"], xlim = c(1, 21), breaks = 15,
     main = expression(tau[0]^2 == 0.2 ~ "," ~ n == 900), xlab = expression(phi))
dev.off()

png("./pics/fig4/i.png",
    width = width, height = height, units = "in", res = 1000)
par(mfrow = c(1, 1))
hist(MLE_0.5_0.2_0.4_n3[, "phi"], xlim = c(1, 21), breaks = 15,
     main = expression(tau[0]^2 == 0.2 ~ "," ~ n == 1600), xlab = expression(phi))
dev.off()

png("./pics/fig4/j.png",
    width = width, height = height, units = "in", res = 1000)
par(mfrow = c(1, 1))
hist(MLE_0.5_0.2_0.4_n1[, "phi"] * MLE_0.5_0.2_0.4_n1[, "sigma2"], 
     xlim = c(1, 20), breaks = 15, xlab = expression(sigma^2 * phi),
     main = expression(tau[0]^2 == 0.2 ~ "," ~ n == 400))
dev.off()

png("./pics/fig4/k.png",
    width = width, height = height, units = "in", res = 1000)
par(mfrow = c(1, 1))
hist(MLE_0.5_0.2_0.4_n2[, "phi"] * MLE_0.5_0.2_0.4_n2[, "sigma2"], 
     xlim = c(1, 20), breaks = 15, xlab = expression(sigma^2 * phi),
     main = expression(tau[0]^2 == 0.2 ~ "," ~ n == 900))
dev.off()

png("./pics/fig4/l.png",
    width = width, height = height, units = "in", res = 1000)
par(mfrow = c(1, 1))
hist(MLE_0.5_0.2_0.4_n3[, "phi"] * MLE_0.5_0.2_0.4_n3[, "sigma2"], 
     xlim = c(1, 20), breaks = 15, xlab = expression(sigma^2 * phi),
     main = expression(tau[0]^2 == 0.2 ~ "," ~ n == 1600))
dev.off()

library(grid)
library(png)
setwd("./pics/fig4")
plots <- lapply(ll <- list.files(patt='*[.]png'),
                function(x){
                    img <- as.raster(readPNG(x))
                    rasterGrob(img, interpolate = FALSE)
                })
library(ggplot2)
library(gridExtra)
ggsave("Figure4.pdf",width=8.5, height=11, dpi = 1000, 
       marrangeGrob(grobs = plots, nrow=4, ncol=3, 
                    layout_matrix = matrix(seq_len(4*3), ncol = 3, byrow = TRUE),
                    top=NULL))



#---- check the percentiles, Bias, SD of tau2, phi, sigma2 and phi*sigma2 -----#
load("./data/sim2_0.5_0_0.15results.RData")
load("./data/sim2_0.5_0_0.4results.RData")
load("./data/sim2_0.5_0_1results.RData")
load("./data/sim2_0.5_0.2_0.15results.RData")
load("./data/sim2_0.5_0.2_0.4results.RData")
load("./data/sim2_0.5_0.2_1results.RData")
load("./data/sim2_0.5_0.8_0.15results.RData")
load("./data/sim2_0.5_0.8_0.4results.RData")
load("./data/sim2_0.5_0.8_1results.RData")

# tau2
round(percentile_BIAS_SD(MLE_0.5_0.2_0.15_n1[, "tau2"], true_value = 0.2), 3)
round(percentile_BIAS_SD(MLE_0.5_0.2_0.15_n2[, "tau2"], true_value = 0.2), 3)
round(percentile_BIAS_SD(MLE_0.5_0.2_0.15_n3[, "tau2"], true_value = 0.2), 3)
round(percentile_BIAS_SD(MLE_0.5_0.2_0.4_n1[, "tau2"], true_value = 0.2), 3)
round(percentile_BIAS_SD(MLE_0.5_0.2_0.4_n2[, "tau2"], true_value = 0.2), 3)
round(percentile_BIAS_SD(MLE_0.5_0.2_0.4_n3[, "tau2"], true_value = 0.2), 3)
round(percentile_BIAS_SD(MLE_0.5_0.2_1_n1[, "tau2"], true_value = 0.2), 3)
round(percentile_BIAS_SD(MLE_0.5_0.2_1_n2[, "tau2"], true_value = 0.2), 3)
round(percentile_BIAS_SD(MLE_0.5_0.2_1_n3[, "tau2"], true_value = 0.2), 3)

round(percentile_BIAS_SD(MLE_0.5_0.8_0.15_n1[, "tau2"], true_value = 0.8), 3)
round(percentile_BIAS_SD(MLE_0.5_0.8_0.15_n2[, "tau2"], true_value = 0.8), 3)
round(percentile_BIAS_SD(MLE_0.5_0.8_0.15_n3[, "tau2"], true_value = 0.8), 3)
round(percentile_BIAS_SD(MLE_0.5_0.8_0.4_n1[, "tau2"], true_value = 0.8), 3)
round(percentile_BIAS_SD(MLE_0.5_0.8_0.4_n2[, "tau2"], true_value = 0.8), 3)
round(percentile_BIAS_SD(MLE_0.5_0.8_0.4_n3[, "tau2"], true_value = 0.8), 3)
round(percentile_BIAS_SD(MLE_0.5_0.8_1_n1[, "tau2"], true_value = 0.8), 3)
round(percentile_BIAS_SD(MLE_0.5_0.8_1_n2[, "tau2"], true_value = 0.8), 3)
round(percentile_BIAS_SD(MLE_0.5_0.8_1_n3[, "tau2"], true_value = 0.8), 3)

# sigma2
percentile_BIAS_SD(MLE_0.5_0_0.15_n1[, "sigma2"], true_value = 1)
percentile_BIAS_SD(MLE_0.5_0_0.15_n2[, "sigma2"], true_value = 1)
percentile_BIAS_SD(MLE_0.5_0_0.15_n3[, "sigma2"], true_value = 1)
percentile_BIAS_SD(MLE_0.5_0_0.4_n1[, "sigma2"], true_value = 1)
percentile_BIAS_SD(MLE_0.5_0_0.4_n2[, "sigma2"], true_value = 1)
percentile_BIAS_SD(MLE_0.5_0_0.4_n3[, "sigma2"], true_value = 1)
percentile_BIAS_SD(MLE_0.5_0_1_n1[, "sigma2"], true_value = 1)
percentile_BIAS_SD(MLE_0.5_0_1_n2[, "sigma2"], true_value = 1)
percentile_BIAS_SD(MLE_0.5_0_1_n3[, "sigma2"], true_value = 1)

percentile_BIAS_SD(MLE_0.5_0.2_0.15_n1[, "sigma2"], true_value = 1)
percentile_BIAS_SD(MLE_0.5_0.2_0.15_n2[, "sigma2"], true_value = 1)
percentile_BIAS_SD(MLE_0.5_0.2_0.15_n3[, "sigma2"], true_value = 1)
percentile_BIAS_SD(MLE_0.5_0.2_0.4_n1[, "sigma2"], true_value = 1)
percentile_BIAS_SD(MLE_0.5_0.2_0.4_n2[, "sigma2"], true_value = 1)
percentile_BIAS_SD(MLE_0.5_0.2_0.4_n3[, "sigma2"], true_value = 1)
percentile_BIAS_SD(MLE_0.5_0.2_1_n1[, "sigma2"], true_value = 1)
percentile_BIAS_SD(MLE_0.5_0.2_1_n2[, "sigma2"], true_value = 1)
percentile_BIAS_SD(MLE_0.5_0.2_1_n3[, "sigma2"], true_value = 1)

percentile_BIAS_SD(MLE_0.5_0.8_0.15_n1[, "sigma2"], true_value = 1)
percentile_BIAS_SD(MLE_0.5_0.8_0.15_n2[, "sigma2"], true_value = 1)
percentile_BIAS_SD(MLE_0.5_0.8_0.15_n3[, "sigma2"], true_value = 1)
percentile_BIAS_SD(MLE_0.5_0.8_0.4_n1[, "sigma2"], true_value = 1)
percentile_BIAS_SD(MLE_0.5_0.8_0.4_n2[, "sigma2"], true_value = 1)
percentile_BIAS_SD(MLE_0.5_0.8_0.4_n3[, "sigma2"], true_value = 1)
percentile_BIAS_SD(MLE_0.5_0.8_1_n1[, "sigma2"], true_value = 1)
percentile_BIAS_SD(MLE_0.5_0.8_1_n2[, "sigma2"], true_value = 1)
percentile_BIAS_SD(MLE_0.5_0.8_1_n3[, "sigma2"], true_value = 1)

# phi
percentile_BIAS_SD(MLE_0.5_0_0.15_n1[, "phi"], true_value = phi0.5_1)
percentile_BIAS_SD(MLE_0.5_0_0.15_n2[, "phi"], true_value = phi0.5_1)
percentile_BIAS_SD(MLE_0.5_0_0.15_n3[, "phi"], true_value = phi0.5_1)
percentile_BIAS_SD(MLE_0.5_0_0.4_n1[, "phi"], true_value = phi0.5_2)
percentile_BIAS_SD(MLE_0.5_0_0.4_n2[, "phi"], true_value = phi0.5_2)
percentile_BIAS_SD(MLE_0.5_0_0.4_n3[, "phi"], true_value = phi0.5_2)
percentile_BIAS_SD(MLE_0.5_0_1_n1[, "phi"], true_value = phi0.5_3)
percentile_BIAS_SD(MLE_0.5_0_1_n2[, "phi"], true_value = phi0.5_3)
percentile_BIAS_SD(MLE_0.5_0_1_n3[, "phi"], true_value = phi0.5_3)

percentile_BIAS_SD(MLE_0.5_0.2_0.15_n1[, "phi"], true_value = phi0.5_1)
percentile_BIAS_SD(MLE_0.5_0.2_0.15_n2[, "phi"], true_value = phi0.5_1)
percentile_BIAS_SD(MLE_0.5_0.2_0.15_n3[, "phi"], true_value = phi0.5_1)
percentile_BIAS_SD(MLE_0.5_0.2_0.4_n1[, "phi"], true_value = phi0.5_2)
percentile_BIAS_SD(MLE_0.5_0.2_0.4_n2[, "phi"], true_value = phi0.5_2)
percentile_BIAS_SD(MLE_0.5_0.2_0.4_n3[, "phi"], true_value = phi0.5_2)
percentile_BIAS_SD(MLE_0.5_0.2_1_n1[, "phi"], true_value = phi0.5_3)
percentile_BIAS_SD(MLE_0.5_0.2_1_n2[, "phi"], true_value = phi0.5_3)
percentile_BIAS_SD(MLE_0.5_0.2_1_n3[, "phi"], true_value = phi0.5_3)

percentile_BIAS_SD(MLE_0.5_0.8_0.15_n1[, "phi"], true_value = phi0.5_1)
percentile_BIAS_SD(MLE_0.5_0.8_0.15_n2[, "phi"], true_value = phi0.5_1)
percentile_BIAS_SD(MLE_0.5_0.8_0.15_n3[, "phi"], true_value = phi0.5_1)
percentile_BIAS_SD(MLE_0.5_0.8_0.4_n1[, "phi"], true_value = phi0.5_2)
percentile_BIAS_SD(MLE_0.5_0.8_0.4_n2[, "phi"], true_value = phi0.5_2)
percentile_BIAS_SD(MLE_0.5_0.8_0.4_n3[, "phi"], true_value = phi0.5_2)
percentile_BIAS_SD(MLE_0.5_0.8_1_n1[, "phi"], true_value = phi0.5_3)
percentile_BIAS_SD(MLE_0.5_0.8_1_n2[, "phi"], true_value = phi0.5_3)
percentile_BIAS_SD(MLE_0.5_0.8_1_n3[, "phi"], true_value = phi0.5_3)


# phi
percentile_BIAS_SD(MLE_0.5_0_0.15_n1[, "phi"] * MLE_0.5_0_0.15_n1[, "sigma2"], true_value = phi0.5_1)
percentile_BIAS_SD(MLE_0.5_0_0.15_n2[, "phi"] * MLE_0.5_0_0.15_n2[, "sigma2"], true_value = phi0.5_1)
percentile_BIAS_SD(MLE_0.5_0_0.15_n3[, "phi"] * MLE_0.5_0_0.15_n3[, "sigma2"], true_value = phi0.5_1)
percentile_BIAS_SD(MLE_0.5_0_0.4_n1[, "phi"] * MLE_0.5_0_0.4_n1[, "sigma2"], true_value = phi0.5_2)
percentile_BIAS_SD(MLE_0.5_0_0.4_n2[, "phi"] * MLE_0.5_0_0.4_n2[, "sigma2"], true_value = phi0.5_2)
percentile_BIAS_SD(MLE_0.5_0_0.4_n3[, "phi"] * MLE_0.5_0_0.4_n3[, "sigma2"], true_value = phi0.5_2)
percentile_BIAS_SD(MLE_0.5_0_1_n1[, "phi"] * MLE_0.5_0_1_n1[, "sigma2"], true_value = phi0.5_3)
percentile_BIAS_SD(MLE_0.5_0_1_n2[, "phi"] * MLE_0.5_0_1_n2[, "sigma2"], true_value = phi0.5_3)
percentile_BIAS_SD(MLE_0.5_0_1_n3[, "phi"] * MLE_0.5_0_1_n3[, "sigma2"], true_value = phi0.5_3)

percentile_BIAS_SD(MLE_0.5_0.2_0.15_n1[, "phi"] * MLE_0.5_0.2_0.15_n1[, "sigma2"], true_value = phi0.5_1)
percentile_BIAS_SD(MLE_0.5_0.2_0.15_n2[, "phi"] * MLE_0.5_0.2_0.15_n2[, "sigma2"], true_value = phi0.5_1)
percentile_BIAS_SD(MLE_0.5_0.2_0.15_n3[, "phi"] * MLE_0.5_0.2_0.15_n3[, "sigma2"], true_value = phi0.5_1)
percentile_BIAS_SD(MLE_0.5_0.2_0.4_n1[, "phi"] * MLE_0.5_0.2_0.4_n1[, "sigma2"], true_value = phi0.5_2)
percentile_BIAS_SD(MLE_0.5_0.2_0.4_n2[, "phi"] * MLE_0.5_0.2_0.4_n2[, "sigma2"], true_value = phi0.5_2)
percentile_BIAS_SD(MLE_0.5_0.2_0.4_n3[, "phi"] * MLE_0.5_0.2_0.4_n3[, "sigma2"], true_value = phi0.5_2)
percentile_BIAS_SD(MLE_0.5_0.2_1_n1[, "phi"] * MLE_0.5_0.2_1_n1[, "sigma2"], true_value = phi0.5_3)
percentile_BIAS_SD(MLE_0.5_0.2_1_n2[, "phi"] * MLE_0.5_0.2_1_n2[, "sigma2"], true_value = phi0.5_3)
percentile_BIAS_SD(MLE_0.5_0.2_1_n3[, "phi"] * MLE_0.5_0.2_1_n3[, "sigma2"], true_value = phi0.5_3)

percentile_BIAS_SD(MLE_0.5_0.8_0.15_n1[, "phi"] * MLE_0.5_0.8_0.15_n1[, "sigma2"], true_value = phi0.5_1)
percentile_BIAS_SD(MLE_0.5_0.8_0.15_n2[, "phi"] * MLE_0.5_0.8_0.15_n2[, "sigma2"], true_value = phi0.5_1)
percentile_BIAS_SD(MLE_0.5_0.8_0.15_n3[, "phi"] * MLE_0.5_0.8_0.15_n3[, "sigma2"], true_value = phi0.5_1)
percentile_BIAS_SD(MLE_0.5_0.8_0.4_n1[, "phi"] * MLE_0.5_0.8_0.4_n1[, "sigma2"], true_value = phi0.5_2)
percentile_BIAS_SD(MLE_0.5_0.8_0.4_n2[, "phi"] * MLE_0.5_0.8_0.4_n2[, "sigma2"], true_value = phi0.5_2)
percentile_BIAS_SD(MLE_0.5_0.8_0.4_n3[, "phi"] * MLE_0.5_0.8_0.4_n3[, "sigma2"], true_value = phi0.5_2)
percentile_BIAS_SD(MLE_0.5_0.8_1_n1[, "phi"] * MLE_0.5_0.8_1_n1[, "sigma2"], true_value = phi0.5_3)
percentile_BIAS_SD(MLE_0.5_0.8_1_n2[, "phi"] * MLE_0.5_0.8_1_n2[, "sigma2"], true_value = phi0.5_3)
percentile_BIAS_SD(MLE_0.5_0.8_1_n3[, "phi"] * MLE_0.5_0.8_1_n3[, "sigma2"], true_value = phi0.5_3)


# check the MLE when fixing phi
load("./data/sim3_MLE_0.5_0.2_0.4_fix0.8phiresults.RData")
load("./data/sim3_MLE_0.5_0.2_0.4_fix1.0phiresults.RData")
load("./data/sim3_MLE_0.5_0.2_0.4_fix1.2phiresults.RData")

Mean_MLE_0.5_0.2_0.4_fix0.8phi <- matrix(NA, 8, 2)

for (i in 1:8){
    for (j in 1:2){
        Mean_MLE_0.5_0.2_0.4_fix0.8phi[i, j] <- mean(MLE_0.5_0.2_0.4_fix0.8phi[, j, i])
    }
}






