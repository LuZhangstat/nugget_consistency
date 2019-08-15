setwd("/Users/luzhang/Documents/github/nugget_consistency")
rm(list = ls())
load("./data/simdata.RData")
source("./projects/utils.R")
library(fields)
D_obs = D[1:n2, 1:n2]
tau2_grid = seq(0, 1, length.out = 30)
sigma2_grid = seq(0.2, 4.2, length.out = 30)
phi_0.5_grid = seq(1 / Matern.cor.to.range(d = 1.2, nu = 0.5, cor.target = 0.05), 
                   1 / Matern.cor.to.range(d = 0.1, nu = 0.5, cor.target = 0.05),
                     length.out = 30)

L_0.5_0_0.4 = matrix(NA, nrow = length(tau2_grid), 
                            ncol = length(phi_0.5_grid))
L_0.5_0_0.4_fix = matrix(NA, nrow = length(tau2_grid), 
                     ncol = length(phi_0.5_grid))
L_0.5_0_0.4_nug = matrix(NA, nrow = length(phi_0.5_grid), 
                         ncol = length(phi_0.5_grid))
L_0.5_0_0.4_wrong_nug = matrix(NA, nrow = length(sigma2_grid), 
                         ncol = length(phi_0.5_grid))

L_0.5_0.2_0.4 = matrix(NA, nrow = length(tau2_grid), 
                       ncol = length(phi_0.5_grid))
L_0.5_0.2_0.4_fix = matrix(NA, nrow = length(tau2_grid), 
                           ncol = length(phi_0.5_grid))
L_0.5_0.2_0.4_nug = matrix(NA, nrow = length(sigma2_grid), 
                         ncol = length(phi_0.5_grid))
L_0.5_0.2_0.4_wrong_nug = matrix(NA, nrow = length(sigma2_grid), 
                           ncol = length(phi_0.5_grid))

L_0.5_0.8_0.4 = matrix(NA, nrow = length(tau2_grid), 
                       ncol = length(phi_0.5_grid))
L_0.5_0.8_0.4_fix = matrix(NA, nrow = length(tau2_grid), 
                           ncol = length(phi_0.5_grid))
L_0.5_0.8_0.4_nug = matrix(NA, nrow = length(sigma2_grid), 
                         ncol = length(phi_0.5_grid))
L_0.5_0.8_0.4_wrong_nug = matrix(NA, nrow = length(sigma2_grid), 
                           ncol = length(phi_0.5_grid))


# for fixed sigma2 = 1
t <- proc.time()
for(i in 1:length(tau2_grid)){
  for(j in 1:length(phi_0.5_grid)){
    K = exp(- phi_0.5_grid[j] * D_obs) + tau2_grid[i] * diag(n2) 
    chol_K = chol(K)
    #L_0.5_0_0.4_fix[i, j] = 
    #  sum(apply(Sim_data_0.5["0", "0.5_0.4", 1:5, 1:n2],
    #            MARGIN = 1, loglik_chol, chol_COV = chol_K, sigma2 = 1))
    L_0.5_0_0.4_fix[i, j] = loglik_chol(Sim_data_0.5["0", "0.5_0.4", 1, 1:n2],
                                        chol_K, sigma2 = 1)
  }
}
proc.time() - t

# for fixing sigma^2 phi^{2nv} = phi_0
t <- proc.time()
for(i in 1:length(tau2_grid)){
  for(j in 1:length(phi_0.5_grid)){
    sigma2_pick = phi0.5_2 / phi_0.5_grid[j]
    K = exp(- phi_0.5_grid[j] * D_obs) + 
      (tau2_grid[i] / sigma2_pick) * diag(n2) 
    chol_K = chol(K)
    #L_0.5_0_0.4[i, j] = 
    #  sum(apply(Sim_data_0.5["0", "0.5_0.4", 1:5, 1:n2],
    #            MARGIN = 1, loglik_chol, chol_COV = chol_K, 
    #            sigma2 = sigma2_pick))
    L_0.5_0_0.4[i, j] = loglik_chol(Sim_data_0.5["0", "0.5_0.4", 1, 1:n2],
                                    chol_COV = chol_K, sigma2 = sigma2_pick)
  }
}
proc.time() - t

# check log likelihood with fixed tau2 = tau21 = 0
t <- proc.time()
for(i in 1:length(sigma2_grid)){
  for(j in 1:length(phi_0.5_grid)){
    
    K = exp(- (phi_0.5_grid[j]) * D_obs) + (tau21 / sigma2_grid[i]) * diag(n2) 
    chol_K = chol(K)
    #L_0.5_0_0.4_nug[i, j] = 
    #  sum(apply(Sim_data_0.5["0", "0.5_0.4", 1:5, 1:n2],
    #            MARGIN = 1, loglik_chol, chol_COV = chol_K, 
    #            sigma2 = sigma2_grid[i]))
    L_0.5_0_0.4_nug[i, j] = 
      loglik_chol(Sim_data_0.5["0", "0.5_0.4", 1, 1:n2], 
                  chol_COV = chol_K, sigma2 = sigma2_grid[i])
  }
}
proc.time() - t

# check log likelihood with fixed tau2 = 0.5
t <- proc.time()
for(i in 1:length(sigma2_grid)){
  for(j in 1:length(phi_0.5_grid)){
    
    K = exp(- (phi_0.5_grid[j]) * D_obs) + (0.5 / sigma2_grid[i]) * diag(n2) 
    chol_K = chol(K)
    #L_0.5_0_0.4_wrong_nug[i, j] = 
    #  sum(apply(Sim_data_0.5["0", "0.5_0.4", 1:5, 1:n2],
    #            MARGIN = 1, loglik_chol, chol_COV = chol_K, 
    #            sigma2 = sigma2_grid[i]))
    L_0.5_0_0.4_wrong_nug[i, j] = 
      loglik_chol(Sim_data_0.5["0", "0.5_0.4", 1, 1:n2], 
                  chol_COV = chol_K, sigma2 = sigma2_grid[i])
  }
}
proc.time() - t


# for fixed sigma2 = 1, tau2 = 0.2
t <- proc.time()
for(i in 1:length(tau2_grid)){
  for(j in 1:length(phi_0.5_grid)){
    K = exp(- phi_0.5_grid[j] * D_obs) + tau2_grid[i] * diag(n2) 
    chol_K = chol(K)
    #L_0.5_0.2_0.4_fix[i, j] = 
    #  sum(apply(Sim_data_0.5["0.2", "0.5_0.4", 1:5, 1:n2],
    #            MARGIN = 1, loglik_chol, chol_COV = chol_K, sigma2 = 1))
    L_0.5_0.2_0.4_fix[i, j] = 
      loglik_chol(Sim_data_0.5["0.2", "0.5_0.4", 1, 1:n2],
                  chol_COV = chol_K, sigma2 = 1)
  }
}
proc.time() - t

# for fixing sigma^2 phi^{2nv} = phi_0, tau20 = 0.2
t <- proc.time()
for(i in 1:length(tau2_grid)){
  for(j in 1:length(phi_0.5_grid)){
    sigma2_pick = phi0.5_2 / phi_0.5_grid[j]
    K = exp(- phi_0.5_grid[j] * D_obs) + 
      (tau2_grid[i] / sigma2_pick) * diag(n2) 
    chol_K = chol(K)
    #L_0.5_0.2_0.4[i, j] = 
    #  sum(apply(Sim_data_0.5["0.2", "0.5_0.4", 1:5, 1:n2],
    #            MARGIN = 1, loglik_chol, chol_COV = chol_K, 
    #            sigma2 = sigma2_pick))
    L_0.5_0.2_0.4[i, j] = 
      loglik_chol(Sim_data_0.5["0.2", "0.5_0.4", 1, 1:n2], 
                  chol_COV = chol_K, sigma2 = sigma2_pick)
  }
}
proc.time() - t

# check log likelihood with fixed tau2 = tau22 = 0.2
t <- proc.time()
for(i in 1:length(sigma2_grid)){
  for(j in 1:length(phi_0.5_grid)){
    
    K = exp(- (phi_0.5_grid[j]) * D_obs) + (tau22 / sigma2_grid[i]) * diag(n2) 
    chol_K = chol(K)
    #L_0.5_0.2_0.4_nug[i, j] = 
    #  sum(apply(Sim_data_0.5["0.2", "0.5_0.4", 1:5, 1:n2],
    #            MARGIN = 1, loglik_chol, chol_COV = chol_K, 
    #            sigma2 = sigma2_grid[i]))
    L_0.5_0.2_0.4_nug[i, j] = 
      loglik_chol(Sim_data_0.5["0.2", "0.5_0.4", 1, 1:n2], 
                  chol_COV = chol_K, sigma2 = sigma2_grid[i])
  }
}
proc.time() - t

# check log likelihood with fixed tau2 = 0.5
t <- proc.time()
for(i in 1:length(sigma2_grid)){
  for(j in 1:length(phi_0.5_grid)){
    
    K = exp(- (phi_0.5_grid[j]) * D_obs) + (0.5 / sigma2_grid[i]) * diag(n2) 
    chol_K = chol(K)
    #L_0.5_0.2_0.4_wrong_nug[i, j] = 
    #  sum(apply(Sim_data_0.5["0.2", "0.5_0.4", 1:5, 1:n2],
    #            MARGIN = 1, loglik_chol, chol_COV = chol_K, 
    #            sigma2 = sigma2_grid[i]))
    L_0.5_0.2_0.4_wrong_nug[i, j] = 
      loglik_chol(Sim_data_0.5["0.2", "0.5_0.4", 1, 1:n2], 
                  chol_COV = chol_K, sigma2 = sigma2_grid[i])
  }
}
proc.time() - t


# for fixed sigma2 = 1, tau2 = 0.8
t <- proc.time()
for(i in 1:length(tau2_grid)){
  for(j in 1:length(phi_0.5_grid)){
    K = exp(- phi_0.5_grid[j] * D_obs) + tau2_grid[i] * diag(n2) 
    chol_K = chol(K)
    #L_0.5_0.8_0.4_fix[i, j] = 
    #  sum(apply(Sim_data_0.5["0.8", "0.5_0.4", 1:5, 1:n2],
    #            MARGIN = 1, loglik_chol, chol_COV = chol_K, sigma2 = 1))
    L_0.5_0.8_0.4_fix[i, j] = 
      loglik_chol(Sim_data_0.5["0.8", "0.5_0.4", 1, 1:n2], 
                  chol_COV = chol_K, sigma2 = 1)
  }
}
proc.time() - t

# for fixing sigma^2 phi^{2nv} = phi_0, tau20 = 0.8
t <- proc.time()
for(i in 1:length(tau2_grid)){
  for(j in 1:length(phi_0.5_grid)){
    sigma2_pick = phi0.5_2 / phi_0.5_grid[j]
    K = exp(- phi_0.5_grid[j] * D_obs) + 
      (tau2_grid[i] / sigma2_pick) * diag(n2) 
    chol_K = chol(K)
    #L_0.5_0.8_0.4[i, j] = 
    #  sum(apply(Sim_data_0.5["0.8", "0.5_0.4", 1:5, 1:n2],
    #            MARGIN = 1, loglik_chol, chol_COV = chol_K, 
    #            sigma2 = sigma2_pick))
    L_0.5_0.8_0.4[i, j] = 
      loglik_chol(Sim_data_0.5["0.8", "0.5_0.4", 1, 1:n2],
                  chol_COV = chol_K, sigma2 = sigma2_pick)
  }
}
proc.time() - t

# check log likelihood with fixed tau2 = tau23 = 0.8
t <- proc.time()
for(i in 1:length(sigma2_grid)){
  for(j in 1:length(phi_0.5_grid)){
    
    K = exp(- (phi_0.5_grid[j]) * D_obs) + (tau23 / sigma2_grid[i]) * diag(n2) 
    chol_K = chol(K)
    #L_0.5_0.8_0.4_nug[i, j] = 
    #  sum(apply(Sim_data_0.5["0.8", "0.5_0.4", 1:5, 1:n2],
    #            MARGIN = 1, loglik_chol, chol_COV = chol_K, 
    #            sigma2 = sigma2_grid[i]))
    L_0.5_0.8_0.4_nug[i, j] = 
      loglik_chol(Sim_data_0.5["0.8", "0.5_0.4", 1, 1:n2],
                  chol_COV = chol_K, sigma2 = sigma2_grid[i])
  }
}
proc.time() - t

# check log likelihood with fixed tau2 = 0.5
t <- proc.time()
for(i in 1:length(sigma2_grid)){
  for(j in 1:length(phi_0.5_grid)){
    
    K = exp(- (phi_0.5_grid[j]) * D_obs) + (0.5 / sigma2_grid[i]) * diag(n2) 
    chol_K = chol(K)
    #L_0.5_0.8_0.4_wrong_nug[i, j] = 
    #  sum(apply(Sim_data_0.5["0.8", "0.5_0.4", 1:5, 1:n2],
    #            MARGIN = 1, loglik_chol, chol_COV = chol_K, 
    #            sigma2 = sigma2_grid[i]))
    L_0.5_0.8_0.4_wrong_nug[i, j] = 
      loglik_chol(Sim_data_0.5["0.8", "0.5_0.4", 1, 1:n2],
                  chol_COV = chol_K, sigma2 = sigma2_grid[i])
  }
}
proc.time() - t

# Check the interpolated map
#load("./data/sim1results.RData")
library(coda)
library(spBayes)
library(MBA)
library(classInt)
library(RColorBrewer)
library(sp)
library(graphics)

h <- 12
grid_0.5 = expand.grid(tau2_grid, phi_0.5_grid)
grid_0.5_nug = expand.grid(sigma2_grid, phi_0.5_grid)
surf.L_0.5_0_0.4 <- mba.surf(cbind(grid_0.5, c(L_0.5_0_0.4)), no.X = 300, 
                             no.Y = 300, exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.L_0.5_0_0.4_fix <- mba.surf(cbind(grid_0.5, c(L_0.5_0_0.4_fix)), 
                                 no.X = 300, no.Y = 300, 
                                 exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.L_0.5_0_0.4_nug <- mba.surf(cbind(grid_0.5_nug, c(L_0.5_0_0.4_nug)), 
                                 no.X = 300, no.Y = 300, 
                                 exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.L_0.5_0_0.4_wrong_nug <- 
  mba.surf(cbind(grid_0.5_nug, c(L_0.5_0_0.4_wrong_nug)), 
           no.X = 300, no.Y = 300, exten = TRUE, sp = TRUE, h = h)$xyz.est

surf.L_0.5_0.2_0.4 <- 
  mba.surf(cbind(grid_0.5, c(L_0.5_0.2_0.4)), no.X = 300, 
           no.Y = 300, exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.L_0.5_0.2_0.4_fix <- 
  mba.surf(cbind(grid_0.5, c(L_0.5_0.2_0.4_fix)), no.X = 300, no.Y = 300, 
           exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.L_0.5_0.2_0.4_nug <- 
  mba.surf(cbind(grid_0.5_nug, c(L_0.5_0.2_0.4_nug)), no.X = 300, no.Y = 300, 
           exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.L_0.5_0.2_0.4_wrong_nug <- 
  mba.surf(cbind(grid_0.5_nug, c(L_0.5_0.2_0.4_wrong_nug)), 
           no.X = 300, no.Y = 300, exten = TRUE, sp = TRUE, h = h)$xyz.est

surf.L_0.5_0.8_0.4 <- 
  mba.surf(cbind(grid_0.5, c(L_0.5_0.8_0.4)), no.X = 300, 
           no.Y = 300, exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.L_0.5_0.8_0.4_fix <- 
  mba.surf(cbind(grid_0.5, c(L_0.5_0.8_0.4_fix)), no.X = 300, no.Y = 300, 
           exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.L_0.5_0.8_0.4_nug <- 
  mba.surf(cbind(grid_0.5_nug, c(L_0.5_0.8_0.4_nug)), no.X = 300, no.Y = 300, 
           exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.L_0.5_0.8_0.4_wrong_nug <- 
  mba.surf(cbind(grid_0.5_nug, c(L_0.5_0.8_0.4_wrong_nug)), 
           no.X = 300, no.Y = 300, exten = TRUE, sp = TRUE, h = h)$xyz.est

# surf.brks.L_0.5_0_0.4 <- 
#   classIntervals(c(surf.L_0.5_0_0.4$z, surf.L_0.5_0_0.4_wrong_nug$z, 1012.719),
#                  500, 'pretty')$brks
# surf.brks.L_0.5_0.2_0.4 <- 
#   classIntervals(c(surf.L_0.5_0.2_0.4$z, surf.L_0.5_0.2_0.4_nug$z), 
#                  500, 'pretty')$brks
# surf.brks.L_0.5_0.8_0.4 <- 
#   classIntervals(c(surf.L_0.5_0.8_0.4_wrong_nug$z, surf.L_0.5_0.8_0.4_nug$z, 
#                    -5500), 500, 'pretty')$brks
surf.brks <- 
  classIntervals(c(surf.L_0.5_0_0.4$z, 
                   surf.L_0.5_0.2_0.4$z,
                   surf.L_0.5_0.8_0.4_wrong_nug$z, 
                   -1000, 220), 500, 'pretty')$brks

col.pal <- colorRampPalette(brewer.pal(9,'Greys')[9:1])
xlim <- c(0.0, 1.3) #c(0.39, 1.13) # c(0.0, 1.7)
xlim_nug <- c(0.2, 5.0)
# zlim.L_0.5_0_0.4 <- range(c(surf.L_0.5_0_0.4[["z"]], 
#                             surf.L_0.5_0_0.4_wrong_nug[["z"]], 1012.719))
# zlim.L_0.5_0.2_0.4 <- range(c(surf.L_0.5_0.2_0.4[["z"]], 
#                               surf.L_0.5_0.2_0.4_nug[["z"]]))
# zlim.L_0.5_0.8_0.4 <- range(c(surf.L_0.5_0.8_0.4_wrong_nug[["z"]], 
#                               surf.L_0.5_0.8_0.4_nug[["z"]], -5500))
zlim = c(-1000.00, 220)

width <- 5.0
height <- 5.0
pointsize <- 16

setEPS()
postscript("./pics/map-L_05_0_04.eps", width = width, height = height)
#png(paste("./pics/map-L_05_0_04.png", sep = ""), 
#    width = width, height = height, pointsize = pointsize, family = "Courier")
par(mfrow = c(1, 1))
i <- as.image.SpatialGridDataFrame(surf.L_0.5_0_0.4)
plot(grid_0.5, typ = "n", cex = 0.5, xlim = xlim, axes = FALSE, 
     ylab = expression(phi), xlab = expression(tau^2)) 
axis(2, las=1)
axis(1)
image.plot(i, add = TRUE, col = rev(col.pal(length(surf.brks) - 1)), 
           zlim = zlim)
contour(i$x, i$y, i$z, col = "black", add = TRUE, nlevels = 25, 
        zlim = zlim)
dev.off()


setEPS()
postscript("./pics/map-L_05_0_04_fix.eps", width = width, height = height)
#png(paste("./pics/map-L_05_0_04_fix.png", sep = ""), 
#    width = width, height = height, pointsize = pointsize, family = "Courier")
par(mfrow = c(1, 1))
i <- as.image.SpatialGridDataFrame(surf.L_0.5_0_0.4_fix)
plot(grid_0.5, typ = "n", cex = 0.5, xlim = xlim, axes = FALSE, 
     ylab = expression(phi), xlab = expression(tau^2)) 
axis(2, las=1)
axis(1)
image.plot(i, add = TRUE, 
           col = rev(col.pal(length(surf.brks) - 1)), 
           zlim = zlim)
contour(i$x, i$y, i$z, col = "black", add = TRUE, nlevels = 25, 
        zlim = zlim)
dev.off()

setEPS()
postscript("./pics/map-L_05_0_04_nug.eps", width = width, height = height)
#png(paste("./pics/map-L_05_0_04_nug.png", sep = ""), 
#    width = width, height = height, pointsize = pointsize, family = "Courier")
par(mfrow = c(1, 1))
i <- as.image.SpatialGridDataFrame(surf.L_0.5_0_0.4_nug)
plot(grid_0.5_nug, typ = "n", cex = 0.5, xlim = xlim_nug, axes = FALSE, 
     ylab = expression(phi), xlab = expression(sigma^2)) 
axis(2, las=1)
axis(1)
image.plot(i, add = TRUE, 
           col = rev(col.pal(length(surf.brks) - 1)), 
           zlim = zlim)
contour(i$x, i$y, i$z, col = "black", add = TRUE, nlevels = 25, 
        zlim = zlim)
dev.off()

setEPS()
postscript("./pics/map-L_05_0_04_wrong_nug.eps", width = width, height = height)
#png(paste("./pics/map-L_05_0_04_wrong_nug.png", sep = ""), 
#    width = width, height = height, pointsize = pointsize, family = "Courier")
par(mfrow = c(1, 1))
i <- as.image.SpatialGridDataFrame(surf.L_0.5_0_0.4_wrong_nug)
plot(grid_0.5_nug, typ = "n", cex = 0.5, xlim = xlim_nug, axes = FALSE, 
     ylab = expression(phi), xlab = expression(sigma^2)) 
axis(2, las=1)
axis(1)
image.plot(i, add = TRUE, 
           col = rev(col.pal(length(surf.brks) - 1)), 
           zlim = zlim)
contour(i$x, i$y, i$z, col = "black", add = TRUE, nlevels = 25,
        zlim = zlim)
dev.off()

setEPS()
postscript("./pics/map-L_05_02_04.eps", width = width, height = height)
#png(paste("./pics/map-L_05_02_04.png", sep = ""), 
#    width = width, height = height, pointsize = pointsize, family = "Courier")
par(mfrow = c(1, 1))
i <- as.image.SpatialGridDataFrame(surf.L_0.5_0.2_0.4)
plot(grid_0.5, typ = "n", cex = 0.5, xlim = xlim, axes = FALSE, 
     ylab = expression(phi), xlab = expression(tau^2)) 
axis(2, las=1)
axis(1)
image.plot(i, add = TRUE, 
           col = rev(col.pal(length(surf.brks) - 1)), 
           zlim = zlim)
contour(i$x, i$y, i$z, col = "black", add = TRUE, nlevels = 25, 
        zlim = zlim)
dev.off()

setEPS()
postscript("./pics/map-L_05_02_04_fix.eps", width = width, height = height)
#png(paste("./pics/map-L_05_02_04_fix.png", sep = ""), 
#    width = width, height = height, pointsize = pointsize, family = "Courier")
#par(mfrow = c(1, 1))
i <- as.image.SpatialGridDataFrame(surf.L_0.5_0.2_0.4_fix)
plot(grid_0.5, typ = "n", cex = 0.5, xlim = xlim, axes = FALSE, 
     ylab = expression(phi), xlab = expression(tau^2)) 
axis(2, las=1)
axis(1)
image.plot(i, add = TRUE, 
           col = rev(col.pal(length(surf.brks) - 1)), 
           zlim = zlim)
contour(i$x, i$y, i$z, col = "black", add = TRUE, nlevels = 25, 
        zlim = zlim)
dev.off()

setEPS()
postscript("./pics/map-L_05_02_04_nug.eps", width = width, height = height)
par(mfrow = c(1, 1))
i <- as.image.SpatialGridDataFrame(surf.L_0.5_0.2_0.4_nug)
plot(grid_0.5_nug, typ = "n", cex = 0.5, xlim = xlim_nug, axes = FALSE, 
     ylab = expression(phi), xlab = expression(sigma^2)) 
axis(2, las=1)
axis(1)
image.plot(i, add = TRUE, 
           col = rev(col.pal(length(surf.brks) - 1)), 
           zlim = zlim)
contour(i$x, i$y, i$z, col = "black", add = TRUE, nlevels = 25, 
        zlim = zlim)
dev.off()

setEPS()
postscript("./pics/map-L_05_02_04_wrong_nug.eps", width = width, height = height)
par(mfrow = c(1, 1))
i <- as.image.SpatialGridDataFrame(surf.L_0.5_0.2_0.4_wrong_nug)
plot(grid_0.5_nug, typ = "n", cex = 0.5, xlim = xlim_nug, axes = FALSE, 
     ylab = expression(phi), xlab = expression(sigma^2)) 
axis(2, las=1)
axis(1)
image.plot(i, add = TRUE, 
           col = rev(col.pal(length(surf.brks) - 1)), 
           zlim = zlim)
contour(i$x, i$y, i$z, col = "black", add = TRUE, nlevels = 25, 
        zlim = zlim)
dev.off()


setEPS()
postscript("./pics/map-L_05_08_04.eps", width = width, height = height)
#png(paste("./pics/map-L_05_08_04.png", sep = ""), 
#    width = width, height = height, pointsize = pointsize, family = "Courier")
par(mfrow = c(1, 1))
i <- as.image.SpatialGridDataFrame(surf.L_0.5_0.8_0.4)
plot(grid_0.5, typ = "n", cex = 0.5, xlim = xlim, axes = FALSE, 
     ylab = expression(phi), xlab = expression(tau^2)) 
axis(2, las=1)
axis(1)
image.plot(i, add = TRUE, 
           col = rev(col.pal(length(surf.brks) - 1)), 
           zlim = zlim)
contour(i$x, i$y, i$z, col = "black", add = TRUE, nlevels = 25, 
        zlim = zlim)
dev.off()

setEPS()
postscript("./pics/map-L_05_08_04_fix.eps", width = width, height = height)
#png(paste("./pics/map-L_05_08_04_fix.png", sep = ""), 
#    width = width, height = height, pointsize = pointsize, family = "Courier")
par(mfrow = c(1, 1))
i <- as.image.SpatialGridDataFrame(surf.L_0.5_0.8_0.4_fix)
plot(grid_0.5, typ = "n", cex = 0.5, xlim = xlim, axes = FALSE, 
     ylab = expression(phi), xlab = expression(tau^2)) 
axis(2, las=1)
axis(1)
image.plot(i, add = TRUE, 
           col = rev(col.pal(length(surf.brks) - 1)), 
           zlim = zlim)
contour(i$x, i$y, i$z, col = "black", add = TRUE, nlevels = 25, 
        zlim = zlim)
dev.off()

setEPS()
postscript("./pics/map-L_05_08_04_nug.eps", width = width, height = height)
#png(paste("./pics/map-L_05_08_04_nug.png", sep = ""), 
#    width = width, height = height, pointsize = pointsize, family = "Courier")
par(mfrow = c(1, 1))
i <- as.image.SpatialGridDataFrame(surf.L_0.5_0.8_0.4_nug)
plot(grid_0.5_nug, typ = "n", cex = 0.5, xlim = xlim_nug, axes = FALSE, 
     ylab = expression(phi), xlab = expression(sigma^2)) 
axis(2, las=1)
axis(1)
image.plot(i, add = TRUE, 
           col = rev(col.pal(length(surf.brks) - 1)), 
           zlim = zlim)
contour(i$x, i$y, i$z, col = "black", add = TRUE, nlevels = 25, 
        zlim = zlim)
dev.off()

setEPS()
postscript("./pics/map-L_05_08_04_wrong_nug.eps", width = width, height = height)
#png(paste("./pics/map-L_05_08_04_wrong_nug.png", sep = ""), 
#    width = width, height = height, pointsize = pointsize, family = "Courier")
par(mfrow = c(1, 1))
i <- as.image.SpatialGridDataFrame(surf.L_0.5_0.8_0.4_wrong_nug)
plot(grid_0.5_nug, typ = "n", cex = 0.5, xlim = xlim_nug, axes = FALSE, 
     ylab = expression(phi), xlab = expression(sigma^2)) 
axis(2, las=1)
axis(1)
image.plot(i, add = TRUE, 
           col = rev(col.pal(length(surf.brks) - 1)), 
           zlim = zlim)
contour(i$x, i$y, i$z, col = "black", add = TRUE, nlevels = 25, 
        zlim = zlim)
dev.off()


# comparision #
png(paste("./pics/map-sigma2phi=const_compar.png", sep = ""), 
    width = width * 3 / 2, height = height /2, pointsize = pointsize, 
    family = "Courier")
par(mfrow = c(1, 3))

dev.off()

png(paste("./pics/map-tau2fixtrue_compar.png", sep = ""), 
    width = width * 3 / 2, height = height /2, pointsize = pointsize, 
    family = "Courier")
par(mfrow = c(1, 3))

dev.off()





save(list = c("L_0.5_0_0.4", "L_0.5_0_0.4_fix", "L_0.5_0_0.4_nug", 
              "L_0.5_0_0.4_wrong_nug",
              "L_0.5_0.2_0.4", "L_0.5_0.2_0.4_fix", "L_0.5_0.2_0.4_nug",
              "L_0.5_0.2_0.4_wrong_nug",
              "L_0.5_0.8_0.4", "L_0.5_0.8_0.4_fix", "L_0.5_0.8_0.4_nug",
              "L_0.5_0.8_0.4_wrong_nug"), 
     file = "./data/sim1results.RData", 
     envir = .GlobalEnv)

