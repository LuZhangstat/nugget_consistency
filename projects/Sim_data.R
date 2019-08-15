setwd("/Users/luzhang/Documents/github/nugget_consistency")
rm(list = ls())
library(fields)
library(MASS)

# simulation setting #

nu1 = 0.5; nu2 = 1.5 
n1 = 400; n2 = 900; n3 = 1600
sigma2 = 1
tau21 = 0; tau22 = 0.2; tau23 = 0.8
tau2_opts = c(tau21, tau22, tau23)
phi0.5_1 = - log(0.05) / 0.15
phi0.5_2 = - log(0.05) / 0.4
phi0.5_3 = - log(0.05) / 1
phi0.5_opts = c(phi0.5_1, phi0.5_2, phi0.5_3)
phi0.5names  = c("0.5_0.15", "0.5_0.4", "0.5_1")

phi1.5_1 = 1 / Matern.cor.to.range(d = 0.15, nu = 1.5, cor.target = 0.05)
phi1.5_2 = 1 / Matern.cor.to.range(d = 0.4, nu = 1.5, cor.target = 0.05)
phi1.5_3 = 1 / Matern.cor.to.range(d = 1.0, nu = 1.5, cor.target = 0.05)
phi1.5_opts = c(phi1.5_1, phi1.5_2, phi1.5_3)
phi1.5names  = c("1.5_0.15", "1.5_0.4", "1.5_1")

# obvervation and prediction locations #

set.seed(123)
grid_edge = seq(from = 0.005, to = 0.995, by = 0.015)
grid = as.matrix(expand.grid(grid_edge, grid_edge))
rand_grid = grid + 
  matrix(runif(2 * (length(grid_edge)^2), min = -0.005, max = 0.005), ncol = 2)
pick_ind3 = sample.int((length(grid_edge)^2), n3, replace = FALSE)
obs_loc3 = rand_grid[pick_ind3, ]
#pick_ind2 = sample.int(n3, n2, replace = FALSE)
#obs_loc2 = obs_loc3[pick_ind2, ]
#pick_ind1 = sample.int(n2, n1, replace = FALSE)
#obs_loc1 = obs_loc2[pick_ind1, ]

grid_edge_pred = seq(from = 0.01, to = 0.99, by = 0.02)
pred_loc = as.matrix(expand.grid(grid_edge_pred, grid_edge_pred)) # loc for prediction

setEPS()
postscript("./pics/obs_locations.eps")
plot(rand_grid, typ="n", cex = 0.1) # plot observed location sets
points(obs_loc3, cex = 0.5, col = 2, pch = 16) # red
points(obs_loc3[1:n2, ], cex = 0.5, col = 3, pch = 16) # green
points(obs_loc3[1:n1, ], cex = 0.5, col = 4, pch = 16) # blue
points(pred_loc, cex = 0.3, col = 1, pch = 4) # blue
dev.off()

# Generate simulations #
coords_total = rbind(obs_loc3, pred_loc) # 1600 + 2500 = 4100 locations
N_coords = nrow(coords_total)
N_sim = 1000   # simulate 1000 realization for each setting
D = rdist(coords_total)
Sim_data_0.5 = array(data = NA, dim = c(length(tau2_opts), length(phi0.5_opts), 
                                        N_sim, N_coords), 
                     dimnames = list(tau2_opts = tau2_opts,
                                     phi0.5_opts = phi0.5names,
                                     No_sim = 1:N_sim,
                                     index = 1:N_coords))
for(i in 1:length(tau2_opts)){
  for(j in 1:length(phi0.5_opts)){
    COV = sigma2 * exp(- phi0.5_opts[j] * D) + tau2_opts[i] * diag(N_coords)
    Sim_data_0.5[i, j, , ] = 
      mvrnorm(n = N_sim, mu = rep(0, N_coords), Sigma = COV)
  }
}

Sim_data_1.5 = array(data = NA, dim = c(length(tau2_opts), length(phi1.5_opts), 
                                        N_sim, N_coords), 
                     dimnames = list(tau2_opts = tau2_opts,
                                     phi1.5_opts = phi1.5names,
                                     No_sim = 1:N_sim,
                                     index = 1:N_coords))
for(i in 1:length(tau2_opts)){
  for(j in 1:length(phi1.5_opts)){
    COV = sigma2 * Matern(D, alpha = phi1.5_opts[j], smoothness = 1.5) + 
      tau2_opts[i] * diag(N_coords)
    Sim_data_1.5[i, j, , ] = 
      mvrnorm(n = N_sim, mu = rep(0, N_coords), Sigma = COV)
  }
}

save(list = ls(all.names = TRUE), 
     file = "./data/simdata.RData", 
     envir = .GlobalEnv)


# plot interpolated surface of the first simulation 18 settings#


