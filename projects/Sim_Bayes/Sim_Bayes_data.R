#setwd("/Users/luzhang/Documents/github/nugget_consistency")
rm(list = ls())
library(fields)
library(MASS)

# simulation setting #

nu = 0.5
n1 = 400; n2 = 900; n3 = 1600
sigmasq = 1
tausq = 0.5; 
phi = - log(0.05) / 0.3

# obvervation and prediction locations #

set.seed(1234)
grid_edge = seq(from = 0.005, to = 0.995, by = 0.015)
grid = as.matrix(expand.grid(grid_edge, grid_edge))
rand_grid = grid + 
  matrix(runif(2 * (length(grid_edge)^2), min = -0.005, max = 0.005), ncol = 2)
pick_ind = sample.int((length(grid_edge)^2), n3, replace = FALSE)
obs_loc = rand_grid[pick_ind, ]

grid_edge_pred = seq(from = 0.01, to = 0.99, by = 0.02)
pred_loc = as.matrix(expand.grid(grid_edge_pred, grid_edge_pred)) # loc for prediction

setEPS()
postscript("./pics/obs_locations_sim_bayes.eps")
plot(rand_grid, typ="n", cex = 0.1) # plot observed location sets
points(obs_loc, cex = 0.5, col = 2, pch = 16) # red
points(obs_loc[1:n2, ], cex = 0.5, col = 3, pch = 16) # green
points(obs_loc[1:n1, ], cex = 0.5, col = 4, pch = 16) # blue
dev.off()

# Generate simulations #
coords_total = obs_loc # 900 locations
N_coords = nrow(coords_total)
D = rdist(coords_total)
COV = sigmasq * exp(- phi * D) + tausq * diag(N_coords)
Sim_Bayes_y = mvrnorm(n = 1, mu = rep(0, N_coords), Sigma = COV)

save(list = ls(all.names = TRUE), 
     file = "./data/sim_bayes_data.RData", 
     envir = .GlobalEnv)


# plot interpolated surface of the first simulation 18 settings#


