#  simultation for prediction #
#setwd("/Users/luzhang/Documents/github/nugget_consistency")
rm(list = ls())
library(fields)
library(MASS)
#load("./data/simdata.RData")
source("./projects/utils_2.R")

# randomly generate 12,000 samples in [0,1]
set.seed(1234)
n = 12000
obs_loc <- runif(n)
pred_loc <- c(0.25, 0.5, 0.75)
coords_total = c(obs_loc, pred_loc) # 1600 + 2500 = 4100 locations
N_coords = length(coords_total)
D = rdist(coords_total)
phi = - log(0.05) / 0.4
sigma2 = 1.0
tau2 = 0.2

# select subset range from 200 to 1600 
Sam_N_list_d1 = c(500, seq(1000, 12000, by= 1000))
D_pred = D[12001:12003, 12001:12003]

# preallocation 
lower_l1 = rep(0, length(Sam_N_list_d1))

# mean squared prediction error based on fixed pamameters #
t <- proc.time()
for (i in 1:length(Sam_N_list_d1)){
  cat(i, "\n");
  pick_ind_ls = 1:Sam_N_list_d1[i]
  D_obs = D[pick_ind_ls, pick_ind_ls]
  D_obs_pred = D[pick_ind_ls, 12001:12003]
  D_pred_obs= D[12001:12003, pick_ind_ls]
  #(a)#
  lower_l1[i] = MSPE_theta0(N_obs = Sam_N_list_d1[i], N_pred = 3, phi, sigma2, 
                           tau2, D_obs, D_obs_pred, D_pred)
  cat(lower_l1[i], "\t")
  cat(proc.time() - t, "\n")
}

proc.time() - t

save(list = c("lower_l1", "Sam_N_list_d1"), 
     file = "./data/sim3_d1_a_results.RData", 
     envir = .GlobalEnv)



