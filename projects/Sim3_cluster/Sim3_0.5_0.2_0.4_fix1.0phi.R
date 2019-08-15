setwd("/u/home/l/luzhangs/project-sudipto/nugget_consist")
rm(list = ls())
load("./data/simdata.RData")
source("./projects/utils.R")
# need to compare the MLE of parameter sigma^2 tau^2 phi

# First obtain the MLE of phi and delta = tau^2 / sigma^2 through gradient decent
# of the profile likelihood, then calculate the sigma^2
library("optimx")
# tt <- optimx(c(phi0.5_1, tau22), neg_loglik_0.5_profile,
#              X = Sim_data_0.5["0.2", "0.5_0.15", 1, 1:n2], D = D[1:n2, 1:n2],
#              n = n2, method = c('L-BFGS-B'), lower = c(0, 0))
# # 7.784 seconds

Sam_N_list = c(200, 400, 600, 800, 1000, 1200, 1400, 1600)

phi_fix = phi0.5_2

MLE_0.5_0.2_0.4_fix1.0phi = array(NA, c(N_sim, 2, 8))

t <- proc.time()
for (j in 1:length(Sam_N_list)){
  cat(j, "\n")
  for (i in 1:N_sim){
    X_data = Sim_data_0.5["0.2", "0.5_0.4", i, 1:Sam_N_list[j]]
    optL <- optimx(tau22, neg_loglik_0.5_fixphi, 
                   phi_fix = phi_fix, X = X_data, 
                   D = D[1:Sam_N_list[j], 1:Sam_N_list[j]], 
                   n = Sam_N_list[j], method = c('L-BFGS-B'), lower = 0)
    sigma2_MLE <- sigma2_MLE_cal(phi_fix, optL$p1,  X = X_data, 
                                 D = D[1:Sam_N_list[j], 1:Sam_N_list[j]], 
                                 n = Sam_N_list[j])
    MLE_0.5_0.2_0.4_fix1.0phi[i, , j]<- c(sigma2_MLE * optL$p1, sigma2_MLE)
  }
  cat(proc.time() - t, "\n")
}
proc.time() - t
# 5 for 6.905 seconds. Estimated time: 1381

save(list = c("MLE_0.5_0.2_0.4_fix1.0phi"), 
     file = "./data/sim3_MLE_0.5_0.2_0.4_fix1.0phiresults.RData", 
     envir = .GlobalEnv)

# need more than one dau, around 28.4 hours, so maybe require 2 days

