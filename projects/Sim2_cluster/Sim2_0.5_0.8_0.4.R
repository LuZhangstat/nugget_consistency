setwd("/u/home/l/luzhangs/project-sudipto/nugget_consist")
rm(list = ls())
load("./data/simdata.RData")
source("./projects/utils.R")
# need to compare the MLE of parameter sigma^2 tau^2 phi

# First obtain the MLE of phi and delta = tau^2 / sigma^2 through gradient decent
# of the profile likelihood, then calculate the sigma^2
library("optimx")
# tt <- optimx(c(phi0.5_1, tau22), neg_loglik_0.5_profile,
#              X = Sim_data_0.5["0.2", "0.5_0.15", 1, 1:n3], D = D[1:n3, 1:n3],
#              n = n3, method = c('L-BFGS-B'), lower = c(0, 0))
# 54.01 seconds


MLE_0.5_0.8_0.4_n1 = matrix(NA, nrow = N_sim, ncol = 3)
MLE_0.5_0.8_0.4_n2 = matrix(NA, nrow = N_sim, ncol = 3)
MLE_0.5_0.8_0.4_n3 = matrix(NA, nrow = N_sim, ncol = 3)
colnames(MLE_0.5_0.8_0.4_n1) = c("tau2", "sigma2", "phi")
colnames(MLE_0.5_0.8_0.4_n2) = c("tau2", "sigma2", "phi")
colnames(MLE_0.5_0.8_0.4_n3) = c("tau2", "sigma2", "phi")

t <- proc.time()
for (i in 1:N_sim){
  X_data = Sim_data_0.5["0.8", "0.5_0.4", i, 1:n1]
  optL <- optimx(c(phi0.5_2, tau23), neg_loglik_0.5_profile, 
                 X = X_data, D = D[1:n1, 1:n1], 
                 n = n1, method = c('L-BFGS-B'), lower = 0)
  sigma2_MLE <- sigma2_MLE_cal(optL$p1, optL$p2,  X = X_data, 
                               D = D[1:n1, 1:n1], n = n1)
  MLE_0.5_0.8_0.4_n1[i, ]<- c(sigma2_MLE * optL$p2, sigma2_MLE, optL$p1)
}
proc.time() - t
# 5 for 8.769 seconds, estimated time: 1753.8

t <- proc.time()
for (i in 1:N_sim){
  X_data = Sim_data_0.5["0.8", "0.5_0.4", i, 1:n2]
  optL <- optimx(c(phi0.5_2, tau23), neg_loglik_0.5_profile, 
                 X = X_data, D = D[1:n2, 1:n2], 
                 n = n2, method = c('L-BFGS-B'), lower = 0)
  sigma2_MLE <- sigma2_MLE_cal(optL$p1, optL$p2,  X = X_data, 
                               D = D[1:n2, 1:n2], n = n2)
  MLE_0.5_0.8_0.4_n2[i, ]<- c(sigma2_MLE * optL$p2, sigma2_MLE, optL$p1)
}
proc.time() - t
# 5 for 68.624 seconds, estimated time: 13724.8

t <- proc.time()
for (i in 1:N_sim){
  X_data = Sim_data_0.5["0.8", "0.5_0.4", i, 1:n3]
  optL <- optimx(c(phi0.5_2, tau23), neg_loglik_0.5_profile, 
                 X = X_data, D = D[1:n3, 1:n3], 
                 n = n3, method = c('L-BFGS-B'), lower = 0)
  sigma2_MLE <- sigma2_MLE_cal(optL$p1, optL$p2,  X = X_data, 
                               D = D[1:n3, 1:n3], n = n3)
  MLE_0.5_0.8_0.4_n3[i, ]<- c(sigma2_MLE * optL$p2, sigma2_MLE, optL$p1)
}
proc.time() - t
# 5 for 384.184 seconds, estimated time: 76836.8


save(list = c("MLE_0.5_0.8_0.4_n1", "MLE_0.5_0.8_0.4_n2", "MLE_0.5_0.8_0.4_n3"), 
     file = "./data/sim2_0.5_0.8_0.4results.RData", 
     envir = .GlobalEnv)


# the total estimated time is around 25.6 hours
