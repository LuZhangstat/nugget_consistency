 #  simultation for prediction #
setwd("/Users/luzhang/Documents/github/nugget_consistency")
rm(list = ls())
load("./data/simdata.RData")
source("./projects/utils_2.R")

# select subset range from 200 to 1600 
set.seed(123)
Sam_N_list = c(200, 400, 600, 800, 1000, 1200, 1400, 1600)
#sample_ind = sample.int(1600, 1600, replace = FALSE)

# choose phi = phi0.5_2 sigma2 = sigma2, tau2 = tau22
phi = phi0.5_2; tau2 = tau22
D_pred = D[1601:4100, 1601:4100]

 ##########################
#    simultion 3 part A    #
 ##########################

# preallocation 
data_A_ratio1 = rep(0, length(Sam_N_list))
data_A_ratio2 = rep(0, length(Sam_N_list))
data_A_ratio3 = rep(0, length(Sam_N_list))
data_A_ratio4 = rep(0, length(Sam_N_list))
data_A_ratio5 = rep(0, length(Sam_N_list))
data_A_ratio6 = rep(0, length(Sam_N_list))

upper1_l = rep(0, length(Sam_N_list))
upper2_l = rep(0, length(Sam_N_list))
upper3_l = rep(0, length(Sam_N_list))
upper4_l = rep(0, length(Sam_N_list))
upper5_l = rep(0, length(Sam_N_list))
upper6_l = rep(0, length(Sam_N_list))
lower_l = rep(0, length(Sam_N_list))

# mean squared prediction error based on fixed pamameters #
t <- proc.time()
for (i in 1:length(Sam_N_list)){
  cat(i, "\n");
  D_obs = D[1:Sam_N_list[i], 1:Sam_N_list[i]]
  D_obs_pred = D[1:Sam_N_list[i], 1601:4100]
  D_pred_obs= D[1601:4100, 1:Sam_N_list[i]]
  #(a)#
  lower_l[i] = MSPE_theta0(N_obs = Sam_N_list[i], N_pred = 2500, phi, sigma2, 
                      tau2, D_obs, D_obs_pred, D_pred)
  # let phi1 = 2*phi 
  upper1_l[i] = MSPE_theta01(N_obs = Sam_N_list[i], N_pred = 2500, phi, sigma2, tau2, 
                        phi1 = phi, sigma21 = sigma2, tau21 = 0.5 * tau2, 
                        D_obs, D_obs_pred, D_pred_obs, D_pred)
  
  upper2_l[i] = MSPE_theta01(N_obs = Sam_N_list[i], N_pred = 2500, phi, sigma2, tau2, 
                        phi1 = phi, sigma21 = sigma2, tau21 = 2 * tau2, 
                        D_obs, D_obs_pred, D_pred_obs, D_pred)
  
  upper3_l[i] = MSPE_theta01(N_obs = Sam_N_list[i], N_pred = 2500, phi, sigma2, tau2, 
                        phi1 = 2 * phi, sigma21 = sigma2, tau21 = tau2, 
                        D_obs, D_obs_pred, D_pred_obs, D_pred)
  
  upper4_l[i] = MSPE_theta01(N_obs = Sam_N_list[i], N_pred = 2500, phi, sigma2, tau2, 
                        phi1 = phi, sigma21 = 2 * sigma2, tau21 = tau2, 
                        D_obs, D_obs_pred, D_pred_obs, D_pred)
  
  upper5_l[i] = MSPE_theta01(N_obs = Sam_N_list[i], N_pred = 2500, phi, sigma2, tau2, 
                        phi1 = 2.0 * phi, sigma21 = 0.5 * sigma2, tau21 = tau2, 
                        D_obs, D_obs_pred, D_pred_obs, D_pred)
  
  upper6_l[i] = MSPE_theta01(N_obs = Sam_N_list[i], N_pred = 2500, phi, sigma2, tau2, 
                        phi1 = 0.5 * phi, sigma21 = 2 * sigma2, tau21 = tau2, 
                        D_obs, D_obs_pred, D_pred_obs, D_pred)
  
  data_A_ratio1[i] = upper1_l[i] / lower_l[i]
  data_A_ratio2[i] = upper2_l[i] / lower_l[i]
  data_A_ratio3[i] = upper3_l[i] / lower_l[i]
  data_A_ratio4[i] = upper4_l[i] / lower_l[i]
  data_A_ratio5[i] = upper5_l[i] / lower_l[i]
  data_A_ratio6[i] = upper6_l[i] / lower_l[i]
  cat(proc.time() - t, "\n")
}
proc.time() - t
data_A_ratio1
data_A_ratio2
data_A_ratio3
data_A_ratio4
data_A_ratio5
data_A_ratio6
# plug-in predictions #

save(list = c("lower_l", "upper1_l", "upper2_l", "upper3_l", "upper4_l",
              "upper5_l", "upper6_l",
              "data_A_ratio1", "data_A_ratio2", "data_A_ratio3", 
              "data_A_ratio4", "data_A_ratio5", "data_A_ratio6"), 
     file = "./data/sim3_2_a_results.RData", 
     envir = .GlobalEnv)

Sam_N_list2 = c(2500, 3200, 4100);
lower_l2 = rep(0, 3);

for (i in 1:length(Sam_N_list2)){
  set.seed(123)
  cat(i, "\n");
  pick_ind_large = sample.int(4100, Sam_N_list2[i])
  D_obs = D[pick_ind_large, pick_ind_large]
  D_obs_pred = D[pick_ind_large, 1601:4100]
  D_pred_obs= D[1601:4100, pick_ind_large]
  #(a)#
  lower_l2[i] = MSPE_theta0(N_obs = Sam_N_list2[i], N_pred = 2500, phi, sigma2, 
                           tau2, D_obs, D_obs_pred, D_pred)
  cat(lower_l2[i], "\n")
  }

 ##########################
#    simultion 3 part B    #
 ##########################

# preallocation 
data_B_ratio1 = rep(0, length(Sam_N_list))
data_B_ratio2 = rep(0, length(Sam_N_list))
data_B_ratio3 = rep(0, length(Sam_N_list))
data_B_ratio4 = rep(0, length(Sam_N_list))
data_B_ratio5 = rep(0, length(Sam_N_list))
data_B_ratio6 = rep(0, length(Sam_N_list))

# mean squared prediction error based on fixed pamameters #
t <- proc.time()
for (i in 1:length(Sam_N_list)){
  cat(i, "\n");
  D_obs = D[1:Sam_N_list[i], 1:Sam_N_list[i]]
  D_obs_pred = D[1:Sam_N_list[i], 1601:4100]
  D_pred_obs= D[1601:4100, 1:Sam_N_list[i]]
  #(b)#
  
  # 1st phi1*sigma21 = phi*sigma
  upper1 = MSPE_theta0(N_obs = Sam_N_list[i], N_pred = 2500, phi = phi * 2,
                       sigma2 = sigma2 * 0.5, tau2, D_obs, D_obs_pred, D_pred)
  lower1 = MSPE_theta01(N_obs = Sam_N_list[i], N_pred = 2500, phi, sigma2, tau2,
                        phi1 = 2 * phi, sigma21 = 0.5 * sigma2, tau21 = tau2,
                        D_obs, D_obs_pred, D_pred_obs, D_pred)

  upper2 = MSPE_theta0(N_obs = Sam_N_list[i], N_pred = 2500, phi = phi * 0.5,
                       sigma2 = sigma2 * 2, tau2, D_obs, D_obs_pred, D_pred)
  lower2 = MSPE_theta01(N_obs = Sam_N_list[i], N_pred = 2500, phi, sigma2, tau2,
                        phi1 = 0.5 * phi, sigma21 = 2 * sigma2, tau21 = tau2,
                        D_obs, D_obs_pred, D_pred_obs, D_pred)
  # 2nd phi1*sigma21 > phi*sigma
  upper3 = MSPE_theta0(N_obs = Sam_N_list[i], N_pred = 2500, phi = 1.1 * phi,
                       sigma2 = 1.1 * sigma2, tau2, D_obs, D_obs_pred, D_pred)
  lower3 = MSPE_theta01(N_obs = Sam_N_list[i], N_pred = 2500, phi, sigma2, tau2,
                        phi1 = 1.1 * phi, sigma21 = 1.1 * sigma2, tau21 = tau2,
                        D_obs, D_obs_pred, D_pred_obs, D_pred)

  # 3rd phi1*sigma21 < phi*sigma
  upper4 = MSPE_theta0(N_obs = Sam_N_list[i], N_pred = 2500, phi = phi * 0.9,
                       sigma2 = 0.9 * sigma2, tau2, D_obs, D_obs_pred, D_pred)
  lower4 = MSPE_theta01(N_obs = Sam_N_list[i], N_pred = 2500, phi, sigma2, tau2,
                        phi1 = 0.9 * phi, sigma21 = 0.9 * sigma2, tau21 = tau2,
                        D_obs, D_obs_pred, D_pred_obs, D_pred)

  # 4nd check the impact of nugget
  upper5 = MSPE_theta0(N_obs = Sam_N_list[i], N_pred = 2500, phi = phi, 
                       sigma2 = sigma2, tau2 = 2.0 * tau2, D_obs, D_obs_pred, D_pred)
  lower5 = MSPE_theta01(N_obs = Sam_N_list[i], N_pred = 2500, phi, sigma2, tau2, 
                        phi1 = phi, sigma21 = sigma2, tau21 = 1.2 * tau2, 
                        D_obs, D_obs_pred, D_pred_obs, D_pred)
  
  upper6 = MSPE_theta0(N_obs = Sam_N_list[i], N_pred = 2500, phi = phi, 
                       sigma2 = sigma2, tau2 = 0.5 * tau2, D_obs, D_obs_pred, D_pred)
  lower6 = MSPE_theta01(N_obs = Sam_N_list[i], N_pred = 2500, phi, sigma2, tau2, 
                        phi1 = phi, sigma21 = sigma2, tau21 = 0.8 * tau2, 
                        D_obs, D_obs_pred, D_pred_obs, D_pred)
  
  data_B_ratio1[i] = upper1 / lower1
  data_B_ratio2[i] = upper2 / lower2
  data_B_ratio3[i] = upper3 / lower3
  data_B_ratio4[i] = upper4 / lower4
  data_B_ratio5[i] = upper5 / lower5
  data_B_ratio6[i] = upper6 / lower6
  cat(proc.time() - t, "\n")
}
proc.time() - t

data_B_ratio1
data_B_ratio2
data_B_ratio3
data_B_ratio4
data_B_ratio5
data_B_ratio6

save(list = c("data_B_ratio1", "data_B_ratio2", "data_B_ratio3", 
              "data_B_ratio4", "data_B_ratio5", "data_B_ratio6"), 
     file = "./data/sim3_2_b_results.RData", 
     envir = .GlobalEnv)



 ##########################
#    simultion 3 part C    #
 ##########################
load("./data/sim3_MLE_0.5_0.2_0.4_fix0.8phiresults.RData")
load("./data/sim3_MLE_0.5_0.2_0.4_fix1.0phiresults.RData")
load("./data/sim3_MLE_0.5_0.2_0.4_fix1.2phiresults.RData")

Mean_MLE_0.5_0.2_0.4_fix0.8phi <- matrix(NA, 8, 2)
Mean_MLE_0.5_0.2_0.4_fix1.0phi <- matrix(NA, 8, 2)
Mean_MLE_0.5_0.2_0.4_fix1.2phi <- matrix(NA, 8, 2)

for (i in 1:8){
  for (j in 1:2){
    Mean_MLE_0.5_0.2_0.4_fix0.8phi[i, j] <- 
      mean(MLE_0.5_0.2_0.4_fix0.8phi[, j, i])
    Mean_MLE_0.5_0.2_0.4_fix1.0phi[i, j] <- 
      mean(MLE_0.5_0.2_0.4_fix1.0phi[, j, i])
    Mean_MLE_0.5_0.2_0.4_fix1.2phi[i, j] <- 
      mean(MLE_0.5_0.2_0.4_fix1.2phi[, j, i])
  }
}

data_C_ratio1 = rep(0, length(Sam_N_list))
data_C_ratio2 = rep(0, length(Sam_N_list))
data_C_ratio3 = rep(0, length(Sam_N_list))

# mean squared prediction error based on fixed pamameters #
t <- proc.time()
for (i in 1:length(Sam_N_list)){
  cat(i, "\n");
  D_obs = D[1:Sam_N_list[i], 1:Sam_N_list[i]]
  D_obs_pred = D[1:Sam_N_list[i], 1601:4100]
  D_pred_obs= D[1601:4100, 1:Sam_N_list[i]]
  #(C)#
  upper1 = MSPE_theta0(N_obs = Sam_N_list[i], N_pred = 2500, phi = 0.8 * phi, 
                       sigma2 = Mean_MLE_0.5_0.2_0.4_fix0.8phi[i, 2], 
                       tau2 = Mean_MLE_0.5_0.2_0.4_fix0.8phi[i, 1], 
                       D_obs, D_obs_pred, D_pred)
  lower1 = MSPE_theta01(N_obs = Sam_N_list[i], N_pred = 2500, phi, sigma2, tau2, 
                        phi1 = 0.8 * phi, 
                        sigma21 = Mean_MLE_0.5_0.2_0.4_fix0.8phi[i, 2], 
                        tau21 = Mean_MLE_0.5_0.2_0.4_fix0.8phi[i, 1], 
                        D_obs, D_obs_pred, D_pred_obs, D_pred)
  
  upper2 = MSPE_theta0(N_obs = Sam_N_list[i], N_pred = 2500, phi = phi, 
                       sigma2 = Mean_MLE_0.5_0.2_0.4_fix1.0phi[i, 2], 
                       tau2 = Mean_MLE_0.5_0.2_0.4_fix1.0phi[i, 1], 
                       D_obs, D_obs_pred, D_pred)
  lower2 = MSPE_theta01(N_obs = Sam_N_list[i], N_pred = 2500, phi, sigma2, tau2, 
                        phi1 = phi, 
                        sigma21 = Mean_MLE_0.5_0.2_0.4_fix1.0phi[i, 2], 
                        tau21 = Mean_MLE_0.5_0.2_0.4_fix1.0phi[i, 1], 
                        D_obs, D_obs_pred, D_pred_obs, D_pred)

  upper3 = MSPE_theta0(N_obs = Sam_N_list[i], N_pred = 2500, phi = 1.2 * phi, 
                       sigma2 = Mean_MLE_0.5_0.2_0.4_fix1.2phi[i, 2], 
                       tau2 = Mean_MLE_0.5_0.2_0.4_fix1.2phi[i, 1], 
                       D_obs, D_obs_pred, D_pred)
  lower3 = MSPE_theta01(N_obs = Sam_N_list[i], N_pred = 2500, phi, sigma2, tau2, 
                        phi1 = 1.2 * phi, 
                        sigma21 = Mean_MLE_0.5_0.2_0.4_fix1.2phi[i, 2], 
                        tau21 = Mean_MLE_0.5_0.2_0.4_fix1.2phi[i, 1], 
                        D_obs, D_obs_pred, D_pred_obs, D_pred)
  

  data_C_ratio1[i] = upper1 / lower1
  data_C_ratio2[i] = upper2 / lower2
  data_C_ratio3[i] = upper3 / lower3
  cat(proc.time() - t, "\n")
}
proc.time() - t

data_C_ratio1
data_C_ratio2
data_C_ratio3

save(list = c("data_C_ratio1", "data_C_ratio2", "data_C_ratio3"), 
     file = "./data/sim3_2_c_results.RData", 
     envir = .GlobalEnv)

##
# The summary plots 
##

width <- 8.0
height <- 3.0
pointsize <- 16


load("./data/sim3_2_a_results.RData")
load("./data/sim3_2_b_results.RData")
load("./data/sim3_2_c_results.RData")
library(ggplot2)

dfa1 <- data.frame(
  sample_size = rep(Sam_N_list[-1], 6),
  MSPE = c(lower_l[-1], upper2_l[-1], upper3_l[-1], 
           upper4_l[-1], upper5_l[-1], upper6_l[-1]),
  group = rep(c("tau^2, sigma^2, phi", 
                "2tau^2, sigma^2, phi", "tau^2, sigma^2, 2phi", 
                "tau^2, 2sigma^2, phi", "tau^2, 0.5sigma^2, 2phi", 
                "tau^2, 2sigma^2, 0.5phi"), each = length(Sam_N_list[-1])))

setEPS()
postscript("./pics/sim3a2_consist.eps", width = width, height = height)
par(mfrow = c(1, 1))
ggplot(data = dfa1, aes(x = sample_size, y = MSPE, group = group)) +
  geom_line(aes(linetype = group))+
  geom_point(aes(shape = group)) + 
  scale_shape_discrete(
    name  ="",
    breaks=c("tau^2, sigma^2, phi", "2tau^2, sigma^2, phi", 
             "tau^2, sigma^2, 2phi", "tau^2, 2sigma^2, phi", 
             "tau^2, 0.5sigma^2, 2phi", "tau^2, 2sigma^2, 0.5phi"),
    labels=c(expression(tau[0]^2 ~ "," ~ sigma[0]^2 ~ "," ~ phi[0]),
             expression(2*tau[0]^2 ~ "," ~ sigma[0]^2 ~ "," ~ phi[0]),
             expression(tau[0]^2 ~ "," ~ sigma[0]^2 ~ "," ~ 2*phi[0]), 
             expression(tau[0]^2 ~ "," ~2*sigma[0]^2 ~ "," ~ phi[0]),
             expression(tau[0]^2 ~ "," ~0.5*sigma[0]^2 ~ "," ~2*phi[0]), 
             expression(tau[0]^2 ~ "," ~2*sigma[0]^2 ~ "," ~0.5*phi[0]))) +
  scale_linetype_discrete(
    name  ="",
    breaks=c("tau^2, sigma^2, phi", "2tau^2, sigma^2, phi", 
             "tau^2, sigma^2, 2phi", "tau^2, 2sigma^2, phi", 
             "tau^2, 0.5sigma^2, 2phi", "tau^2, 2sigma^2, 0.5phi"),
    labels=c(expression(tau[0]^2 ~ "," ~ sigma[0]^2 ~ "," ~ phi[0]),
             expression(2*tau[0]^2 ~ "," ~ sigma[0]^2 ~ "," ~ phi[0]),
             expression(tau[0]^2 ~ "," ~ sigma[0]^2 ~ "," ~ 2*phi[0]), 
             expression(tau[0]^2 ~ "," ~2*sigma[0]^2 ~ "," ~ phi[0]),
             expression(tau[0]^2 ~ "," ~0.5*sigma[0]^2 ~ "," ~2*phi[0]), 
             expression(tau[0]^2 ~ "," ~2*sigma[0]^2 ~ "," ~0.5*phi[0]))) +
  scale_x_continuous(name = "sample size") + 
  theme_bw()
dev.off()

dfa <- data.frame(
  sample_size = rep(Sam_N_list[-1], 5),
  ratio = c(data_A_ratio2[-1], data_A_ratio3[-1], 
            data_A_ratio4[-1], data_A_ratio5[-1], data_A_ratio6[-1]),
  group = rep(c("2tau^2, sigma^2, phi", 
                "tau^2, sigma^2, 2phi", "tau^2, 2sigma^2, phi", 
                "tau^2, 0.5sigma^2, 2phi", "tau^2, 2sigma^2, 0.5phi"), 
              each = length(Sam_N_list[-1])))

setEPS()
postscript("./pics/sim3a2_ratio.eps", width = width, height = height)
par(mfrow = c(1, 1))
ggplot(data = dfa, aes(x = sample_size, y = ratio, group = group)) +
  geom_line(aes(linetype = group))+
  geom_point(aes(shape = group)) + 
  scale_shape_discrete(
    name  ="",
    breaks=c("2tau^2, sigma^2, phi", 
             "tau^2, sigma^2, 2phi", "tau^2, 2sigma^2, phi", 
             "tau^2, 0.5sigma^2, 2phi", "tau^2, 2sigma^2, 0.5phi"),
    labels=c(expression(2*tau[0]^2 ~ "," ~ sigma[0]^2 ~ "," ~ phi[0]),
             expression(tau[0]^2 ~ "," ~ sigma[0]^2 ~ "," ~ 2*phi[0]), 
             expression(tau[0]^2 ~ "," ~2*sigma[0]^2 ~ "," ~ phi[0]),
             expression(tau[0]^2 ~ "," ~0.5*sigma[0]^2 ~ "," ~2*phi[0]), 
             expression(tau[0]^2 ~ "," ~2*sigma[0]^2 ~ "," ~0.5*phi[0]))) +
  scale_linetype_discrete(
    name  ="",
    breaks=c("2tau^2, sigma^2, phi", 
             "tau^2, sigma^2, 2phi", "tau^2, 2sigma^2, phi", 
             "tau^2, 0.5sigma^2, 2phi", "tau^2, 2sigma^2, 0.5phi"),
    labels=c(expression(2*tau[0]^2 ~ "," ~ sigma[0]^2 ~ "," ~ phi[0]),
             expression(tau[0]^2 ~ "," ~ sigma[0]^2 ~ "," ~ 2*phi[0]), 
             expression(tau[0]^2 ~ "," ~2*sigma[0]^2 ~ "," ~ phi[0]),
             expression(tau[0]^2 ~ "," ~0.5*sigma[0]^2 ~ "," ~2*phi[0]), 
             expression(tau[0]^2 ~ "," ~2*sigma[0]^2 ~ "," ~0.5*phi[0]))) +
  scale_x_continuous(name = "sample size") + 
  theme_bw()
dev.off()



dfb <- data.frame(sample_size = rep(Sam_N_list[-1], 6),
  ratio = c(data_B_ratio1[-1], data_B_ratio2[-1], data_B_ratio3[-1], 
            data_B_ratio4[-1], data_B_ratio5[-1], data_B_ratio6[-1]),
                 group = rep(c("tau, 2phi, 0.5sigma^2", "tau, 0.5phi, 2sigma^2", 
                               "tau, 1.1phi, 1.1sigma^2", "tau, 0.9phi, 0.9sigma^2", 
                               "2tau, phi, sigma^2", "0.5tau, phi, sigma^2"), 
                             each = length(Sam_N_list[-1])))

setEPS()
postscript("./pics/sim3b2_ratio.eps", width = width, height = height)
par(mfrow = c(1, 1))
ggplot(data = dfb, aes(x = sample_size, y = ratio, group = group)) +
  geom_line(aes(linetype = group))+
  geom_point(aes(shape = group)) + 
  scale_shape_discrete(
    name  ="",
    breaks=c("tau, 2phi, 0.5sigma^2", "tau, 0.5phi, 2sigma^2", 
             "tau, 1.1phi, 1.1sigma^2", "tau, 0.9phi, 0.9sigma^2", 
             "2tau, phi, sigma^2", "0.5tau, phi, sigma^2"),
    labels=c(expression(tau[0]^2 ~ "," ~ 0.5*sigma[0]^2 ~ "," ~ 2*phi[0]),
             expression(tau[0]^2 ~ "," ~ 2*sigma[0]^2 ~ "," ~ 0.5*phi[0]), 
             expression(tau[0]^2 ~ "," ~1.1*sigma[0]^2 ~ "," ~ 1.1*phi[0]),
             expression(tau[0]^2 ~ "," ~0.9*sigma[0]^2 ~ "," ~0.9*phi[0]), 
             expression(2*tau[0]^2 ~ "," ~sigma[0]^2 ~ "," ~phi[0]),
             expression(0.5*tau[0]^2 ~ "," ~sigma[0]^2 ~ "," ~phi[0]))) +
  scale_linetype_discrete(
    name  ="",
    breaks=c("tau, 2phi, 0.5sigma^2", "tau, 0.5phi, 2sigma^2", 
             "tau, 1.1phi, 1.1sigma^2", "tau, 0.9phi, 0.9sigma^2", 
             "2tau, phi, sigma^2", "0.5tau, phi, sigma^2"),
    labels=c(expression(tau[0]^2 ~ "," ~ 0.5*sigma[0]^2 ~ "," ~ 2*phi[0]),
             expression(tau[0]^2 ~ "," ~ 2*sigma[0]^2 ~ "," ~ 0.5*phi[0]), 
             expression(tau[0]^2 ~ "," ~1.1*sigma[0]^2 ~ "," ~ 1.1*phi[0]),
             expression(tau[0]^2 ~ "," ~0.9*sigma[0]^2 ~ "," ~0.9*phi[0]), 
             expression(2*tau[0]^2 ~ "," ~sigma[0]^2 ~ "," ~phi[0]),
             expression(0.5*tau[0]^2 ~ "," ~sigma[0]^2 ~ "," ~phi[0]))) +
  scale_x_continuous(name = "sample size") + 
  theme_bw()
dev.off()



dfc <- data.frame(sample_size = rep(Sam_N_list, 3),
                  ratio = c(data_C_ratio1, data_C_ratio2, data_C_ratio3),
                  group = rep(c("0.8phi", "phi", "1.2phi"), 
                              each = length(Sam_N_list)))


setEPS()
postscript("./pics/sim3c2_ratio.eps", width = width, height = height)
par(mfrow = c(1, 1))
ggplot(data = dfc, aes(x = sample_size, y = ratio, group = group)) +
  geom_line(aes(linetype = group))+
  geom_point(aes(shape = group)) + 
  scale_shape_discrete(
    name  ="",
    breaks=c("0.8phi", "phi", "1.2phi"),
    labels=c(expression(0.8*phi[0]), expression(phi[0]), 
             expression(1.2*phi[0]))) +
  scale_linetype_discrete(
    name  ="",
    breaks=c("0.8phi", "phi", "1.2phi"),
    labels=c(expression(0.8*phi[0]), expression(phi[0]), 
             expression(1.2*phi[0]))) +
  scale_x_continuous(name = "sample size") + 
  theme_bw()
dev.off()




png(paste("./pics/sim3c2_ratio.png", sep=""), 
    width = width, height = height, pointsize = pointsize, family="Courier")
ggplot(data = dfc, aes(x = sample_size, y = ratio, group = group)) +
  geom_line(aes(linetype = group, color = group))+
  geom_point(aes(shape = group, color = group))
dev.off()


#########
# extra #
#########

x = grid_edge_pred
pred_ind = (1600 + 25*50 + 1):(1600 + 26*50)
y = Sim_data_0.5["0.2", "0.5_0.4", 1, pred_ind]
tau2 = tau22
plot(x, y_pred1, type = "l", ylim = c(-0.8, 0.8))
#points(x, y_pred1, col = "green")
points(x, y_pred2, col = "blue")
points(x, y_pred3, col = "blue")
points(x, y_pred4, col = "red")
points(x, y_pred5, col = "red")
points(x, y_pred6, col = "purple")
points(x, y_pred7, col = "purple")

## use true value ##
y_pred1 = sigma2 * exp(-phi * D[pred_ind, 1:1600]) %*% 
  chol2inv(chol(sigma2 * exp(-phi * D[1:1600, 1:1600]) + 
                  tau2 * diag(1:1600))) %*% 
  Sim_data_0.5["0.2", "0.5_0.4", 1, 1:1600]

## keep phi * sigma2 and tau2 ##
# phi = 2 * phi, sigma2 = 0.5 * sigma2, tau2 = tau2
y_pred2 = 0.5 * sigma2 * exp(- 2 * phi * D[pred_ind, 1:1600]) %*% 
  chol2inv(chol(0.5 * sigma2 * exp(- 2 * phi * D[1:1600, 1:1600]) + 
                  tau2 * diag(1:1600))) %*% 
  Sim_data_0.5["0.2", "0.5_0.4", 1, 1:1600]

# phi = 0.5 * phi, sigma2 = 2 * sigma2, tau2 = tau2
y_pred3 = 2 * sigma2 * exp(- 0.5 * phi * D[pred_ind, 1:1600]) %*% 
  chol2inv(chol(2 * sigma2 * exp(- 0.5 * phi * D[1:1600, 1:1600]) + 
                  tau2 * diag(1:1600))) %*% 
  Sim_data_0.5["0.2", "0.5_0.4", 1, 1:1600]


## keep tau2, double phi * sigma2 ##
# phi = phi, sigma2 = 2 * sigma2, tau2 = tau2
y_pred4 = 2 * sigma2 * exp(- phi * D[pred_ind, 1:1600]) %*% 
  chol2inv(chol(2 * sigma2 * exp(- phi * D[1:1600, 1:1600]) + 
                  tau2 * diag(1:1600))) %*% 
  Sim_data_0.5["0.2", "0.5_0.4", 1, 1:1600]

# phi = 2 * phi, sigma2 = sigma2, tau2 = tau2
y_pred5 = sigma2 * exp(- 2 * phi * D[pred_ind, 1:1600]) %*% 
  chol2inv(chol(sigma2 * exp(- 2 * phi * D[1:1600, 1:1600]) + 
                  tau2 * diag(1:1600))) %*% 
  Sim_data_0.5["0.2", "0.5_0.4", 1, 1:1600]

## double and half tau2, keep phi * sigma2 ##
# phi = phi, sigma2 = sigma2, tau2 = 2*tau2
y_pred6 = sigma2 * exp(- phi * D[pred_ind, 1:1600]) %*% 
  chol2inv(chol(sigma2 * exp(- phi * D[1:1600, 1:1600]) + 
                  2 * tau2 * diag(1:1600))) %*% 
  Sim_data_0.5["0.2", "0.5_0.4", 1, 1:1600]

# phi = phi, sigma2 = sigma2, tau2 = 0.5*tau2
y_pred7 = sigma2 * exp(- phi * D[pred_ind, 1:1600]) %*% 
  chol2inv(chol(sigma2 * exp(- phi * D[1:1600, 1:1600]) + 
                  0.5 * tau2 * diag(1:1600))) %*% 
  Sim_data_0.5["0.2", "0.5_0.4", 1, 1:1600]

tau2 = 0.0
x = grid_edge_pred
pred_ind = (1600 + 25*50 + 1):(1600 + 26*50)
y_tau0 = Sim_data_0.5["0", "0.5_0.4", 1, pred_ind]

plot(x, y_pred1_tau0, type = "l")
#points(x, y_pred1, col = "green")
points(x, y_pred2_tau0, col = "blue")
points(x, y_pred3_tau0, col = "blue")
points(x, y_pred4_tau0, col = "red")
points(x, y_pred5_tau0, col = "red")
points(x, y_pred6_tau0, col = "purple")
points(x, y_pred7_tau0, col = "purple")
points(x, y_tau0)

## use true value ##
y_pred1_tau0 = sigma2 * exp(-phi * D[pred_ind, 1:1600]) %*% 
  chol2inv(chol(sigma2 * exp(-phi * D[1:1600, 1:1600]) + 
                  tau2 * diag(1:1600))) %*% 
  Sim_data_0.5["0", "0.5_0.4", 1, 1:1600]

## keep phi * sigma2 and tau2 ##
# phi = 2 * phi, sigma2 = 0.5 * sigma2, tau2 = tau2
y_pred2_tau0 = 0.5 * sigma2 * exp(- 2 * phi * D[pred_ind, 1:1600]) %*% 
  chol2inv(chol(0.5 * sigma2 * exp(- 2 * phi * D[1:1600, 1:1600]) + 
                  tau2 * diag(1:1600))) %*% 
  Sim_data_0.5["0", "0.5_0.4", 1, 1:1600]

# phi = 0.5 * phi, sigma2 = 2 * sigma2, tau2 = tau2
y_pred3_tau0 = 2 * sigma2 * exp(- 0.5 * phi * D[pred_ind, 1:1600]) %*% 
  chol2inv(chol(2 * sigma2 * exp(- 0.5 * phi * D[1:1600, 1:1600]) + 
                  tau2 * diag(1:1600))) %*% 
  Sim_data_0.5["0", "0.5_0.4", 1, 1:1600]


## keep tau2, double phi * sigma2 ##
# phi = phi, sigma2 = 2 * sigma2, tau2 = tau2
y_pred4_tau0 = 2 * sigma2 * exp(- phi * D[pred_ind, 1:1600]) %*% 
  chol2inv(chol(2 * sigma2 * exp(- phi * D[1:1600, 1:1600]) + 
                  tau2 * diag(1:1600))) %*% 
  Sim_data_0.5["0", "0.5_0.4", 1, 1:1600]

# phi = 2 * phi, sigma2 = sigma2, tau2 = tau2
y_pred5_tau0 = sigma2 * exp(- 2 * phi * D[pred_ind, 1:1600]) %*% 
  chol2inv(chol(sigma2 * exp(- 2 * phi * D[1:1600, 1:1600]) + 
                  tau2 * diag(1:1600))) %*% 
  Sim_data_0.5["0", "0.5_0.4", 1, 1:1600]

## double and half tau2, keep phi * sigma2 ##
# phi = phi, sigma2 = sigma2, tau2 = 2*tau2
y_pred6_tau0 = sigma2 * exp(- phi * D[pred_ind, 1:1600]) %*% 
  chol2inv(chol(sigma2 * exp(- phi * D[1:1600, 1:1600]) + 
                  2 * tau2 * diag(1:1600))) %*% 
  Sim_data_0.5["0", "0.5_0.4", 1, 1:1600]

# phi = phi, sigma2 = sigma2, tau2 = 0.5*tau2
y_pred7_tau0 = sigma2 * exp(- phi * D[pred_ind, 1:1600]) %*% 
  chol2inv(chol(sigma2 * exp(- phi * D[1:1600, 1:1600]) + 
                  0.5 * tau2 * diag(1:1600))) %*% 
  Sim_data_0.5["0", "0.5_0.4", 1, 1:1600]








plot(x, y_pred1, type = "l", ylim = c(-0.8, 0.8))
points(x, y_pred2_delta, col = "green")
points(x, y_pred3_delta, col = "blue")
points(x, y_pred4_delta, col = "blue")
points(x, y_pred5_delta, col = "red")
points(x, y_pred6_delta, col = "red")
#points(x, y_pred7, col = "purple")

# fix phi and delta = tau2 / sigma2
y_pred2_delta <- 2 * sigma2 * exp(-phi * D[pred_ind, 1:1600]) %*% 
  chol2inv(chol(2 * sigma2 * exp(-phi * D[1:1600, 1:1600]) + 
                  2 * tau2 * diag(1:1600))) %*% 
  Sim_data_0.5["0.2", "0.5_0.4", 1, 1:1600]

# fix phi*sigma2 and delta = tau2 / sigma2
y_pred3_delta <- 2 * sigma2 * exp(- 0.5 * phi * D[pred_ind, 1:1600]) %*% 
  chol2inv(chol(2 * sigma2 * exp(- 0.5 * phi * D[1:1600, 1:1600]) + 
                  2 * tau2 * diag(1:1600))) %*% 
  Sim_data_0.5["0.2", "0.5_0.4", 1, 1:1600]

y_pred4_delta <- 0.5 * sigma2 * exp(- 2 *phi * D[pred_ind, 1:1600]) %*% 
  chol2inv(chol(0.5 * sigma2 * exp(-2 * phi * D[1:1600, 1:1600]) + 
                  0.5 * tau2 * diag(1:1600))) %*% 
  Sim_data_0.5["0.2", "0.5_0.4", 1, 1:1600]

# fix tau2 and phi*sigma2
y_pred5_delta <- 0.5 * sigma2 * exp(- 2 * phi * D[pred_ind, 1:1600]) %*% 
  chol2inv(chol(0.5 * sigma2 * exp(- 2 * phi * D[1:1600, 1:1600]) + 
                  tau2 * diag(1:1600))) %*% 
  Sim_data_0.5["0.2", "0.5_0.4", 1, 1:1600]


y_pred6_delta <- 2.0 * sigma2 * exp(- 0.5 * phi * D[pred_ind, 1:1600]) %*% 
  chol2inv(chol(2.0 * sigma2 * exp(- 0.5 * phi * D[1:1600, 1:1600]) + 
                  tau2 * diag(1:1600))) %*% 
  Sim_data_0.5["0.2", "0.5_0.4", 1, 1:1600]


###########
# extra 2 #
###########
library(fields)
library(MASS)
sigma2 = 1; tau2 = tau22;
phi = phi0.5_2

set.seed(1234)
pred.x = seq(from = 0.2, to = 0.4, length.out = 68);
pred_loc2 = cbind(rep(0.4, 68), pred.x)
coords_total2 = rbind(obs_loc3, pred_loc2)
D2 = rdist(coords_total2)
COV2 = sigma2 * exp(- phi * D2) + tau2 * diag(nrow(coords_total2))
pred_data <- mvrnorm(n = 1, mu = rep(0, nrow(coords_total2)), Sigma = COV2)

plot(pred.x, pred_data[(n3+1):(n3 + 68)], type = "l")
points(pred.x, y_pred1_test, col = "green")
points(pred.x, y_pred2_test, col = "blue")
points(pred.x, y_pred3_test, col = "blue")
points(pred.x, y_pred4_test, col = "red")
points(pred.x, y_pred5_test, col = "red")
## use true value ##
y_pred1_test = sigma2 * exp(-phi * D2[1601:1668, 1:1600]) %*% 
  chol2inv(chol(sigma2 * exp(-phi * D2[1:1600, 1:1600]) + 
                  tau2 * diag(1:1600))) %*% pred_data[1:1600]

## keep phi * sigma2 and tau2 ##
# phi = 2 * phi, sigma2 = 0.5 * sigma2, tau2 = tau2
y_pred2_test = 0.5 * sigma2 * exp(- 2 * phi * D2[1601:1668, 1:1600]) %*% 
  chol2inv(chol(0.5 * sigma2 * exp(- 2 * phi * D2[1:1600, 1:1600]) + 
                  tau2 * diag(1:1600))) %*%  pred_data[1:1600]
y_pred3_test = 2 * sigma2 * exp(- 0.5 * phi * D2[1601:1668, 1:1600]) %*% 
  chol2inv(chol(2 * sigma2 * exp(- 0.5 * phi * D2[1:1600, 1:1600]) + 
                  tau2 * diag(1:1600))) %*%  pred_data[1:1600]


## keep tau2, double phi * sigma2 ##
# phi = phi, sigma2 = 2 * sigma2, tau2 = tau2
y_pred4_test = 2 * sigma2 * exp(- phi * D2[1601:1668, 1:1600]) %*% 
  chol2inv(chol(2 * sigma2 * exp(- phi * D2[1:1600, 1:1600]) + 
                  tau2 * diag(1:1600))) %*%  pred_data[1:1600]
y_pred5_test = sigma2 * exp(- 2 * phi * D2[1601:1668, 1:1600]) %*% 
  chol2inv(chol(sigma2 * exp(- 2 * phi * D2[1:1600, 1:1600]) + 
                  tau2 * diag(1:1600))) %*%  pred_data[1:1600]


