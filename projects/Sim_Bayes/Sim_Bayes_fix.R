rm(list = ls())
library(spBayes)
library(classInt)
library(RColorBrewer)
library(MBA)
library(fields)
library(geoR)
library(bayesplot)

# read in data #
load("./data/sim_bayes_data.RData")

x.res <- 100
y.res <- 100
col.br <- colorRampPalette(c("blue", "cyan", "yellow", "red"))
surf <- mba.surf(cbind(coords_total, Sim_Bayes_y), no.X = x.res, no.Y=y.res, h = 5, m = 2,
                 extend = FALSE)$xyz.est
image.plot(surf, xaxs = "r", yaxs = "r", xlab = "Easting (m)",
           ylab = "Northing (m)", col = col.br(25))
contour(surf, add = T)

## fit model with cmdstanr, fix phi at a fixed point ##
library(cmdstanr)
N_warmup = 500
N_sample = 500
N_chain = 4
fit_flag_fix = FALSE

# preallocate fitting results
fit_results_fixphi <- array(data = NA, dim = c(5, 3, 3, N_sample * N_chain), 
                            dimnames = list(
                              label = c("phi", "0.2phi", "0.5phi", "2phi",
                                        "5phi"),
                              pars = c("sigmasq", "tausq", "kappa"),
                              No_location = c("400", "900", "1600"),
                              index = 1:(N_sample * N_chain)))


file_fix <- file.path(getwd(), "projects", "Sim_Bayes", "GP_fix.stan")
mod_fix <- cmdstan_model(file_fix)

## n = 400, phi = phi0 ##
data1 <- list(N = n1,
              D = as.matrix(dist(obs_loc[1:n1, ])), # distance matrix
              coords = obs_loc[1:n1, ],
              y = Sim_Bayes_y[1:n1],
              phi = phi, 
              sb = 1.0,
              tb = 0.5)

if(fit_flag_fix){
  fit1 <- mod_fix$sample(
    data = data1,
    seed = 123,
    init = 2,
    chains = N_chain,
    parallel_chains = N_chain,
    save_warmup = TRUE,
    iter_warmup = N_warmup,
    iter_sampling = N_sample,
    refresh = 50,
    sig_figs = 16
  )
  
  fit1$cmdstan_diagnose()
  fit1$print()
  fit1$save_object(file = "./results/sim_bayes_phi1.RDS")
} else{
  fit1 <- readRDS(file = "./results/sim_bayes_phi1.RDS")
}

fit_results_fixphi["phi", "sigmasq", "400", ] <- 
  fit1$draws("sigmasq", inc_warmup = FALSE)
fit_results_fixphi["phi", "tausq", "400", ] <- 
  fit1$draws("tausq", inc_warmup = FALSE)
fit_results_fixphi["phi", "kappa", "400", ] <- 
  fit1$draws("kappa", inc_warmup = FALSE)


## n = 900, phi = phi0 ##
data2 <- list(N = n2,
              D = as.matrix(dist(obs_loc[1:n2, ])), # distance matrix
              coords = obs_loc[1:n2, ],
              y = Sim_Bayes_y[1:n2],
              phi = phi, 
              sb = 1.0,
              tb = 0.5)

if(fit_flag_fix){
  fit2 <- mod_fix$sample(
    data = data2,
    seed = 123,
    init = 2,
    chains = N_chain,
    parallel_chains = N_chain,
    save_warmup = TRUE,
    iter_warmup = N_warmup,
    iter_sampling = N_sample,
    refresh = 50,
    sig_figs = 16
  )
  
  fit2$cmdstan_diagnose()
  fit2$print()
  fit2$save_object(file = "./results/sim_bayes_phi2.RDS")
} else{
  fit2 <- readRDS(file = "./results/sim_bayes_phi2.RDS")
}

p1_2 <- mcmc_trace(fit2$draws("lp__", inc_warmup = TRUE)[-(1:20), ,])#[60:80, ,]) #c("lp__", "phi[1]", "lambda[1]", "theta1[1]")
p1_2
mcmc_hist(fit2$draws("tausq", inc_warmup = FALSE))
mcmc_hist(fit2$draws("sigmasq", inc_warmup = FALSE))
mcmc_hist(fit2$draws("kappa", inc_warmup = FALSE))

fit_results_fixphi["phi", "sigmasq", "900", ] <- 
  fit2$draws("sigmasq", inc_warmup = FALSE)
fit_results_fixphi["phi", "tausq", "900", ] <- 
  fit2$draws("tausq", inc_warmup = FALSE)
fit_results_fixphi["phi", "kappa", "900", ] <- 
  fit2$draws("kappa", inc_warmup = FALSE)


## n = 1600, phi = phi0 ##
data3 <- list(N = n3,
              D = as.matrix(dist(obs_loc[1:n3, ])), # distance matrix
              coords = obs_loc[1:n3, ],
              y = Sim_Bayes_y[1:n3],
              phi = phi, 
              sb = 1.0,
              tb = 0.5)

if(fit_flag_fix){
  fit3 <- mod_fix$sample(
    data = data3,
    seed = 123,
    init = 2,
    chains = N_chain,
    parallel_chains = N_chain,
    save_warmup = TRUE,
    iter_warmup = N_warmup,
    iter_sampling = N_sample,
    refresh = 50,
    sig_figs = 16
  )
  
  fit3$cmdstan_diagnose()
  fit3$print()
  fit3$save_object(file = "./results/sim_bayes_phi3.RDS")
}else{
  fit3 <- readRDS(file = "./results/sim_bayes_phi3.RDS")
}  

p1_3 <- mcmc_trace(fit3$draws("lp__", inc_warmup = TRUE)[-(1:20), ,])#[60:80, ,]) #c("lp__", "phi[1]", "lambda[1]", "theta1[1]")
p1_3
mcmc_hist(fit3$draws("tausq", inc_warmup = FALSE))
mcmc_hist(fit3$draws("sigmasq", inc_warmup = FALSE))
mcmc_hist(fit3$draws("kappa", inc_warmup = FALSE))

fit_results_fixphi["phi", "sigmasq", "1600", ] <- 
  fit3$draws("sigmasq", inc_warmup = FALSE)
fit_results_fixphi["phi", "tausq", "1600", ] <- 
  fit3$draws("tausq", inc_warmup = FALSE)
fit_results_fixphi["phi", "kappa", "1600", ] <- 
  fit3$draws("kappa", inc_warmup = FALSE)

## fix phi at 0.2phi ##
data1$phi = 0.2 * phi

if(fit_flag_fix){
  fit1_02 <- mod_fix$sample(
    data = data1,
    seed = 123,
    init = 2,
    chains = N_chain,
    parallel_chains = N_chain,
    save_warmup = TRUE,
    iter_warmup = N_warmup,
    iter_sampling = N_sample,
    refresh = 50,
    sig_figs = 16)
  
  fit1_02$print()
  fit1_02$save_object(file = "./results/sim_bayes_02phi1.RDS")
}else{
  fit1_02 <- readRDS(file = "./results/sim_bayes_02phi1.RDS")
}


fit_results_fixphi["0.2phi", "sigmasq", "400", ] <- 
  fit1_02$draws("sigmasq", inc_warmup = FALSE)
fit_results_fixphi["0.2phi", "tausq", "400", ] <- 
  fit1_02$draws("tausq", inc_warmup = FALSE)
fit_results_fixphi["0.2phi", "kappa", "400", ] <- 
  fit1_02$draws("kappa", inc_warmup = FALSE)

#
data2$phi = 0.2 * phi

if(fit_flag_fix){
  fit2_02 <- mod_fix$sample(
    data = data2,
    seed = 123,
    init = 2,
    chains = N_chain,
    parallel_chains = N_chain,
    save_warmup = TRUE,
    iter_warmup = N_warmup,
    iter_sampling = N_sample,
    refresh = 50,
    sig_figs = 16
  )
  fit2_02$print()
  fit2_02$save_object(file = "./results/sim_bayes_02phi2.RDS")
}else{
  fit2_02 <- readRDS(file = "./results/sim_bayes_02phi2.RDS")
}

fit_results_fixphi["0.2phi", "sigmasq", "900", ] <- 
  fit2_02$draws("sigmasq", inc_warmup = FALSE)
fit_results_fixphi["0.2phi", "tausq", "900", ] <- 
  fit2_02$draws("tausq", inc_warmup = FALSE)
fit_results_fixphi["0.2phi", "kappa", "900", ] <- 
  fit2_02$draws("kappa", inc_warmup = FALSE)

# 
data3$phi = 0.2 * phi

if(fit_flag_fix){
  fit3_02 <- mod_fix$sample(
    data = data3,
    seed = 123,
    init = 2,
    chains = N_chain,
    parallel_chains = N_chain,
    save_warmup = TRUE,
    iter_warmup = N_warmup,
    iter_sampling = N_sample,
    refresh = 50,
    sig_figs = 16
  )
  fit3_02$print()
  fit3_02$save_object(file = "./results/sim_bayes_02phi3.RDS")
}else{
  fit3_02 <- readRDS(file = "./results/sim_bayes_02phi3.RDS")
}

fit_results_fixphi["0.2phi", "sigmasq", "1600", ] <- 
  fit3_02$draws("sigmasq", inc_warmup = FALSE)
fit_results_fixphi["0.2phi", "tausq", "1600", ] <- 
  fit3_02$draws("tausq", inc_warmup = FALSE)
fit_results_fixphi["0.2phi", "kappa", "1600", ] <- 
  fit3_02$draws("kappa", inc_warmup = FALSE)


## fix phi at 0.5phi ##
data1$phi = 0.5 * phi

if(fit_flag_fix){
  fit1_05 <- mod_fix$sample(
    data = data1,
    seed = 123,
    init = 2,
    chains = N_chain,
    parallel_chains = N_chain,
    save_warmup = TRUE,
    iter_warmup = N_warmup,
    iter_sampling = N_sample,
    refresh = 50,
    sig_figs = 16
  )
  
  fit1_05$print()
  fit1_05$save_object(file = "./results/sim_bayes_05phi1.RDS")
}else{
  fit1_05 <- readRDS(file = "./results/sim_bayes_05phi1.RDS")
}


fit_results_fixphi["0.5phi", "sigmasq", "400", ] <- 
  fit1_05$draws("sigmasq", inc_warmup = FALSE)
fit_results_fixphi["0.5phi", "tausq", "400", ] <- 
  fit1_05$draws("tausq", inc_warmup = FALSE)
fit_results_fixphi["0.5phi", "kappa", "400", ] <- 
  fit1_05$draws("kappa", inc_warmup = FALSE)

#
data2$phi = 0.5 * phi

if(fit_flag_fix){
  fit2_05 <- mod_fix$sample(
    data = data2,
    seed = 123,
    init = 2,
    chains = N_chain,
    parallel_chains = N_chain,
    save_warmup = TRUE,
    iter_warmup = N_warmup,
    iter_sampling = N_sample,
    refresh = 50,
    sig_figs = 16
  )
  
  fit2_05$print()
  fit2_05$save_object(file = "./results/sim_bayes_05phi2.RDS")
} else{
  fit2_05 <- readRDS(file = "./results/sim_bayes_05phi2.RDS")
}


fit_results_fixphi["0.5phi", "sigmasq", "900", ] <- 
  fit2_05$draws("sigmasq", inc_warmup = FALSE)
fit_results_fixphi["0.5phi", "tausq", "900", ] <- 
  fit2_05$draws("tausq", inc_warmup = FALSE)
fit_results_fixphi["0.5phi", "kappa", "900", ] <- 
  fit2_05$draws("kappa", inc_warmup = FALSE)

# 
data3$phi = 0.5 * phi

if(fit_flag_fix){
  fit3_05 <- mod_fix$sample(
    data = data3,
    seed = 123,
    init = 2,
    chains = N_chain,
    parallel_chains = N_chain,
    save_warmup = TRUE,
    iter_warmup = N_warmup,
    iter_sampling = N_sample,
    refresh = 50,
    sig_figs = 16
  )
  fit3_05$print()
  fit3_05$save_object(file = "./results/sim_bayes_05phi3.RDS")
}else{
  fit3_05 <- readRDS(file = "./results/sim_bayes_05phi3.RDS")
}

fit_results_fixphi["0.5phi", "sigmasq", "1600", ] <- 
  fit3_05$draws("sigmasq", inc_warmup = FALSE)
fit_results_fixphi["0.5phi", "tausq", "1600", ] <- 
  fit3_05$draws("tausq", inc_warmup = FALSE)
fit_results_fixphi["0.5phi", "kappa", "1600", ] <- 
  fit3_05$draws("kappa", inc_warmup = FALSE)

## fix phi at 2phi ##
data1$phi = 2 * phi

if(fit_flag_fix){
  fit1_2 <- mod_fix$sample(
    data = data1,
    seed = 123,
    init = 2,
    chains = N_chain,
    parallel_chains = N_chain,
    save_warmup = TRUE,
    iter_warmup = N_warmup,
    iter_sampling = N_sample,
    refresh = 50,
    sig_figs = 16
  )
  fit1_2$print()
  fit1_2$save_object(file = "./results/sim_bayes_2phi1.RDS")
} else {
  fit1_2 <- readRDS(file = "./results/sim_bayes_2phi1.RDS")
}


fit_results_fixphi["2phi", "sigmasq", "400", ] <- 
  fit1_2$draws("sigmasq", inc_warmup = FALSE)
fit_results_fixphi["2phi", "tausq", "400", ] <- 
  fit1_2$draws("tausq", inc_warmup = FALSE)
fit_results_fixphi["2phi", "kappa", "400", ] <- 
  fit1_2$draws("kappa", inc_warmup = FALSE)

#
data2$phi = 2 * phi

if(fit_flag_fix){
  fit2_2 <- mod_fix$sample(
    data = data2,
    seed = 123,
    init = 2,
    chains = N_chain,
    parallel_chains = N_chain,
    save_warmup = TRUE,
    iter_warmup = N_warmup,
    iter_sampling = N_sample,
    refresh = 50,
    sig_figs = 16
  )
  fit2_2$print()
  fit2_2$save_object(file = "./results/sim_bayes_2phi2.RDS")
} else {
  fit2_2 <- readRDS(file = "./results/sim_bayes_2phi2.RDS")
}

fit_results_fixphi["2phi", "sigmasq", "900", ] <- 
  fit2_2$draws("sigmasq", inc_warmup = FALSE)
fit_results_fixphi["2phi", "tausq", "900", ] <- 
  fit2_2$draws("tausq", inc_warmup = FALSE)
fit_results_fixphi["2phi", "kappa", "900", ] <- 
  fit2_2$draws("kappa", inc_warmup = FALSE)

# 
data3$phi = 2 * phi

if(fit_flag_fix){
  fit3_2 <- mod_fix$sample(
    data = data3,
    seed = 123,
    init = 2,
    chains = N_chain,
    parallel_chains = N_chain,
    save_warmup = TRUE,
    iter_warmup = N_warmup,
    iter_sampling = N_sample,
    refresh = 50,
    sig_figs = 16
  )
  fit3_2$print()
  fit3_2$save_object(file = "./results/sim_bayes_2phi3.RDS")
} else {
  fit3_2 <- readRDS(file = "./results/sim_bayes_2phi3.RDS")
}

fit_results_fixphi["2phi", "sigmasq", "1600", ] <- 
  fit3_2$draws("sigmasq", inc_warmup = FALSE)
fit_results_fixphi["2phi", "tausq", "1600", ] <- 
  fit3_2$draws("tausq", inc_warmup = FALSE)
fit_results_fixphi["2phi", "kappa", "1600", ] <- 
  fit3_2$draws("kappa", inc_warmup = FALSE)



## fix phi at 5 * phi ##
data1$phi = 5 * phi

if(fit_flag_fix){
  fit1_5 <- mod_fix$sample(
    data = data1,
    seed = 123,
    init = 2,
    chains = N_chain,
    parallel_chains = N_chain,
    save_warmup = TRUE,
    iter_warmup = N_warmup,
    iter_sampling = N_sample,
    refresh = 50,
    sig_figs = 16
  )
  fit1_5$print()
  fit1_5$save_object(file = "./results/sim_bayes_5phi1.RDS")
} else {
  fit1_5 <- readRDS(file = "./results/sim_bayes_5phi1.RDS")
}

fit_results_fixphi["5phi", "sigmasq", "400", ] <- 
  fit1_5$draws("sigmasq", inc_warmup = FALSE)
fit_results_fixphi["5phi", "tausq", "400", ] <- 
  fit1_5$draws("tausq", inc_warmup = FALSE)
fit_results_fixphi["5phi", "kappa", "400", ] <- 
  fit1_5$draws("kappa", inc_warmup = FALSE)

#
data2$phi = 5 * phi

if(fit_flag_fix){
  fit2_5 <- mod_fix$sample(
    data = data2,
    seed = 123,
    init = 2,
    chains = N_chain,
    parallel_chains = N_chain,
    save_warmup = TRUE,
    iter_warmup = N_warmup,
    iter_sampling = N_sample,
    refresh = 50,
    sig_figs = 16
  )
  fit2_5$print()
  fit2_5$save_object(file = "./results/sim_bayes_5phi2.RDS")
}else{
  fit2_5 <- readRDS(file = "./results/sim_bayes_5phi2.RDS")
}

fit_results_fixphi["5phi", "sigmasq", "900", ] <- 
  fit2_5$draws("sigmasq", inc_warmup = FALSE)
fit_results_fixphi["5phi", "tausq", "900", ] <- 
  fit2_5$draws("tausq", inc_warmup = FALSE)
fit_results_fixphi["5phi", "kappa", "900", ] <- 
  fit2_5$draws("kappa", inc_warmup = FALSE)

# 
data3$phi = 5 * phi

if(fit_flag_fix){
  fit3_5 <- mod_fix$sample(
    data = data3,
    seed = 123,
    init = 2,
    chains = N_chain,
    parallel_chains = N_chain,
    save_warmup = TRUE,
    iter_warmup = N_warmup,
    iter_sampling = N_sample,
    refresh = 50,
    sig_figs = 16
  )
  fit3_5$print()
  fit3_5$save_object(file = "./results/sim_bayes_5phi3.RDS")
}else{
  fit3_5 <- readRDS(file = "./results/sim_bayes_5phi3.RDS")
}

fit_results_fixphi["5phi", "sigmasq", "1600", ] <- 
  fit3_5$draws("sigmasq", inc_warmup = FALSE)
fit_results_fixphi["5phi", "tausq", "1600", ] <- 
  fit3_5$draws("tausq", inc_warmup = FALSE)
fit_results_fixphi["5phi", "kappa", "1600", ] <- 
  fit3_5$draws("kappa", inc_warmup = FALSE)


save(list = c("fit_results_fixphi"), 
     file = "./data/sim_bayes_fit_fix_results.RData", 
     envir = .GlobalEnv)





