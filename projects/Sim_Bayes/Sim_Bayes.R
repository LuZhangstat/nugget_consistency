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

## fit model with cmdstanr, assume phi is unknown##
library(cmdstanr)
N_warmup = 500
N_sample = 500
N_chain = 4

fit_flag = FALSE
# preallocate fitting results
fit_results <- array(data = NA, dim = c(4, 3, N_sample * N_chain), 
                     dimnames = list(pars = c("sigmasq", "tausq", "phi", "kappa"),
                                     No_location = c("400", "900", "1600"),
                                     index = 1:(N_sample * N_chain)))


file <- file.path(getwd(), "projects", "Sim_Bayes", "GP.stan")
mod <- cmdstan_model(file)

## sample size = 400 ##
data1 <- list(N = n1,
              D = as.matrix(dist(obs_loc[1:n1, ])), # distance matrix
              coords = obs_loc[1:n1, ],
              y = Sim_Bayes_y[1:n1],
              pa = 2.0,
              pb = 2 / phi, 
              sb = 1.0,
              tb = 0.5)
if(fit_flag){
  fit1 <- mod$sample(
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
  fit1$save_object(file = "./results/sim_bayes_full1.RDS")
}else{
  fit1 <- readRDS(file = "./results/sim_bayes_full1.RDS")
}

p1 <- mcmc_trace(fit1$draws("lp__", inc_warmup = TRUE)[-(1:20), ,])#[60:80, ,]) #c("lp__", "phi[1]", "lambda[1]", "theta1[1]")
p1
p2 <- mcmc_trace(fit1$draws("phi", inc_warmup = TRUE)[-(1:20), ,]) #c("lp__", "phi[1]", "lambda[1]", "theta1[1]")
p2
mcmc_hist(fit1$draws("phi", inc_warmup = FALSE))
mcmc_hist(fit1$draws("tausq", inc_warmup = FALSE))
mcmc_hist(fit1$draws("sigmasq", inc_warmup = FALSE))
mcmc_hist(fit1$draws("kappa", inc_warmup = FALSE))

fit_results["sigmasq", "400", ] <- fit1$draws("sigmasq", inc_warmup = FALSE)
fit_results["tausq", "400", ] <- fit1$draws("tausq", inc_warmup = FALSE)
fit_results["phi", "400", ] <- fit1$draws("phi", inc_warmup = FALSE)
fit_results["kappa", "400", ] <- fit1$draws("kappa", inc_warmup = FALSE)


## sample size = 900 ##
data2 <- list(N = n2,
              D = as.matrix(dist(obs_loc[1:n2, ])), # distance matrix
              coords = obs_loc[1:n2, ],
              y = Sim_Bayes_y[1:n2],
              pa = 2.0,
              pb = 2 / phi, 
              sb = 1.0,
              tb = 0.5)
if(fit_flag){
fit2 <- mod$sample(
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

fit2$save_object(file = "./results/sim_bayes_full2.RDS")
}else{
fit2 <- readRDS(file = "./results/sim_bayes_full2.RDS")
}

p1_2 <- mcmc_trace(fit2$draws("lp__", inc_warmup = TRUE)[-(1:20), ,])#[60:80, ,]) #c("lp__", "phi[1]", "lambda[1]", "theta1[1]")
p1_2
p2_2 <- mcmc_trace(fit2$draws("phi", inc_warmup = TRUE)[-(1:20), ,]) #c("lp__", "phi[1]", "lambda[1]", "theta1[1]")
p2_2
mcmc_hist(fit2$draws("phi", inc_warmup = FALSE))
mcmc_hist(fit2$draws("tausq", inc_warmup = FALSE))
mcmc_hist(fit2$draws("sigmasq", inc_warmup = FALSE))
mcmc_hist(fit2$draws("kappa", inc_warmup = FALSE))

fit_results["sigmasq", "900", ] <- fit2$draws("sigmasq", inc_warmup = FALSE)
fit_results["tausq", "900", ] <- fit2$draws("tausq", inc_warmup = FALSE)
fit_results["phi", "900", ] <- fit2$draws("phi", inc_warmup = FALSE)
fit_results["kappa", "900", ] <- fit2$draws("kappa", inc_warmup = FALSE)


## sample size = 1600 ##
data3 <- list(N = n3,
              D = as.matrix(dist(obs_loc[1:n3, ])), # distance matrix
              coords = obs_loc[1:n3, ],
              y = Sim_Bayes_y[1:n3],
              pa = 2.0,
              pb = 2 / phi, 
              sb = 1.0,
              tb = 0.5)

if(fit_flag){
fit3 <- mod$sample(
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
fit3$save_object(file = "./results/sim_bayes_full3.RDS")
}else{
  fit3 <- readRDS(file = "./results/sim_bayes_full3.RDS")
}

p1_3 <- mcmc_trace(fit3$draws("lp__", inc_warmup = TRUE)[-(1:20), ,])#[60:80, ,]) #c("lp__", "phi[1]", "lambda[1]", "theta1[1]")
p1_3
p2_3 <- mcmc_trace(fit3$draws("phi", inc_warmup = TRUE)[-(1:20), ,]) #c("lp__", "phi[1]", "lambda[1]", "theta1[1]")
p2_3
mcmc_hist(fit3$draws("phi", inc_warmup = FALSE))
mcmc_hist(fit3$draws("tausq", inc_warmup = FALSE))
mcmc_hist(fit3$draws("sigmasq", inc_warmup = FALSE))
mcmc_hist(fit3$draws("kappa", inc_warmup = FALSE))

fit_results["sigmasq", "1600", ] <- fit3$draws("sigmasq", inc_warmup = FALSE)
fit_results["tausq", "1600", ] <- fit3$draws("tausq", inc_warmup = FALSE)
fit_results["phi", "1600", ] <- fit3$draws("phi", inc_warmup = FALSE)
fit_results["kappa", "1600", ] <- fit3$draws("kappa", inc_warmup = FALSE)


save(list = c("fit_results"), 
     file = "./data/sim_bayes_fit_results.RData", 
     envir = .GlobalEnv)

