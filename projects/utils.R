## utils function ##

loglik_chol<- function(X, chol_COV, sigma2){
  # X: data
  # chol_COV: cholesky decomposition of the covariance matrix / sigma2
  # sigma2: partial sill
  u = forwardsolve(chol_COV, X, transpose = T, upper.tri = T)
  n = length(X)
  return(- sum(log(diag(chol_COV))) - 0.5 * n * log(sigma2) - 
           0.5 * sum(u^2) / sigma2)
}


neg_loglik_0.5_profile<- function(phi_delta, X, D, n){
  # phi_delta: decay and nugget / psill
  # D: distance matrix
  # X: data
  
  chol_COV = chol(exp(- phi_delta[1] * D) + phi_delta[2] * diag(n)) 
  u = forwardsolve(chol_COV, X, transpose = T, upper.tri = T)
  return(sum(log(diag(chol_COV))) +  0.5 * n * log(sum(u^2) / n))
}

neg_loglik_0.5_profile_X<- function(phi_delta_beta, X, Y, D, n){
  # phi_delta: decay and nugget / psill
  # D: distance matrix
  # X: data
  
  chol_COV = chol(exp(- phi_delta_beta[1] * D) + phi_delta_beta[2] * diag(n)) 
  mu = Y - X %*% phi_delta_beta[c(-1, -2)]
  u = forwardsolve(chol_COV, mu, transpose = T, upper.tri = T)
  return(sum(log(diag(chol_COV))) +  0.5 * n * log(sum(u^2) / n))
}

neg_loglik_0.5_profile_0nug<- function(phi, X, D, n){
  # phi_delta: decay and nugget / psill
  # D: distance matrix
  # X: data
  
  chol_COV = chol(exp(- phi * D)) 
  u = forwardsolve(chol_COV, X, transpose = T, upper.tri = T)
  return(sum(log(diag(chol_COV))) +  0.5 * n * log(sum(u^2) / n))
}

neg_loglik_0.5_fixphi<- function(delta, phi_fix, X, D, n){
  # phi_delta: decay and nugget / psill
  # D: distance matrix
  # X: data
  
  chol_COV = chol(exp(- phi_fix * D) + delta * diag(n)) 
  u = forwardsolve(chol_COV, X, transpose = T, upper.tri = T)
  return(sum(log(diag(chol_COV))) +  0.5 * n * log(sum(u^2) / n))
}

sigma2_MLE_cal <- function(phi, delta, X, D, n){
  chol_COV = chol(exp(- phi * D) + delta * diag(n))
  u = forwardsolve(chol_COV, X, transpose = T, upper.tri = T)
  return(sum(u^2) / n)
}

sigma2_MLE_cal_0nug <- function(phi, X, D, n){
  chol_COV = chol(exp(- phi * D))
  u = forwardsolve(chol_COV, X, transpose = T, upper.tri = T)
  return(sum(u^2) / n)
}

percentile_BIAS_SD <- function(x, true_value, 
                               probs = c(0.05, 0.25, 0.5, 0.75, 0.95), 
                               round = 3){
  return(paste0(round(c(quantile(x, probs = probs), mean(x) - true_value, sd(x)), 
                      round), collapse = " & "))
}

MSPE_theta01 <- function(N_obs, N_pred, phi, sigma2, tau2, phi1, sigma21, tau21, 
                         D_obs, D_obs_pred, D_pred_obs, D_pred){
  invKSn = chol2inv(chol(sigma21 * exp(- phi1 * D_obs) + 
                           diag(rep(1, N_obs)) * tau21))
  K1USinvKK0SU = sigma21 * exp(- phi1 * D_pred_obs) %*% invKSn %*% 
    (sigma2 * exp(- phi * D_obs_pred))
  VAR = sigma2 * exp(- phi * D_pred) + diag(rep(1, N_pred)) * tau2 +
   sigma21 * exp(- phi1 * D_pred_obs) %*% invKSn %*%
   (sigma2 * exp(- phi * D_obs) + diag(rep(1, N_obs)) * tau2) %*% invKSn %*%
   (sigma21 * exp(- phi1 * D_obs_pred)) - K1USinvKK0SU - t(K1USinvKK0SU)
  # VAR = sigma2 * exp(- phi * D_pred) + 
  #   sigma21 * exp(- phi1 * D_pred_obs) %*% invKSn %*% 
  #   (sigma2 * exp(- phi * D_obs) + diag(rep(1, N_obs)) * tau2) %*% invKSn %*%
  #   (sigma21 * exp(- phi1 * D_obs_pred)) - K1USinvKK0SU - t(K1USinvKK0SU)
  return(mean(diag(VAR)))
}

MSPE_theta0 <- function(N_obs, N_pred, phi, sigma2, tau2, D_obs, D_obs_pred, 
                        D_pred){
  CholL0Sn = t(chol(sigma2 * exp(- phi * D_obs) + 
                      diag(rep(1, N_obs)) * tau2))
  vUSn = forwardsolve(CholL0Sn, sigma2 * exp(- phi * D_obs_pred))
  VAR = sigma2 * exp(- phi * D_pred) + diag(rep(1, N_pred)) * tau2 - 
    crossprod(vUSn)
  #VAR = sigma2 * exp(- phi * D_pred) - crossprod(vUSn)
  return(mean(diag(VAR)))
}
