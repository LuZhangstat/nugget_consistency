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
  # phi_delta: decay and psill / nugget
  # D: distance matrix
  # X: data
  
  chol_COV = chol(exp(- phi_delta[1] * D) + phi_delta[2] * diag(n)) 
  u = forwardsolve(chol_COV, X, transpose = T, upper.tri = T)
  return(sum(log(diag(chol_COV))) +  0.5 * n * log(sum(u^2) / n))
}

neg_loglik_0.5_profile_0nug<- function(phi, X, D, n){
  # phi_delta: decay and psill / nugget
  # D: distance matrix
  # X: data
  
  chol_COV = chol(exp(- phi * D)) 
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



