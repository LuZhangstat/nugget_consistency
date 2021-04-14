data {
  int<lower=1> N;
  matrix[N,N] D;
  matrix[N, 2] coords;
  vector[N] y;       // response
  real pa;           // shape parameter in prior of phi
  real pb;           // rate parameter in prior of phi
  real sb;           // scale parameter in prior of signmasq
  real tb;           // scale parameter in prior of tausq
}

parameters {
  real<lower=0> tausq;    
  real<lower=0> phi;    // decay  
  real<lower=0> sigmasq;  // psill
}

transformed parameters {
	real kappa = sigmasq * phi;    // kappa
}

model {
  matrix[N, N] cov =  sigmasq * exp(- D * phi) 
                     + diag_matrix(rep_vector(tausq, N));
  matrix[N, N] L_cov = cholesky_decompose(cov);
  
  //phi ~ inv_gamma(2, pb);
  phi ~ gamma(pa, pb);
  sigmasq ~ inv_gamma(2, sb);
  tausq ~ inv_gamma(2, tb);

  y ~ multi_normal_cholesky(rep_vector(0, N), L_cov);
}
