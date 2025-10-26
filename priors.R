library(mvtnorm)
source("simul data.R")
# https://statswithr.github.io/book/introduction-to-bayesian-regression.html
# Le prior pour lambda est p(lambda*)
# lambda* = (lambda + 1)/2
# lambda* ~ B(a, b)
prior_lambda_beta <- function(a, b) {
  # a : (hyperparamètre)
  # b : (hyperparamètre)
  lambda_star = rbeta(n = 1,
                      shape1 = a,
                      shape2 = b)
  lambda_prior = 2*lambda_star -1
  
  return(lambda_prior)
}

# Le prior pour sigma2 est p(sigma2|lambda)
# p(sigma2|lambda) ~ Gamma-1(a, b)
prior_sigma2_ig <- function(a, b) {
  # a : (hyperparamètre)
  # b : (hyperparamètre)
  sigma2_prior = 1/rgamma(n = 1,
                          shape = a,
                          scale = b)
  
  return(sigma2_prior)
}

prior_beta_normal <- function(m, sigma2, V = diag(length(m))) {
  # Prior de beta|sigma2
  # m : espérance (hyperparamètre)
  # sigma2 : variance provenant du prior de sigma2
  # V : matrice de précision (pour ajuster variance et covariance) (hyperparamètre)
  beta_prior = mvtnorm::rmvnorm(n = 1,
                                mean = m,
                                sigma = sigma2*V)
  
  return(beta_prior)
}


vrais_prior_lambda_beta <- function(lambda, a, b) {
  d
}
