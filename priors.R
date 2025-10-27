library(mvtnorm)
source("simul data.R")
source("outils.R")
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


d_prior_lambda <- function(lambda, a, b) {
  dl = dbeta((lambda +1)/2, 
             shape1 = a, 
             shape2 = b,
             log = TRUE)
  return(dl)
}

d_prior_sigma2 <- function(sigma2, a, b) {
  ds = dgamma(1/sigma2,
              shape = a,
              scale = b,
              log = TRUE)
  return(ds)
}

d_prior_beta <- function(beta, sigma2, V = diag(length(m))) {
  db = mvtnorm::dmvnorm(beta,
                        mean = m,
                        sigma = sigma2*V,
                        log = TRUE)
  return(db)
}
