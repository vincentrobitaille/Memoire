library(mvtnorm)

source("priors.R")
# Hyperparamètres
# beta
m <- c(0, 0)
k <- length(m)
V <- diag(k)
a_beta <- 1
b_beta <- 1
a_gamma <- 1
b_gamma <- 1

n <- 500
# Matrice pour garder les paramètres échantillonnés

N <- 5000
param_hist <- data.frame(
  beta0 = rep(NA, N),
  beta1 = rep(NA, N),
  lambda = rep(NA, N),
  sigma2 = rep(NA, N),
  rho = rep(NA, N),
  accept = rep(NA, N)
)

# Valeurs de départ pour les paramètres

param_hist[1, ] <- c(0, 0, 0.5, 1, NA, NA)


# Outils

tau = 5*diag(k)
phi = 2

rapport_q <- function(beta, sigma2, lambda, df, i, W, y, X) {
  dbeta_cand = dmvnorm(c(df[i -1,]$beta0, 
                    df[i -1,]$beta1), 
                  mean = beta,
                  sigma = tau, 
                  log = TRUE) - 
    dmvnorm(beta, 
            mean = c(df[i -1,]$beta0, 
                     df[i -1,]$beta1),
            sigma = tau, 
            log = TRUE)
  
  dsigma2_cand = dnorm(sqrt(df[i -1,]$sigma2),
                  mean = sqrt(sigma),
                  sd = phi,
                  log = TRUE) - 
    dnorm(sqrt(sigma),
          mean = sqrt(df[i -1,]$sigma2),
          sd = phi,
          log = TRUE)
  
  Sn_inv = solve(diag(n) - lambda*W)
  
  Sn_inv2 = solve(diag(n) - df[i -1,]$lambda*W)
  
  dvrais_cand = dnorm(y,
                 mean = Sn_inv %*% X %*% t(beta),
                 sd = sqrt(sigma2 * Sn_inv %*% t(Sn_inv)),
                 log = TRUE) - 
    dnorm(y,
          mean = Sn_inv2 %*% X %*% t(df[i -1, c("beta0", "beta1")]),
          sd = sqrt(df[i -1,]$sigma2 * Sn_inv2 %*% t(Sn_inv2)),
          log = TRUE)
  
  rho = exp(dbeta_cand + dsigma2_cand + dvrais_cand)
  
  return(rho)
}


# Algo de MH
for (i in 2:N) {
  # Probabilité de ref pour acceptation
  u = runif(n = 1,
            min = 0,
            max = 1)
  # Candidats pour betas
  beta = rmvnorm(1, 
                 mean = c(param_hist[i -1,]$beta0, param_hist[i -1,]$beta1),
                 sigma = tau)
  # Candidat pour sigma2
  sigma2 = rnorm(1,
                  mean = sqrt(param_hist[i -1,]$sigma2),
                  sd = phi)^2
  lambda = runif(1, 
                 min = 0, 
                 max = 1)
  
  # Calcul probabilité d'acceptation
  rho = rapport_q(beta = beta,
                  sigma2 = sigma2,
                  lambda = lambda,
                  df = param_hist,
                  i = i,
                  W = W,
                  y = y,
                  X = X)
  
  # Acceptation/rejet
  ar = 1*(rho > u)
  
  param_hist[i, ] = c(beta[1,1], beta[1,2], lambda, sigma2, rho, ar)
} 