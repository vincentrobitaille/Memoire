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

# Outils

tau = 5*diag(k)
phi = 2


# Algo de MH
MH_spatial <- function(y, X, h_beta = c(a = 1, b = 1), h_gamma = c(a = 1, b = 1),
                       m, V, W, N = 5000, init_param, tau = 5, phi = 2) {
  # Notes à faire: Généraliser pour nombre arbitraire de beta
  #                Intégrer hyperparamètres
  
  # y : variable endo
  # W : matrice spatiale (exogène)
  # X : matrice [1, x1, ... xk-1]
  # h_beta : vecteur nommé c(a_beta = a, b_beta = b) (prior lambda)
  # h_gamma : vecteur nommé c(a_gamma = a, h_gamma = b) (prior lambda)
  # m : vecteur nommé c(beta0 = b0, ..., betak-1 = bk-1) (prior beta)
  # V : matrice de précision (pour ajuster variance et covariance) (prior beta)
  # init_param : vecteur nommé des valeur initiales pour les paramètres
  # N : nombre d'échantillons dans le M-H
  # tau : variance pour le noyau de beta
  # phi : variance pour le noyau de sigma2
  k = length(m)
  jk = which(substr(names(init_param), 1, 1) == "b")
  js = which(substr(names(init_param), 1, 1) == "s")
  jl = which(substr(names(init_param), 1, 1) == "l")
  n = nrow(X)
  tau = tau * diag(k)
  
  
  # Valeurs de départ pour les paramètres
  param_hist[1, ] = c(init_param, rho = NA,  accept = NA)
  
  jra = which(substr(names(init_param), 1, 1) %in% c("r", "a"))
  
  param_hist = dplyr::bind_rows(
    data.frame(as.list(init_param)),
    param_hist
  )
  
  for (i in 2:N) {
    # Probabilité de ref pour acceptation
    u = runif(n = 1,
              min = 0,
              max = 1)
    # Candidats pour betas
    beta_candid = rmvnorm(1, 
                   mean = unlist(param_hist[i -1, jk]),
                   sigma = tau)
    # Candidat pour sigma2
    sigma2_candid = rnorm(1,
                   mean = sqrt(param_hist[i -1,]$sigma2),
                   sd = phi)^2
    lambda_candid = runif(1, 
                   min = 0, 
                   max = 1)
    
    # Rapport des noyaux (log)
    rq = rapport_q_SAR(beta = param_hist[i -1, jk],
                       beta_candid = beta_candid,
                       sigma2 = param_hist[i -1, js],
                       sigma2_candid = sigma2_candid,
                       lambda = param_hist[i -1, jl],
                       lambda_candid = lambda_candid,
                       tau = tau,
                       phi = phi,
                       log = TRUE)
    
    rf = rapport_f_SAR(beta = param_hist[i -1, jk],
                       beta_candid = beta_candid,
                       sigma2 = param_hist[i -1, js],
                       sigma2_candid = sigma2_candid,
                       lambda = param_hist[i -1, jl],
                       lambda_candid = lambda_candid,
                       h_beta = h_beta,
                       h_gamma = h_gamma,
                       y = y, X = X, W = W, n = n, V = V)
    
    rho = exp(rq + rf)
    
    # Acceptation/rejet
    ar = 1*(rho > u)
    
    param_hist[i, jk] = unlist(beta)
    param_hist[i, js] = sigma2
    param_hist[i, jl] = lambda
    param_hist[i, jra] = c(rho, ar)
  } 
}
