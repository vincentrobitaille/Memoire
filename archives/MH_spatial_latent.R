library(mvtnorm)

source("priors.R")
source("simul data.R")
source("outils.R")
# Hyperparamètres
# beta
# m <- c(0, 0)
# k <- length(m)
# V <- diag(k)
# a_beta <- 1
# b_beta <- 1
# a_gamma <- 1
# b_gamma <- 1
# 
# # Outils
# 
# tau = 5*diag(k)
# phi = 2


# Algo de MH
MH_spatial <- function(y, X, h_beta = c(a = 1, b = 1), h_gamma = c(a = 1, b = 1),
                       m, V, W, N = 5000, init_param, tau = 5, phi = 2, c_lambda = 1) {
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
  
  beta_acc = matrix(init_param[jk], ncol = 2)
  colnames(beta_acc) = names(init_param[jk])
  sigma2_acc = init_param[js]
  lambda_acc = init_param[jl]
  
  # Valeurs de départ pour les paramètres
  param_hist = data.frame(as.list(init_param), 
                          rho = NA,
                          accept = NA)
  # param_hist[1, ] = c(init_param, rho = NA,  accept = NA)
  
  jra = which(substr(names(init_param), 1, 1) %in% c("r", "a"))
  
  p2 = matrix(nrow = N, ncol = k+4)
  colnames(p2) = colnames(param_hist)
  
  param_hist = dplyr::bind_rows(
    param_hist,
    p2 |> as.data.frame()
  )
  
  rm(p2)
  
  pb_mh = txtProgressBar(min = 2, max = N+1)
  
  ar_hist = 0
  
  for (i in 2:(N+1)) {
    # Probabilité de ref pour acceptation
    u = runif(n = 1,
              min = 0,
              max = 1)
    # Candidats pour betas
    beta_candid = beta_acc + 
      rmvnorm(1, 
              mean = rep(0, k),
              sigma = tau)
    # Candidat pour sigma2
    sigma2_candid = (sqrt(sigma2_acc) + 
                       rnorm(1, mean = 0, sd = phi))^2
    # sigma2_candid = rnorm(1,
    #                mean = sqrt(sigma2_acc),
    #                sd = phi)^2
    lambda_candid = lambda_acc + c_lambda*rnorm(1,
                                                mean = 0,
                                                sd = 1)
    # lambda_candid = runif(1, 
    #                min = -1, 
    #                max = 1)
    
    # Rapport des noyaux (log)
    rq = rapport_q_SAR(beta = beta_acc,
                       beta_candid = beta_candid,
                       sigma2 = sigma2_acc,
                       sigma2_candid = sigma2_candid,
                       lambda = lambda_acc,
                       lambda_candid = lambda_candid,
                       tau = tau,
                       phi = phi,
                       log = TRUE)
    
    rf = rapport_f_SAR(beta = beta_acc,
                       beta_candid = beta_candid,
                       sigma2 = sigma2_acc,
                       sigma2_candid = sigma2_candid,
                       lambda = lambda_acc,
                       lambda_candid = lambda_candid,
                       h_beta = h_beta,
                       h_gamma = h_gamma,
                       y = y, X = X, W = W, n = n, V = V, m = m,
                       log = TRUE)
    
    rho = min(c(exp(rq + rf), 1))
    
    # Acceptation/rejet
    ar = 1*(rho > u)
    ar_hist = ar_hist + ar
    ar_m = ar_hist/(i-1)
    
    # if (ar_m < 0.4) {
    #   c_lambda = c_lambda/1.1
    # } else if (ar_m > 0.6) {
    #   c_lambda = c_lambda*1.1
    # }
    
    if (ar == 1) {
      beta_acc = beta_candid
      sigma2_acc = sigma2_candid
      lambda_acc = lambda_candid
    }
    
    param_hist[i, jk] = beta_candid[1,]
    param_hist[i, js] = sigma2_candid
    param_hist[i, jl] = lambda_candid
    param_hist[i, "rho"] = rho
    param_hist[i, "accept"] = ar
    
    # Sys.sleep(0.0001)
    setTxtProgressBar(pb_mh, i)
  }
  
  return(param_hist)
}
