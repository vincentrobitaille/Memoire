library(dplyr)
library(tidyr)
library(tmvtnorm)
library(msm)

n = 10

tobit_latent_post <- function(Y, index_censure, W, rho, lambda, beta, sigma2, n) {
  In = diag(n)
  
  log_pi = -n*log(2*pi*sqrt(sigma2))
  
  log_det = log(det(In - lambda %*% W))
  
  sq_i = t((In - lambda %*% W)%*%Y - X%*%beta)%*%((In - lambda %*% W)%*%Y - X%*%beta)
  
  log_vrais = log_pi + log_det - 1/(2*sigma2) * sq_i
  
  return(log_vrais)
}

sample_pdf_latent <- function(Y, W, lambda, beta, sigma2, index_censure, m) {
  
  # Méthode d'échantillonnage de Gibbs (Geweke) pour l'estimation des y latents
  
  In = diag(n)
  Sn = In - lambda * W
  Sn_inv = solve(Sn)
  y_tilde = Sn_inv %*% X %*% beta
  Xbeta = X %*% beta
  
  z_latent_t = as.vector(rep(0, length(index_censure)))
  
  Y_latent = Y
  
  for (m_pass in (1:m)) {
    for (mn in (1:length(index_censure))) {
      epsilon = Sn %*% Y_latent - Xbeta
      omega = Sn_inv %*% epsilon
      i = index_censure[mn]
      z_i = msm::rtnorm(n = 1,
                        mean = y_tilde[i],
                        sd = sqrt(sigma2*(omega[i])^2),
                        upper = 0)
      
      z_latent_t[mn] = z_i
      
      Y_latent[i] = z_i
    }
    }
  return(Y_latent)
  }

