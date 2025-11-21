library(mvtnorm)
library(truncnorm)

source("outils.R")

MH_Gibbs_latent <- function(y, X, W, c_beta = c(0, 0), c_lambda = 1, h_T = 1e+10*diag(2), N = 5000,
                            h_igamma = c(a = 0, b = 0), censure = NULL, m_step = 1,
                            theta_0 = c(beta0 = 0, beta1 = 0, sigma2 = 2, lambda = 0.1),
                            tuning_lambda = FALSE) {
  # y : vecteur de la variable endogène
  # X : matrice [1, x1, ... xk-1] des variables exogènes
  # W : matrice spatiale
  # c_lambda : Paramètre initial le tuning (écart-type de la proposition de lambda)
  # h_T : Hyperparamètre pour les betas
  # h_gamma : hyperparamètres (a, b) pour le prior gamma de sigma2
  # N : Taille d'échantillonnage
  # theta_0 : Vecteur nommé de paramètres initiaux (pour beta, sigma2 et lambda)
  # censure : si vrai, y fourni est censuré, faux si pas censuré (variable latente)
  # tuning_lambda : Si vrai, fait un tuning pour la proposition de lambda
  
  # Diverses préparations pour faciliter le processus
  # Calculs à l'avance
  # Matrice pour stocker échantillons, etc
  {
    param_hist = data.frame(
      as.list(theta_0),
      c_lambda = c_lambda,
      rho = NA,
      accept = NA
    )
    
    ph2 = as.data.frame(matrix(nrow = N, ncol = ncol(param_hist)))
    colnames(ph2) = colnames(param_hist)
    
    param_hist = rbind(
      param_hist,
      ph2
    )
    rm(ph2)
    
    n_beta = names(theta_0)[substr(names(theta_0), 1, 4) == "beta"]
    
    k = ncol(X)
    n = nrow(X)
    y = matrix(y, ncol = 1)
    c_beta = matrix(c_beta, ncol = 1)
    a_star = h_igamma["a"] + n/2
    h_T = h_T*diag(k)
    h_T_inv = solve(h_T)
    In = diag(n)
  }
  
  # Rapport des vraisemblances pour p(Y|lambda',...)/p(Y|lambda, ...)
  lambda_cond_prop = function(y, Xb, Sn, Sn_candid, sigma2) {
    Sny_Xb_candid = Sn_candid %*% y - Xb
    Sny_Xb = Sn %*% y - Xb
    
    d_lambda = exp(
      log(det(Sn_candid)) - 
        log(det(Sn)) - 
        (1/(2*sigma2))*t(Sny_Xb_candid)%*%Sny_Xb_candid +
        (1/(2*sigma2))*t(Sny_Xb)%*%Sny_Xb
      )
    return(d_lambda)
  }
  
  sum_accept = 0
  
  lambda = theta_0["lambda"]
  sigma2 = theta_0["sigma2"]
  beta = matrix(theta_0[substr(names(theta_0), 1, 4) == "beta"], ncol = 2)
  
  # Préparation si utilisation de Geweke
  if (!is.null(censure)) {
    y_obs = y
    y_lat_ind = which(y_obs == censure)
    y_latent_noms <- paste0("y", y_lat_ind)
    param_hist[, y_latent_noms] = 0
  }
  
  # Chaine de MCMC MH et Gibbs
  pb_mh = txtProgressBar(min = 1, max = N-1)
  for (i in 1:N) {
    Sn = In - lambda*W
    Sn_inv = solve(Sn)
    Xb = X%*%t(beta)
    
    # Échantillonneur de Gibbs-Geweke pour les variables latentes
    if (!is.null(censure)) {
      for (m_s in 1:m_step) {
        for (j in y_lat_ind) {
          # Échantillonne d'une normale multivariée tronquée
          y[j] = rtruncnorm(n = 1, mean = ((lambda*W%*%y)[j] + Xb[j]), sd = sqrt(sigma2), b = 0)
        }
      }
      param_hist[i, y_latent_noms] = y[y_lat_ind]
    }
    
    # Échantillonnage posterior conditionnel beta (Gibbs)
    # p(beta|sigma2, lambda, y)
    h_T_star = solve(t(X) %*% X + h_T_inv)
    c_star = h_T_star %*% (t(X) %*% Sn %*% y + h_T_inv %*% c_beta)
    
    beta = rmvnorm(1, c_star, sigma2*h_T_star)

    Xb = X %*% t(beta)    
    # Échantillonnage posterior conditionnel sigma2 (Gibbs)
    # p(sigma2|beta, lambda, y)
    Sny_Xbeta = Sn %*% y - Xb
    b_star = h_igamma["b"] + t(Sny_Xbeta) %*% (Sny_Xbeta/2)
    
    sigma2 = 1/rgamma(1, shape = a_star, rate = b_star)
    
    # Échantillonnage posterior conditionnel lambda (MH)
    # p(sigma2|beta, lambda, y)
    # Ancienne méthode
    # lambda_candid = lambda + c_lambda*rnorm(1)
    lambda_transf_candid = log((1+lambda)/
                                 (1-lambda)) + c_lambda*rnorm(n = 1, 
                                                            mean = 0, 
                                                            sd = 1)
    e_lambda_candid = exp(lambda_transf_candid) 
    lambda_candid = (e_lambda_candid - 1)/(e_lambda_candid + 1)
    
    Sn_candid = In - lambda_candid*W
    
    # Avec propositions et priors symétriques, rho est le rapport des vraisemblances
    rho = lambda_cond_prop(y = y, Xb = Xb, Sn = Sn, Sn_candid = Sn_candid,
                           sigma2 = sigma2)
    
    accept = 1*(rho > runif(1))
    
    sum_accept = sum_accept + accept
    
    mean_accept = sum_accept/i
    
    # Tuning pour la variance de la marche aléatoire de la proposition de lambda
    if(tuning_lambda) {
      c_lambda = max(c_lambda + (mean_accept-0.45)/(i^0.6), 1e-8) 
    }
    
    # Actualise la valeur de lambda si le candidat est accepté
    if ( accept == 1 ) {
      lambda = lambda_candid
    }
    
    param_hist[i+1, 
               c(n_beta, "sigma2", "lambda", 
                 "rho", "accept", "c_lambda")] = c(beta[1,], sigma2, lambda, 
                                                   rho, accept, c_lambda)
    
    setTxtProgressBar(pb_mh, i)
  }
  
  return(param_hist)
}