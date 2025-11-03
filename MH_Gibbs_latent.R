library(mvtnorm)

source("priors.R")
source("simul data.R")
source("outils.R")

MH_Gibbs_latent <- function(y, X, W, c_beta = c(0, 0), c_lambda = 0.1, h_T = 1e+10*diag(2), N = 5000,
                            h_igamma = c(a = 0, b = 0), censure = NULL, m_step = 1,
                            theta_0 = c(beta0 = 0, beta1 = 0, sigma2 = 2, lambda = 0.1)) {
  # y : vecteur de la variable endogène
  # X : matrice [1, x1, ... xk-1] des variables exogènes
  # W : matrice spatiale
  # c : Paramètre de tuning ?
  # h_T : Hyperparamètre pour les betas
  # h_gamma : hyperparamètres (a, b) pour le prior gamma de sigma2
  # N : Taille d'échantillonnage
  # theta_0 : Vecteur nommé de paramètres initiaux (pour beta, sigma2 et lambda)
  # censure : si vrai, y fourni est censuré, faux si pas censuré (variable latente)
  
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
  
  SAR_vrais_gibbs = function(y, lambda, W, Sn_inv, X, beta, sigma2, n, log) {
    # Vraisemblance pour un modèle SAR (normale)
    dvrais_cand = dmvnorm(t(y),
                          mean = Sn_inv %*% (X %*% t(beta)),
                          sigma = sigma2 * (Sn_inv %*% t(Sn_inv)),
                          log = TRUE)
    return(dvrais_cand)
  }
  
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
    }
  
  # MCMC
  pb_mh = txtProgressBar(min = 1, max = N-1)
  for (i in 1:N) {
    Sn = diag(n) - lambda*W
    Sn_inv = solve(Sn)
    Xb = X%*%t(beta)
    
    # Échantillonneur de Gibbs-Geweke pour les variables latentes
    if (!is.null(censure)) {
      for (m_s in 1:m_step) {
        for (j in y_lat_ind) {
          y[j] = lambda*W[j]%*%y + Xb + rnorm(1, mean = 0, sd = sqrt(sigma2))
        }
      }
    }
    
    # Échantillonnage posterior conditionnel beta (Gibbs)
    # p(beta|sigma2, lambda, y)
    h_T_star = solve(t(X) %*% X + h_T_inv)
    c_star = h_T_star %*% (t(X) %*% Sn %*% y + h_T_inv %*% c_beta)
    
    beta = rmvnorm(1, c_star, sigma2*h_T_star)
    
    # Échantillonnage posterior conditionnel sigma2 (Gibbs)
    # p(sigma2|beta, lambda, y)
    Sny_Xbeta = Sn %*% y - X %*% t(beta)
    b_star = h_igamma["b"] + t(Sny_Xbeta) %*% (Sny_Xbeta/2)
    
    sigma2 = 1/rgamma(1, shape = a_star, rate = b_star)
    
    u = runif(1)
    
    # Échantillonnage posterior conditionnel lambda (MH)
    # p(sigma2|beta, lambda, y)
    lambda_candid = lambda + c_lambda*rnorm(1)
    
    Sn_candid = diag(n) - lambda_candid*W
    
    # On rejette les candidats de lambda qui ne sont pas dans (-1, 1)
    if (lambda_candid < 1 & 
        lambda_candid > -1) {
      rho = lambda_cond_prop(y = y, Xb = Xb, Sn = Sn, Sn_candid = Sn_candid,
                             sigma2 = sigma2)
    } else {
      rho = 0
    }
    
    accept = 1*(rho > u)
    
    sum_accept = sum_accept + accept
    
    mean_accept = sum_accept/i
    
    # Tuning pour la variance de la marche aléatoire de la proposition de lambda
    if ( mean_accept > 0.6) {
      c_lambda = c_lambda*1.1
    } else if ( mean_accept < 0.4) {
      c_lambda = c_lambda/1.1
    }
    
    if ( accept == 1 ) {
      lambda = lambda_candid
      # À Vérifier pour beta et sigma2 ?
      # beta = beta_candid
    }
    
    param_hist[i+1, 
               c(n_beta, "sigma2", "lambda", 
                 "rho", "accept", "c_lambda")] = c(beta[1,], sigma2, lambda, 
                                                   rho, accept, c_lambda)
    
    setTxtProgressBar(pb_mh, i)
  }
  
  return(param_hist)
}