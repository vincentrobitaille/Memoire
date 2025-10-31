library(mvtnorm)
source("priors.R")
SAR_vrais <- function(y, lambda, W, X, beta, sigma2, n, log) {
  # Vraisemblance pour un modèle SAR (normale)
  Sn_inv = solve(diag(n) - lambda*W)
  dvrais_cand = dmvnorm(y,
                        mean = Sn_inv %*% X %*% t(beta),
                        sigma = sigma2 * Sn_inv %*% t(Sn_inv),
                        log = TRUE)
  return(dvrais_cand)
}

# rapport_q_SAR <- function(beta, sigma2, lambda, df, i, W, y, X, n) {
#   dbeta_cand = dmvnorm(c(df[i -1,]$beta0, 
#                          df[i -1,]$beta1), 
#                        mean = beta,
#                        sigma = tau, 
#                        log = TRUE) - 
#     dmvnorm(beta, 
#             mean = c(df[i -1,]$beta0, 
#                      df[i -1,]$beta1),
#             sigma = tau, 
#             log = TRUE)
#   
#   dsigma2_cand = dnorm(sqrt(df[i -1,]$sigma2),
#                        mean = sqrt(sigma),
#                        sd = phi,
#                        log = TRUE) - 
#     dnorm(sqrt(sigma),
#           mean = sqrt(df[i -1,]$sigma2),
#           sd = phi,
#           log = TRUE)
#   
#   Sn_inv = solve(diag(n) - lambda*W)
#   
#   Sn_inv2 = solve(diag(n) - df[i -1,]$lambda*W)
#   
#   dvrais_cand = dnorm(y,
#                       mean = Sn_inv %*% X %*% t(beta),
#                       sd = sqrt(sigma2 * Sn_inv %*% t(Sn_inv)),
#                       log = TRUE) - 
#     dnorm(y,
#           mean = Sn_inv2 %*% X %*% t(df[i -1, c("beta0", "beta1")]),
#           sd = sqrt(df[i -1,]$sigma2 * Sn_inv2 %*% t(Sn_inv2)),
#           log = TRUE)
#   
#   rho = exp(dbeta_cand + dsigma2_cand + dvrais_cand)
#   
#   return(rho)
# }

rapport_q_SAR <- function(beta, beta_candid, sigma2, sigma2_candid, lambda, 
                          lambda_candid, tau, phi, log = TRUE) {
  dbeta_cand = dmvnorm(unlist(beta), 
                       mean = beta_candid,
                       sigma = tau, 
                       log = TRUE) - 
    dmvnorm(beta_candid, 
            mean = unlist(beta),
            sigma = tau, 
            log = TRUE)
  
  dsigma2_cand = dnorm(sqrt(sigma2),
                       mean = sqrt(sigma2_candid),
                       sd = phi,
                       log = TRUE) - 
    dnorm(sqrt(sigma2_candid),
          mean = sqrt(sigma2),
          sd = phi,
          log = TRUE)
  
  if (log) {
    rap_q = dbeta_cand + dsigma2_cand
  } else {
    rap_q = exp(dbeta_cand + dsigma2_cand)
  }
  
  return(rap_q)
}

rapport_f_SAR <- function(beta, beta_candid, m, sigma2, sigma2_candid, 
                          lambda, lambda_candid, h_beta, h_gamma, 
                          y, X, W, n, V, log = TRUE) {
  dbeta_prior = d_prior_beta(beta = beta_candid,
                             m = m,
                             sigma2 = sigma2_candid, 
                             V = V) - 
    d_prior_beta(beta = unlist(beta),
                 m = m,
                 sigma2 = sigma2, 
                 V = V)
  dsigma2_prior = d_prior_sigma2(sigma2 = sigma2_candid,
                                 a = h_gamma["a"], 
                                 b = h_gamma["b"]) -
    d_prior_sigma2(sigma2 = sigma2,
                   a = h_gamma["a"], 
                   b = h_gamma["b"])
  dlambda_prior = d_prior_lambda(lambda = lambda_candid,
                                 a = h_beta["a"],
                                 b = h_beta["b"]) -
    d_prior_lambda(lambda = lambda,
                   a = h_beta["a"],
                   b = h_beta["b"])
  dSAR_vrais = SAR_vrais(y = y, lambda = lambda_candid, W = W, X = X, 
                         beta = beta_candid, sigma2 = sigma2_candid, 
                         n, log = TRUE) - 
    SAR_vrais(y = y, lambda = lambda, W = W, X = X, beta = beta, 
              sigma2 = sigma2, n, log = TRUE)
  
  if (log) {
    rap_f = dbeta_prior + dsigma2_prior + dlambda_prior + dSAR_vrais
  } else {
    rap_f = exp(dbeta_prior + dsigma2_prior + dlambda_prior + dSAR_vrais)
  }
  
  return(rap_f)
}
