SAR_vrais <- function(y, lambda, W, X, beta, sigma2, n, log) {
  # Vraisemblance pour un modèle SAR (normale)
  Sn_inv = solve(diag(n) - lambda*W)
  dvrais_cand = dnorm(y,
                      mean = Sn_inv %*% X %*% t(beta),
                      sd = sqrt(sigma2 * Sn_inv %*% t(Sn_inv)),
                      log = TRUE)
  return(dvrais_cand)
}

rapport_q_SAR <- function(beta, sigma2, lambda, df, i, W, y, X, n) {
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
