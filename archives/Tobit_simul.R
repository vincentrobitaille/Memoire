library(dplyr)
library(tidyr)


W <- matrix()
lambda <- c()
sigma2 <- c()
beta <- vector()
yX <- data.frame(y = NA, x1 = 1, x2 = NA)
n = nrow(yX)

# Certains éléments définis et calculés à l'avance car fixes
ind1 <- which(yX$y == 0)
m = length(ind1)
ind2 <- which(yX$y > 0)
Inm <- diag(n-m)
W22 <- W[ind2, ind2]
W1 <- W[ind1,]
W2 <- W[ind2,]
X1 <- yX[ind1, -1]
X2 <- yX[ind2, -1]
Y = yX[, 1]
Y1 <- yX[ind1, 1]
Y2 <- yX[ind2, 1]

tobit_simul_vrais <- function(Y, Y1, Y2, ind0, ind1, W1, W2, W22, sigma2, beta, lambda) {
  #  À améliorer pour prendre en compte la séparation des données, simplification pour l'instant
  # Précalcul certains morceaux pour simplification du code
  pnorm_i = pnorm((lambda * W[i,] %*% Y + X[i,] %*% beta)/sqrt(sigma2))
  sq_i = (Y[i] - lambda %*% W[i,] %*% Y - X[i,] %*% beta)^2
  log_det = log(det(Inm - lambda * W22))
  
  log_vrais = sum(log(1 - pnorm_i)) - 
    sum(0.5*(log(2*pi*sigma2) + sq_i/sigma2)) +
    log_det
  
  return(log_vrais)
}


