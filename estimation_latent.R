
# Hyperparamètres
# beta
m <- c(1, 1)
k <- length(m)
V <- diag(k)
a_beta <- 1
b_beta <- 1
a_gamma <- 1
b_gamma <- 1

# Matrice pour garder les paramètres échantillonnés

N <- 5000
param_hist <- data.frame(
  beta0 = rep(NA, N),
  beta1 = rep(NA, N),
  lambda = rep(NA, N),
  sigma2 = rep(NA, N)
)

# Valeurs de départ pour les paramètres

param_hist[1, ] <- c(0, 0, 0.5, 1)

# Algo de MH
for (i in 1:N) {
 u = runif(n = 1,
           min = 0,
           max = 1)
 
 
} 