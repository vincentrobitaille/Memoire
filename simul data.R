library(spdep)
library(Matrix)
library(dplyr)
library(tidyr)
# library(ggplot2)

set.seed(6390)

n = 50

In <- diag(n)

W <- matrix(rbinom(n = n^2, size = 1, prob = 0.4), 
            nrow = n, ncol = n)
diag(W) <- 0
  
W <- W |> forceSymmetric()

W_r <- W |> rowSums()

W <- W / W_r

beta <- matrix(c(1, 1), ncol = 1)

lambda <- 0.6

sigma2 <- 3

epsilon <- rnorm(n, mean = 0, sd = sqrt(sigma2))

X <- matrix(data = c(rep(1, n), rnorm(n = n, mean = 0, sd = 3)), ncol = 2)

y_latent <- solve((In - lambda*W)) %*%( X %*% beta + epsilon)

# mean(y_latent[,1] < 0)
# sum(y_latent[,1] < 0)

y_obs <- matrix(ifelse(y_latent[,1] < 0, 0, y_latent[,1]), ncol = 1)

df_censure <- data.frame(
  y_latent = as.vector(y_latent), 
  y_obs = as.vector(y_obs), 
  X0 = as.vector(X[,1]), 
  X1 = as.vector(X[,2])
)

# 
# plot(X[,2], y_latent[,1])
# abline(h= 0)
# 
# df <- data.frame(y_obs = y_obs[,1], y_latent = y_latent[,1], x2 = X[,2])
# 
# # Reg lin classique
# fit1 <- lm(y_obs ~ x2, data = df)
# 
# summary(fit1)
# 
# df |> 
#   ggplot(aes(x = x2, y = y_obs)) + 
#   geom_point() +
#   geom_smooth(method = "lm", colour = "blue", formula = y ~ x, se = FALSE)



