// The input data is a vector 'y' of length 'N'.
data {
  int N;
  real Y[N];
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real mu;
  real<lower=0> sigma;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model { 
  for(i in 1:N)
    Y[i] ~ normal(mu, sigma);
  mu ~ normal(1.7, 0.3);
  sigma ~ cauchy(0, 1);
}

