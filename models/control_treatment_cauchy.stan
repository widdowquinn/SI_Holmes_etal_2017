data {
  int<lower=0> M;  # number of training datapoints
  int<lower=0> N;  # number of test datapoints
  int<lower=0> J;  # number of unique probes
  int<lower=1, upper=J> tagidx_train[M];  # probe indices (training)
  int<lower=1, upper=J> tagidx_test[N];   # probe indices (test)
  vector[M] t_train;
  vector[N] t_test;
  vector[M] x_train;  # training input datapoints
  vector[N] x_test;   # test input datapoints
  vector[M] y;
}
parameters {
  vector[J] a;
  vector[J] b;
  vector[J] g;
  vector[J] d;
  real mu_a;
  real mu_b;  
  real mu_g;
  real mu_d;  
  real<lower=0> sigma_y;
  real<lower=0,upper=100> sigma_a;
  real<lower=0,upper=100> sigma_b;
  real<lower=0,upper=100> sigma_g;  
  real<lower=0,upper=100> sigma_d;  
}
transformed parameters{
  vector[M] y_hat;
  vector[N] mu_pred;  

  for (i in 1:M)
    y_hat[i] = a[tagidx_train[i]] + b[tagidx_train[i]] * x_train[i] +
               g[tagidx_train[i]] * t_train[i] + d[tagidx_train[i]] * t_train[i] * x_train[i];
    
  for (j in 1:N)
    mu_pred[j] = a[tagidx_test[j]] + b[tagidx_test[j]] * x_test[j] +
                 g[tagidx_test[j]] * t_test[j] + d[tagidx_test[j]] * t_test[j] * x_test[j];
}
model {
  sigma_a ~ uniform(0, 100);
  a ~ cauchy(mu_a, sigma_a);

  sigma_b ~ uniform(0, 100);
  b ~ cauchy(mu_b, sigma_b);

  sigma_g ~ uniform(0, 100);
  g ~ cauchy(mu_g, sigma_g);

  sigma_d ~ uniform(0, 100);
  d ~ cauchy(mu_d, sigma_d);

  y ~ normal(y_hat, sigma_y);
}
generated quantities {
  vector[N] y_pred;
  
  for (i in 1:N)
    y_pred[i] = normal_rng(mu_pred[i], sigma_y);
}
