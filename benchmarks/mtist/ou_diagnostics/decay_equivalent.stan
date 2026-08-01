data {
  int<lower=1> N;
  vector[N] y;
  vector[N] xi;
  vector[N] xj;
  int<lower=1> S;
  array[N] int<lower=1, upper=S> sid;
  array[N] int<lower=0, upper=N> prev;
  vector[N] dt;
}
transformed data {
  real dt_unit;
  real sum_dt = 0;
  int cnt = 0;
  for (n in 1:N) if (prev[n] != 0 && dt[n] > 0) {
    sum_dt += dt[n]; cnt += 1;
  }
  dt_unit = cnt > 0 ? sum_dt / cnt : 1;
}
parameters {
  real r0;
  vector[S] r0_raw;
  real<lower=0> sd_r0;
  real a_ii;
  real a_ij;
  real<lower=0> sigma;
  real<lower=0> sd_ou;
  real<lower=log(-log(0.99) / dt_unit),
       upper=log(-log(0.01) / dt_unit)> log_lambda;
  real log_nu_minus_two;
  vector[N] z_e;
}
transformed parameters {
  real<lower=0> lambda = exp(log_lambda);
  real<lower=0.01, upper=0.99> phi = exp(-lambda * dt_unit);
  vector[S] r0_sub = sd_r0 * r0_raw;
  vector[N] mu = r0 + r0_sub[sid] + a_ii .* xi + a_ij .* xj;
  real<lower=2> nu = 2 + exp(log_nu_minus_two);
  vector[N] e;
  for (n in 1:N) {
    if (prev[n] == 0) e[n] = sd_ou * z_e[n];
    else {
      real rho = exp(-lambda * dt[n]);
      e[n] = rho * e[prev[n]] + sd_ou * sqrt(1 - square(rho)) * z_e[n];
    }
  }
}
model {
  r0 ~ normal(0, 1);
  r0_raw ~ normal(0, 1);
  sd_r0 ~ normal(0, 0.5);
  a_ii ~ normal(0, 0.7);
  a_ij ~ normal(0, 0.7);
  sigma ~ normal(0, 0.5);
  sd_ou ~ normal(0, 1);
  target += beta_lpdf(phi | 8, 2) + log(dt_unit) + log_lambda + log(phi);
  log_nu_minus_two ~ normal(log(3), 0.75);
  z_e ~ normal(0, 1);
  y ~ student_t(nu, mu + e, sigma);
}
generated quantities {
  vector[N] log_lik;
  real sigma_pred = sqrt(square(sigma) + square(sd_ou));
  for (n in 1:N) log_lik[n] = student_t_lpdf(y[n] | nu, mu[n] + e[n], sigma);
}
