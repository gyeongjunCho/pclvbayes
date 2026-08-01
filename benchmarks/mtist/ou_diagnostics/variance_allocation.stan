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
  real<lower=0> s_total;
  real<lower=0, upper=1> omega;
  real<lower=log(-log(0.99) / dt_unit),
       upper=log(-log(0.01) / dt_unit)> log_lambda;
  real log_nu_minus_two;
  vector[N] z_e;
}
transformed parameters {
  real<lower=0> sigma = s_total * sqrt(1 - omega);
  real<lower=0> sd_ou = s_total * sqrt(omega);
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

  // Exact transformation of independent canonical half-normal priors:
  // sigma ~ half-normal(0, 0.5), sd_ou ~ half-normal(0, 1).
  target += normal_lpdf(sigma | 0, 0.5);
  target += normal_lpdf(sd_ou | 0, 1);
  target += log(s_total) - log(2) - 0.5 * log(omega)
            - 0.5 * log1m(omega);

  target += beta_lpdf(phi | 8, 2) + log(dt_unit) + log_lambda + log(phi);
  log_nu_minus_two ~ normal(log(3), 0.75);
  z_e ~ normal(0, 1);
  y ~ student_t(nu, mu + e, sigma);
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N) log_lik[n] = student_t_lpdf(y[n] | nu, mu[n] + e[n], sigma);
}
