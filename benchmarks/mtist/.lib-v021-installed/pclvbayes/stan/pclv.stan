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
  real sum_dt = 0.0;
  int cnt = 0;
  vector[N] dt_ratio;

  for (n in 1:N) {
    if (prev[n] != 0 && dt[n] > 0) {
      sum_dt += dt[n];
      cnt += 1;
    }
  }

  dt_unit = (cnt > 0) ? sum_dt / cnt : 1.0;

  for (n in 1:N) {
    dt_ratio[n] = (prev[n] == 0) ? 0.0 : dt[n] / dt_unit;
  }
}

parameters {
  real r0;

  vector[S] r0_raw;
  real<lower=0> sd_r0;

  real a_ii;
  real a_ij;

  real<lower=0> sigma;

  real<lower=0> sd_ou;
  real<lower=0.01, upper=0.99> phi;

  real log_nu_minus_two;

  vector[N] z_e;
}

transformed parameters {
  real<lower=2> nu = 2 + exp(log_nu_minus_two);
}

model {
  vector[S] r0_sub;
  vector[N] mu;
  vector[N] e;
  real log_phi = log(phi);

  r0_sub = sd_r0 * r0_raw;

  mu =
    r0
    + r0_sub[sid]
    + a_ii .* xi
    + a_ij .* xj;

  for (n in 1:N) {
    if (prev[n] == 0) {
      e[n] = sd_ou * z_e[n];
    } else {
      real rho = exp(log_phi * dt_ratio[n]);
      real sd_innov = sd_ou * sqrt(1.0 - square(rho));

      e[n] =
        rho * e[prev[n]]
        + sd_innov * z_e[n];
    }
  }

  r0 ~ normal(0, 1);
  r0_raw ~ normal(0, 1);
  sd_r0 ~ normal(0, 0.5);

  a_ii ~ normal(0, 0.7);
  a_ij ~ normal(0, 0.7);

  sigma ~ normal(0, 0.5);

  sd_ou ~ normal(0, 1);
  phi ~ beta(8, 2);

  log_nu_minus_two ~ normal(log(3), 0.75);

  z_e ~ normal(0, 1);

  y ~ student_t(nu, mu + e, sigma);
}

generated quantities {
  real lambda = -log(phi) / dt_unit;
}
