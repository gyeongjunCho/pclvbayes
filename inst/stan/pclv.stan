// glv_pairwise_ar1_orthogonalized.stan
data {
  int<lower=1> N;
  vector[N] y;          // standardized ΔALR_i/Δt
  vector[N] xi;         // standardized ALR_i (lagged)
  vector[N] xj;         // standardized ALR_j (lagged)
  int<lower=1> S;
  array[N] int<lower=1, upper=S> sid;
  array[N] int<lower=0, upper=N> prev;
  vector[N] dt;         // 0 if first within subject, else > 0 (preprocessed)
}

transformed data {
  real dt_unit;
  real sum_dt = 0.0;
  int cnt = 0;
  for (n in 1:N) {
    if (prev[n] != 0 && dt[n] > 0) {
      sum_dt = sum_dt + dt[n];
      cnt = cnt + 1;
    }
  }
  dt_unit = (cnt > 0) ? sum_dt / cnt : 1.0;
}

parameters {
  // global intercept
  real r0;

  // subject random intercepts (non-centered)
  vector[S] r0_raw;
  real<lower=0> sd_r0;

  // coefficients
  real a_ii;
  real a_ij;

  // noise scales
  real<lower=0> sigma;       // white noise (obs)
  real<lower=0> sd_ou;       // OU stationary SD  ⟵  (CHANGED)

  // AR persistence at unit dt
  real<lower=0.01, upper=0.99> phi;

  // shared Student-t degrees of freedom
  real log_nu_minus_two;

  // whitened innovations
  vector[N] z_e;
}

transformed parameters {
  // vector[N] mu = r0 + a_ii .* xi + a_ij .* xj;
  vector[S] r0_sub = sd_r0 * r0_raw;
  vector[N] mu = r0 + r0_sub[sid] + a_ii .* xi + a_ij .* xj;
  real<lower=2> nu = 2 + exp(log_nu_minus_two);

  // build e via whitening in stationary-SD parameterization
  vector[N] e;
  for (n in 1:N) {
    if (prev[n] == 0) {
      e[n] = sd_ou * z_e[n];
    } else {
      real rho = pow(phi, dt[n] / dt_unit);
      real sd_innov = sd_ou * sqrt(1.0 - rho * rho);
      e[n] = rho * e[prev[n]] + sd_innov * z_e[n];
    }
  }
}

model {
  r0        ~ normal(0, 1);

  r0_raw    ~ normal(0, 1);
  sd_r0     ~ normal(0, 0.5);

  a_ii      ~ normal(0, 0.7);
  a_ij      ~ normal(0, 0.7);

  sigma     ~ normal(0, 0.5);
  sd_ou     ~ normal(0, 1.0);
  phi       ~ beta(8, 2); // mean ~0.8
  log_nu_minus_two ~ normal(log(3), 0.75);

  z_e       ~ normal(0, 1);

  // likelihood
  y ~ student_t(nu, mu + e, sigma);
}

generated quantities {
  vector[N] log_lik;
  vector[N] y_rep;

  // 예전 파이프라인 호환: lambda, k, sigma_ou 복원
  real lambda = (-log(phi)) / dt_unit;
  real k_gq   = fmax(2.0 * lambda, 1e-12);   //
  real sigma_ou = sd_ou * sqrt(k_gq);        //
  real tau_r    = 0;

  // 새 정의에 맞는 예측 표준편차
  real sigma_pred = sqrt(square(sigma) + square(sd_ou));

  for (n in 1:N) {
    log_lik[n] = student_t_lpdf(y[n] | nu, mu[n] + e[n], sigma);
    y_rep[n] = student_t_rng(nu, mu[n] + e[n], sigma);
  }
}

