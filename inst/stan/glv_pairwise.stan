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
  int<lower=0,upper=1> use_student_t;
  real<lower=2> nu_fixed;     // ← 고정 자유도 (예: 4 또는 7)
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
  // random intercepts (non-centered)
  real r0;

  // coefficients
  real a_ii;
  real a_ij;

  // noise scales
  real<lower=0> sigma;       // white noise (obs)
  real<lower=0> sd_ou;       // OU stationary SD  ⟵  (CHANGED)

  // AR persistence at unit dt
  real<lower=0.01, upper=0.99> phi;

  // whitened innovations
  vector[N] z_e;
}

transformed parameters {
  vector[N] mu = r0 + a_ii .* xi + a_ij .* xj;

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

  a_ii      ~ normal(0, 0.7);
  a_ij      ~ normal(0, 0.7);

  sigma     ~ normal(0, 0.5);
  sd_ou     ~ normal(0, 1.0);
  phi       ~ beta(8, 2); // mean ~0.8

  z_e       ~ normal(0, 1);

  // likelihood
  if (use_student_t == 1)
    y ~ student_t(nu_fixed, mu + e, sigma);
  else
    y ~ normal(mu + e, sigma);
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
    log_lik[n] = normal_lpdf(y[n] | mu[n] + e[n], sigma);
    y_rep[n]   = normal_rng(mu[n] + e[n], sigma);
  }
}

