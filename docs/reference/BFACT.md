# Bayesian Latent Factor Model for Climate Projections (BFACT)

Runs a Gibbs sampler to estimate latent factor trajectories U, spline
coefficients gamma, loadings beta, observation/model precisions phiep,
and factor precisions phiet.

## Usage

``` r
BFACT(
  Y,
  z,
  T,
  T1,
  T2,
  H,
  iseed,
  J,
  nsim,
  mga = 0.1,
  mbe = 0.1,
  aep = 0.1,
  bep = 0.1,
  aet = 0.1,
  bet = 0.1,
  a1 = 3,
  a2 = 3
)
```

## Arguments

- Y:

  Numeric matrix T x m of model outputs (each column a model series).

- z:

  Numeric vector of observed anomalies for years T1:T2.

- T:

  Integer total number of time points.

- T1, T2:

  Integers giving the observation window indices (1 ≤ T1 ≤ T2 ≤ T).

- H:

  Integer number of latent factors.

- iseed:

  Integer RNG seed for reproducibility.

- J:

  Integer number of spline basis functions (natural splines via
  splines::ns).

- nsim:

  Integer number of posterior draws.

- mga:

  Prior precision (ridge) for spline coefficients gamma. Default 0.1.

- mbe:

  Prior precision (ridge) for loadings beta. Default 0.1.

- aep:

  Gamma shape for observation/model precisions phiep. Default 0.1.

- bep:

  Gamma rate for observation/model precisions phiep. Default 0.1.

- aet:

  Gamma shape for factor precisions phiet. Default 0.1.

- bet:

  Gamma rate for factor precisions phiet. Default 0.1.

- a1:

  Gamma shape for shrinkage delta_1. Default 3.

- a2:

  Gamma shape for shrinkage delta\_\>=2. Default 3.

## Value

A list with elements: X, ysim, U_samples, beta_samples, gamma_samples
