# Fit BFACT Model for a Single H Value

Fits a Bayesian Flexible Adaptation to Climate Trends (BFACT) latent
factor model to data. This is the core fitting function that runs MCMC
to estimate posterior distributions of latent factors and model weights.

## Usage

``` r
fit_bfact_model(sim, H, nsim = 10000, iseed = 123)
```

## Arguments

- sim:

  A list containing:

  - Y: matrix of climate model outputs (n_time × K)

  - z: vector of observations (n_time)

  - years: vector of years

- H:

  Integer. Number of latent factors to fit.

- nsim:

  Integer. Number of MCMC iterations (default 10000).

- iseed:

  Integer. Random seed for reproducibility.

## Value

A list containing:

- H: fitted H value

- posterior: list with samples for beta, U, phi, sigma2

- metrics: list with RMSE and other diagnostics

- settings: list of input parameters

## Details

The BFACT model represents observations as: z_t = mu + U %\*% beta_t +
epsilon_t

where U is a basis matrix (n_time × H), beta_t are time-varying weights,
and epsilon_t is Gaussian error.

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- simulate_data(H_true = 2, K = 20)
fit <- fit_bfact_model(sim, H = 2, nsim = 5000)
} # }
```
