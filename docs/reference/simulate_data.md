# Simulate Climate Data

Generate synthetic climate model outputs and observations for testing.
Creates a latent factor model with H latent factors, K climate models,
and realistic temporal structure.

## Usage

``` r
simulate_data(
  H_true = 2,
  K = 20,
  years_start = 1850,
  years_end = 2100,
  obs_end = 2022,
  seed = NULL
)
```

## Arguments

- H_true:

  Integer. True number of latent factors.

- K:

  Integer. Number of climate models (default 20).

- years_start:

  Integer. Start year (default 1850).

- years_end:

  Integer. End year (default 2100).

- obs_end:

  Integer. Last year of observations (default 2022).

- seed:

  Integer. Random seed for reproducibility.

## Value

A list containing:

- Y: matrix of climate model outputs (n_time × K)

- z: vector of observations (n_time)

- years: vector of years

- y_obs_full: full observations including future

- K: number of models

- H_true: true H value

## Details

The generative model: Y_kt ~ N(mu_k + U %*% beta_k(t), tau2_k) z_t ~
N(mu + U %*% beta(t), sigma2)

where U is a basis matrix (n_time × H) with smooth temporal patterns.

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- simulate_data(H_true = 2, K = 20, seed = 123)
dim(sim$Y) # [1] 251 20 (years 1850-2100)
length(sim$z) # [1] 173 (years 1850-2022)
} # }
```
