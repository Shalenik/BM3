# BFACT

Bayesian Factor Analysis for Climate Trajectories (BFACT) is an R
package for modeling climate data using a Bayesian implementation of a
latent factor model.

## Installation

You can install the development version of BFACT from GitHub:

``` r

# install.packages("devtools")
devtools::install_github("Shalenik/BM3")
```

## Features

- **Bayesian Factor Analysis**: Fit BFACT models to climate data with
  customizable latent factors

- **Model Selection**: Tools for selecting optimal hyperparameters

- **Posterior Simulation**: Generate and analyze posterior samples

## Getting Started

See the vignettes for detailed examples:

- [Fitting BFACT to Real Climate
  Data](https://shalenik.github.io/BM3/articles/01-fit-real-data.md) -
  Learn how to fit the BFACT model to real temperature anomaly data
- [Generate Posterior
  Simulations](https://shalenik.github.io/BM3/articles/02-posterior-simulation.md) -
  Work with posterior samples and synthetic data generation
- [Model Selection Across Multiple
  Replications](https://shalenik.github.io/BM3/articles/03-multiple-replicates.md) -
  Evaluate model selection consistency across replications

## Quick Example

``` r

library(BFACT)

# Simulate data with H_true = 2 latent factors
sim <- simulate_data(H_true = 2, K = 20, T_obs = 50)

# Fit BFACT model with H = 2
fit <- fit_bfact_model(sim, H = 2, nsim = 5000)

# Extract posterior mean curve
post_mean <- posterior_mean_curve(fit)

# Plot results
plot_data_with_posterior(sim$Y, sim$z, 
                        posterior_samples = sample_posterior(fit, ncurves = 50),
                        title = "BFACT Fit")
```

## Main Functions

- [`BFACT()`](https://shalenik.github.io/BM3/reference/BFACT.md) - Core
  Bayesian factor analysis function
- [`fit_bfact_model()`](https://shalenik.github.io/BM3/reference/fit_bfact_model.md) -
  Convenient wrapper for model fitting
- [`consolidate_results()`](https://shalenik.github.io/BM3/reference/consolidate_results.md) -
  Aggregate results across multiple H values
- [`find_elbow()`](https://shalenik.github.io/BM3/reference/find_elbow.md) -
  Identify optimal model complexity
- [`detect_outliers()`](https://shalenik.github.io/BM3/reference/detect_outliers.md) -
  Remove outlier climate models
- [`simulate_data()`](https://shalenik.github.io/BM3/reference/simulate_data.md) -
  Generate synthetic data for testing

## Citation

If you use BFACT in your research, please cite:

    [Citation information to be added]

## License

\[License information to be added\]
