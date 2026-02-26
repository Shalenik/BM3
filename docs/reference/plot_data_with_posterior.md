# Plot Data with Posterior Samples or True Trend

Create a plot of climate models with observations and either posterior
samples or true trend.

## Usage

``` r
plot_data_with_posterior(
  Y,
  z,
  posterior_samples = NULL,
  true_trend = NULL,
  title = NULL,
  years = NULL,
  obs_years = 1:length(z)
)
```

## Arguments

- Y:

  Matrix of climate model outputs (time x models).

- z:

  Vector of observations.

- posterior_samples:

  Optional matrix of posterior samples (samples x time).

- true_trend:

  Optional vector of true trend (e.g., mu from simulation).

- title:

  Optional plot title (default NULL).

- years:

  Optional vector of years for x-axis (default: 1:nrow(Y)).

- obs_years:

  Indices for where observations are available (default: 1:length(z)).

## Value

ggplot2 object.

## Examples

``` r
if (FALSE) { # \dontrun{
# With posterior samples
p <- plot_data_with_posterior(Y, z, posterior_samples, title = "Real Data Fit")

# With true trend
p <- plot_data_with_posterior(Y, z, true_trend = mu, title = "Synthetic Data")
} # }
```
