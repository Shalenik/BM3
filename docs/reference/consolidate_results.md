# Consolidate Fitting Results

Create summary statistics from fitted BFACT models across multiple H
values. Combines individual fit objects into posterior summaries and
model comparison metrics.

## Usage

``` r
consolidate_results(fits, sim, by_replicate = FALSE)
```

## Arguments

- fits:

  List of fit objects (output from
  [`fit_bfact_model()`](https://shalenik.github.io/BM3/reference/fit_bfact_model.md)),
  or a list of such lists for multiple replications. See details.

- sim:

  List containing simulation data (Y, z, years), or list of such lists
  for multiple replications.

- by_replicate:

  Logical. If TRUE, compute per-replicate statistics (default FALSE). If
  fits and sim are lists of lists, this will aggregate across
  replications.

## Value

A list containing:

- timewise_summaries: tibble with posterior means/CIs per H

- model_comparison: tibble comparing metrics across H values

- per_replicate_stats: per-replicate RMSE if by_replicate=TRUE

- replicate_comparison: tibble with mean/sd RMSE by H across
  replications (if multiple reps)

## Details

Consolidates multiple H fits into a single summary table with columns:

- H_fit: fitted H value

- time: time index

- post_mean: posterior mean prediction

- post_ci_lower, post_ci_upper: 95% credible interval

- se: squared error at each time point

- rmse: overall RMSE for that H

For multiple replications, pass:

- fits: list of lists, where each inner list is fits for one replicate

- sim: list of simulation lists, one per replicate

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- simulate_data(H_true = 2, K = 20)
fit_H2 <- fit_bfact_model(sim, H = 2)
fit_H3 <- fit_bfact_model(sim, H = 3)

results <- consolidate_results(
    fits = list(fit_H2, fit_H3),
    sim = sim
)
head(results$timewise_summaries)
} # }
```
