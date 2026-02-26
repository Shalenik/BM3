# Plot Variance of Posterior Means Across Replicates

Creates a plot showing how the variance of posterior means across
replicates changes with H. Lower variance indicates more consistent
predictions.

## Usage

``` r
plot_posterior_variance(results_multi, title = NULL)
```

## Arguments

- results_multi:

  List from
  [`consolidate_results()`](https://shalenik.github.io/BM3/reference/consolidate_results.md)
  with `by_replicate=TRUE`.

- title:

  Optional plot title (default NULL).

## Value

ggplot2 object.

## Details

For each H value, computes the variance of posterior means across
replicates at each time point, then averages these variances over time.
Lower variance indicates that the fitted H produces more consistent
predictions across different random realizations of the data.

## Examples

``` r
if (FALSE) { # \dontrun{
results <- consolidate_results(fits_list, sims_list, by_replicate = TRUE)
p <- plot_posterior_variance(results)
print(p)
} # }
```
