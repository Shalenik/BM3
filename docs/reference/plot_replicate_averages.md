# Plot Average Posterior Means Across Replicates with True Trend

Creates a plot showing the average posterior mean for each fitted H
value across multiple replications, overlaid with the true trend (mu)
and simulated data.

## Usage

``` r
plot_replicate_averages(results_multi, sims_list, title = NULL)
```

## Arguments

- results_multi:

  List from
  [`consolidate_results()`](https://shalenik.github.io/BM3/reference/consolidate_results.md)
  with `by_replicate=TRUE`.

- sims_list:

  List of simulation objects (one per replicate), each containing Y, z,
  mu, and years.

- title:

  Optional plot title (default NULL).

## Value

ggplot2 object.

## Details

The function:

1.  Computes the average posterior mean for each H across all replicates

2.  Plots these averages as colored lines

3.  Overlays the true trend (mu) as a red line

4.  Shows the averaged Y (climate models) from each replicate as grey
    lines

## Examples

``` r
if (FALSE) { # \dontrun{
results <- consolidate_results(fits_list, sims_list, by_replicate = TRUE)
p <- plot_replicate_averages(results, sims_list)
print(p)
} # }
```
