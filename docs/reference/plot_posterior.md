# Plot Posterior Distribution

Create a plot of posterior means and credible intervals.

## Usage

``` r
plot_posterior(results, H, data = NULL, title = NULL)
```

## Arguments

- results:

  List from
  [`consolidate_results()`](https://shalenik.github.io/BM3/reference/consolidate_results.md).

- H:

  Integer. Which H value to plot.

- data:

  Optional. Observations to overlay (default NULL).

- title:

  Optional plot title (default NULL).

## Value

ggplot2 object.

## Examples

``` r
if (FALSE) { # \dontrun{
plot_posterior(results, H = 2, title = "My Title")
} # }
```
